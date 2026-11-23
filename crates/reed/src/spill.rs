//! Exclusive temporary spill files (doc_v20 SEC-1, doc_v25 V25-09).
//!
//! Every feature that spills to disk (stdin replay, convert's rendered rows)
//! creates its file here. The name stays predictable (`pid` + a process
//! counter) so operators can recognise spill files, but creation uses
//! `create_new` (O_EXCL): an existing file *or symlink* at the candidate
//! name fails the open instead of being followed or truncated. On a shared
//! host a pre-planted symlink therefore cannot redirect the spill into a
//! victim file. A collision retries the next counter value, bounded.

use std::fs::{File, OpenOptions};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

/// Process-wide spill-name counter: one sequence per process, starting at 0.
static SPILL_COUNTER: AtomicU64 = AtomicU64::new(0);

/// How many candidate names one spill may burn through before giving up.
const SPILL_ATTEMPTS: usize = 16;

/// How old a leftover spill must be before the startup sweep may remove it.
///
/// The same window as the output temporaries: young files may belong to a
/// writer that has not touched them yet, and age is the only portable signal.
const STALE_SPILL_SECS: u64 = 60;

/// Directories already swept in this process. The sweep is a **startup**
/// action (doc_v30 R10): it runs once, before the first spill is created in a
/// directory, and never again — a per-create walk of `/tmp` would cost more
/// than the leaks it prevents.
static SWEPT_DIRS: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());

/// Where a spill's candidate names come from.
///
/// Production uses one process-wide counter; a test constructs its own so it
/// never predicts another test's next value (doc_v26 V26-12: the SEC-1 tests
/// read `SPILL_COUNTER` and planted the next 16 names, but other tests
/// incremented the same counter between the read and the call, so the planted
/// names were stale and the call succeeded instead of failing).
pub(crate) trait NameSource {
    /// The next candidate counter value.
    fn next(&self) -> u64;
}

/// The process-wide source (production).
pub(crate) struct GlobalCounter;

impl NameSource for GlobalCounter {
    fn next(&self) -> u64 {
        SPILL_COUNTER.fetch_add(1, Ordering::Relaxed)
    }
}

/// A source a single test owns outright: a private starting point and its own
/// counter, so nothing else can advance it.
#[cfg(test)]
pub(crate) struct TestCounter {
    next: AtomicU64,
}

#[cfg(test)]
impl TestCounter {
    /// Starts at `first`.
    pub(crate) fn starting_at(first: u64) -> Self {
        Self {
            next: AtomicU64::new(first),
        }
    }
}

#[cfg(test)]
impl NameSource for TestCounter {
    fn next(&self) -> u64 {
        self.next.fetch_add(1, Ordering::Relaxed)
    }
}

/// The candidate spill path `counter` generations on in `dir`.
fn spill_path(dir: &Path, tag: &str, ext: &str, counter: u64) -> PathBuf {
    dir.join(format!("reed-{tag}-{}-{counter}.{ext}", std::process::id()))
}

/// Creates the next spill file in `dir`, exclusively.
///
/// # Errors
///
/// The underlying I/O error, or `AlreadyExists` after every candidate name
/// was taken (nothing is clobbered in that case).
pub(crate) fn create_spill_file(
    dir: &Path,
    tag: &str,
    ext: &str,
) -> std::io::Result<(File, PathBuf)> {
    create_spill_file_from(dir, tag, ext, &GlobalCounter)
}

/// [`create_spill_file`] against an explicit name source (tests own theirs,
/// doc_v26 V26-12).
///
/// # Errors
///
/// As [`create_spill_file`].
pub(crate) fn create_spill_file_from(
    dir: &Path,
    tag: &str,
    ext: &str,
    names: &dyn NameSource,
) -> std::io::Result<(File, PathBuf)> {
    for _ in 0..SPILL_ATTEMPTS {
        let path = spill_path(dir, tag, ext, names.next());
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((file, path)),
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    Err(std::io::Error::new(
        std::io::ErrorKind::AlreadyExists,
        format!(
            "cannot create a spill file in {}: all {} candidate names are taken",
            dir.display(),
            SPILL_ATTEMPTS
        ),
    ))
}

/// The owner pid encoded in a spill name, if the name has the exact
/// `reed-<tag>-<pid>-<counter>.<ext>` shape.
///
/// Anchored on purpose: `reed-notes.txt`, a missing counter, a non-numeric
/// field or a name without an extension are user files and are never swept.
/// The tag may contain hyphens; the pid and counter are parsed from the
/// right, so a future tag with a dash does not change ownership.
fn spill_owner(name: &str) -> Option<u32> {
    let rest = name.strip_prefix("reed-")?;
    let (rest, ext) = rest.rsplit_once('.')?;
    if ext.is_empty() || !ext.bytes().all(|b| b.is_ascii_alphanumeric()) {
        return None;
    }
    let (rest, counter) = rest.rsplit_once('-')?;
    if counter.is_empty() || !counter.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let (tag, pid) = rest.rsplit_once('-')?;
    if tag.is_empty() || pid.is_empty() || !pid.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    pid.parse().ok()
}

/// Whether `pid` still names a process. Only Linux exposes a race-free-enough
/// answer through `/proc`; elsewhere the conservative "alive" keeps the file.
#[cfg(target_os = "linux")]
fn process_is_alive(pid: u32) -> bool {
    Path::new("/proc").join(pid.to_string()).exists()
}

#[cfg(not(target_os = "linux"))]
fn process_is_alive(_pid: u32) -> bool {
    true
}

/// Removes spill files this process can prove orphaned, keeping every name in
/// the shared live table.
///
/// A spill is swept only when it is older than [`STALE_SPILL_SECS`] **and**
/// either names this pid (a crashed earlier run of this process, since a live
/// one is registered) or names a pid that no longer exists. Another live
/// process's spill is reported, never deleted: a reused pid cannot be
/// distinguished from a live writer from the name alone (the same rule
/// [`reed_io::write::atomic`] applies to output temporaries).
fn sweep_spill_dir(dir: &Path, live: &[PathBuf]) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    let now = std::time::SystemTime::now();
    let own_pid = std::process::id();
    for entry in entries.flatten() {
        let name = entry.file_name();
        let Some(pid) = spill_owner(&name.to_string_lossy()) else {
            continue;
        };
        let path = entry.path();
        // A registered spill is active, however quiet it has been.
        if live.contains(&path) {
            continue;
        }
        let old = entry
            .metadata()
            .and_then(|m| m.modified())
            .ok()
            .and_then(|m| now.duration_since(m).ok())
            .is_some_and(|age| age.as_secs() >= STALE_SPILL_SECS);
        if !old {
            continue;
        }
        if pid == own_pid || !process_is_alive(pid) {
            // unlink(2) never follows a symlink; a matching user file is
            // protected by the anchored name shape above.
            let _ = std::fs::remove_file(&path);
        } else {
            tracing::warn!(
                path = %path.display(),
                pid,
                "stale spill file remains; remove it manually"
            );
        }
    }
}

/// The startup sweep: once per directory, before the first spill there.
fn sweep_stale_spills(dir: &Path) {
    {
        let Ok(mut swept) = SWEPT_DIRS.lock() else {
            return;
        };
        if swept.iter().any(|d| d == dir) {
            return;
        }
        swept.push(dir.to_path_buf());
    }
    sweep_spill_dir(dir, &reed_io::write::atomic::live_temp_paths());
}

/// A spill file this process created; removed on drop.
///
/// The guard owns an open handle as well as the name (doc_v30 R10): the
/// inode created with `O_EXCL` stays pinned for the guard's lifetime, so a
/// later writer reopens *this* file through [`Self::try_clone`] instead of
/// trusting the path. Readers still open by path (they need an independent
/// offset); the O_EXCL creation and unlink-on-drop cover the symlink case.
#[derive(Debug)]
pub(crate) struct TempSpill {
    /// A handle on the exclusively-created inode.
    file: File,
    path: PathBuf,
}

impl TempSpill {
    /// Creates an exclusively-named spill file (`None` = the platform temp
    /// directory). The caller writes through the returned handle; the
    /// [`TempSpill`] owns the path and removes it on drop.
    ///
    /// # Errors
    ///
    /// [`create_spill_file`]'s errors.
    pub(crate) fn create(
        dir: Option<&Path>,
        tag: &str,
        ext: &str,
    ) -> std::io::Result<(File, Self)> {
        let dir = dir.map_or_else(std::env::temp_dir, Path::to_path_buf);
        sweep_stale_spills(&dir);
        let (file, path) = create_spill_file(&dir, tag, ext)?;
        // The shared live table is `AtomicFile`'s (doc_v30 R10): a sweep in
        // either subsystem can now see that this spill is still being written.
        reed_io::write::atomic::register_live_temp(&path);
        // The clone shares the file offset with the returned handle, which is
        // what append-only spill writers want (`try_clone` + `write_all`).
        let owned = match file.try_clone() {
            Ok(owned) => owned,
            Err(e) => {
                // No guard exists yet to clean up on drop.
                reed_io::write::atomic::unregister_live_temp(&path);
                let _ = std::fs::remove_file(&path);
                return Err(e);
            }
        };
        Ok((file, Self { file: owned, path }))
    }

    /// A new handle on the pinned inode, for append writers.
    ///
    /// # Errors
    ///
    /// The platform's descriptor-limit error.
    pub(crate) fn try_clone(&self) -> std::io::Result<File> {
        self.file.try_clone()
    }

    /// The file's path (for reopen/append/read).
    #[must_use]
    pub(crate) fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for TempSpill {
    fn drop(&mut self) {
        // No longer active: release the shared-table entry before unlinking,
        // so a concurrent sweep never mistakes the name for a live spill.
        reed_io::write::atomic::unregister_live_temp(&self.path);
        // `remove_file` is unlink(2): it never follows a symlink, so even a
        // name swapped after creation removes the link entry, not a victim.
        let _ = std::fs::remove_file(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sec1_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("reed-sec1-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: scratch dir");
        dir
    }

    /// Ages `path` past the sweep window.
    fn age(path: &Path, secs: u64) {
        let file = OpenOptions::new()
            .write(true)
            .open(path)
            .expect("invariant: open fixture");
        file.set_modified(std::time::SystemTime::now() - std::time::Duration::from_secs(secs))
            .expect("invariant: age fixture");
    }

    fn live_contains(path: &Path) -> bool {
        reed_io::write::atomic::live_temp_paths().contains(&path.to_path_buf())
    }

    /// SEC-1: a pre-created file at the next name must not be truncated.
    ///
    /// doc_v26 V26-12: the test drives a private [`TestCounter`] instead of
    /// the process-wide counter, so no other test can advance the name
    /// between the plant and the call. The old version read
    /// `SPILL_COUNTER` under a test-only lock that the remaining spill tests
    /// did not all take, which made this test order-dependent.
    #[test]
    fn precreated_file_is_not_truncated() {
        let dir = sec1_dir("file");
        let names = TestCounter::starting_at(0);
        let victim = spill_path(&dir, "stdin", "vcf", 0);
        std::fs::write(&victim, b"victim bytes").expect("invariant: plant victim");
        let (_spill, path) =
            create_spill_file_from(&dir, "stdin", "vcf", &names).expect("invariant: create");
        assert_eq!(
            std::fs::read(&victim).expect("invariant: re-read victim"),
            b"victim bytes".as_slice()
        );
        assert_ne!(path, victim, "the spill must pick another name");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// SEC-1: a pre-planted symlink must not be followed.
    #[cfg(unix)]
    #[test]
    fn precreated_symlink_is_not_followed() {
        let dir = sec1_dir("symlink");
        let names = TestCounter::starting_at(0);
        let victim = dir.join("victim.vcf");
        std::fs::write(&victim, b"victim bytes").expect("invariant: plant victim");
        let predicted = spill_path(&dir, "convert", "spill", 0);
        std::os::unix::fs::symlink(&victim, &predicted).expect("invariant: plant symlink");
        let (_spill, path) =
            create_spill_file_from(&dir, "convert", "spill", &names).expect("invariant: create");
        assert_eq!(
            std::fs::read(&victim).expect("invariant: re-read victim"),
            b"victim bytes".as_slice(),
            "the symlink's target must be untouched"
        );
        assert_ne!(path, predicted, "the spill must not write through the link");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// SEC-1: exhausted names error without clobbering any planted file.
    ///
    /// doc_v26 V26-12 reproduced this as flaky under the default parallel
    /// harness: it predicted the process-wide counter's next 16 values, and
    /// [`Self::temp_spill_removes_its_file_on_drop`] advanced the counter
    /// between the plant and the call. It now owns its name source, so the
    /// prediction is exact and the test needs no lock.
    #[test]
    fn exhausted_candidates_error_without_clobbering() {
        let dir = sec1_dir("exhausted");
        let first = 0;
        let names = TestCounter::starting_at(first);
        for n in first..first + SPILL_ATTEMPTS as u64 {
            std::fs::write(spill_path(&dir, "stdin", "vcf", n), b"victim")
                .expect("invariant: plant candidate");
        }
        let err =
            create_spill_file_from(&dir, "stdin", "vcf", &names).expect_err("every name is taken");
        assert_eq!(err.kind(), std::io::ErrorKind::AlreadyExists);
        for n in first..first + SPILL_ATTEMPTS as u64 {
            assert_eq!(
                std::fs::read(spill_path(&dir, "stdin", "vcf", n)).expect("invariant: re-read"),
                b"victim".as_slice(),
                "candidate {n} must not have been touched"
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The lock that used to serialize these tests is gone: each owns its
    /// names, so the whole module runs in parallel with anything else.
    #[test]
    fn the_name_source_is_private_per_test() {
        let a = TestCounter::starting_at(0);
        let b = TestCounter::starting_at(0);
        assert_eq!(a.next(), b.next());
        assert_eq!(a.next(), 1, "each source advances independently");
        assert_eq!(b.next(), 1);
    }

    /// The owner removes the file on drop, success and failure alike.
    #[test]
    fn temp_spill_removes_its_file_on_drop() {
        let dir = sec1_dir("drop");
        let (mut file, spill) =
            TempSpill::create(Some(&dir), "convert", "spill").expect("invariant: create");
        std::io::Write::write_all(&mut file, b"rows").expect("invariant: write");
        let path = spill.path().to_path_buf();
        assert!(path.exists());
        drop(spill);
        assert!(!path.exists(), "the spill file must be removed on drop");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// doc_v30 R10: spills and output temporaries share one live table, so a
    /// sweep in either subsystem sees the other's active files.
    #[test]
    fn a_spill_is_registered_in_the_shared_live_table() {
        let dir = sec1_dir("shared-table");
        let (file, spill) =
            TempSpill::create(Some(&dir), "convert", "spill").expect("invariant: spill");
        let spill_path = spill.path().to_path_buf();
        let (atomic, atomic_file) =
            reed_io::write::atomic::AtomicFile::create(&dir.join("out.txt"))
                .expect("invariant: atomic");
        let atomic_path = atomic.temp_path().to_path_buf();

        let live = reed_io::write::atomic::live_temp_paths();
        assert!(
            live.contains(&spill_path),
            "the spill must be in AtomicFile's live table"
        );
        assert!(live.contains(&atomic_path));

        drop(spill);
        assert!(!live_contains(&spill_path), "drop deregisters the spill");
        assert!(!spill_path.exists());
        drop(atomic);
        drop((file, atomic_file));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A registered spill survives a sweep even when its mtime is hours old:
    /// age cannot tell "quiet writer" from "crashed writer", ownership can.
    #[test]
    fn a_live_spill_survives_the_sweep() {
        let dir = sec1_dir("live-sweep");
        let (mut file, spill) =
            TempSpill::create(Some(&dir), "convert", "spill").expect("invariant: create");
        std::io::Write::write_all(&mut file, b"rows").expect("invariant: write");
        let path = spill.path().to_path_buf();
        age(&path, 3600);

        sweep_spill_dir(&dir, &reed_io::write::atomic::live_temp_paths());
        assert!(
            path.exists(),
            "a registered spill is active, however old its mtime"
        );
        drop(spill);
        assert!(!path.exists());
        drop(file);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// An unregistered, stale, same-pid name is a crashed run's leftover.
    #[test]
    fn an_unregistered_stale_spill_is_swept() {
        let dir = sec1_dir("stale-sweep");
        let (file, path) = create_spill_file(&dir, "convert", "spill").expect("invariant: create");
        drop(file);
        age(&path, 3600);

        sweep_spill_dir(&dir, &[]);
        assert!(!path.exists(), "an orphaned same-pid spill is removed");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The startup sweep runs once per directory, before the first spill
    /// there — not on every create (a per-create `/tmp` walk is not worth it).
    #[test]
    fn the_startup_sweep_runs_once_per_directory() {
        let dir = sec1_dir("startup-once");
        let path = dir.join(format!("reed-convert-{}-0.spill", std::process::id()));
        std::fs::write(&path, b"x").expect("invariant: write");
        age(&path, 3600);
        sweep_stale_spills(&dir);
        assert!(!path.exists(), "the first sweep removes the leftover");

        std::fs::write(&path, b"x").expect("invariant: write again");
        age(&path, 3600);
        sweep_stale_spills(&dir);
        assert!(
            path.exists(),
            "the second call is a no-op: the sweep is a startup action"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A dead owner's spill can never be written again: it is orphaned.
    #[cfg(target_os = "linux")]
    #[test]
    fn a_dead_owners_spill_is_swept() {
        let dir = sec1_dir("dead-pid");
        let mut child = std::process::Command::new("true")
            .spawn()
            .expect("invariant: spawn true");
        let pid = child.id();
        assert!(child.wait().expect("invariant: wait").success());

        let path = dir.join(format!("reed-convert-{pid}-0.spill"));
        std::fs::write(&path, b"x").expect("invariant: write");
        age(&path, 3600);
        sweep_spill_dir(&dir, &[]);
        assert!(!path.exists(), "a dead pid's spill is orphaned");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A live owner's spill is never deleted, even when it looks stale: the
    /// pid may be alive and merely between writes.
    #[cfg(unix)]
    #[test]
    fn a_live_owners_spill_is_kept() {
        let dir = sec1_dir("live-pid");
        let mut child = std::process::Command::new("sleep")
            .arg("30")
            .spawn()
            .expect("invariant: spawn sleep");
        let pid = child.id();

        let path = dir.join(format!("reed-convert-{pid}-0.spill"));
        std::fs::write(&path, b"x").expect("invariant: write");
        age(&path, 3600);
        sweep_spill_dir(&dir, &[]);
        assert!(path.exists(), "a live process's spill is never deleted");
        let _ = child.kill();
        let _ = child.wait();
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Only the exact `reed-<tag>-<pid>-<counter>.<ext>` shape is a spill; a
    /// user file that happens to share the prefix is never touched.
    #[test]
    fn spill_like_user_files_are_never_swept() {
        let dir = sec1_dir("user-files");
        let names = [
            "reed-notes.txt",
            "reed-convert-abc-0.spill",
            "reed-convert-1.spill",
            "reed-convert-1-x.spill",
            "notreed-convert-1-0.spill",
            "reed-convert-1-0",
        ];
        let paths: Vec<PathBuf> = names.iter().map(|name| dir.join(name)).collect();
        for path in &paths {
            std::fs::write(path, b"x").expect("invariant: write");
            age(path, 3600);
        }
        sweep_spill_dir(&dir, &[]);
        for path in &paths {
            assert!(path.exists(), "{} must be kept", path.display());
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}

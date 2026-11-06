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
use std::sync::atomic::{AtomicU64, Ordering};

/// Process-wide spill-name counter: one sequence per process, starting at 0.
static SPILL_COUNTER: AtomicU64 = AtomicU64::new(0);

/// How many candidate names one spill may burn through before giving up.
const SPILL_ATTEMPTS: usize = 16;

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
    dir.join(format!(
        "reed-{tag}-{}-{counter}.{ext}",
        std::process::id()
    ))
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

/// A spill file this process created; removed on drop.
#[derive(Debug)]
pub(crate) struct TempSpill {
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
        let (file, path) = create_spill_file(&dir, tag, ext)?;
        Ok((file, Self { path }))
    }

    /// The file's path (for reopen/append/read).
    #[must_use]
    pub(crate) fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for TempSpill {
    fn drop(&mut self) {
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
        let err = create_spill_file_from(&dir, "stdin", "vcf", &names)
            .expect_err("every name is taken");
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
}

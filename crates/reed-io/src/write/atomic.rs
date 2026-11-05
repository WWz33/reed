//! Atomic output files (doc_v26 V26-08 / V26-17).
//!
//! Every command used to `File::create` its output directly, truncating the
//! user's previous result before a single byte of the new one was written. A
//! mid-run failure, a process interrupt, or a later output failing left a
//! truncated or half-written file behind — and for multi-file reports, a mix
//! of old and new.
//!
//! [`AtomicFile`] writes to a sibling temporary name in the same directory
//! (created with `create_new`, so a pre-planted file or symlink is never
//! followed or truncated) and only on [`AtomicFile::commit`] syncs and
//! `rename`s it over the target. A rename within one directory is atomic on
//! every supported filesystem, so a reader of the final path sees either the
//! complete old file or the complete new one, never a partial write.
//!
//! # Contract
//!
//! The caller owns the [`File`] returned by [`AtomicFile::create`] and must
//! flush it (directly or by dropping its `BufWriter`) before calling
//! [`AtomicFile::commit`]; `commit` then syncs the data and renames. Nothing
//! is published until `commit`, and a dropped value removes the temporary.
//!
//! # Promise
//!
//! This is **per-file** atomicity. A command that writes several files
//! (stats reports, `popg`'s windowed tables, PCA's three outputs, DAPC's four)
//! commits them one by one, so a failure in the middle leaves the earlier
//! files new and the later ones old. That is stated in the CLI documentation;
//! a group-wide transaction would need a staging directory and rollback and is
//! deliberately not implemented.

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// Process-wide counter for temporary names (one sequence per process).
static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// How many candidate names one output may burn through before giving up.
const TEMP_ATTEMPTS: usize = 16;

/// Removes temporaries left for `target` by a **crashed run of this same
/// process** — the one case `Drop` cannot handle.
///
/// The previous version swept every file matching `<name>.reed-tmp-*` older
/// than 60 s, which was unsafe in two ways (review item 3):
///
/// * it did not verify the trailing `-<pid>-<counter>` shape, so a user file
///   that merely shared the prefix would be deleted;
/// * it judged by age, not by liveness, so a *concurrent* run of the same
///   command that had been quiet for a minute would have its temporary
///   removed mid-write, and its commit would then fail.
///
/// It now removes only names whose pid field names this very process
/// (`<name>.reed-tmp-<our pid>-<digits>`) and that are older than
/// [`STALE_TEMP_SECS`]. Another process's temporary is never touched — a
/// crashed *other* process's leftovers stay, which is the conservative
/// choice: stale bytes are recoverable, deleted data is not. There is no
/// portable, race-free way to prove another pid is dead from a filename, so
/// this implementation does not pretend to.
fn sweep_stale_temps(target: &Path) {
    let Some(parent) = target.parent() else { return };
    let name = target
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    let own_pid = std::process::id();
    let Ok(entries) = std::fs::read_dir(if parent.as_os_str().is_empty() {
        Path::new(".")
    } else {
        parent
    }) else {
        return;
    };
    let now = std::time::SystemTime::now();
    for entry in entries.flatten() {
        let file_name = entry.file_name().to_string_lossy().into_owned();
        // `<target name>.reed-tmp-<pid>-<counter>`, fully anchored: a suffix
        // that is not exactly digits is a user file and is left alone.
        let Some(suffix) = file_name
            .strip_prefix(&format!("{name}.reed-tmp-"))
        else {
            continue;
        };
        let Some((pid_text, counter_text)) = suffix.split_once('-') else {
            continue;
        };
        if pid_text.parse::<u32>().ok() != Some(own_pid)
            || counter_text.parse::<u64>().is_err()
        {
            continue;
        }
        let old = entry
            .metadata()
            .and_then(|m| m.modified())
            .ok()
            .and_then(|m| now.duration_since(m).ok())
            .is_some_and(|age| age.as_secs() >= STALE_TEMP_SECS);
        if old {
            let _ = std::fs::remove_file(entry.path());
        }
    }
}

/// How old a leftover temporary must be before this process removes its own.
const STALE_TEMP_SECS: u64 = 60;

/// The candidate temporary path for `final_path`, `counter` generations on.
///
/// The name is a plain (non-dot) sibling: `.<file>` names are avoided because
/// a CLI test asserts nothing dot-prefixed is left in the working directory,
/// and because a leading dot would hide a long run's temporary in `ls`. It
/// does not look like the spill naming (`reed-<tag>-...`) either.
#[must_use]
pub fn temp_path(final_path: &Path, counter: u64) -> PathBuf {
    let name = final_path
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    final_path.with_file_name(format!(
        "{name}.reed-tmp-{}-{counter}",
        std::process::id()
    ))
}

/// An output file that becomes visible only on commit.
///
/// # Examples
///
/// ```no_run
/// # use std::io::Write as _;
/// # use reed_io::write::atomic::AtomicFile;
/// # fn demo() -> std::io::Result<()> {
/// let (out, mut file) = AtomicFile::create(std::path::Path::new("out.txt"))?;
/// file.write_all(b"contents")?;
/// file.flush()?;
/// drop(file);
/// out.commit()?;          // out.txt now holds "contents"
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct AtomicFile {
    temp: PathBuf,
    target: PathBuf,
    committed: bool,
}

impl AtomicFile {
    /// Creates the temporary sibling of `target` and hands out its handle.
    ///
    /// The parent directory is created when missing (matching the old
    /// `create_dir_all` calls for report prefixes); the temporary file itself
    /// uses `create_new` with bounded retries.
    ///
    /// # Errors
    ///
    /// I/O failures creating the directory or the temporary file, and
    /// `AlreadyExists` when every candidate name is taken (nothing was
    /// clobbered in that case).
    pub fn create(target: &Path) -> std::io::Result<(Self, File)> {
        Self::create_with_parents(target, false)
    }

    /// [`AtomicFile::create`], optionally creating missing parent
    /// directories.
    ///
    /// The default does **not** create them: a mistyped `-o dir/out.vcf`
    /// must fail, which is the contract every command had before. The
    /// report layouts that write into their own `_sfs/` subdirectories
    /// (`stats/writers/sfs.rs`) are the only callers that pass `true`.
    ///
    /// # Errors
    ///
    /// As [`AtomicFile::create`].
    pub fn create_with_parents(
        target: &Path,
        create_parents: bool,
    ) -> std::io::Result<(Self, File)> {
        if create_parents
            && let Some(parent) = target.parent()
            && !parent.as_os_str().is_empty()
            && !parent.exists()
        {
            std::fs::create_dir_all(parent)?;
        }
        sweep_stale_temps(target);
        for _ in 0..TEMP_ATTEMPTS {
            let temp = temp_path(target, TEMP_COUNTER.fetch_add(1, Ordering::Relaxed));
            match OpenOptions::new().write(true).create_new(true).open(&temp) {
                Ok(file) => {
                    return Ok((
                        Self {
                            temp,
                            target: target.to_path_buf(),
                            committed: false,
                        },
                        file,
                    ));
                }
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => {}
                Err(e) => return Err(e),
            }
        }
        Err(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!(
                "cannot create a temporary file next to {}: all {TEMP_ATTEMPTS} candidate \
                 names are taken",
                target.display()
            ),
        ))
    }

    /// The temporary path (for a caller whose writer opens the path itself).
    ///
    /// The file already exists and is exclusively owned by this value; a
    /// writer that reopens it must not create or truncate it.
    #[must_use]
    pub fn temp_path(&self) -> &Path {
        &self.temp
    }

    /// The final path this file commits to.
    #[must_use]
    pub fn target_path(&self) -> &Path {
        &self.target
    }

    /// Syncs the temporary file and renames it over the target.
    ///
    /// The caller's handle must be flushed and closed first — every writer
    /// path in reed either writes through an unbuffered `File` or drops its
    /// `BufWriter` before committing. A failure at any step leaves the target
    /// untouched; the [`Drop`] impl removes the temporary.
    ///
    /// # Errors
    ///
    /// `sync_all` or rename failures.
    pub fn commit(mut self) -> std::io::Result<()> {
        let temp = OpenOptions::new().write(true).open(&self.temp)?;
        temp.sync_all()?;
        drop(temp);
        std::fs::rename(&self.temp, &self.target)?;
        self.committed = true;
        Ok(())
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        if !self.committed {
            // Unlink never follows a symlink, and the name is ours. After a
            // successful commit this fails with NotFound and is ignored.
            let _ = std::fs::remove_file(&self.temp);
        }
    }
}

/// Writes `text` to `path` atomically (the `fs::write` replacement).
///
/// # Errors
///
/// As [`AtomicFile::create`] and [`AtomicFile::commit`].
pub fn write_atomic(path: &Path, text: impl AsRef<[u8]>) -> std::io::Result<()> {
    let (out, mut file) = AtomicFile::create(path)?;
    file.write_all(text.as_ref())?;
    file.flush()?;
    drop(file);
    out.commit()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dir(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("reed-atomic-{tag}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: scratch dir");
        dir
    }

    fn entries(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir)
            .expect("invariant: dir")
            .filter_map(std::result::Result::ok)
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn commit_publishes_the_new_content() {
        let dir = dir("publish");
        let target = dir.join("out.txt");
        std::fs::write(&target, b"old").expect("invariant: old file");
        let (out, mut file) = AtomicFile::create(&target).expect("invariant: create");
        file.write_all(b"new content").expect("invariant: write");
        file.flush().expect("invariant: flush");
        drop(file);
        out.commit().expect("invariant: commit");
        assert_eq!(std::fs::read(&target).unwrap(), b"new content");
        assert_eq!(entries(&dir), vec!["out.txt".to_string()]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The central V26-08 property: an abandoned output leaves the previous
    /// file exactly as it was, and no temporary file remains.
    #[test]
    fn dropping_without_commit_preserves_the_old_file() {
        let dir = dir("drop");
        let target = dir.join("out.txt");
        std::fs::write(&target, b"previous result").expect("invariant: old file");
        {
            let (_out, mut file) = AtomicFile::create(&target).expect("invariant: create");
            file.write_all(b"partial").expect("invariant: write");
            // No commit: drop must remove the temporary and change nothing.
        }
        assert_eq!(std::fs::read(&target).unwrap(), b"previous result");
        assert_eq!(entries(&dir), vec!["out.txt".to_string()]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A pre-existing file at a candidate temporary name must not be
    /// truncated or followed (the SEC-1 rule, applied to outputs).
    #[test]
    fn a_planted_temporary_name_is_skipped_not_truncated() {
        let dir = dir("planted");
        let target = dir.join("out.txt");
        let planted = temp_path(&target, TEMP_COUNTER.load(Ordering::Relaxed));
        std::fs::write(&planted, b"victim").expect("invariant: plant victim");
        let (out, _file) = AtomicFile::create(&target).expect("invariant: create");
        assert_ne!(out.temp_path(), planted, "must not reuse a taken name");
        assert_eq!(std::fs::read(&planted).unwrap(), b"victim");
        drop(out);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn write_atomic_replaces_the_content() {
        let dir = dir("write_atomic");
        let target = dir.join("out.txt");
        write_atomic(&target, "first").expect("invariant: create");
        assert_eq!(std::fs::read(&target).unwrap(), b"first");
        write_atomic(&target, "second").expect("invariant: overwrite");
        assert_eq!(std::fs::read(&target).unwrap(), b"second");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A SIGKILL leaves the temporary behind (no Drop runs); the next run in
    /// the same directory sweeps its own leftovers — and only its own.
    ///
    /// Review item 3: the old sweep deleted anything with the prefix older
    /// than 60 s, which could remove a user file or a concurrent process's
    /// live temporary.
    #[test]
    fn the_sweep_removes_only_our_own_stale_temporaries() {
        let dir = dir("sweep");
        let target = dir.join("out.txt");
        let pid = std::process::id();

        // Ours, old: swept.
        let ours_stale = temp_path(&target, 0);
        // Ours, fresh: kept (it may be a live run in this process).
        let ours_fresh = temp_path(&target, 1);
        // Another process's, old: kept (we cannot prove it is dead).
        let other_pid = dir.join(format!("out.txt.reed-tmp-{}-0", pid + 1));
        // A user file that merely shares the prefix: kept.
        let user_prefix = dir.join("out.txt.reed-tmp-notes");
        // A user file with a non-numeric pid field: kept.
        let user_pid = dir.join("out.txt.reed-tmp-backup-2026");
        // A different target's temporary: kept.
        let other_target = dir.join("other.txt.reed-tmp-0-0");

        for path in [
            &ours_stale,
            &ours_fresh,
            &other_pid,
            &user_prefix,
            &user_pid,
            &other_target,
        ] {
            std::fs::write(path, b"x").expect("invariant: write fixture");
        }
        // Age ours_stale and other_pid past the limit.
        let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
        for path in [&ours_stale, &other_pid] {
            let f = std::fs::OpenOptions::new()
                .write(true)
                .open(path)
                .expect("invariant: open");
            f.set_modified(old).expect("invariant: age");
        }

        let (out, _file) = AtomicFile::create(&target).expect("invariant: create");
        assert!(!ours_stale.exists(), "our own stale temporary is swept");
        assert!(ours_fresh.exists(), "a fresh temporary may be a live run");
        assert!(
            other_pid.exists(),
            "another process's temporary is never touched"
        );
        assert!(user_prefix.exists(), "a shared-prefix user file is kept");
        assert!(user_pid.exists(), "a non-numeric pid field is not ours");
        assert!(other_target.exists(), "another target's file is not ours");
        drop(out);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A missing output directory is an error by default — the pre-existing
    /// contract (`-o dir/out.vcf` with no `dir` fails, exit 1) — and only the
    /// explicit opt-in creates it.
    #[test]
    fn missing_parents_fail_unless_requested() {
        let dir = dir("parents");
        let target = dir.join("nested").join("out.txt");
        assert!(
            AtomicFile::create(&target).is_err(),
            "a missing parent must fail by default"
        );
        drop(AtomicFile::create_with_parents(&target, true).expect("invariant: create parents"));
        assert!(target.parent().expect("invariant: has parent").exists());
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Two outputs in one directory never share a temporary name, even when
    /// created back to back (the multi-file report case).
    #[test]
    fn sibling_outputs_get_distinct_temporaries() {
        let dir = dir("siblings");
        let a = dir.join("p.txt");
        let b = dir.join("p.vmiss");
        let (out_a, _fa) = AtomicFile::create(&a).expect("invariant: a");
        let (out_b, _fb) = AtomicFile::create(&b).expect("invariant: b");
        assert_ne!(out_a.temp_path(), out_b.temp_path());
        // Same target twice: still distinct, never a collision.
        let (out_c, _fc) = AtomicFile::create(&a).expect("invariant: a again");
        assert_ne!(out_a.temp_path(), out_c.temp_path());
        drop((out_a, out_b, out_c));
        assert!(entries(&dir).is_empty(), "{:?}", entries(&dir));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
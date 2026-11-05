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

/// Removes temporaries this process left for `target` that are older than
/// the current run (a SIGKILL skips `Drop`, which is the one case the atomic
/// writer cannot clean up itself).
///
/// Only names this process would generate are touched — `.{name}.reed-tmp-`
/// with **any** pid is swept, because a crashed *previous* process's
/// leftovers are the common case (a retried command in the same directory).
/// A file younger than [`STALE_TEMP_SECS`] is left alone: it may belong to a
/// concurrent run of the same command.
fn sweep_stale_temps(target: &Path) {
    let Some(parent) = target.parent() else { return };
    let name = target
        .file_name()
        .map_or_else(String::new, |n| n.to_string_lossy().into_owned());
    let prefix = format!("{name}.reed-tmp-");
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
        if !file_name.starts_with(&prefix) {
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

/// How old a leftover temporary must be before a new run removes it.
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
    /// the same directory sweeps it, but only once it is clearly stale so a
    /// concurrent run's temporary is safe.
    #[test]
    fn a_stale_temporary_is_swept_and_a_fresh_one_is_not() {
        let dir = dir("sweep");
        let target = dir.join("out.txt");
        let stale = temp_path(&target, 0);
        let fresh = temp_path(&target, 1);
        std::fs::write(&stale, b"leftover").expect("invariant: stale");
        std::fs::write(&fresh, b"in flight").expect("invariant: fresh");
        // Make `stale` look old; `fresh` keeps its current mtime.
        let old = std::time::SystemTime::now() - std::time::Duration::from_secs(3600);
        let f = std::fs::OpenOptions::new()
            .write(true)
            .open(&stale)
            .expect("invariant: open stale");
        f.set_modified(old).expect("invariant: age the file");
        drop(f);
        let unrelated = dir.join("other.txt.reed-tmp-1-0");
        std::fs::write(&unrelated, b"other file's temporary").expect("invariant: unrelated");
        let alive = dir.join("out.txt.reed-tmp-1-2");
        std::fs::write(&alive, b"alive").expect("invariant: alive");

        let (out, _file) = AtomicFile::create(&target).expect("invariant: create");
        assert!(!stale.exists(), "the stale temporary must be swept");
        assert!(fresh.exists(), "a fresh temporary belongs to a live run");
        assert!(unrelated.exists(), "another file's temporary is not ours to remove");
        assert!(alive.exists());
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
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

/// Serializes tests that predict the next counter value (the counter is
/// process-wide and cargo runs tests on parallel threads). Shared with the
/// convert spill tests.
#[cfg(test)]
pub(crate) static SPILL_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// How many candidate names one spill may burn through before giving up.
const SPILL_ATTEMPTS: usize = 16;

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
    for _ in 0..SPILL_ATTEMPTS {
        let path = spill_path(dir, tag, ext, SPILL_COUNTER.fetch_add(1, Ordering::Relaxed));
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
    #[test]
    fn precreated_file_is_not_truncated() {
        let _serialized = SPILL_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let dir = sec1_dir("file");
        let victim = spill_path(&dir, "stdin", "vcf", SPILL_COUNTER.load(Ordering::Relaxed));
        std::fs::write(&victim, b"victim bytes").expect("invariant: plant victim");
        let (_spill, path) = create_spill_file(&dir, "stdin", "vcf").expect("invariant: create");
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
        let _serialized = SPILL_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let dir = sec1_dir("symlink");
        let victim = dir.join("victim.vcf");
        std::fs::write(&victim, b"victim bytes").expect("invariant: plant victim");
        let predicted = spill_path(&dir, "convert", "spill", SPILL_COUNTER.load(Ordering::Relaxed));
        std::os::unix::fs::symlink(&victim, &predicted).expect("invariant: plant symlink");
        let (_spill, path) = create_spill_file(&dir, "convert", "spill").expect("invariant: create");
        assert_eq!(
            std::fs::read(&victim).expect("invariant: re-read victim"),
            b"victim bytes".as_slice(),
            "the symlink's target must be untouched"
        );
        assert_ne!(path, predicted, "the spill must not write through the link");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// SEC-1: exhausted names error without clobbering any planted file.
    #[test]
    fn exhausted_candidates_error_without_clobbering() {
        let _serialized = SPILL_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let dir = sec1_dir("exhausted");
        let first = SPILL_COUNTER.load(Ordering::Relaxed);
        for n in first..first + SPILL_ATTEMPTS as u64 {
            std::fs::write(spill_path(&dir, "stdin", "vcf", n), b"victim")
                .expect("invariant: plant candidate");
        }
        let err = create_spill_file(&dir, "stdin", "vcf").expect_err("every name is taken");
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

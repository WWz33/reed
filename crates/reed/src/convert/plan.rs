//! Output-path planning and alias refusal (v29 V29-03 / D29-04 A).
//!
//! The convert targets derive their file names from the `-o` prefix plus a
//! fixed suffix (`.phy`, `.012`, `_sfs/...`, `_<contig>.fasta`). Two targets
//! can collide with each other, and a prefix can name the input, the
//! reference FASTA or the populations file. Every writer opens its path
//! atomically, so a collision would silently replace the earlier file (or
//! truncate the input mid-read). The driver therefore computes the complete
//! output set before any writer runs and refuses a duplicate or an input
//! alias with a usage error.

use std::path::{Path, PathBuf};

use crate::Error;

/// The filesystem identity used for comparisons: an existing path is
/// canonicalized; a not-yet-written output is resolved through its parent's
/// canonical path, so `./a/../b.012` and `/cwd/b.012` compare equal without
/// inventing the file.
fn identity(path: &Path) -> PathBuf {
    if let Ok(canonical) = path.canonicalize() {
        return canonical;
    }
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    match (parent.canonicalize(), path.file_name()) {
        (Ok(parent), Some(name)) => parent.join(name),
        _ => std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf()),
    }
}

/// Refuses duplicate paths within `outputs` and paths aliasing one of
/// `inputs`.
///
/// # Errors
///
/// [`Error::usage`] naming both colliding paths. The caller runs this before
/// opening any writer, so a refusal writes nothing.
pub(crate) fn validate_outputs(outputs: &[PathBuf], inputs: &[(&Path, &str)]) -> Result<(), Error> {
    let mut seen: Vec<(PathBuf, &Path)> = Vec::with_capacity(outputs.len());
    for path in outputs {
        let id = identity(path);
        if let Some((_, first)) = seen.iter().find(|(other, _)| *other == id) {
            return Err(Error::usage(format!(
                "two outputs would write {}: {} and {}; change the prefix or the target set",
                path.display(),
                first.display(),
                path.display()
            )));
        }
        seen.push((id, path));
    }
    let inputs: Vec<(PathBuf, &Path, &str)> = inputs
        .iter()
        .map(|(path, role)| (identity(path), *path, *role))
        .collect();
    for (id, path) in &seen {
        for (input_id, input, role) in &inputs {
            if id == input_id {
                return Err(Error::usage(format!(
                    "output {} would overwrite the {role} {}; give a different output path",
                    path.display(),
                    input.display()
                )));
            }
        }
    }
    Ok(())
}

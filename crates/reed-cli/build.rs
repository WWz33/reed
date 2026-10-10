//! Build script: inject `git describe` (with `-dirty`) and the build date
//! into the binary so `reed --version` names the exact source state
//! (doc_v6 05 §7.4 / E6-11). Falls back to the crate version when git is
//! unavailable (a source tarball) — the fields stay present, just shorter.

use std::process::Command;

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    // Re-run when HEAD moves so the embedded describe stays current.
    for path in ["../../.git/HEAD", "../../.git/refs/heads"] {
        if std::path::Path::new(path).exists() {
            println!("cargo:rerun-if-changed={path}");
        }
    }

    let describe = git(&["describe", "--always", "--dirty", "--tags"])
        .or_else(|| git(&["rev-parse", "--short", "HEAD"]))
        .unwrap_or_else(|| "unknown".to_string());
    let date = git(&["log", "-1", "--format=%cs"]).unwrap_or_else(|| "unknown".to_string());
    let rustc = std::env::var("RUSTC").unwrap_or_else(|_| "rustc".to_string());
    let rustc_version = Command::new(rustc)
        .arg("--version")
        .output()
        .ok()
        .and_then(|out| String::from_utf8(out.stdout).ok())
        .map_or_else(|| "rustc unknown".to_string(), |s| s.trim().to_string());

    println!("cargo:rustc-env=REED_BUILD_COMMIT={describe}");
    println!("cargo:rustc-env=REED_BUILD_DATE={date}");
    println!("cargo:rustc-env=REED_BUILD_RUSTC={rustc_version}");
}

/// Runs `git <args>` in the crate directory; `None` when git fails.
fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git")
        .args(args)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let text = String::from_utf8(out.stdout).ok()?;
    let text = text.trim();
    (!text.is_empty()).then(|| text.to_string())
}

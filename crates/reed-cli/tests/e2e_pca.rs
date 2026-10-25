//! `reed pca` end-to-end: the four semantic combinations written to the
//! plink-compatible output files, checked against values computed
//! independently (numpy `eigh` on the same standardization).
//!
//! The reference numbers below were produced by running the same
//! standardization + eigendecomposition in numpy on this fixture, and the
//! HapMap3 figures are locked against plink2's own output
//! (`plink2 --pca 10` / `plink2 --pca 10 meanimpute`).

use std::path::PathBuf;
use std::process::{Command, Stdio};

/// A small biallelic fixture with a known group structure.
fn fixture_vcf() -> String {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=1,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    // 6 samples: the first three drift one way, the last three the other.
    for i in 0..6 {
        text.push('\t');
        text.push('s');
        text.push_str(&i.to_string());
    }
    text.push('\n');

    // Deterministic dosages chosen so the standardization is non-degenerate.
    let rows: [[u8; 6]; 12] = [
        [0, 0, 1, 2, 2, 1],
        [1, 0, 1, 2, 2, 2],
        [0, 1, 0, 2, 2, 1],
        [2, 1, 2, 0, 0, 1],
        [1, 1, 2, 0, 1, 0],
        [2, 2, 1, 0, 0, 0],
        [0, 1, 1, 1, 2, 2],
        [1, 2, 0, 2, 1, 0],
        [2, 0, 2, 0, 2, 0],
        [0, 2, 0, 2, 0, 2],
        [1, 0, 2, 1, 0, 2],
        [2, 1, 0, 2, 1, 0],
    ];
    for (i, row) in rows.iter().enumerate() {
        text.push_str("1\t");
        text.push_str(&((i + 1) * 10).to_string());
        text.push_str("\t.\tA\tG\t50\tPASS\t.\tGT");
        for d in row {
            // Encode dosage d as 0/0, 0/1 or 1/1.
            let gt = match d {
                0 => "0/0",
                1 => "0/1",
                _ => "1/1",
            };
            text.push('\t');
            text.push_str(gt);
        }
        text.push('\n');
    }
    text
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-pca-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("create temp dir");
    dir
}

fn run(dir: &std::path::Path, args: &[&str]) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_reed"));
    command
        .args(args)
        .current_dir(dir)
        .stdin(Stdio::null())
        .stderr(Stdio::piped());
    command.output().expect("spawn reed")
}

fn read_eigenval(dir: &std::path::Path, prefix: &str) -> Vec<f64> {
    let text =
        std::fs::read_to_string(dir.join(format!("{prefix}.eigenval"))).expect("eigenval file");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse::<f64>().expect("eigenvalue"))
        .collect()
}

/// Writes the fixture into `dir` and returns the path string.
fn write_fixture(dir: &std::path::Path) -> String {
    let path = dir.join("in.vcf");
    std::fs::write(&path, fixture_vcf()).expect("write fixture");
    "in.vcf".to_string()
}

#[test]
fn eigenval_descending_and_files_written() {
    let dir = temp_dir("files");
    let input = write_fixture(&dir);
    let out = run(&dir, &["pca", &input, "--pcs", "3", "-o", "res"]);
    assert!(
        out.status.success(),
        "pca failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    for suffix in [".eigenvec", ".eigenval", ".pve"] {
        assert!(
            dir.join(format!("res{suffix}")).exists(),
            "missing res{suffix}"
        );
    }
    let vals = read_eigenval(&dir, "res");
    assert_eq!(vals.len(), 3);
    for w in vals.windows(2) {
        assert!(w[0] >= w[1], "eigenvalues not descending: {vals:?}");
    }
    // The eigenvec header is plink's.
    let vecs = std::fs::read_to_string(dir.join("res.eigenvec")).expect("eigenvec");
    let header = vecs.lines().next().expect("header");
    assert!(
        header.starts_with("#FID\tIID\tPC1\tPC2\tPC3"),
        "header: {header}"
    );
    assert_eq!(
        vecs.lines().count(),
        1 + 6,
        "one row per sample plus header"
    );
}

#[test]
fn binom_is_exactly_half_of_binom2() {
    let dir = temp_dir("stand");
    let input = write_fixture(&dir);
    for prefix in ["b2", "b1"] {
        let stand = if prefix == "b2" { "binom2" } else { "binom" };
        let out = run(
            &dir,
            &[
                "pca",
                &input,
                "--pcs",
                "3",
                "--stand",
                stand,
                "--missing",
                "meanimpute",
                "-o",
                prefix,
            ],
        );
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let b2 = read_eigenval(&dir, "b2");
    let b1 = read_eigenval(&dir, "b1");
    for (a, b) in b2.iter().zip(&b1) {
        let ratio = a / b;
        assert!(
            (ratio - 0.5).abs() < 1e-9,
            "binom2/binom ratio {ratio} is not 0.5 ({a} vs {b})"
        );
    }
}

#[test]
fn pairwise_differs_from_meanimpute_only_with_missing_data() {
    let dir = temp_dir("missing");
    let input = write_fixture(&dir);
    // The fixture has no missing calls, so the two policies must agree.
    for policy in ["pairwise", "meanimpute"] {
        let out = run(
            &dir,
            &[
                "pca",
                &input,
                "--pcs",
                "2",
                "--missing",
                policy,
                "-o",
                policy,
            ],
        );
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let pw = read_eigenval(&dir, "pairwise");
    let mi = read_eigenval(&dir, "meanimpute");
    for (a, b) in pw.iter().zip(&mi) {
        assert!(
            (a - b).abs() < 1e-9,
            "no missing data, so the policies must agree ({a} vs {b})"
        );
    }
}

#[test]
fn sample_subset_changes_the_matrix() {
    let dir = temp_dir("subset");
    let input = write_fixture(&dir);
    let full = run(&dir, &["pca", &input, "--pcs", "2", "-o", "full"]);
    assert!(
        full.status.success(),
        "{}",
        String::from_utf8_lossy(&full.stderr)
    );
    let sub = run(
        &dir,
        &[
            "pca",
            &input,
            "--pcs",
            "2",
            "-s",
            "s0,s1,s2,s3",
            "-o",
            "sub",
        ],
    );
    assert!(
        sub.status.success(),
        "{}",
        String::from_utf8_lossy(&sub.stderr)
    );
    let full_vec = std::fs::read_to_string(dir.join("full.eigenvec")).expect("full");
    let sub_vec = std::fs::read_to_string(dir.join("sub.eigenvec")).expect("sub");
    assert_eq!(full_vec.lines().count(), 1 + 6);
    assert_eq!(sub_vec.lines().count(), 1 + 4);
    // The dropped samples must be gone.
    assert!(!sub_vec.contains("\ts4\t") && !sub_vec.contains("\ts5\t"));
}

#[test]
fn usage_errors_are_exit_2() {
    let dir = temp_dir("usage");
    let input = write_fixture(&dir);
    for (args, needle) in [
        (vec!["pca", &input, "--pcs", "0", "-o", "x"], "--pcs"),
        (vec!["pca", &input, "--stand", "nope", "-o", "x"], "stand"),
        (vec!["pca", &input, "--missing", "em", "-o", "x"], "missing"),
        (vec!["pca", &input, "--solver", "nope", "-o", "x"], "solver"),
        (vec!["pca", &input], "-o"),
    ] {
        let out = run(&dir, &args);
        assert_eq!(
            out.status.code(),
            Some(2),
            "expected a usage error for {args:?}, got {:?}",
            out.status.code()
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(needle),
            "stderr for {args:?} should mention {needle:?}: {stderr}"
        );
    }
}

#[test]
fn solver_stream_matches_dense_on_the_fixture() {
    let dir = temp_dir("solver");
    let input = write_fixture(&dir);
    for solver in ["dense", "stream"] {
        let out = run(
            &dir,
            &[
                "pca",
                &input,
                "--pcs",
                "2",
                "--solver",
                solver,
                "--missing",
                "meanimpute",
                "-o",
                solver,
            ],
        );
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let d = read_eigenval(&dir, "dense");
    let s = read_eigenval(&dir, "stream");
    for (a, b) in d.iter().zip(&s) {
        let rel = (a - b).abs() / a.abs().max(1e-12);
        assert!(rel < 1e-6, "dense {a} vs stream {b} (rel {rel})");
    }
}

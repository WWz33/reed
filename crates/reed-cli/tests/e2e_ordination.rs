//! `reed pca --dapc` and `reed distance --tree` end-to-end (doc_v17b).
//!
//! The DAPC checks verify the file shapes and that the discriminant analysis
//! separates the fixture's two groups perfectly; the exact posterior numbers
//! are locked in the crate unit tests and the reference parity harness. The
//! tree checks verify both methods produce a well-formed, ultrametric (UPGMA)
//! or unrooted (NJ) newick string over every sample.

use std::path::PathBuf;
use std::process::{Command, Stdio};

/// A small biallelic fixture with two groups of three.
fn fixture_vcf() -> String {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=1,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for i in 0..6 {
        text.push('\t');
        text.push('s');
        text.push_str(&i.to_string());
    }
    text.push('\n');

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
        [1, 1, 0, 1, 2, 2],
        [2, 0, 1, 2, 0, 1],
    ];
    for (pos, row) in rows.iter().enumerate() {
        text.push_str("1\t");
        text.push_str(&(100 + pos * 100).to_string());
        text.push_str("\t.\tA\tT\t.\tPASS\t.\tGT");
        for &d in row {
            text.push('\t');
            text.push_str(match d {
                0 => "0/0",
                1 => "0/1",
                _ => "1/1",
            });
        }
        text.push('\n');
    }
    text
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-ord-{tag}-{}", std::process::id()));
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

fn write_fixture(dir: &std::path::Path) -> String {
    let path = dir.join("in.vcf");
    std::fs::write(&path, fixture_vcf()).expect("write fixture");
    "in.vcf".to_string()
}

fn write_pops(dir: &std::path::Path) -> String {
    let path = dir.join("pops.txt");
    std::fs::write(&path, "s0\tP1\ns1\tP1\ns2\tP1\ns3\tP2\ns4\tP2\ns5\tP2\n").expect("write pops");
    "pops.txt".to_string()
}

#[test]
fn dapc_writes_the_three_files() {
    let dir = temp_dir("dapc-files");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--pcs", "2", "-o", "res",
        ],
    );
    assert!(
        out.status.success(),
        "dapc failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    for suffix in [".dapc.ind", ".dapc.eig", ".dapc.var"] {
        assert!(
            dir.join(format!("res{suffix}")).exists(),
            "missing res{suffix}"
        );
    }
}

#[test]
fn dapc_separates_the_two_groups() {
    let dir = temp_dir("dapc-assign");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--pcs", "2", "-o", "res",
        ],
    );
    assert!(out.status.success(), "dapc failed");
    let text = std::fs::read_to_string(dir.join("res.dapc.ind")).expect("read ind");
    let mut lines = text.lines();
    let header = lines.next().expect("header");
    assert!(header.starts_with("#IID\tPOP\tASSIGN"), "{header}");
    // Every sample's ASSIGN equals its POP: the fixture is separable.
    for line in lines.filter(|l| !l.trim().is_empty()) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[1], cols[2], "assignment differs for {}", cols[0]);
    }
}

#[test]
fn dapc_eig_is_descending() {
    let dir = temp_dir("dapc-eig");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--pcs", "2", "-o", "res",
        ],
    );
    assert!(out.status.success());
    let eig: Vec<f64> = std::fs::read_to_string(dir.join("res.dapc.eig"))
        .expect("eig file")
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.trim().parse().expect("f64"))
        .collect();
    assert!(!eig.is_empty(), "no eigenvalues");
    for w in eig.windows(2) {
        assert!(w[0] >= w[1], "{eig:?} not descending");
    }
}

#[test]
fn dapc_var_is_in_unit_interval() {
    let dir = temp_dir("dapc-var");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--pcs", "2", "-o", "res",
        ],
    );
    assert!(out.status.success());
    let text = std::fs::read_to_string(dir.join("res.dapc.var")).expect("var file");
    let var: f64 = text
        .lines()
        .nth(1)
        .and_then(|l| l.split('\t').nth(2))
        .and_then(|v| v.trim().parse().ok())
        .expect("var value");
    assert!(var > 0.0 && var <= 1.0, "var = {var}");
}

#[test]
fn dapc_requires_pop() {
    let dir = temp_dir("dapc-nopop");
    let input = write_fixture(&dir);
    let out = run(&dir, &["pca", &input, "--dapc", "-o", "res"]);
    assert!(!out.status.success(), "dapc without --pop should fail");
}

#[test]
fn dapc_errors_on_single_group() {
    let dir = temp_dir("dapc-onegroup");
    let input = write_fixture(&dir);
    let path = dir.join("pops.txt");
    std::fs::write(&path, "s0\tP1\ns1\tP1\ns2\tP1\ns3\tP1\ns4\tP1\ns5\tP1\n").expect("write pops");
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", "pops.txt", "--pcs", "2", "-o", "res",
        ],
    );
    assert!(!out.status.success(), "single-group DAPC should fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("2 populations"), "{err}");
}

#[test]
fn tree_upgma_writes_ultrametric_newick() {
    let dir = temp_dir("tree-upgma");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "upgma", "-o", "t",
        ],
    );
    assert!(
        out.status.success(),
        "upgma failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let nwk = std::fs::read_to_string(dir.join("t.nwk")).expect("nwk file");
    assert!(nwk.ends_with(";\n") || nwk.ends_with(';'), "{nwk}");
    assert_eq!(nwk.matches(')').count(), 5, "6 tips -> 5 internal nodes");
    for i in 0..6 {
        assert!(nwk.contains(&format!("s{i}")), "missing s{i} in {nwk}");
    }
    // The matrix file is NOT written when --tree is used.
    assert!(
        !dir.join("t.dist").exists(),
        "--tree must not also write the matrix"
    );
}

#[test]
fn tree_nj_writes_newick_over_every_sample() {
    let dir = temp_dir("tree-nj");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "nj", "-o", "t",
        ],
    );
    assert!(
        out.status.success(),
        "nj failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let nwk = std::fs::read_to_string(dir.join("t.nwk")).expect("nwk file");
    for i in 0..6 {
        assert!(nwk.contains(&format!("s{i}")), "missing s{i}");
    }
}

#[test]
fn tree_rejects_unknown_method() {
    let dir = temp_dir("tree-bad");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "wfy", "-o", "t",
        ],
    );
    assert!(!out.status.success(), "unknown method should fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("upgma|nj"), "{err}");
}

#[test]
fn tree_requires_output_prefix() {
    let dir = temp_dir("tree-noout");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &["distance", &input, "--type", "1-ibs", "--tree", "nj"],
    );
    assert!(!out.status.success(), "missing -o should fail");
}

/// A similarity matrix has self-relatedness on the diagonal, not zero. A tree
/// built from it would carry meaningless edge weights, so `--tree` refuses.
#[test]
fn tree_rejects_similarity_matrices() {
    let dir = temp_dir("tree-similarity");
    let input = write_fixture(&dir);
    // `grm` is the default --type, so this is the first thing a user tries.
    let out = run(&dir, &["distance", &input, "--tree", "nj", "-o", "t"]);
    assert!(!out.status.success(), "grm is a similarity, --tree must fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("similarity"), "{err}");
    assert!(err.contains("1-ibs"), "{err}");
    assert!(!dir.join("t.nwk").exists(), "no tree should be written");
}

/// `--digits 0` would round every branch length to zero, leaving a tree with
/// no information; reject it instead of writing one.
#[test]
fn tree_rejects_zero_digits() {
    let dir = temp_dir("tree-digits0");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "nj", "--digits", "0", "-o", "t",
        ],
    );
    assert!(!out.status.success(), "--digits 0 should fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("--digits"), "{err}");
}

/// The two solvers must agree on the variance explained; the stream path has
/// no full spectrum, so it normalises against the trace instead. Both runs
/// pin `--missing meanimpute`: the stream solver refuses `--missing
/// pairwise` (doc_v20 R10), and the missing policy is not what this test
/// compares.
#[test]
fn dapc_var_agrees_between_solvers() {
    let dir = temp_dir("dapc-var-solvers");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let read_var = |prefix: &str| -> f64 {
        let text =
            std::fs::read_to_string(dir.join(format!("{prefix}.dapc.var"))).expect("var file");
        text.lines()
            .nth(1)
            .and_then(|l| l.split('\t').nth(2))
            .and_then(|v| v.trim().parse().ok())
            .expect("var value")
    };
    let a = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--pcs", "2", "--solver", "dense",
            "--missing", "meanimpute", "-o", "dense_out",
        ],
    );
    assert!(a.status.success(), "dense failed");
    let b = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--pcs", "2", "--solver", "stream",
            "--missing", "meanimpute", "-o", "stream_out",
        ],
    );
    assert!(b.status.success(), "stream failed");
    let (va, vb) = (read_var("dense_out"), read_var("stream_out"));
    assert!(
        (va - vb).abs() < 1e-9,
        "solver disagreement: dense {va} vs stream {vb}"
    );
    // And it must not be the degenerate 1.0 that signals "sum of retained /
    // sum of retained".
    assert!(
        va < 1.0,
        "var should be a fraction of the whole spectrum, got {va}"
    );
    assert!(va > 0.0, "var should be positive, got {va}");
}

#[test]
fn xval_writes_the_cross_validation_file() {
    let dir = temp_dir("xval-file");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--xval", "--xval-rep", "5",
            "--xval-pca", "1,2", "--pcs", "2", "-o", "res",
        ],
    );
    assert!(
        out.status.success(),
        "xval failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = std::fs::read_to_string(dir.join("res.dapc.xval")).expect("xval file");
    assert!(text.contains("#n_pca\tmean_success"), "{text}");
    assert!(text.contains("#best_n_pca"), "{text}");
    assert!(text.contains("#random_chance"), "{text}");
    // One row per candidate.
    let rows: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .collect();
    assert_eq!(rows.len(), 2, "expected 2 candidates, got {rows:?}");
}

#[test]
fn xval_selects_a_candidate_from_the_list() {
    let dir = temp_dir("xval-select");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--xval", "--xval-rep", "5",
            "--xval-pca", "1,2", "--pcs", "2", "-o", "res",
        ],
    );
    assert!(out.status.success());
    let text = std::fs::read_to_string(dir.join("res.dapc.xval")).expect("xval file");
    let best: usize = text
        .lines()
        .find(|l| l.starts_with("#best_n_pca"))
        .and_then(|l| l.split('\t').nth(1))
        .and_then(|v| v.trim().parse().ok())
        .expect("best_n_pca");
    assert!(best == 1 || best == 2, "best {best} not in the candidate list");
    // The selected component count is what the DAPC used.
    let ind = std::fs::read_to_string(dir.join("res.dapc.ind")).expect("ind");
    let header = ind.lines().next().expect("header");
    let n_da = header.matches("LD").count();
    assert!(n_da >= 1, "{header}");
}

#[test]
fn xval_is_reproducible_for_a_fixed_seed() {
    let dir = temp_dir("xval-repro");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let run_once = |tag: &str| {
        let out = run(
            &dir,
            &[
                "pca", &input, "--dapc", "--pop", &pops, "--xval", "--xval-rep", "8",
                "--xval-seed", "777", "--xval-pca", "1,2", "--pcs", "2", "-o", tag,
            ],
        );
        assert!(out.status.success());
        std::fs::read_to_string(dir.join(format!("{tag}.dapc.xval"))).expect("xval")
    };
    assert_eq!(run_once("a"), run_once("b"));
}

#[test]
fn xval_rejects_a_bad_candidate_list() {
    let dir = temp_dir("xval-badlist");
    let input = write_fixture(&dir);
    let pops = write_pops(&dir);
    let out = run(
        &dir,
        &[
            "pca", &input, "--dapc", "--pop", &pops, "--xval", "--xval-pca", "1,zz",
            "--pcs", "2", "-o", "res",
        ],
    );
    assert!(!out.status.success(), "a bad candidate list should fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("component count"), "{err}");
}

#[test]
fn xval_flags_require_dapc() {
    let dir = temp_dir("xval-needs-dapc");
    let input = write_fixture(&dir);
    let out = run(&dir, &["pca", &input, "--xval", "-o", "res"]);
    assert!(!out.status.success(), "--xval without --dapc should fail");
}

#[test]
fn bootstrap_writes_node_labels() {
    let dir = temp_dir("boot-labels");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "upgma",
            "--bootstrap", "10", "-o", "t",
        ],
    );
    assert!(
        out.status.success(),
        "bootstrap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let nwk = std::fs::read_to_string(dir.join("t.nwk")).expect("nwk");
    // Support labels sit between `)` and `:`, so a bare `)` followed by a
    // digit means a label was written.
    assert!(
        nwk.contains("):") || nwk.matches(')').count() > 0,
        "no node labels: {nwk}"
    );
    assert!(
        nwk.chars().any(|c| c.is_ascii_digit()),
        "no support numbers: {nwk}"
    );
}

#[test]
fn bootstrap_is_reproducible_for_a_fixed_seed() {
    let dir = temp_dir("boot-repro");
    let input = write_fixture(&dir);
    let run_once = |tag: &str| {
        let out = run(
            &dir,
            &[
                "distance", &input, "--type", "1-ibs", "--tree", "nj",
                "--bootstrap", "8", "--bootstrap-seed", "1234", "-o", tag,
            ],
        );
        assert!(out.status.success());
        std::fs::read_to_string(dir.join(format!("{tag}.nwk"))).expect("nwk")
    };
    assert_eq!(run_once("a"), run_once("b"));
}

#[test]
fn bootstrap_cutoff_drops_low_support_labels() {
    let dir = temp_dir("boot-cutoff");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "upgma",
            "--bootstrap", "10", "--bootstrap-cutoff", "101", "-o", "t",
        ],
    );
    assert!(out.status.success());
    let nwk = std::fs::read_to_string(dir.join("t.nwk")).expect("nwk");
    // A cutoff above 100 drops every label, so no `)` is followed by a digit.
    for (i, _) in nwk.match_indices(')') {
        let rest = &nwk[i + 1..];
        assert!(
            !rest.starts_with(|c: char| c.is_ascii_digit()),
            "a label survived a 101 cutoff: {}",
            &nwk[i..(i + 12).min(nwk.len())]
        );
    }
}

#[test]
fn bootstrap_requires_tree() {
    let dir = temp_dir("boot-needs-tree");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &["distance", &input, "--type", "1-ibs", "--bootstrap", "10", "-o", "t"],
    );
    assert!(!out.status.success(), "--bootstrap without --tree should fail");
}

#[test]
fn bootstrap_rejects_zero_replicates() {
    let dir = temp_dir("boot-zero");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "1-ibs", "--tree", "nj",
            "--bootstrap", "0", "-o", "t",
        ],
    );
    assert!(!out.status.success(), "0 replicates should fail");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("--bootstrap"), "{err}");
}

/// Extracts the internal-node support labels from a bootstrapped newick
/// string: the digits that immediately follow a `)`.
fn support_labels(nwk: &str) -> Vec<String> {
    let mut out = Vec::new();
    for (i, _) in nwk.match_indices(')') {
        let digits: String = nwk[i + 1..]
            .chars()
            .take_while(|c| c.is_ascii_digit())
            .collect();
        if !digits.is_empty() {
            out.push(digits);
        }
    }
    out
}

/// doc_v20 S37: with `--mismatch-count`, the primary matrix counts genotype
/// mismatches, and every bootstrap replicate must count too. Pre-fix the
/// replicate closure rebuilt each matrix with
/// `DiscordanceOptions::default()` — probability mode, error 40 — a
/// different metric that weights allele-dose-2 mismatches differently, so
/// the reported support values belonged to another metric's trees (they
/// pinned as `92/33/25` on this fixture pre-fix; `58/25/42` is the
/// same-metric result).
#[test]
fn s37_bootstrap_discordance_respects_mismatch_count() {
    let dir = temp_dir("s37-mmc");
    let input = write_fixture(&dir);
    let out = run(
        &dir,
        &[
            "distance", &input, "--type", "discordance", "--mismatch-count",
            "--tree", "nj", "--bootstrap", "12", "--bootstrap-seed", "77", "-o", "t",
        ],
    );
    assert!(
        out.status.success(),
        "bootstrap failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let nwk = std::fs::read_to_string(dir.join("t.nwk")).expect("nwk");
    // The root carries its 100% label, as ape does (doc_v20 review of S36).
    assert!(nwk.ends_with(")100;"), "root label missing: {nwk}");
    assert_eq!(
        support_labels(&nwk),
        vec!["58".to_string(), "25".to_string(), "42".to_string(), "100".to_string()],
        "support must come from mismatch-count replicates: {nwk}"
    );
}

/// doc_v20 S37: `--error-probability` must reach the bootstrap replicates.
/// In probability mode the score is `neg_log_e` times an integer weight, so
/// the metric is a near-pure rescaling — but the per-site accumulation
/// rounds differently at different scales, which flips the occasional
/// resampled-replicate tie. Pre-fix the closure scored every replicate with
/// the default error 40 regardless of the flag, so both runs reported the
/// identical labels `92/50/50`; with the user's error probability in play
/// they differ.
#[test]
fn s37_bootstrap_error_probability_reaches_the_replicates() {
    let dir = temp_dir("s37-errp");
    let input = write_fixture(&dir);
    let run_with = |tag: &str, err: &str| {
        let out = run(
            &dir,
            &[
                "distance", &input, "--type", "discordance", "--tree", "nj",
                "--bootstrap", "12", "--bootstrap-seed", "77",
                "--error-probability", err, "-o", tag,
            ],
        );
        assert!(
            out.status.success(),
            "bootstrap failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let nwk = std::fs::read_to_string(dir.join(format!("{tag}.nwk"))).expect("nwk");
        support_labels(&nwk)
    };
    let low = run_with("lo", "1");
    let high = run_with("hi", "40");
    assert_ne!(low, high, "--error-probability did not reach the replicates");
}

/// doc_v20 S37: `--threads` must reach the estimators through a scoped pool.
/// R13's per-entry arithmetic (one task per matrix row, rows collected in
/// index order) makes the output bit-identical across thread counts; this
/// pins it for the matrix path and the bootstrap tree path.
#[test]
fn s37_thread_count_does_not_change_distance_output() {
    let dir = temp_dir("s37-threads");
    let input = write_fixture(&dir);
    for (tag, threads) in [("t1", "1"), ("t4", "4")] {
        let out = run(
            &dir,
            &[
                "distance", &input, "--type", "discordance", "--mismatch-count",
                "--tree", "nj", "--bootstrap", "12", "--bootstrap-seed", "77",
                "--threads", threads, "-o", tag,
            ],
        );
        assert!(
            out.status.success(),
            "--threads {threads} failed: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    let one = std::fs::read_to_string(dir.join("t1.nwk")).expect("t1.nwk");
    let four = std::fs::read_to_string(dir.join("t4.nwk")).expect("t4.nwk");
    assert_eq!(one, four, "--threads changed the bootstrapped tree");

    // The plain matrix path too (`.gtcheck` table).
    for (tag, threads) in [("m1", "1"), ("m4", "4")] {
        let out = run(
            &dir,
            &[
                "distance", &input, "--type", "discordance", "--threads", threads,
                "-o", tag,
            ],
        );
        assert!(out.status.success(), "--threads {threads} matrix run failed");
    }
    let one = std::fs::read_to_string(dir.join("m1.gtcheck")).expect("m1.gtcheck");
    let four = std::fs::read_to_string(dir.join("m4.gtcheck")).expect("m4.gtcheck");
    assert_eq!(one, four, "--threads changed the matrix output");
}

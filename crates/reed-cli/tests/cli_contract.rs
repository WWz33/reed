//! CLI contract gates (doc_v4 06-roadmap W1): the help text, the verbosity
//! contract and the file-only `index` input are user-visible promises, so
//! they are locked here instead of drifting silently (01 F4-3, F4-4, F4-9,
//! F4-13).

use std::process::{Command, Stdio};

const HEADER: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
";

fn fixture() -> Vec<u8> {
    let mut text = String::from(HEADER);
    for i in 0..10_u32 {
        text.push_str(&format!(
            "Chr01\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
            i + 1
        ));
    }
    let mut deflater = reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater.deflate_block(block, &mut out).expect("invariant: block");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    out
}

fn run(args: &[&str]) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_reed"));
    command.args(args).stdin(Stdio::null()).stderr(Stdio::piped());
    command.output().expect("invariant: spawn reed")
}

fn help(args: &[&str]) -> String {
    let out = run(args);
    assert!(out.status.success(), "--help must succeed");
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn temp_dir(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-cli-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    dir
}

/// 01 F4-13: the `<INPUT>` help is one constant for `filter` and `stats`.
#[test]
fn input_help_is_shared() {
    const LINE: &str = "Input VCF/VCF.gz/BCF, or `-` for stdin";
    for sub in ["filter", "stats"] {
        let text = help(&[sub, "--help"]);
        assert!(text.contains(LINE), "{sub}: {text}");
    }
}

/// W0 + F4-9: the new threshold spellings are the only ones on the surface;
/// the drop-direction flags are gone and `-g` no longer abbreviates `--pop`.
#[test]
fn help_lists_the_renamed_flags_only() {
    let filter = help(&["filter", "--help"]);
    for flag in ["--geno", "--mind"] {
        assert!(filter.contains(flag), "filter help lacks {flag}");
    }
    assert!(!filter.contains("--max-missing"), "the withdrawn --max-missing is still listed");

    let stats = help(&["stats", "--help"]);
    assert!(stats.contains("--geno"), "stats shares the site flags");
    assert!(!stats.contains("--max-missing"), "stats must not list --max-missing");
    assert!(!stats.contains("-g, --pop"), "-g must not abbreviate --pop");
    assert!(stats.contains("--pop"), "the long spelling stays");
    assert!(stats.contains("--populations"), "the visible alias stays");
}

/// 01 F4-4: `index` needs a seekable file; stdin fails with that reason.
#[test]
fn index_rejects_stdin_with_a_clear_error() {
    let index_help = help(&["index", "--help"]);
    assert!(index_help.contains("file only"), "{index_help}");
    assert!(!index_help.contains("for stdin"), "{index_help}");

    let out = run(&["index", "-"]);
    assert_eq!(out.status.code(), Some(1), "driver errors exit 1");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("stdin is not supported"), "{stderr}");
    assert!(!stderr.contains("No such file"), "{stderr}");
}

/// 01 F4-3 / 03 §3.3: one verbosity contract. 0 is silent, 1 prints the
/// summary, 2 adds the timing block, 3 is accepted, 4 is rejected.
#[test]
fn verbosity_levels_are_uniform() {
    let dir = temp_dir("verbosity");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let out_path = dir.join("out.vcf.gz");
    let out_arg = out_path.display().to_string();
    let base: Vec<String> = vec![
        "filter".into(),
        input.display().to_string(),
        "-Oz".into(),
        "-o".into(),
        out_arg.clone(),
    ];
    let run_with = |level: &str| {
        let mut argv: Vec<String> = base.clone();
        argv.push("--verbosity".into());
        argv.push(level.into());
        let refs: Vec<&str> = argv.iter().map(String::as_str).collect();
        run(&refs)
    };
    let silent = run_with("0");
    assert!(silent.status.success());
    assert!(silent.stderr.is_empty(), "level 0 must be silent");

    let summary = run_with("1");
    let text = String::from_utf8_lossy(&summary.stderr);
    assert!(text.contains("reed filter:"), "{text}");
    assert!(text.contains("funnel:"), "level 1 lists the funnel: {text}");
    assert!(!text.contains("total_ms"), "level 1 has no timing: {text}");

    let timing = run_with("2");
    let text = String::from_utf8_lossy(&timing.stderr);
    assert!(text.contains("total_ms"), "level 2 prints the timing: {text}");
    assert!(text.contains("phase.open_ms"), "{text}");

    assert!(run_with("3").status.success(), "level 3 is accepted");
    let bad = run_with("4");
    assert_eq!(bad.status.code(), Some(2), "unimplemented levels are rejected");
    assert!(
        String::from_utf8_lossy(&bad.stderr).contains("unsupported --verbosity"),
        "{}",
        String::from_utf8_lossy(&bad.stderr)
    );

    // `stats` follows the same contract: silent at 0, timing at 2.
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopB\n").expect("invariant: groups");
    let stats = |level: &str, prefix: &str| -> String {
        let out = Command::new(env!("CARGO_BIN_EXE_reed"))
            .args([
                "stats",
                input.to_str().unwrap(),
                "--pop",
                pops.to_str().unwrap(),
                "--pi",
                "--variant-only",
                "-o",
                prefix,
                "--verbosity",
                level,
            ])
            .stdin(Stdio::null())
            .stderr(Stdio::piped())
            .output()
            .expect("invariant: spawn reed");
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        String::from_utf8_lossy(&out.stderr).into_owned()
    };
    let quiet = stats("0", &dir.join("s0").display().to_string());
    assert!(quiet.is_empty(), "stats level 0 must be silent: {quiet}");
    let timed = stats("2", &dir.join("s2").display().to_string());
    assert!(timed.contains("phase.run_ms"), "{timed}");
    let _ = std::fs::remove_dir_all(&dir);
}

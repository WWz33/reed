//! CLI contract gates (doc_v4 06-roadmap W1): the help text, the verbosity
//! contract and the file-only `index` input are user-visible promises, so
//! they are locked here instead of drifting silently (01 F4-3, F4-4, F4-9,
//! F4-13).

use std::fmt::Write as _;
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
        let _ = writeln!(text, "Chr01\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0", i + 1);
    }
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block");
    }
    out.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    out
}

fn run(args: &[&str]) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_reed"));
    command
        .args(args)
        .stdin(Stdio::null())
        .stderr(Stdio::piped());
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
    assert!(
        !filter.contains("--max-missing"),
        "the withdrawn --max-missing is still listed"
    );

    let stats = help(&["stats", "--help"]);
    assert!(stats.contains("--geno"), "stats shares the site flags");
    assert!(
        !stats.contains("--max-missing"),
        "stats must not list --max-missing"
    );
    assert!(!stats.contains("-g, --pop"), "-g must not abbreviate --pop");
    assert!(stats.contains("--pop"), "the long spelling stays");
    assert!(stats.contains("--populations"), "the visible alias stays");
}

/// X4: `--sfs`/`--fst-global` are on the surface; `--fst-global` needs
/// `--fst` (clap enforces it) and `--sfs` needs `--proj`.
#[test]
fn help_lists_the_new_statistics_flags() {
    let popg = help(&["popg", "--help"]);
    for flag in ["--sfs", "--proj", "--sfs-int", "--fst-global"] {
        assert!(popg.contains(flag), "popg help lacks {flag}");
    }

    // `--fst-global` without `--fst`: clap's `requires` fails at parse time.
    let out = run(&["popg", "x.vcf", "--fst-global"]);
    assert!(!out.status.success(), "--fst-global alone must fail");

    // `--sfs` without `--proj`: the translator names the missing flag; a
    // missing flag the operator can add is a usage error (cli.md §9 → 2).
    let dir = temp_dir("sfs-proj");
    let input = dir.join("in.vcf");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let out = run(&[
        "popg",
        input.to_str().unwrap(),
        "--sfs",
        "-o",
        dir.join("p").to_str().unwrap(),
    ]);
    assert_eq!(out.status.code(), Some(2), "usage errors exit 2");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("--proj"), "{stderr}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// 01 F4-4: `index` needs a seekable file; stdin fails with that reason.
#[test]
fn index_rejects_stdin_with_a_clear_error() {
    let index_help = help(&["index", "--help"]);
    assert!(index_help.contains("file only"), "{index_help}");
    assert!(!index_help.contains("for stdin"), "{index_help}");

    let out = run(&["index", "-"]);
    assert_eq!(out.status.code(), Some(2), "stdin is a usage error");
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("stdin is not supported"), "{stderr}");
    assert!(!stderr.contains("No such file"), "{stderr}");
}

/// cli.md §9 through Y6 batch D: `reed::Error::kind()` decides the exit
/// code. Usage (2): a windowed statistic without --pop, --fst-global with
/// hudson, --sfs without --proj. Runtime (1): the input file does not exist.
#[test]
fn exit_codes_split_usage_from_runtime() {
    let dir = temp_dir("exit-codes");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let inp = input.to_str().unwrap();
    let pops = dir.join("pops.tsv");
    std::fs::write(&pops, "s0\tpop1\ns1\tpop2\n").expect("invariant: pops");
    let pop_arg = pops.to_str().unwrap();

    // Usage 2: a windowed statistic needs --pop (population missing, §9).
    let out = run(&["popg", inp, "--pi"]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "missing population: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Usage 2: hudson has no global form.
    let out = run(&[
        "popg",
        inp,
        "--fst-global",
        "--fst",
        "hudson",
        "--pop",
        pop_arg,
    ]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "fst-global+hudson: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Usage 2: --sfs needs --proj.
    let out = run(&["popg", inp, "--sfs", "--pop", pop_arg]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "sfs without proj: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Runtime 1: the input itself is missing (a QC flag keeps the request
    // valid so the failure is the missing file, not the request).
    let out = run(&[
        "stats",
        dir.join("nope.vcf.gz").to_str().unwrap(),
        "--freq",
        "-o",
        dir.join("p").to_str().unwrap(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(1),
        "missing input: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let _ = std::fs::remove_dir_all(&dir);
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
    assert!(
        text.contains("total_ms"),
        "level 2 prints the timing: {text}"
    );
    assert!(text.contains("phase.open_ms"), "{text}");

    assert!(run_with("3").status.success(), "level 3 is accepted");
    let bad = run_with("4");
    assert_eq!(
        bad.status.code(),
        Some(2),
        "unimplemented levels are rejected"
    );
    assert!(
        String::from_utf8_lossy(&bad.stderr).contains("unsupported --verbosity"),
        "{}",
        String::from_utf8_lossy(&bad.stderr)
    );

    // `popg` follows the same contract: silent at 0, timing at 2.
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopB\n").expect("invariant: groups");
    let stats = |level: &str, prefix: &str| -> String {
        let out = Command::new(env!("CARGO_BIN_EXE_reed"))
            .args([
                "popg",
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
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8_lossy(&out.stderr).into_owned()
    };
    let quiet = stats("0", &dir.join("s0").display().to_string());
    assert!(quiet.is_empty(), "popg level 0 must be silent: {quiet}");
    let timed = stats("2", &dir.join("s2").display().to_string());
    assert!(timed.contains("phase.run_ms"), "{timed}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// 01 F4-8: the body is the deterministic artefact; the recorded argv
/// (binary path, `--threads`, …) is environment information and can be
/// suppressed with `--no-command`.
#[test]
fn no_command_drops_only_the_argv_line() {
    let dir = temp_dir("nocommand");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let out = dir.join("out.vcf");
    let run_with = |extra: &[&str]| {
        let mut argv = vec![
            "filter",
            input.to_str().unwrap(),
            "-Ov",
            "-o",
            out.to_str().unwrap(),
            "--verbosity",
            "0",
        ];
        argv.extend_from_slice(extra);
        let o = run(&argv);
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        std::fs::read_to_string(&out).expect("invariant: read")
    };
    assert!(
        run_with(&[]).contains("##reed_command="),
        "argv is recorded by default"
    );
    let text = run_with(&["--no-command"]);
    assert!(
        !text.contains("##reed_command="),
        "--no-command drops the argv line"
    );
    assert!(
        text.contains("##reed_version="),
        "the other provenance lines stay"
    );
    assert!(
        !run_with(&["--no-version"]).contains("##reed_"),
        "--no-version drops all"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// 01 F4-10: report mode accepts `-t/-T` — the region predicate gates the
/// records in the same pass (it used to error "not supported yet").
#[test]
fn report_mode_accepts_targets() {
    let dir = temp_dir("report-targets");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let prefix = dir.join("p").display().to_string();
    let o = run(&[
        "stats",
        input.to_str().unwrap(),
        "--freq",
        "--summary",
        "-t",
        "Chr01:1-5",
        "-o",
        &prefix,
        "--verbosity",
        "0",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let afreq = std::fs::read_to_string(format!("{prefix}.afreq")).expect("invariant: afreq");
    assert_eq!(
        afreq.lines().filter(|l| !l.starts_with('#')).count(),
        5,
        "only the target range's records: {afreq}"
    );
    let stats = std::fs::read_to_string(format!("{prefix}.stats")).expect("invariant: stats");
    assert!(
        stats.contains("number of records:\t5"),
        "the summary sees the same five records"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// Y6 batch G: `reed index` really builds an index next to the input, and
/// `reed popg` really writes the window-layout product — the minimal
/// end-to-end fixture runs for the two commands the round-trip tests cannot
/// reach (binary-level only here).
#[test]
fn index_and_stats_run_on_the_minimal_fixture() {
    let dir = temp_dir("rt-fixture");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let inp = input.to_str().unwrap();

    let out = run(&["index", inp]);
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let csi = dir.join("in.vcf.gz.csi");
    assert!(csi.is_file(), "index must write in.vcf.gz.csi");
    assert!(csi.metadata().expect("invariant: stat").len() > 0);

    let pops = dir.join("pops.tsv");
    std::fs::write(&pops, "s0\tA\ns1\tA\n").expect("invariant: pops");
    let prefix = dir.join("out");
    let out = run(&[
        "popg",
        inp,
        "--pop",
        pops.to_str().unwrap(),
        "--pi",
        "--window",
        "5",
        "--step",
        "5",
        "--variant-only",
        "-o",
        prefix.to_str().unwrap(),
    ]);
    assert!(
        out.status.success(),
        "stats failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let pi = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi file");
    let header = pi.lines().next().expect("invariant: header");
    assert!(
        header.starts_with("pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi"),
        "{header}"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The hapmap fixture for the `convert --vcf` gates: one-letter IUPAC cells
/// over two samples (doc_v13 §3.12).
const HAPMAP: &str = "\
##fileformat=hapmap
rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	s0	s1
chr1:10	A/G	chr1	10	+	NA	NA	NA	NA	NA	NA	A	R
chr1:20	C/T	chr1	20	+	NA	NA	NA	NA	NA	NA	T	N
";

/// v13 §3.12: `convert --vcf` translates a hapmap back to VCF — the records
/// carry the declared REF/ALT and the IUPAC-expanded GT, the header has the
/// provenance lines and no `##contig`, and `-Oz` writes BGZF.
#[test]
fn convert_hapmap_to_vcf_outputs_records_and_provenance() {
    let dir = temp_dir("convert-vcf");
    let input = dir.join("in.hmp.txt");
    std::fs::write(&input, HAPMAP).expect("invariant: fixture");
    let out = dir.join("out.vcf");
    let o = run(&[
        "convert",
        input.to_str().unwrap(),
        "--vcf",
        "-o",
        out.to_str().unwrap(),
        "--verbosity",
        "0",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let text = std::fs::read_to_string(&out).expect("invariant: read");
    assert!(text.starts_with("##fileformat=VCFv4.3\n"), "{text}");
    assert!(text.contains("##reed_version="), "{text}");
    assert!(text.contains("##reed_command="), "{text}");
    assert!(!text.contains("##contig"), "{text}");
    assert!(
        text.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n"),
        "{text}"
    );
    assert!(
        text.contains("chr1\t10\tchr1:10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"),
        "{text}"
    );
    assert!(
        text.contains("chr1\t20\tchr1:20\tC\tT\t.\tPASS\t.\tGT\t1/1\t./.\n"),
        "{text}"
    );

    // --no-version drops the provenance block; the records stay.
    let bare = dir.join("bare.vcf");
    let o = run(&[
        "convert",
        input.to_str().unwrap(),
        "--vcf",
        "--no-version",
        "-o",
        bare.to_str().unwrap(),
        "--verbosity",
        "0",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let bare_text = std::fs::read_to_string(&bare).expect("invariant: read");
    assert!(!bare_text.contains("##reed_"), "{bare_text}");
    assert!(bare_text.contains("\tGT\t0/0\t0/1\n"), "{bare_text}");

    // -Oz writes BGZF (gzip magic).
    let gz = dir.join("out.vcf.gz");
    let o = run(&[
        "convert",
        input.to_str().unwrap(),
        "--vcf",
        "-Oz",
        "-o",
        gz.to_str().unwrap(),
        "--verbosity",
        "0",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let bytes = std::fs::read(&gz).expect("invariant: read");
    assert_eq!(&bytes[..2], &[0x1f, 0x8b], "-Oz must be BGZF");

    let _ = std::fs::remove_dir_all(&dir);
}

/// v13 §1.1 rule 6: `--vcf` reads a hapmap, so it is exclusive with every
/// VCF-input target, takes no VCF-input knobs, and rejects the BCF output
/// spellings.
#[test]
fn convert_vcf_rejects_the_vcf_input_surface() {
    let dir = temp_dir("convert-vcf-gates");
    let input = dir.join("in.hmp.txt");
    std::fs::write(&input, HAPMAP).expect("invariant: fixture");
    let path = input.to_str().unwrap();

    let cases: &[&[&str]] = &[
        &["convert", path, "--vcf", "--phylip"],
        &["convert", path, "--vcf", "--012"],
        &["convert", path, "--vcf", "--sfs", "unfolded"],
        &["convert", path, "--vcf", "--outgroup", "s0"],
        &["convert", path, "--vcf", "-s", "s0"],
        &["convert", path, "--vcf", "-Ob"],
        &["convert", path, "--vcf", "-Ou"],
        &["convert", path, "--phylip", "-Ov"],
        &["convert", path, "--phylip", "--no-version"],
    ];
    for args in cases {
        let o = run(args);
        assert_eq!(
            o.status.code(),
            Some(2),
            "{args:?} must be a usage error: {}",
            String::from_utf8_lossy(&o.stderr)
        );
    }

    // The two mixed-target messages name the reason.
    let o = run(&["convert", path, "--vcf", "--phylip"]);
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(stderr.contains("cannot mix --vcf"), "{stderr}");
    let o = run(&["convert", path, "--vcf", "-Ob"]);
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(stderr.contains("-O v or -O z"), "{stderr}");

    // A non-hapmap input is a data error (exit 1), not a usage error — the
    // same split `filter`/`stats` make when handed a non-VCF file.
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n").expect("fixture");
    let o = run(&[
        "convert",
        vcf.to_str().unwrap(),
        "--vcf",
        "--verbosity",
        "0",
    ]);
    assert_eq!(o.status.code(), Some(1), "wrong-format input is runtime");
    let stderr = String::from_utf8_lossy(&o.stderr);
    assert!(stderr.contains("not a hapmap file"), "{stderr}");

    let _ = std::fs::remove_dir_all(&dir);
}

/// v13 §3.12: `reed convert --vcf` streams the VCF to stdout when no `-o`
/// is given (bcftools' `-Ov` shape), so a hapmap pipes into a VCF pipeline.
#[test]
fn convert_vcf_streams_to_stdout_without_o() {
    let dir = temp_dir("convert-vcf-stdout");
    let input = dir.join("in.hmp.txt");
    std::fs::write(&input, HAPMAP).expect("invariant: fixture");
    let o = run(&[
        "convert",
        input.to_str().unwrap(),
        "--vcf",
        "--verbosity",
        "0",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let text = String::from_utf8_lossy(&o.stdout);
    assert!(text.starts_with("##fileformat=VCFv4.3\n"), "{text}");
    assert!(
        text.contains("chr1\t10\tchr1:10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n"),
        "{text}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 usage gates: `--levels` with no level, `--step` beside
/// `--bed-file`, `--idspace-to` beside `--set-missing-var-ids`, a repeated
/// `--pbs` group and a contig-rename collision are refused up front, and a
/// `POS=0` record is skipped as malformed instead of crashing the debug
/// build (doc_v15 D15-21/22/25/27/30/33).
#[test]
fn v15_usage_gates() {
    let dir = temp_dir("v15-usage");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##contig=<ID=chr2,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr1\t0\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr2\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    )
    .expect("invariant: write");
    let pops = dir.join("pops.tsv");
    std::fs::write(&pops, "s0\tA\n").expect("invariant: write");
    let bed = dir.join("w.bed");
    std::fs::write(&bed, "chr1\t0\t100\n").expect("invariant: write");
    let chrs = dir.join("chrs.txt");
    std::fs::write(&chrs, "chr1\tX\nchr2\tX\n").expect("invariant: write");
    let out = dir.join("out");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, p, b, c, o) = (sv(&vcf), sv(&pops), sv(&bed), sv(&chrs), sv(&out));

    let case = |args: Vec<&str>, needle: &str| {
        let result = run(&args);
        assert_eq!(
            result.status.code(),
            Some(2),
            "usage error expected: {args:?}"
        );
        let err = String::from_utf8_lossy(&result.stderr);
        assert!(err.contains(needle), "want {needle:?} in: {err}");
    };
    case(
        vec!["popg", &v, "--pop", &p, "--pi", "--levels", "", "-o", &o],
        "--levels lists no level",
    );
    case(
        vec![
            "popg",
            &v,
            "--pop",
            &p,
            "--pi",
            "--bed-file",
            &b,
            "--step",
            "5",
            "-o",
            &o,
        ],
        "cannot be used with",
    );
    case(
        vec![
            "rename",
            &v,
            "--set-missing-var-ids",
            "a:b",
            "--idspace-to",
            "_",
        ],
        "--idspace-to cannot be combined",
    );
    case(
        vec!["popg", &v, "--pop", &p, "--pbs", "A,A,A", "-o", &o],
        "repeats (want three distinct populations)",
    );
    case(vec!["rename", &v, "--rename-chrs", &c], "map both");

    // The POS=0 record must be skipped with a counter, never a panic.
    let result = run(&["convert", &v, "--012", "-o", &o]);
    assert!(result.status.success(), "convert runs: {result:?}");
    let pos = std::fs::read_to_string(format!("{o}.012.pos")).expect("invariant: pos file");
    assert!(!pos.contains("\t0\n"), "POS=0 must be skipped: {pos}");
    assert!(pos.contains("chr1\t100"));

    // `--sfs --proj` above a group's chromosome count is refused up front,
    // through the shared plan validation (doc_v15 D15-9; popg and convert).
    let result = run(&[
        "convert", &v, "--pop", &p, "--sfs", "--proj", "40", "-o", &o,
    ]);
    assert_eq!(result.status.code(), Some(2));
    let err = String::from_utf8_lossy(&result.stderr);
    assert!(err.contains("exceeds population"), "{err}");

    // One sample assigned to two populations is a usage error, not a silent
    // last-group-wins (doc_v15 D15-10).
    let dup = dir.join("dup.tsv");
    std::fs::write(&dup, "s0\tA\ns0\tB\n").expect("invariant: write");
    let result = run(&["popg", &v, "--pop", &sv(&dup), "--pi", "-o", &o]);
    assert_eq!(result.status.code(), Some(2));
    let err = String::from_utf8_lossy(&result.stderr);
    assert!(err.contains("assigned to both"), "{err}");
}

/// doc_v15: `popg` validates `-o` before the input is opened (D15-13), and
/// the bare route honors the invariant-site contract (D15-15).
#[test]
fn popg_validates_before_running() {
    let dir = temp_dir("v15-popg-order");
    let missing = dir.join("does-not-exist.vcf");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();

    // `-o` is checked before the input is even opened: the missing file
    // never gets the chance to produce its own error.
    let result = run(&["popg", &sv(&missing), "--pi"]);
    assert_eq!(result.status.code(), Some(2));
    let err = String::from_utf8_lossy(&result.stderr);
    assert!(err.contains("pass -o PREFIX"), "{err}");

    // Bare popg on an all-variant file refuses without `--variant-only`,
    // exactly like the windowed route (D15-15).
    let vcf = dir.join("allvar.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    )
    .expect("invariant: write");
    let pops = dir.join("p.tsv");
    std::fs::write(&pops, "s0\tA\n").expect("invariant: write");
    let result = run(&["popg", &sv(&vcf), "--pop", &sv(&pops)]);
    assert_eq!(result.status.code(), Some(2), "{result:?}");
    let err = String::from_utf8_lossy(&result.stderr);
    assert!(err.contains("no invariant sites"), "{err}");
    let ok = run(&["popg", &sv(&vcf), "--pop", &sv(&pops), "--variant-only"]);
    assert!(ok.status.success(), "{ok:?}");
}

/// doc_v15 D15-14: with the pipe's read end closed, the bare `popg` TSV
/// write must fail as an ordinary error, not panic (exit 101).
#[test]
fn popg_broken_pipe_is_not_a_panic() {
    let dir = temp_dir("v15-epipe");
    let mut text = String::from("##fileformat=VCFv4.2\n");
    for c in 0..20 {
        let _ = writeln!(text, "##contig=<ID=ctg{c:03},length=1000>");
    }
    text.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n");
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n");
    for c in 0..20 {
        let _ = writeln!(
            text,
            "ctg{c:03}\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/1"
        );
    }
    let vcf = dir.join("many.vcf");
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("p.tsv");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: write");

    let mut child = Command::new(env!("CARGO_BIN_EXE_reed"))
        .args([
            "popg",
            vcf.to_str().expect("invariant: utf8"),
            "--pop",
            pops.to_str().expect("invariant: utf8"),
            "--variant-only",
        ])
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("invariant: spawn reed");
    drop(child.stdout.take());
    let out = child.wait_with_output().expect("invariant: wait");
    assert_ne!(out.status.code(), Some(101), "must not panic: {out:?}");
    assert_eq!(out.status.code(), Some(1), "{out:?}");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("Broken pipe"), "{err}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-29: a gzipped hapmap piped through stdin sniffs and decodes
/// like the same bytes given as a path (the stdin branch used to hand the
/// raw gzip stream to the text parser).
#[test]
fn gzipped_hapmap_on_stdin_decodes() {
    use std::io::Write as _;

    let dir = temp_dir("v15-hmp-gz");
    let hmp = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\ts0\ts1\n\
               rs1\tA/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNN\tAA\tAT\n\
               rs2\tC/G\tchr1\t200\t+\tNA\tNA\tNA\tNA\tNA\tNN\tCC\tCG\n";
    // BGZF is gzip-compatible: the stdin sniff only looks for the magic.
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut compressed = Vec::new();
    for block in hmp.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut compressed)
            .expect("invariant: block");
    }
    compressed.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);

    let out_path = dir.join("out.vcf");
    let mut child = Command::new(env!("CARGO_BIN_EXE_reed"))
        .args([
            "convert",
            "-",
            "--vcf",
            "-Ov",
            "-o",
            out_path.to_str().expect("invariant: utf8"),
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("invariant: spawn reed");
    child
        .stdin
        .as_mut()
        .expect("invariant: stdin pipe")
        .write_all(&compressed)
        .expect("invariant: write stdin");
    drop(child.stdin.take());
    let out = child.wait_with_output().expect("invariant: wait");
    assert!(out.status.success(), "{out:?}");
    let text = std::fs::read_to_string(&out_path).expect("invariant: read");
    let records = text.lines().filter(|l| !l.starts_with('#')).count();
    assert_eq!(records, 2, "{text}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-18: adding `--f2` must not silently switch the π/dxy/θ
/// convention to variant-only — the π rows stay byte-identical — and an
/// all-variant file stays refused with or without `--f2`.
#[test]
fn f2_does_not_switch_the_pi_convention() {
    let dir = temp_dir("v15-f2");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();

    let header = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=1000>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n";
    let variants = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n\
        chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\n";
    let invariant = "chr1\t300\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n";
    let allsites = dir.join("allsites.vcf");
    std::fs::write(&allsites, format!("{header}{variants}{invariant}")).expect("invariant: write");
    let all_variant = dir.join("allvar.vcf");
    std::fs::write(&all_variant, format!("{header}{variants}")).expect("invariant: write");
    let pops = dir.join("p.tsv");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: write");
    let (v, var, p) = (sv(&allsites), sv(&all_variant), sv(&pops));

    // Identical π rows with and without `--f2` on a file with invariant
    // sites (the all-sites check passes either way).
    let out_plain = sv(&dir.join("plain"));
    let out_f2 = sv(&dir.join("withf2"));
    let r = run(&["popg", &v, "--pop", &p, "--pi", "-o", &out_plain]);
    assert!(r.status.success(), "{r:?}");
    let r = run(&[
        "popg", &v, "--pop", &p, "--pi", "--f2", "A,B", "-o", &out_f2,
    ]);
    assert!(r.status.success(), "{r:?}");
    let plain = std::fs::read_to_string(dir.join("plain_pi.txt")).expect("invariant: read");
    let with_f2 = std::fs::read_to_string(dir.join("withf2_pi.txt")).expect("invariant: read");
    assert_eq!(plain, with_f2, "pi must not change when --f2 is added");

    // The all-variant file is refused by the invariant-site check in both
    // shapes: `--f2` must not silently flip the convention.
    for (i, extra) in [Vec::<&str>::new(), vec!["--f2", "A,B"]]
        .into_iter()
        .enumerate()
    {
        let out = sv(&dir.join(format!("allvar-out{i}")));
        let mut args = vec!["popg", var.as_str(), "--pop", p.as_str(), "--pi"];
        args.extend(extra);
        args.extend(["-o", out.as_str()]);
        let r = run(&args);
        assert_eq!(r.status.code(), Some(2), "{r:?}");
        let err = String::from_utf8_lossy(&r.stderr);
        assert!(err.contains("no invariant sites"), "{err}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-39: `-o -` streams to stdout instead of creating a file
/// literally named `-` in the working directory.
#[test]
fn dash_output_streams_to_stdout() {
    let dir = temp_dir("v15-dash-out");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    )
    .expect("invariant: write");
    let out = Command::new(env!("CARGO_BIN_EXE_reed"))
        .args([
            "filter",
            vcf.to_str().expect("invariant: utf8"),
            "-o",
            "-",
            "-Ov",
        ])
        .current_dir(&dir)
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("invariant: spawn reed");
    assert!(out.status.success(), "{out:?}");
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(
        text.contains("#CHROM"),
        "records must stream to stdout: {text}"
    );
    assert!(
        !dir.join("-").exists(),
        "no file literally named '-' may be created"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v20 R5: `-o` equal to the input is refused before the output is
/// opened — the writers `File::create` the output after the input is
/// already open, which truncated the input mid-read (unrecoverable loss).
/// Every gate also proves the input file is byte-identical afterwards.
#[test]
fn r5_filter_refuses_output_equal_to_input() {
    let dir = temp_dir("r5-filter");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let before = std::fs::read(&input).expect("invariant: read");
    let inp = input.to_str().expect("invariant: utf8");
    let out = run(&["filter", inp, "-o", inp, "-Oz", "--verbosity", "0"]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(2),
        "output == input is a usage error: {stderr}"
    );
    assert!(stderr.contains("would overwrite the input"), "{stderr}");
    let after = std::fs::read(&input).expect("invariant: read");
    assert_eq!(after, before, "the input must survive the refused run");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v20 R5 (rename flavor): same gate for `reed rename -o INPUT` — the
/// plain-VCF rewrite used to truncate the input it was still reading.
#[test]
fn r5_rename_refuses_output_equal_to_input() {
    let dir = temp_dir("r5-rename");
    let input = dir.join("in.vcf");
    let mut text = String::from(HEADER);
    for i in 0..10_u32 {
        let _ = writeln!(text, "Chr01\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0", i + 1);
    }
    std::fs::write(&input, text).expect("invariant: fixture");
    let before = std::fs::read(&input).expect("invariant: read");
    let inp = input.to_str().expect("invariant: utf8");
    let out = run(&["rename", inp, "-o", inp, "--clear-ids", "--verbosity", "0"]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(2),
        "output == input is a usage error: {stderr}"
    );
    assert!(stderr.contains("would overwrite the input"), "{stderr}");
    let after = std::fs::read(&input).expect("invariant: read");
    assert_eq!(after, before, "the input must survive the refused run");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v20 R5 (convert flavor): `convert --vcf` writes `-o` verbatim, so
/// `-o HAPMAP` truncated the hapmap being converted.
#[test]
fn r5_convert_refuses_output_equal_to_input() {
    let dir = temp_dir("r5-convert");
    let input = dir.join("in.hmp.txt");
    std::fs::write(&input, HAPMAP).expect("invariant: fixture");
    let before = std::fs::read(&input).expect("invariant: read");
    let inp = input.to_str().expect("invariant: utf8");
    let out = run(&["convert", inp, "--vcf", "-o", inp, "--verbosity", "0"]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(2),
        "output == input is a usage error: {stderr}"
    );
    assert!(stderr.contains("would overwrite the input"), "{stderr}");
    let after = std::fs::read(&input).expect("invariant: read");
    assert_eq!(after, before, "the input must survive the refused run");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v20 R5 (index flavor): `index -f -o INPUT` overwrote the VCF with
/// its own index — `-f` was the destructive flavor.
#[test]
fn r5_index_refuses_output_equal_to_input() {
    let dir = temp_dir("r5-index");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let before = std::fs::read(&input).expect("invariant: read");
    let inp = input.to_str().expect("invariant: utf8");
    let out = run(&["index", inp, "-f", "-o", inp]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(2),
        "output == input is a usage error: {stderr}"
    );
    assert!(stderr.contains("would overwrite the input"), "{stderr}");
    let after = std::fs::read(&input).expect("invariant: read");
    assert_eq!(after, before, "the input must survive the refused run");
    let _ = std::fs::remove_dir_all(&dir);
}

/// Compresses `text` into a BGZF file at `path`.
fn write_bgzf(path: &std::path::Path, text: &str) {
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block");
    }
    out.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    std::fs::write(path, out).expect("invariant: write");
}

/// Reads a BGZF file back as text.
fn read_bgzf(path: &std::path::Path) -> String {
    let bytes = std::fs::read(path).expect("invariant: read");
    let mut inflater = reed::bgzf::inflate::Inflater::new();
    let mut text = Vec::new();
    let mut rest = &bytes[..];
    while let Some(total) = reed::bgzf::block::parse_block_size(rest) {
        inflater
            .inflate_block(&rest[..total], &mut text)
            .expect("invariant: inflate");
        rest = &rest[total..];
    }
    String::from_utf8(text).expect("invariant: utf8")
}

/// doc_v22: three real-world VCF shapes bcftools accepts but reed rejected.
#[test]
fn v22_accepts_sites_only_haploid_and_trailing_tabs() {
    let dir = temp_dir("v22-shapes");

    // 1. Sites-only VCF (8 columns, `bcftools view -G` output): filter and
    //    index must work and the output must stay 8 columns.
    let sites = dir.join("sites.vcf.gz");
    write_bgzf(
        &sites,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\n\
         chr1\t20\t.\tC\tT\t60\tPASS\t.\n",
    );
    let sites_out = dir.join("sites.out.vcf.gz");
    let out = run(&[
        "filter",
        sites.to_str().expect("invariant: utf8"),
        "-Oz",
        "-o",
        sites_out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = read_bgzf(&sites_out);
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        assert_eq!(line.split('\t').count(), 8, "{line}");
    }
    let indexed = run(&["index", sites.to_str().expect("invariant: utf8")]);
    assert!(
        indexed.status.success(),
        "sites-only index: {}",
        String::from_utf8_lossy(&indexed.stderr)
    );
    // Record-plane edits work on sites-only input too (V22-R8): only the
    // first five columns are read, the rest is copied verbatim.
    let renamed = dir.join("sites.renamed.vcf.gz");
    let rename = run(&[
        "rename",
        sites.to_str().expect("invariant: utf8"),
        "--set-id",
        "%CHROM_%POS",
        "-Oz",
        "-o",
        renamed.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        rename.status.success(),
        "sites-only rename: {}",
        String::from_utf8_lossy(&rename.stderr)
    );
    let renamed_text = read_bgzf(&renamed);
    let ids: Vec<&str> = renamed_text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.split('\t').nth(2).expect("invariant: ID column"))
        .collect();
    assert_eq!(ids, ["chr1_10", "chr1_20"]);

    // 2. Haploid calls (mtDNA / male X): the default policy rejects them,
    //    `--ploidy allow-haploid` accepts them.
    let haploid = dir.join("haploid.vcf.gz");
    write_bgzf(
        &haploid,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1\t0\n",
    );
    let haploid_out = dir.join("haploid.out.vcf.gz");
    let rejected = run(&[
        "filter",
        haploid.to_str().expect("invariant: utf8"),
        "-Oz",
        "-o",
        haploid_out.to_str().expect("invariant: utf8"),
    ]);
    assert_eq!(rejected.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&rejected.stderr).contains("haploid"),
        "{}",
        String::from_utf8_lossy(&rejected.stderr)
    );
    let allowed = run(&[
        "filter",
        haploid.to_str().expect("invariant: utf8"),
        "--ploidy",
        "allow-haploid",
        "-Oz",
        "-o",
        haploid_out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        allowed.status.success(),
        "{}",
        String::from_utf8_lossy(&allowed.stderr)
    );
    let haploid_text = read_bgzf(&haploid_out);
    let haploid_rows: Vec<&str> = haploid_text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(haploid_rows.len(), 1);
    // 8 site columns + FORMAT + 2 samples.
    assert_eq!(haploid_rows[0].split('\t').count(), 11);
    assert!(
        haploid_rows[0].ends_with("\tGT\t1\t0"),
        "{}",
        haploid_rows[0]
    );

    // 3. A record with a trailing empty sample column (row ends with `\t`):
    //    bcftools ignores it; reed must not reject the whole file.
    let trailing = dir.join("trailing.vcf.gz");
    write_bgzf(
        &trailing,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t\n",
    );
    let trailing_out = dir.join("trailing.out.vcf.gz");
    let ok = run(&[
        "filter",
        trailing.to_str().expect("invariant: utf8"),
        "-Oz",
        "-o",
        trailing_out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        ok.status.success(),
        "{}",
        String::from_utf8_lossy(&ok.stderr)
    );
    let trailing_text = read_bgzf(&trailing_out);
    let trailing_rows: Vec<&str> = trailing_text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    assert_eq!(trailing_rows.len(), 1);
    assert_eq!(
        trailing_rows[0].trim_end_matches('\t').split('\t').count(),
        11,
        "the trailing tab must not become a sample column: {}",
        trailing_rows[0]
    );
    assert!(
        trailing_rows[0].contains("\tGT\t0/1\t1/1"),
        "{}",
        trailing_rows[0]
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// The `--min-dp` flag itself (not a hand-built spec) carries the AD
/// fallback through the CLI translation (V22-R9). A regression that dropped
/// `depth: true` would otherwise only show up as the original
/// "FORMAT field DP is not declared" error.
#[test]
fn v22_cli_min_dp_falls_back_to_ad() {
    let dir = temp_dir("v22-min-dp");
    let vcf = dir.join("ad_only.vcf.gz");
    write_bgzf(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:2,3\n\
         chr1\t20\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:10,20\n",
    );
    let out = dir.join("masked.vcf.gz");
    let run_out = run(&[
        "filter",
        vcf.to_str().expect("invariant: utf8"),
        "--min-dp",
        "10",
        "-Oz",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        run_out.status.success(),
        "{}",
        String::from_utf8_lossy(&run_out.stderr)
    );
    let text = read_bgzf(&out);
    let rows: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].rsplit('\t').next(), Some("./.:2,3"), "{}", rows[0]);
    assert_eq!(
        rows[1].rsplit('\t').next(),
        Some("0/1:10,20"),
        "{}",
        rows[1]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// V22-R2-4/R2-6: `--hist dp` accepts an AD-only header through the CLI
/// gate (it used to demand FORMAT/DP even though the engine has the AD
/// fallback), and the depth/GQ consumers reject non-Integer tags with the
/// usage exit instead of silently treating every value as missing.
#[test]
fn v22_r2_hist_dp_ad_only_and_stats_type_checks() {
    let dir = temp_dir("v22-r2-hist");
    let ad_only = dir.join("ad_only.vcf.gz");
    write_bgzf(
        &ad_only,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:10,20\t0/1:5,5\n\
         chr1\t20\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/0:0,0\t0/1:.\n",
    );
    let prefix = dir.join("ad");
    let out = run(&[
        "stats",
        ad_only.to_str().expect("invariant: utf8"),
        "--hist",
        "dp",
        "-o",
        prefix.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        out.status.success(),
        "the AD fallback must satisfy the --hist dp gate: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let hist = std::fs::read_to_string(format!("{}.dp.hist", prefix.display()))
        .expect("invariant: the dp histogram exists");
    assert!(hist.contains("10\t10\t1\t"), "{hist}");
    assert!(hist.contains("30\t30\t1\t"), "{hist}");
    assert!(
        !hist.lines().any(|line| line.starts_with("0\t0\t")),
        "{hist}"
    );

    // Float depth/GQ tags are rejected by type, not silently decoded as
    // missing. The three fixtures differ only in the offending tag.
    let float_dp = dir.join("float_dp.vcf.gz");
    write_bgzf(
        &float_dp,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"d\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:9.5\n",
    );
    let float_ad = dir.join("float_ad.vcf.gz");
    write_bgzf(
        &float_ad,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=AD,Number=R,Type=Float,Description=\"a\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:9.5,0.5\n",
    );
    let float_gq = dir.join("float_gq.vcf.gz");
    write_bgzf(
        &float_gq,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"q\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/1:9.5\n",
    );
    let cases: [(&str, &std::path::Path, &[&str], &str); 3] = [
        ("float-dp", &float_dp, &["--depth"], "DP is Float"),
        ("float-ad", &float_ad, &["--hist", "dp"], "AD is Float"),
        ("float-gq", &float_gq, &["--hist", "gq"], "GQ is Float"),
    ];
    for (tag, vcf, report, expected) in cases {
        let prefix = dir.join(tag);
        let mut args = vec!["stats", vcf.to_str().expect("invariant: utf8"), report[0]];
        args.extend(&report[1..]);
        args.extend(["-o", prefix.to_str().expect("invariant: utf8")]);
        let out = run(&args);
        assert_eq!(
            out.status.code(),
            Some(2),
            "{tag}: a non-Integer depth/GQ tag is a usage error: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains(expected), "{tag}: {stderr}");
    }
    // `--summary` shares the depth type check with `--depth`.
    let prefix = dir.join("float-summary");
    let out = run(&[
        "stats",
        float_dp.to_str().expect("invariant: utf8"),
        "--summary",
        "-o",
        prefix.to_str().expect("invariant: utf8"),
    ]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("DP is Float"), "{stderr}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v23 V23-21 (decision A): every prefix-writing command refuses `-o -`
/// and an empty prefix *before* opening the input. Previously `stats`,
/// `convert`, `distance` and `pca` created literal `-.afreq`/`-.012`/
/// `-.mibs`/`-.eigenvec` files and exited 0, and `popg` only noticed after
/// the whole computation.
#[test]
fn prefix_commands_refuse_dash_and_empty_output() {
    let dir = temp_dir("v23-prefix");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let inp = input.to_str().expect("invariant: utf8");
    let cwd = dir.to_str().expect("invariant: utf8");
    let bare = dir.join("bare");
    let bare = bare.to_str().expect("invariant: utf8");

    // (subcommand arguments, an argument that makes it write a prefixed file)
    let cases: &[(&[&str], &[&str])] = &[
        (&["stats"], &["--freq"]),
        (&["convert"], &["--012"]),
        (&["pca"], &["--pcs", "2"]),
        (&["distance"], &["--type", "ibs"]),
        (&["popg"], &["--pi"]),
    ];
    for (sub, extra) in cases {
        for prefix in ["-", "", "   "] {
            let mut args: Vec<&str> = sub.to_vec();
            args.push(inp);
            args.extend(*extra);
            args.extend(["-o", prefix, "--verbosity", "0"]);
            let out = std::process::Command::new(env!("CARGO_BIN_EXE_reed"))
                .args(&args)
                .current_dir(cwd)
                .stdin(Stdio::null())
                .stderr(Stdio::piped())
                .output()
                .expect("invariant: spawn reed");
            let stderr = String::from_utf8_lossy(&out.stderr);
            assert_eq!(
                out.status.code(),
                Some(2),
                "{sub:?} with -o {prefix:?} must be a usage error: {stderr}"
            );
            assert!(
                stderr.contains("-o") || stderr.contains("prefix"),
                "{sub:?} -o {prefix:?}: {stderr}"
            );
            // No literal `-`-prefixed or hidden output may appear.
            let dash_files: Vec<_> = std::fs::read_dir(cwd)
                .expect("invariant: dir")
                .filter_map(Result::ok)
                .map(|e| e.file_name().to_string_lossy().into_owned())
                .filter(|n| n.starts_with('-') || n.starts_with('.'))
                .collect();
            assert!(
                dash_files.is_empty(),
                "{sub:?} -o {prefix:?} created {dash_files:?}"
            );
        }
    }

    // The control: a real prefix still works for one of them.
    let out = run(&["stats", inp, "--freq", "-o", bare, "--verbosity", "0"]);
    assert!(out.status.success(), "{out:?}");
    assert!(
        dir.join("bare.afreq").exists(),
        "control run must write its file"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v23 V23-49 (decisions B): the two region families AND together, the
/// same family twice is a usage error, and `-T` is an include file (`^`
/// excludes). Values measured against bcftools 1.24-2 on the same fixture.
#[test]
fn region_flag_combinations_match_bcftools() {
    let dir = temp_dir("v23-regions");
    let vcf = dir.join("in.vcf.gz");
    let text = "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##contig=<ID=chr2,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr2\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut compressed = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut compressed)
            .expect("invariant: block");
    }
    compressed.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    std::fs::write(&vcf, compressed).expect("invariant: write");
    let out = run(&["index", vcf.to_str().expect("invariant: utf8")]);
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let inp = vcf.to_str().expect("invariant: utf8");
    let c2only = dir.join("c2only.tab");
    std::fs::write(&c2only, "chr2\t1\t1000\n").expect("invariant: write");
    let c2 = c2only.to_str().expect("invariant: utf8");

    // (args, expected CHROM/POS pairs). The first block is byte-identical
    // to bcftools 1.24-2; the repeated-flag cases follow reed's documented
    // union rule (see the `-r`+`-R` case below).
    type RegionCase<'a> = (&'a [&'a str], &'a [(&'a str, &'a str)]);
    let cases: &[RegionCase<'_>] = &[
        (&["-r", "chr2"], &[("chr2", "10")]),
        // Cross-family intersection: chr2 (from -r) AND chr2:10 (from -t).
        (&["-r", "chr1,chr2", "-t", "chr2:10"], &[("chr2", "10")]),
        // Disjoint families select nothing (bcftools prints no records here).
        (&["-r", "chr2", "-t", "chr1:10"], &[]),
        (&["-t", "^chr1"], &[("chr2", "10")]),
        // `-T` is an include file, not an exclude list.
        (&["-T", c2], &[("chr2", "10")]),
        // `-r` + `-T`: both families again.
        (&["-r", "chr1,chr2", "-T", c2], &[("chr2", "10")]),
        // Repeated `-t`: include and exclude lists union, order-independent.
        // "keep chr1, minus chr1:10" leaves nothing on this fixture.
        (&["-t", "chr1", "-t", "^chr1:10"], &[]),
        (&["-t", "^chr1:10", "-t", "chr1"], &[]),
        // Two includes union: both contigs.
        (
            &["-t", "chr1:10", "-t", "chr2:10"],
            &[("chr1", "10"), ("chr2", "10")],
        ),
    ];
    for (args, want) in cases {
        let mut argv: Vec<&str> = vec!["filter", inp];
        argv.extend(*args);
        argv.extend(["-Ov", "--verbosity", "0"]);
        let out = run(&argv);
        assert!(out.status.success(), "{args:?}: {out:?}");
        let text = String::from_utf8_lossy(&out.stdout);
        let got: Vec<(&str, &str)> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let mut f = l.split('\t');
                (f.next().expect("chrom"), f.next().expect("pos"))
            })
            .collect();
        assert_eq!(got, *want, "{args:?}");
    }

    // `-r` and `-R` together union (both are the indexed family).
    let r12 = dir.join("r12.tab");
    std::fs::write(&r12, "chr2\t1\t1000\n").expect("invariant: write");
    let mut argv: Vec<&str> = vec![
        "filter",
        inp,
        "-r",
        "chr1",
        "-R",
        r12.to_str().expect("utf8"),
    ];
    argv.extend(["-Ov", "--verbosity", "0"]);
    let out = run(&argv);
    assert!(out.status.success(), "{out:?}");
    let text = String::from_utf8_lossy(&out.stdout);
    let got: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').next().expect("chrom"))
        .collect();
    assert_eq!(
        got,
        vec!["chr1", "chr2"],
        "repeated -r/-R union, not last-wins"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v23 V23-32: `--format` is a matrix-only option. The per-pair table
/// kinds used to ignore it silently, so a script could configure a layout
/// that was never written.
#[test]
fn format_is_rejected_for_table_distance_kinds() {
    let dir = temp_dir("v23-format-table");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let inp = input.to_str().expect("invariant: utf8");
    let out_dir = dir.join("out");
    let outp = out_dir.to_str().expect("invariant: utf8");

    let out = run(&[
        "distance",
        inp,
        "--type",
        "pihat",
        "--format",
        "square",
        "-o",
        outp,
        "--verbosity",
        "0",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(2),
        "must be a usage error: {stderr}"
    );
    assert!(stderr.contains("--format"), "{stderr}");

    // The matrix control still accepts it.
    let out = run(&[
        "distance",
        inp,
        "--type",
        "ibs",
        "--format",
        "square",
        "-o",
        outp,
        "--verbosity",
        "0",
    ]);
    assert!(
        out.status.success(),
        "matrix --format must still work: {out:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v23 V23-32: a missing `--contig-lengths` file names the path, not
/// just the OS message.
#[test]
fn contig_lengths_error_names_the_file() {
    let dir = temp_dir("v23-contig-lengths");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let inp = input.to_str().expect("invariant: utf8");
    let out_dir = dir.join("out");
    let outp = out_dir.to_str().expect("invariant: utf8");
    let missing = dir.join("no-such-lengths.tsv");
    let missing = missing.to_str().expect("invariant: utf8");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let pops = pops.to_str().expect("invariant: utf8");

    // `stats` redirects --contig-lengths to `popg`, which is where the
    // windowed modes live.
    let out = run(&[
        "popg",
        inp,
        "--pop",
        pops,
        "--pi",
        "--contig-lengths",
        missing,
        "-o",
        outp,
        "--verbosity",
        "0",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success(), "a missing file must fail: {out:?}");
    assert!(stderr.contains(missing), "the path must be named: {stderr}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v24 V24-04: region parsing uses the header's contig vocabulary. A
/// contig whose name contains `:` is a whole name when declared, and still
/// accepts an interval; the plan predicate and the index seek are derived
/// from the same parse.
#[test]
fn hla_contig_regions_use_the_header_vocabulary() {
    let dir = temp_dir("v24-hla-regions");
    let vcf = dir.join("in.vcf.gz");
    let text = "##fileformat=VCFv4.2\n\
         ##contig=<ID=HLA-A*01:01,length=1000>\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         HLA-A*01:01\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         HLA-A*01:01\t20\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         HLA-A*01:01\t30\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut compressed = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut compressed)
            .expect("invariant: block");
    }
    compressed.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    std::fs::write(&vcf, compressed).expect("invariant: write");
    let inp = vcf.to_str().expect("invariant: utf8");
    let out = run(&["index", inp]);
    assert!(
        out.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let positions = |args: &[&str]| -> Vec<(String, String)> {
        let mut argv: Vec<&str> = vec!["filter", inp];
        argv.extend(args);
        argv.extend(["-Ov", "--verbosity", "0"]);
        let out = run(&argv);
        assert!(out.status.success(), "{args:?}: {out:?}");
        String::from_utf8_lossy(&out.stdout)
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let mut f = l.split('\t');
                (
                    f.next().expect("chrom").to_string(),
                    f.next().expect("pos").to_string(),
                )
            })
            .collect()
    };

    assert_eq!(
        positions(&["-r", "HLA-A*01:01"]),
        vec![
            ("HLA-A*01:01".to_string(), "10".to_string()),
            ("HLA-A*01:01".to_string(), "20".to_string()),
            ("HLA-A*01:01".to_string(), "30".to_string()),
        ],
        "a declared colon contig is one whole region"
    );
    assert_eq!(
        positions(&["-r", "HLA-A*01:01:10-20"]),
        vec![
            ("HLA-A*01:01".to_string(), "10".to_string()),
            ("HLA-A*01:01".to_string(), "20".to_string()),
        ],
        "the last ':' with an interval splits the range"
    );
    assert_eq!(
        positions(&["-t", "HLA-A*01:01:10-20"]),
        vec![
            ("HLA-A*01:01".to_string(), "10".to_string()),
            ("HLA-A*01:01".to_string(), "20".to_string()),
        ],
        "streamed targets use the same vocabulary"
    );
    assert_eq!(
        positions(&["-r", "chr1"]),
        vec![("chr1".to_string(), "10".to_string())]
    );

    // A typo is loud, not a silent empty result.
    let out = run(&[
        "filter",
        inp,
        "-r",
        "HLA-A*01:99",
        "-Ov",
        "--verbosity",
        "0",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success(), "undeclared colon name must fail");
    assert!(stderr.contains("not declared in the header"), "{stderr}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v25 V25-06: parameter domains are usage errors (exit 2), not a panic
/// (`--threads 0` on distance used to mean all cores while every budget path
/// meant sequential) and not a runtime failure (`--memory 0` used to fail at
/// the first chunk). Each case also names the expected message fragment.
#[test]
fn parameter_domains_are_usage_errors() {
    let cases: &[(&[&str], &str)] = &[
        (&["filter", "--threads", "0", "-"], "must be >= 1"),
        (&["stats", "--threads", "0", "-"], "must be >= 1"),
        (&["popg", "--threads", "0", "-"], "must be >= 1"),
        (&["pca", "--threads", "0", "-"], "must be >= 1"),
        (&["distance", "--threads", "0", "-"], "must be >= 1"),
        (&["convert", "--threads", "0", "-"], "must be >= 1"),
        (&["rename", "--threads", "0", "-"], "must be >= 1"),
        (&["index", "--threads", "0", "in.vcf.gz"], "must be >= 1"),
        (&["filter", "--memory", "0", "-"], "must be >= 1 MiB"),
        (&["popg", "--memory", "0", "-"], "must be >= 1 MiB"),
        (&["filter", "-O", "", "in.vcf"], "invalid -O ''"),
        (&["filter", "-O", "z0", "in.vcf"], "want 1-9"),
        (&["filter", "-O", "z10", "in.vcf"], "want 1-9"),
        (&["rename", "-O", "", "in.vcf"], "invalid -O ''"),
        (&["rename", "-O", "b", "in.vcf"], "accepts v or z"),
        (&["popg", "--ld-decay", "0", "in.vcf"], "must be >= 1"),
        (&["filter", "-W", "-o", "-", "in.vcf"], "-W needs"),
    ];
    for (args, want) in cases {
        let out = run(args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert_eq!(
            out.status.code(),
            Some(2),
            "{args:?} must be a usage error; stderr: {stderr}"
        );
        assert!(
            stderr.contains(want),
            "{args:?}: expected {want:?} in stderr, got: {stderr}"
        );
    }
}

/// doc_v26 V26-03: `--threads` has a documented upper bound; anything above
/// it is a usage error on every command, not an attempt to spawn an
/// arbitrary number of threads or allocate `threads + 4` channel slots.
#[test]
fn threads_above_the_supported_maximum_are_usage_errors() {
    // The number is literal here on purpose: the front end depends on the
    // facade only (doc_v5 08 §6.4), and the engine's unit test pins the
    // constant to this value (`MAX_THREADS_documented_in_the_cli`).
    let over = "1025";
    let at_max = "1024";
    for command in [
        "filter", "stats", "popg", "pca", "distance", "convert", "rename",
    ] {
        let out = run(&[command, "--threads", over, "-"]);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert_eq!(
            out.status.code(),
            Some(2),
            "{command} --threads {over} must be a usage error; stderr: {stderr}"
        );
        assert!(
            stderr.contains("exceeds the supported maximum"),
            "{command}: {stderr}"
        );
        // The cap itself is accepted (parsing succeeds; the run then fails on
        // the empty stdin, which is a different, non-usage-code path).
        let out = run(&[command, "--threads", at_max, "-"]);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            !stderr.contains("exceeds the supported maximum"),
            "{command} --threads {at_max} must be accepted; stderr: {stderr}"
        );
    }
}

/// A small plain VCF with the given sample names, one record per POS.
fn plain_vcf(samples: &[&str]) -> String {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=Chr01,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n",
    );
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in samples {
        text.push('\t');
        text.push_str(s);
    }
    text.push('\n');
    for i in 0..8_u32 {
        let gts: Vec<&str> = if samples.len() == 1 {
            vec!["0/1"]
        } else {
            ["0/1", "0/0", "1/1", "0/0"]
                .iter()
                .copied()
                .cycle()
                .skip(i as usize)
                .take(samples.len())
                .collect()
        };
        let _ = writeln!(
            text,
            "Chr01\t{}\t.\tA\tG\t50\tPASS\t.\tGT\t{}",
            i + 1,
            gts.join("\t")
        );
    }
    text
}

/// doc_v25 V25-08 (R14): the piHat table writes the sample name in both the
/// FID and IID columns; a name containing `_` is no longer split at a
/// guessed boundary.
#[test]
fn pihat_table_writes_the_sample_name_in_both_columns() {
    let dir = temp_dir("pihat");
    let input = dir.join("in.vcf");
    std::fs::write(&input, plain_vcf(&["fam_a", "fam_b", "fam_c"])).expect("invariant: fixture");
    let prefix = dir.join("out");
    let out = run(&[
        "distance",
        "--type",
        "pihat",
        "-o",
        prefix.to_str().expect("utf8"),
        input.to_str().expect("utf8"),
    ]);
    assert!(
        out.status.success(),
        "distance failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = std::fs::read_to_string(dir.join("out.genome")).expect("invariant: table");
    let mut checked = 0;
    for line in text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols[0], cols[1], "FID1 != IID1: {line}");
        assert_eq!(cols[2], cols[3], "FID2 != IID2: {line}");
        assert!(cols[0].contains('_'), "the name was split: {line}");
        checked += 1;
    }
    assert_eq!(checked, 3, "three pairs");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v25 V25-08: a request-shaped tree failure is usage (exit 2); the
/// data-shaped ones are runtime (exit 1). Too few taxa is the CLI-fixable
/// side and must not be confused with a missing input file.
#[test]
fn tree_with_too_few_taxa_is_a_usage_error() {
    let dir = temp_dir("taxa");
    let input = dir.join("one.vcf");
    std::fs::write(&input, plain_vcf(&["only"])).expect("invariant: fixture");
    let out = run(&[
        "distance",
        "--tree",
        "nj",
        "-o",
        dir.join("out").to_str().expect("utf8"),
        input.to_str().expect("utf8"),
    ]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("need at least"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v25 V25-08: an unwritable output is a runtime failure (exit 1) on
/// rename, filter and both convert routes — it used to depend on the
/// subcommand (the hapmap→VCF route still exited 2).
#[test]
fn output_create_failure_is_runtime_on_rename_and_filter() {
    let dir = temp_dir("badout");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let bad = dir.join("missing").join("out.vcf");
    for command in ["rename", "filter"] {
        let mut args = vec![command, input.to_str().expect("utf8")];
        if command == "rename" {
            args.push("--clear-ids");
        }
        args.extend(["-o", bad.to_str().expect("utf8")]);
        let out = run(&args);
        assert_eq!(
            out.status.code(),
            Some(1),
            "{command}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
    }
    // The hapmap input route creates its output in a different writer; it
    // used to map `File::create` to a usage error (exit 2).
    let hmp = dir.join("in.hmp.txt");
    std::fs::write(
        &hmp,
        "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\t\
         assayLSID\tpanelLSID\tQCcode\ts0\n\
         rs1\tA/G\t1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\n",
    )
    .expect("invariant: hapmap fixture");
    let out = run(&[
        "convert",
        hmp.to_str().expect("utf8"),
        "--vcf",
        "-o",
        bad.to_str().expect("utf8"),
    ]);
    assert_eq!(
        out.status.code(),
        Some(1),
        "hapmap convert: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v25 V25-08 (R11): report-mode usage validation runs before the input
/// opens, so a missing file cannot turn a bad command line into exit 1.
#[test]
fn stats_usage_errors_precede_a_missing_input() {
    let dir = temp_dir("r11");
    let missing = dir.join("missing.vcf");
    let prefix = dir.join("out");

    // Missing input + invalid sample filter: usage wins.
    let out = run(&[
        "stats",
        missing.to_str().expect("utf8"),
        "--freq",
        "-s",
        "s0",
        "-o",
        prefix.to_str().expect("utf8"),
    ]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Missing input + valid flags: the missing file is runtime.
    let out = run(&[
        "stats",
        missing.to_str().expect("utf8"),
        "--freq",
        "-o",
        prefix.to_str().expect("utf8"),
    ]);
    assert_eq!(
        out.status.code(),
        Some(1),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Existing input + the same invalid sample filter: still usage.
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let out = run(&[
        "stats",
        input.to_str().expect("utf8"),
        "--freq",
        "-s",
        "s0",
        "-o",
        prefix.to_str().expect("utf8"),
    ]);
    assert_eq!(
        out.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v25 V25-08b: every command's `--timing` carries the same engine key
/// set, so one performance-regression parser works across subcommands.
#[test]
fn timing_engine_keys_are_the_same_on_every_command() {
    let dir = temp_dir("timing");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: fixture");
    let input = input.to_str().expect("utf8");
    let prefix = dir.join("out");
    let prefix = prefix.to_str().expect("utf8");

    let pops = dir.join("pops.tsv");
    std::fs::write(&pops, "s0\tpop1\ns1\tpop1\n").expect("invariant: pops");
    let pops = pops.to_str().expect("utf8");
    let filter_out = format!("{prefix}.vcf");
    let index_out = format!("{prefix}.csi");
    let cases: Vec<(&str, Vec<&str>)> = vec![
        (
            "filter",
            vec!["filter", input, "-o", &filter_out, "--verbosity", "2"],
        ),
        (
            "popg",
            vec![
                "popg",
                input,
                "--pi",
                "--pop",
                pops,
                "--variant-only",
                "--window",
                "1000",
                "-o",
                prefix,
                "--verbosity",
                "2",
            ],
        ),
        (
            "stats",
            vec!["stats", input, "--freq", "-o", prefix, "--verbosity", "2"],
        ),
        (
            "convert",
            vec![
                "convert",
                input,
                "--phylip",
                "-o",
                prefix,
                "--verbosity",
                "2",
            ],
        ),
        (
            "index",
            vec![
                "index",
                input,
                "--csi",
                "-o",
                &index_out,
                "--verbosity",
                "2",
            ],
        ),
    ];
    const ENGINE_KEYS: &[&str] = &[
        "chunks_read",
        "worker_busy_ms",
        "reader_ms",
        "reduce_ms",
        "chain_wait_ms",
        "reader_wait_ms",
        "reduce_wait_ms",
        "inflate_ms",
        "job_ms",
        "deflate_ms",
        "in_flight_peak_bytes",
    ];
    for (command, args) in &cases {
        let out = run(args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(out.status.success(), "{command}: {stderr}");
        let keys: Vec<&str> = stderr
            .lines()
            .filter_map(|line| line.split_once('\t').map(|(key, _)| key))
            .collect();
        for key in ENGINE_KEYS {
            assert!(
                keys.contains(key),
                "{command}: missing timing key {key}; keys: {keys:?}"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-01: the bare `popg` route writes its wide TSV to `-o` as a
/// *file path*, and used to do so without the `-o` == input guard every
/// other file-writing command has. `-o` naming the input truncated it
/// mid-read (md5 changed) and still exited 0 — unrecoverable data loss that
/// reported success (doc_v20 R5).
#[test]
fn popg_bare_route_refuses_output_naming_the_input() {
    let dir = temp_dir("v27-popg-overwrite");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n",
    )
    .expect("invariant: write");
    let before = std::fs::read(&vcf).expect("invariant: read");
    let pops = dir.join("p.tsv");
    std::fs::write(&pops, "s0\tA\n").expect("invariant: write");
    let v = sv(&vcf);

    let result = run(&["popg", &v, "--pop", &sv(&pops), "--variant-only", "-o", &v]);
    assert_eq!(result.status.code(), Some(2), "{result:?}");
    let err = String::from_utf8_lossy(&result.stderr);
    assert!(err.contains("would overwrite the input"), "{err}");
    assert_eq!(
        std::fs::read(&vcf).expect("invariant: read"),
        before,
        "the input must be byte-for-byte unchanged"
    );

    // The windowed route is a PREFIX: `-o <input>` only adds suffixes, so it
    // stays legal (and must not touch the input either).
    let ok = run(&[
        "popg",
        &v,
        "--pop",
        &sv(&pops),
        "--pi",
        "--variant-only",
        "--window",
        "1000",
        "-o",
        &v,
    ]);
    assert!(ok.status.success(), "{ok:?}");
    assert_eq!(
        std::fs::read(&vcf).expect("invariant: read"),
        before,
        "the windowed route writes PREFIX_pi.txt, never the input"
    );
    assert!(dir.join("in.vcf_pi.txt").exists(), "the pi file is written");

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-14/V27-15: the three site families -- `-r/-R`, `-t/-T` and
/// `--sites-file` -- **intersect**. They used to share one list, so every
/// pair was a union: `-r chr1:201-400 -t chr1:1-60` counted two records
/// where `reed filter` (and bcftools) count none, and `--sites-file` naming
/// one position plus `-t` naming three more counted all four.
#[test]
fn popg_and_stats_intersect_the_site_families() {
    let dir = temp_dir("v27-region-and");
    let vcf = dir.join("in.vcf.gz");
    // Records 100..500 on chr1, all biallelic so the counts are the record
    // counts. Pop A holds s0/s1 and pop B s2.
    let text = "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
         chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
         chr1\t300\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
         chr1\t400\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
         chr1\t500\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n";
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut compressed = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut compressed)
            .expect("invariant: block");
    }
    compressed.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    std::fs::write(&vcf, compressed).expect("invariant: write");
    let indexed = run(&["index", vcf.to_str().expect("utf8"), "-c"]);
    assert!(
        indexed.status.success(),
        "index failed: {}",
        String::from_utf8_lossy(&indexed.stderr)
    );
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\n").expect("invariant: pops");
    let sites = dir.join("sites.txt");
    std::fs::write(&sites, "chr1\t200\n").expect("invariant: sites");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, p, s) = (sv(&vcf), sv(&pops), sv(&sites));

    // `pi.txt` data rows: the record count is the summary line's `N records`
    // (a zero-record run writes only the header, which is what the
    // differential gate below keys on).
    let popg_records = |name: &str, extra: &[&str]| -> (u64, String) {
        let out = dir.join(name);
        let o = sv(&out);
        let mut argv: Vec<&str> = vec!["popg", &v, "--pop", &p, "--pi", "--variant-only", "-o", &o];
        argv.extend_from_slice(extra);
        let result = run(&argv);
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let stderr = String::from_utf8_lossy(&result.stderr).into_owned();
        let n = stderr
            .split_once(" records,")
            .and_then(|(head, _)| head.rsplit(' ').next())
            .and_then(|n| n.parse::<u64>().ok())
            .unwrap_or_else(|| panic!("{name}: no record count in {stderr}"));
        let body = std::fs::read_to_string(format!("{o}_pi.txt")).expect("invariant: pi file");
        (n, body)
    };
    // The same record counts from `reed filter`, the command that already
    // had the intersection (V23-61) and the reference for these cases.
    let filter_records = |extra: &[&str]| -> u64 {
        let mut argv: Vec<&str> = vec!["filter", &v];
        argv.extend_from_slice(extra);
        argv.extend(["-Ov", "--verbosity", "0"]);
        let result = run(&argv);
        assert!(
            result.status.success(),
            "{extra:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        String::from_utf8_lossy(&result.stdout)
            .lines()
            .rev()
            .take_while(|l| !l.starts_with('#'))
            .filter(|l| !l.is_empty())
            .count() as u64
    };

    // V27-15: cross-family intersection. The two ranges are disjoint, so
    // both sides must select nothing — the union used to count 300 and 400.
    let (n, body) = popg_records("cross", &["-r", "chr1:201-400", "-t", "chr1:1-60"]);
    assert_eq!(
        n,
        filter_records(&["-r", "chr1:201-400", "-t", "chr1:1-60"])
    );
    assert_eq!(n, 0, "disjoint -r and -t must select nothing");
    assert_eq!(body.lines().count(), 1, "only the header row: {body:?}");

    // Each family on its own keeps its own meaning.
    let (n, _) = popg_records("r_only", &["-r", "chr1:201-400"]);
    assert_eq!(n, 2, "-r alone keeps 300/400");
    let (n, _) = popg_records("t_only", &["-t", "chr1:1-60"]);
    assert_eq!(n, 0, "-t alone keeps nothing in 1-60");
    // Overlapping families keep the overlap.
    let (n, _) = popg_records("overlap", &["-r", "chr1:201-400", "-t", "chr1:350-500"]);
    assert_eq!(
        n,
        filter_records(&["-r", "chr1:201-400", "-t", "chr1:350-500"])
    );
    assert_eq!(n, 1, "only 400 is in both ranges");

    // V27-14: `--sites-file` is its own family. Alone it keeps 200; with a
    // target list that excludes 200 the intersection is empty, not the
    // union {200,300,400,500}.
    let (n, _) = popg_records("sites_only", &["--sites-file", &s]);
    assert_eq!(n, 1, "--sites-file alone keeps its one position");
    let (n, body) = popg_records("sites_and_t", &["--sites-file", &s, "-t", "chr1:300-500"]);
    assert_eq!(n, 0, "--sites-file and -t must intersect");
    assert_eq!(body.lines().count(), 1, "only the header row: {body:?}");
    let (n, _) = popg_records(
        "sites_and_t_hit",
        &["--sites-file", &s, "-t", "chr1:100-300"],
    );
    assert_eq!(n, 1, "200 is inside both: the intersection keeps it");

    // The bare route (no windowed kind) counts records the same way and
    // must honour the same intersection: its wide TSV keeps one row per
    // default-kind window, so a non-empty `pi.pair.0` cell is the marker.
    let bare_nonempty_rows = |extra: &[&str]| -> u64 {
        let mut argv: Vec<&str> = vec!["popg", &v, "--pop", &p, "--variant-only"];
        argv.extend_from_slice(extra);
        let result = run(&argv);
        assert!(
            result.status.success(),
            "{extra:?}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        String::from_utf8_lossy(&result.stdout)
            .lines()
            .skip(1)
            .filter(|l| l.split('\t').nth(3).is_some_and(|cell| cell != "."))
            .count() as u64
    };
    assert_eq!(
        bare_nonempty_rows(&["-r", "chr1:201-400", "-t", "chr1:1-60"]),
        0,
        "bare route must intersect -r and -t"
    );
    assert_eq!(
        bare_nonempty_rows(&["--sites-file", &s, "-t", "chr1:300-500"]),
        0,
        "bare route must intersect --sites-file and -t"
    );
    // Both surviving records (300 and 400) share the default 100k window,
    // so `-r` alone leaves exactly one non-empty row.
    assert_eq!(
        bare_nonempty_rows(&["-r", "chr1:201-400"]),
        1,
        "bare route keeps -r alone"
    );

    // `stats --freq` shares the prelude, so it shares the fix: `-r` and `-t`
    // intersect there too, and `--sites-file` no longer widens `-t`.
    let stats_rows = |name: &str, extra: &[&str]| -> u64 {
        let out = dir.join(name);
        let o = sv(&out);
        let mut argv: Vec<&str> = vec!["stats", &v, "--freq", "-o", &o];
        argv.extend_from_slice(extra);
        let result = run(&argv);
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        // The afreq header is one line; every data row is a reported site.
        std::fs::read_to_string(format!("{o}.afreq"))
            .expect("invariant: afreq")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .count() as u64
    };
    assert_eq!(
        stats_rows("s_cross", &["-r", "chr1:201-400", "-t", "chr1:1-60"]),
        0,
        "stats --freq must intersect -r and -t"
    );
    assert_eq!(
        stats_rows("s_r", &["-r", "chr1:201-400"]),
        2,
        "stats --freq keeps -r alone"
    );
    assert_eq!(
        stats_rows("s_sites", &["--sites-file", &s, "-t", "chr1:300-500"]),
        0,
        "stats --freq must intersect --sites-file and -t"
    );
    assert_eq!(
        stats_rows("s_sites_hit", &["--sites-file", &s, "-t", "chr1:100-300"]),
        1,
        "200 is inside both"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-07: SFS, f-statistics and LD are genome-wide computations;
/// with no window-level statistic they must not require contig lengths.
#[test]
fn genome_wide_popg_statistics_do_not_require_contig_lengths() {
    let dir = temp_dir("v27-nolen");
    let vcf = dir.join("in.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n\
        chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n\
        chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n";
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, p) = (sv(&vcf), sv(&pops));
    for (name, extra) in [
        ("sfs", vec!["--sfs", "--proj", "4"]),
        ("f2", vec!["--f2", "A,B"]),
        ("ld", vec!["--ld-decay", "1000"]),
        // The SFS-only route forces genome level even if the operator writes
        // `--levels window`; it still has no trailing Bp window to clip.
        (
            "sfs_window",
            vec!["--sfs", "--proj", "4", "--levels", "window"],
        ),
    ] {
        let out = dir.join(name);
        let o = sv(&out);
        let mut argv = vec!["popg", &v, "--pop", &p, "-o", &o];
        argv.extend(extra);
        let result = run(&argv);
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }
    assert!(dir.join("sfs_sfs/dadi/A-4.sfs").exists());
    assert!(dir.join("sfs_window_sfs/dadi/A-4.sfs").exists());
    assert!(dir.join("f2.fstats").exists());
    assert!(dir.join("ld.ld.decay").exists());
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-05/V27-20: output membership follows the request, not whether
/// a nonempty result happened to arrive. A zero-pair run replaces the old
/// `.ld` with a header, writes 0/NA decay bins, and lists `.ld` in the
/// success summary. The same holds for an empty SFS and f-statistic report.
#[test]
fn empty_popg_outputs_replace_the_previous_run_and_are_listed() {
    let dir = temp_dir("v27-empty-session");
    let vcf = dir.join("in.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=1000>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n\
        chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n\
        chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0\n";
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, p, o) = (sv(&vcf), sv(&pops), sv(&out));
    let run_all = |extra: &[&str]| {
        let mut argv = vec![
            "popg",
            &v,
            "--pop",
            &p,
            "--sfs",
            "--proj",
            "4",
            "--ld-decay",
            "1000",
            "--ld-pairs",
            "--f2",
            "A,B",
            "-o",
            &o,
        ];
        argv.extend_from_slice(extra);
        let result = run(&argv);
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        String::from_utf8_lossy(&result.stderr).into_owned()
    };
    run_all(&[]);
    assert_ne!(
        std::fs::read_to_string(dir.join("out.ld"))
            .expect("nonempty pairs")
            .lines()
            .count(),
        0
    );
    let stderr = run_all(&["-t", "chr1:999-1000"]);
    assert!(
        stderr.contains("out.ld"),
        "the streamed pair path is listed: {stderr}"
    );
    assert_eq!(
        std::fs::read_to_string(dir.join("out.ld")).expect("invariant: empty pairs"),
        "#CHROM_A\tPOS_A\tID_A\tCHROM_B\tPOS_B\tID_B\tUNPHASED_R2\n"
    );
    let decay = std::fs::read_to_string(dir.join("out.ld.decay")).expect("invariant: decay");
    assert!(decay.contains("all\t0\t1000\t0\tNA"), "{decay}");
    let fstats = std::fs::read_to_string(dir.join("out.fstats")).expect("invariant: fstats");
    assert!(fstats.contains("f2\tA,B\t0\tNA\tNA\t0\t0"), "{fstats}");
    let dadi = std::fs::read_to_string(dir.join("out_sfs/dadi/A-4.sfs")).expect("invariant: dadi");
    assert_eq!(dadi, "5 unfolded\n0 0 0 0 0\n1 0 0 0 1\n");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-06: the real reader must accept sorted sites between sparse
/// windows even when a chunk seam falls there. The long, deterministic IDs
/// make the compressed input exceed one chunk (about 1.3 MiB), so this is an
/// end-to-end repro rather than a test of the sink in isolation.
#[test]
fn sparse_windows_accept_sorted_sites_across_a_real_chunk_seam() {
    let dir = temp_dir("v27-sparse-seam");
    let vcf = dir.join("in.vcf.gz");
    let alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=40000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
    );
    let mut state = 27_u64;
    for pos in 1..=40_000_u64 {
        let mut id = String::with_capacity(36);
        for _ in 0..36 {
            // A fixed LCG: no entropy source, identical bytes on every run.
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            id.push(alphabet[usize::try_from(state % 62).expect("invariant: small")] as char);
        }
        writeln!(text, "chr1\t{pos}\t{id}\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0")
            .expect("invariant: write to String");
    }
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut compressed = Vec::new();
    for block in text.as_bytes().chunks(65_000) {
        deflater
            .deflate_block(block, &mut compressed)
            .expect("invariant: block");
    }
    compressed.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    assert!(compressed.len() > 1 << 20, "fixture must span input chunks");
    std::fs::write(&vcf, compressed).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--variant-only",
        "--window",
        "100",
        "--step",
        "1M",
        "--skip-empty",
        "--threads",
        "1",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let body = std::fs::read_to_string(format!("{}_pi.txt", out.display())).expect("pi file");
    let row = body
        .lines()
        .find(|line| line.starts_with("A\t"))
        .expect("pop A row");
    let fields: Vec<&str> = row.split('\t').collect();
    assert_eq!(&fields[1..4], ["chr1", "1", "100"]);
    assert_eq!(
        fields[5], "100",
        "all first-window sites are counted: {body}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-04: the real CLI must keep all 4096 cross-block pairs when
/// the LD partition has a one-site tail. The `.ld` rows and the decay bins
/// must agree, so this also catches a writer that silently loses a tail.
#[test]
fn ld_pairs_keep_the_single_site_tail_block() {
    let dir = temp_dir("v27-ld-tail");
    let vcf = dir.join("in.vcf");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=5000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
    );
    for pos in 1..=4097 {
        writeln!(
            text,
            "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1"
        )
        .expect("invariant: write to String");
    }
    std::fs::write(&vcf, text).expect("invariant: fixture");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--ld-decay",
        "100000",
        "--ld-pairs",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let pairs = std::fs::read_to_string(dir.join("out.ld")).expect("invariant: pair table");
    let pair_count = pairs.lines().count() - 1;
    assert_eq!(pair_count, 4097 * 4096 / 2);
    let decay = std::fs::read_to_string(dir.join("out.ld.decay")).expect("invariant: decay");
    let decay_count: u64 = decay
        .lines()
        .skip(1)
        .map(|line| {
            line.split('\t')
                .nth(3)
                .expect("invariant: n_pairs")
                .parse::<u64>()
                .expect("invariant: integer")
        })
        .sum();
    assert_eq!(
        usize::try_from(decay_count).expect("invariant: count fits"),
        pair_count
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-39: popg is biallelic-only. `-M 3` is a usage error rather
/// than an accepted option whose records are silently discarded later;
/// multiallelic records are filtered before they can affect denominators.
#[test]
fn popg_rejects_multiallelic_mode_and_filters_multiallelic_records() {
    let dir = temp_dir("v27-biallelic");
    let vcf = dir.join("in.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=200>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
        chr1\t100\t.\tA\tT,G\t50\tPASS\t.\tGT\t0/1\t0/0\n\
        chr1\t200\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/0\n";
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, p) = (sv(&vcf), sv(&pops));
    let bad = run(&[
        "popg",
        &v,
        "--pop",
        &p,
        "--pi",
        "-M",
        "3",
        "-o",
        &sv(&dir.join("bad")),
    ]);
    assert_eq!(bad.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&bad.stderr).contains("does not support multiallelic"));

    let out = dir.join("out");
    let ok = run(&[
        "popg",
        &v,
        "--pop",
        &p,
        "--pi",
        "--variant-only",
        "-M",
        "2",
        "-o",
        &sv(&out),
    ]);
    assert!(
        ok.status.success(),
        "{}",
        String::from_utf8_lossy(&ok.stderr)
    );
    let body = std::fs::read_to_string(format!("{}_pi.txt", out.display())).expect("invariant: pi");
    assert!(
        body.contains("A\tchr1\t1\t200\t0.005\t1\t1\t1\t0"),
        "{body}"
    );

    let stderr = String::from_utf8_lossy(&ok.stderr);
    assert!(stderr.contains("1 multiallelic skipped"), "{stderr}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-03/V27-12/V27-13/V27-27: this small fixture is also a pixy
/// acceptance case. A's site 1 has only one called allele in B, so it is
/// excluded from Hudson's numerator AND denominator; π uses per-population
/// occupied sites and dxy their intersection. `no_snps` counts variants.
#[test]
fn popg_matches_pixy_on_low_coverage_and_occupied_site_counts() {
    let dir = temp_dir("v27-pixy-counts");
    let vcf = dir.join("in.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=100>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n\
        chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/.\t./.\n\
        chr1\t20\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1\n";
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--dxy",
        "--fst",
        "hudson",
        "--bypass-invariant-check",
        "--window",
        "100",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let pi = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    assert!(pi.contains("A\tchr1\t1\t100\t0.25\t2\t3\t12\t0"), "{pi}");
    assert!(
        pi.contains("B\tchr1\t1\t100\t0.6666666666666666\t2\t4\t6\t6"),
        "{pi}"
    );
    let dxy = std::fs::read_to_string(dir.join("out_dxy.txt")).expect("invariant: dxy");
    assert!(
        dxy.contains("A\tB\tchr1\t1\t100\t0.45\t2\t9\t20\t12"),
        "{dxy}"
    );
    let fst = std::fs::read_to_string(dir.join("out_fst.txt")).expect("invariant: fst");
    assert!(
        fst.contains("A\tB\tchr1\t1\t100\t0.33333333333333337\t2"),
        "{fst}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-12/V27-13: a record with no calls in A is absent from A's
/// |S|, and dxy counts it only when both populations have a call. Missing
/// audit counts are still preserved exactly (pixy's 6 and 16 here).
#[test]
fn popg_matches_pixy_on_pairwise_occupancy_and_missing_audits() {
    let dir = temp_dir("v27-pixy-missing");
    let vcf = dir.join("in.vcf");
    let text = "##fileformat=VCFv4.2\n\
        ##contig=<ID=chr1,length=100>\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n\
        chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t./.\t./.\t0/1\t0/1\n\
        chr1\t20\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\n";
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--dxy",
        "--bypass-invariant-check",
        "--window",
        "100",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let pi = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    assert!(
        pi.contains("A\tchr1\t1\t100\t0.6666666666666666\t1\t4\t6\t6"),
        "{pi}"
    );
    assert!(
        pi.contains("B\tchr1\t1\t100\t0.6666666666666666\t2\t8\t12\t0"),
        "{pi}"
    );
    let dxy = std::fs::read_to_string(dir.join("out_dxy.txt")).expect("invariant: dxy");
    assert!(
        dxy.contains("A\tB\tchr1\t1\t100\t0.5\t1\t8\t16\t16"),
        "{dxy}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-18: the CLI's window-output session is streaming and atomic.
/// A complete 100-site window must not be confused by the per-row watermark;
/// a later malformed input must not publish its temporary session files.
#[test]
fn streaming_window_output_counts_rows_and_cleans_failed_sessions() {
    let dir = temp_dir("v27-stream-session");
    let vcf = dir.join("in.vcf");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n",
    );
    for pos in 1..=100 {
        writeln!(text, "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0")
            .expect("invariant: write to String");
    }
    std::fs::write(&vcf, text).expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let args = [
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--variant-only",
        "--window",
        "100",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ];
    let ok = run(&args);
    assert!(
        ok.status.success(),
        "{}",
        String::from_utf8_lossy(&ok.stderr)
    );
    let pi = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    assert!(pi.contains("A\tchr1\t1\t100\t1\t100\t100\t100\t0"), "{pi}");
    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .expect("invariant: read dir")
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.file_name().to_string_lossy().contains("reed-tmp"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "temporary session files remain: {leftovers:?}"
    );

    // Reuse the prefix with an input whose malformed last record fails after
    // the session has already written rows. The old good file must survive.
    let before = std::fs::read(dir.join("out_pi.txt")).expect("invariant: old pi");
    std::fs::write(&vcf, "chr1\t1\t.\tA\tG\n").expect("invariant: bad write");
    let failed = run(&args);
    assert!(!failed.status.success());
    assert_eq!(
        std::fs::read(dir.join("out_pi.txt")).expect("invariant: old pi"),
        before
    );
    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .expect("invariant: read dir")
        .filter_map(|entry| entry.ok())
        .filter(|entry| entry.file_name().to_string_lossy().contains("reed-tmp"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "failed session left temporaries: {leftovers:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-21: invalid level names are usage errors on every route, and
/// BED without a window-level statistic is refused with the actual reason.
#[test]
fn popg_validates_levels_and_bed_window_requirements() {
    let dir = temp_dir("v27-levels");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .expect("invariant: write");
    let bed = dir.join("windows.bed");
    std::fs::write(&bed, "chr1\t0\t1000\n").expect("invariant: bed");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, b, p) = (sv(&vcf), sv(&bed), sv(&pops));
    for args in [
        vec!["--sfs", "--proj", "2", "--levels", "bogus"],
        vec!["--ld-decay", "1000", "--levels", "bogus"],
    ] {
        let out = sv(&dir.join("out"));
        let mut argv = vec!["popg", &v, "--pop", &p, "-o", &out];
        argv.extend(args);
        let result = run(&argv);
        assert_eq!(result.status.code(), Some(2), "{argv:?}");
        assert!(String::from_utf8_lossy(&result.stderr).contains("unknown level"));
    }
    let result = run(&[
        "popg",
        &v,
        "--pop",
        &p,
        "--sfs",
        "--proj",
        "2",
        "--bed-file",
        &b,
        "--levels",
        "window",
        "-o",
        &sv(&dir.join("bed")),
    ]);
    assert_eq!(result.status.code(), Some(2));
    let error = String::from_utf8_lossy(&result.stderr);
    assert!(
        error.contains("--bed-file needs a window-level statistic"),
        "{error}"
    );
    assert!(!error.contains("require the window level"), "{error}");
    let result = run(&[
        "popg",
        &v,
        "--pop",
        &p,
        "--f2",
        "A,B",
        "--bed-file",
        &b,
        "-o",
        &sv(&dir.join("f2")),
    ]);
    assert_eq!(result.status.code(), Some(2));
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-22: pair statistics need two populations, and an empty group
/// assignment is an error rather than a successful run with no output.
#[test]
fn popg_rejects_pair_statistics_with_fewer_than_two_groups() {
    let dir = temp_dir("v27-groups");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .expect("invariant: write");
    let one = dir.join("one.txt");
    std::fs::write(&one, "s0\tA\ns1\tA\n").expect("invariant: pops");
    let empty = dir.join("empty.txt");
    std::fs::write(&empty, "").expect("invariant: pops");
    for (pops, args) in [
        (&one, vec!["--fst", "hudson"]),
        (&one, vec!["--dxy"]),
        (&empty, vec!["--pi"]),
    ] {
        let out = dir.join("out");
        let mut argv = vec![
            "popg",
            vcf.to_str().expect("invariant: utf8"),
            "--pop",
            pops.to_str().expect("invariant: utf8"),
            "-o",
            out.to_str().expect("invariant: utf8"),
        ];
        argv.extend(args);
        let result = run(&argv);
        assert_eq!(result.status.code(), Some(2), "{argv:?}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-23: a contig absent from the BED has zero BED windows, not an
/// invalid `window_pos_1 = 1, window_pos_2 = 0` output row.
#[test]
fn popg_bed_uncovered_contigs_do_not_emit_reversed_coordinates() {
    let dir = temp_dir("v27-bed-uncovered");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##contig=<ID=chr2,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n\
         chr2\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .expect("invariant: write");
    let bed = dir.join("one.bed");
    std::fs::write(&bed, "chr1\t0\t400\n").expect("invariant: bed");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--variant-only",
        "--bed-file",
        bed.to_str().expect("invariant: utf8"),
        "--levels",
        "window,contig,genome",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let body = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    assert!(!body.contains("chr2\t1\t0"), "{body}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-24: explicit --contig-lengths values take precedence over the
/// header, rather than being overwritten by it in the driver.
#[test]
fn popg_contig_lengths_file_overrides_header_lengths() {
    let dir = temp_dir("v27-contig-lengths");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t350\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .expect("invariant: write");
    let lengths = dir.join("short.fai");
    std::fs::write(&lengths, "chr1\t400\n").expect("invariant: lengths");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--variant-only",
        "--contig-lengths",
        lengths.to_str().expect("invariant: utf8"),
        "--window",
        "100",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let body = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    assert!(body.contains("chr1\t301\t400"), "{body}");
    assert!(
        !body.contains("401-500") && !body.contains("401\t500"),
        "{body}"
    );
    assert!(!body.contains("501\t600"), "{body}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-26: the run log follows the requested methods and includes
/// SFS discard counts even for a bare SFS run (which has no window rows).
#[test]
fn popg_logs_methods_and_sfs_counts_for_windowed_and_bare_routes() {
    let dir = temp_dir("v27-log-methods");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n\
         chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\n",
    )
    .expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
    let sv = |p: &std::path::Path| p.to_str().expect("invariant: utf8").to_string();
    let (v, p) = (sv(&vcf), sv(&pops));
    let theta_out = sv(&dir.join("theta"));
    let result = run(&[
        "popg",
        &v,
        "--pop",
        &p,
        "--theta",
        "--tajima-d",
        "--fst",
        "hudson",
        "--variant-only",
        "-o",
        &theta_out,
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let log = std::fs::read_to_string(format!("{theta_out}.log")).expect("invariant: log");
    assert!(log.contains("theta"), "{log}");
    assert!(log.contains("tajima_d"), "{log}");
    assert!(log.contains("fst"), "{log}");
    assert!(!log.contains("  pi "), "{log}");
    assert!(!log.contains("  dxy "), "{log}");

    let sfs_out = sv(&dir.join("sfs"));
    let result = run(&[
        "popg", &v, "--pop", &p, "--sfs", "--proj", "4", "-o", &sfs_out,
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let log = std::fs::read_to_string(format!("{sfs_out}.log")).expect("invariant: sfs log");
    assert!(log.contains("sfs"), "{log}");
    assert!(log.contains("sfs sites absorbed: 4"), "{log}");
    assert!(log.contains("sfs skipped (m < n): 0"), "{log}");

    let ld_out = sv(&dir.join("ld"));
    let result = run(&["popg", &v, "--pop", &p, "--ld-decay", "1000", "-o", &ld_out]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let log = std::fs::read_to_string(format!("{ld_out}.log")).expect("invariant: ld log");
    assert!(log.contains("ld"), "{log}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-28: filtering a multiallelic record must not hide a backwards
/// position move. The input contract is checked before the record is skipped.
#[test]
fn multiallelic_records_do_not_bypass_position_order_validation() {
    let dir = temp_dir("v27-multi-order");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n\
         chr1\t100\t.\tA\tT,G\t50\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--variant-only",
        "-o",
        dir.join("out").to_str().expect("invariant: utf8"),
    ]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("strictly ascending"));
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-29: an interrupted process's old temporary is not deleted
/// blindly (PID reuse is possible), but it is surfaced on the next run.
#[test]
fn stale_ld_temporary_is_reported_on_the_next_run() {
    let dir = temp_dir("v27-stale-temp");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n\
         chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\n",
    )
    .expect("invariant: write");
    let out = dir.join("out");
    let stale = dir.join("out.ld.reed-tmp-999999999-0");
    std::fs::write(&stale, "interrupted output").expect("invariant: stale");
    std::fs::File::options()
        .write(true)
        .open(&stale)
        .expect("invariant: open stale")
        .set_modified(std::time::SystemTime::now() - std::time::Duration::from_secs(120))
        .expect("invariant: set old mtime");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--ld-decay",
        "1000",
        "--ld-pairs",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("stale temporary file"));
    assert!(stale.exists(), "another process's temporary is not deleted");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-38: the bare wide table must honor variant-only and keep
/// pair-comparable, physical-bp and raw-theta quantities distinct.
#[test]
fn popg_bare_variant_only_matches_windows_and_preserves_column_units() {
    let dir = temp_dir("v27-bare-variant-only");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=100>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n\
         chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t./.\t0/0\t0/0\n",
    )
    .expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tA\ns2\tB\ns3\tB\n").expect("invariant: pops");
    let input = vcf.to_str().expect("invariant: utf8");
    let pops = pops.to_str().expect("invariant: utf8");
    let bare = run(&["popg", input, "--pop", pops, "--variant-only"]);
    assert!(
        bare.status.success(),
        "{}",
        String::from_utf8_lossy(&bare.stderr)
    );
    let text = String::from_utf8(bare.stdout).expect("invariant: utf8");
    let header = text
        .lines()
        .next()
        .expect("invariant: header")
        .split('\t')
        .collect::<Vec<_>>();
    let cells = text
        .lines()
        .nth(1)
        .expect("invariant: data")
        .split('\t')
        .collect::<Vec<_>>();
    let cell = |key: &str| -> f64 {
        let index = header
            .iter()
            .position(|column| *column == key)
            .unwrap_or_else(|| panic!("missing {key}: {text}"));
        cells[index]
            .parse()
            .unwrap_or_else(|_| panic!("non-numeric {key}: {text}"))
    };
    // A has four cohort chromosomes but two called alleles: diffs=1, comps=1.
    // Pair completion adds 99*C(4,2)=594 pairs, while physical-bp π is 1/100.
    assert!((cell("pi.pair.0") - 1.0 / 595.0).abs() < 1e-12, "{text}");
    assert!((cell("pi.bp.0") - 0.01).abs() < 1e-12, "{text}");
    assert!((cell("dxy.pair.0.1") - 1.0 / 398.0).abs() < 1e-12, "{text}");
    assert!((cell("dxy.bp.0.1") - 0.005).abs() < 1e-12, "{text}");
    assert!((cell("theta_w.total.0") - 1.0).abs() < 1e-12, "{text}");
    assert!((cell("theta_w.bp.0") - 0.01).abs() < 1e-12, "{text}");
    let out = dir.join("out");
    let window = run(&[
        "popg",
        input,
        "--pop",
        pops,
        "--pi",
        "--theta",
        "--variant-only",
        "--window",
        "100",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        window.status.success(),
        "{}",
        String::from_utf8_lossy(&window.stderr)
    );
    let pi = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    let theta =
        std::fs::read_to_string(dir.join("out_watterson_theta.txt")).expect("invariant: theta");
    assert!(pi.contains("0.001680672268907563"), "{pi}");
    assert!(theta.contains("0.01"), "{theta}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v27 V27-10: an empty window is zero under variant-only (the absent
/// positions are monomorphic), while the all-sites/Require route keeps NA.
#[test]
fn variant_only_empty_windows_are_zero_not_na() {
    let dir = temp_dir("v27-empty-window");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=300>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n\
         chr1\t250\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .expect("invariant: write");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tA\ns1\tB\n").expect("invariant: pops");
    let out = dir.join("out");
    let result = run(&[
        "popg",
        vcf.to_str().expect("invariant: utf8"),
        "--pop",
        pops.to_str().expect("invariant: utf8"),
        "--pi",
        "--theta",
        "--variant-only",
        "--window",
        "100",
        "-o",
        out.to_str().expect("invariant: utf8"),
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let pi = std::fs::read_to_string(dir.join("out_pi.txt")).expect("invariant: pi");
    let theta =
        std::fs::read_to_string(dir.join("out_watterson_theta.txt")).expect("invariant: theta");
    for body in [&pi, &theta] {
        assert!(
            body.contains("101\t200") && body.contains("201\t300"),
            "{body}"
        );
        assert!(
            !body.contains("101\t200\tNA") && !body.contains("201\t300\tNA"),
            "{body}"
        );
        assert!(
            body.contains("101\t200\t0") && body.contains("201\t300\t0"),
            "{body}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

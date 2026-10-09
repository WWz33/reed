//! doc_v4 02-cli-semantics.md 回退版 §7.1–7.2 gates: the missing-rate
//! boundaries are exact and agree with the reference tools.
//!
//! Each fixture is a small plain VCF whose per-record missing rate is known
//! by construction (`N = 10` and `N = 7`, so `X·2N` is sometimes an integer
//! and sometimes not; with and without half calls). reed's surviving site set
//! is compared against the set computed independently here; when `vcftools`
//! and `plink2` are present (the gold tool paths), the same fixtures run
//! through them and the sets must agree — `plink2 --geno X` directly, and
//! `vcftools --max-missing 1−X` with the subtraction done in exact decimal.

use std::path::{Path, PathBuf};
use std::process::Command;

use reed::{run_filter, FilterOptions};
use reed_filter::spec::{FilterSpec, MissingRate, TailRules};

/// Writes `text` as a plain VCF under a per-test temp directory.
fn write_vcf(dir: &Path, name: &str, text: &str) -> PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).expect("invariant: fixture write");
    path
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-missing-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    dir
}

/// Runs `reed filter` (library driver) with the given after-mask rules and
/// returns the surviving POS values in file order.
fn reed_positions(input: &Path, out: &Path, tail: TailRules) -> Vec<u64> {
    let spec = FilterSpec { site_after_mask: tail, ..FilterSpec::default() };
    let opts = FilterOptions {
        threads: 2,
        output: Some(out.to_path_buf()),
        write_index: reed::WriteIndex::None,
        compress_level: 6,
        strategy: None,
        command_line: "e2e_missing_boundary".to_string(),
        timing: None,
        report_json: None,
        force_samples: false,
        regions: Vec::new(),
        tmp_dir: None,
        output_compressed: None,
        output_bcf: None,
        mem_bytes: None,
        no_version: false,
        no_command: false,
    };
    run_filter(Some(input), &spec, &opts).expect("invariant: filter run");
    let text = std::fs::read_to_string(out).expect("invariant: read output");
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).expect("POS").parse::<u64>().expect("number"))
        .collect()
}

fn positions(text: &str) -> Vec<u64> {
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(1).expect("POS").parse::<u64>().expect("number"))
        .collect()
}

/// A VCF body of `n` samples where record `pos` has `missing` fully-missing
/// samples (the rest are `0/0`).
fn missing_body(n: usize, records: &[(u64, usize)]) -> String {
    let header = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=Chr01,length=100000>\n\
         ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n",
        (0..n).map(|i| format!("s{i}")).collect::<Vec<_>>().join("\t")
    );
    let mut text = header;
    for &(pos, missing) in records {
        let mut gts = vec!["0/0".to_string(); n - missing];
        gts.extend(std::iter::repeat_n("./.".to_string(), missing));
        text.push_str(&format!(
            "Chr01\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t{}\n",
            gts.join("\t")
        ));
    }
    text
}

/// Resolves an external reference tool: `$<VAR>` when set, else the gold
/// default under `$HOME`.
fn tool(var: &str, default: &str) -> Option<PathBuf> {
    if let Ok(path) = std::env::var(var) {
        let path = PathBuf::from(path);
        return path.exists().then_some(path);
    }
    let home = std::env::var("HOME").ok()?;
    let path = PathBuf::from(home).join(default);
    path.exists().then_some(path)
}

/// `vcftools --vcf in --max-missing X --recode --stdout` → surviving POS.
fn vcftools_positions(vcftools: &Path, dir: &Path, input: &Path, x: &str) -> Option<Vec<u64>> {
    let out = Command::new(vcftools)
        .current_dir(dir)
        .args(["--vcf"])
        .arg(input)
        .args(["--max-missing", x, "--recode", "--stdout"])
        .output()
        .ok()?;
    assert!(
        out.status.success(),
        "vcftools failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    Some(positions(&String::from_utf8_lossy(&out.stdout)))
}

/// `plink2 --vcf in --geno G --export vcf --out PREFIX` → surviving POS.
fn plink2_positions(plink2: &Path, dir: &Path, input: &Path, geno: &str, tag: &str) -> Option<Vec<u64>> {
    // No `with_extension`: tags like `0.9` must not turn into file extensions.
    let prefix = dir.join(format!("plink2-{}", tag.replace('.', "_")));
    let out = Command::new(plink2)
        .current_dir(dir)
        .args(["--vcf"])
        .arg(input)
        .args(["--geno", geno, "--export", "vcf", "--out"])
        .arg(&prefix)
        .output()
        .ok()?;
    assert!(
        out.status.success(),
        "plink2 failed: {}",
        String::from_utf8_lossy(&out.stdout)
    );
    let text = std::fs::read_to_string(format!("{}.vcf", prefix.display())).ok()?;
    Some(positions(&text))
}

/// The 10-sample boundary fixture (no half calls): record `pos` has `pos−1`
/// missing samples except pos 7 (all missing → call rate 0) — call rates
/// 1.0, 0.9, 0.8, 0.7, 0.5, 0.1, 0.0.
fn fixture_10(dir: &Path) -> PathBuf {
    write_vcf(
        dir,
        "ten.vcf",
        &missing_body(10, &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 5), (6, 9), (7, 10)]),
    )
}

#[test]
fn site_call_rate_boundaries_match_reference_sets_and_tools() {
    let dir = temp_dir("site");
    let input = fixture_10(&dir);
    let out = dir.join("out.vcf");
    let vcftools = tool("VCFTOOLS", "biosoft/vcftools");
    let plink2 = tool("PLINK2", "biosoft/plink2");

    // X (--geno: missing rate ≤ X) → expected survivors. The fixture's
    // records have missing rates 0, 0.1, 0.2, 0.3, 0.5, 0.9, 1.0.
    let cases: [(&str, &[u64]); 7] = [
        ("0", &[1]),
        ("0.1", &[1, 2]),
        ("0.2", &[1, 2, 3]),
        ("0.3", &[1, 2, 3, 4]),
        ("0.5", &[1, 2, 3, 4, 5]),
        ("0.9", &[1, 2, 3, 4, 5, 6]),
        ("1", &[1, 2, 3, 4, 5, 6, 7]),
    ];
    for (x, expected) in cases {
        let rate = MissingRate::parse(x).expect("valid");
        let got = reed_positions(
            &input,
            &out,
            TailRules { max_missing: Some(rate), ..TailRules::default() },
        );
        assert_eq!(got, expected, "reed --geno {x}");
        if let Some(vcftools) = &vcftools {
            // vcftools states the keep direction; the conversion is exact
            // decimal (doc_v4 02 回退版 §5.2).
            let keep = rate.complement().render();
            assert_eq!(
                vcftools_positions(vcftools, &dir, &input, &keep).expect("vcftools runs"),
                expected,
                "vcftools --max-missing {keep}"
            );
        }
        // plink2 spells the same threshold the same way; the fixtures have no
        // half calls, so its genotype-missing numerator agrees.
        if let Some(plink2) = &plink2 {
            assert_eq!(
                plink2_positions(plink2, &dir, &input, x, x).expect("plink2 runs"),
                expected,
                "plink2 --geno {x}"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn seven_sample_boundaries_are_exact() {
    // N = 7 → 2N = 14: X·2N is an integer for X = 0.5 but not for 0.8/0.7.
    let dir = temp_dir("seven");
    let input = write_vcf(
        &dir,
        "seven.vcf",
        &missing_body(7, &[(1, 0), (2, 1), (3, 2), (4, 3), (5, 7)]),
    );
    let out = dir.join("out.vcf");
    let vcftools = tool("VCFTOOLS", "biosoft/vcftools");
    let plink2 = tool("PLINK2", "biosoft/plink2");
    // Missing alleles out of 14: 0, 2 (≈0.142857), 4 (≈0.285714),
    // 6 (≈0.428571), 14. The cases straddle the repeating decimals.
    let cases: [(&str, &[u64]); 6] = [
        ("0.1428572", &[1, 2]),
        ("0.1428571", &[1]),
        ("0.2857143", &[1, 2, 3]),
        ("0.28", &[1, 2]),
        ("0.4285715", &[1, 2, 3, 4]),
        ("0.42", &[1, 2, 3]),
    ];
    for (x, expected) in cases {
        let rate = MissingRate::parse(x).expect("valid");
        let got = reed_positions(
            &input,
            &out,
            TailRules { max_missing: Some(rate), ..TailRules::default() },
        );
        assert_eq!(got, expected, "reed --geno {x} (N = 7)");
        if let Some(vcftools) = &vcftools {
            let keep = rate.complement().render();
            assert_eq!(
                vcftools_positions(vcftools, &dir, &input, &keep).expect("vcftools runs"),
                expected,
                "vcftools --max-missing {keep} (N = 7)"
            );
        }
        if let Some(plink2) = &plink2 {
            assert_eq!(
                plink2_positions(plink2, &dir, &input, x, &format!("seven{x}")).expect("plink2 runs"),
                expected,
                "plink2 --geno {x} (N = 7)"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v4 02 §3.1: the numerator is called **alleles**; `0/.` contributes
/// one. vcftools agrees on this fixture (measured 2026-10-09); plink2 is
/// not comparable here (it needs `--vcf-half-call` and counts genotypes).
#[test]
fn half_calls_count_as_one_called_allele() {
    let dir = temp_dir("half");
    let text = "##fileformat=VCFv4.2\n##contig=<ID=Chr01,length=100000>\n\
        ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\
        s0\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\n\
        Chr01\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\n\
        Chr01\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/.\t0/.\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\n\
        Chr01\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/.\t./.\t./.\t./.\t./.\t./.\t./.\t./.\t./.\t./.\n\
        Chr01\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/.\t0/.\t0/.\t0/.\t0/.\t0/.\t0/.\t0/.\t0/.\t0/.\n";
    let input = write_vcf(&dir, "half.vcf", text);
    let out = dir.join("out.vcf");
    // pos 1: 20/20 = 1.0; pos 2: 16/20 = 0.8; pos 3: 1/20 = 0.05;
    // pos 4: 10/20 = 0.5.
    let cases: [(&str, &[u64]); 5] = [
        ("0", &[1]),
        ("0.2", &[1, 2]),
        ("0.5", &[1, 2, 4]),
        ("0.95", &[1, 2, 3, 4]),
        ("1", &[1, 2, 3, 4]),
    ];
    let vcftools = tool("VCFTOOLS", "biosoft/vcftools");
    for (x, expected) in cases {
        let rate = MissingRate::parse(x).expect("valid");
        let got = reed_positions(
            &input,
            &out,
            TailRules { max_missing: Some(rate), ..TailRules::default() },
        );
        assert_eq!(got, expected, "reed --geno {x} (half calls)");
        if let Some(vcftools) = &vcftools {
            let keep = rate.complement().render();
            assert_eq!(
                vcftools_positions(vcftools, &dir, &input, &keep).expect("vcftools runs"),
                expected,
                "vcftools --max-missing {keep} (half calls)"
            );
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v4 02 回退版 §7.2: the sample-level rule (`--min-gq 20 --mind 0.1`)
/// keeps exactly what plink2's `--vcf-min-gq 20 --mind 0.1` keeps. s1 is
/// masked in 2 of 10 records (0.2 → dropped), s2 is missing exactly once
/// (0.1 → kept at the boundary), s3 is missing 3 times (0.3 → dropped).
#[test]
fn sample_call_rate_matches_plink2_mind() {
    let dir = temp_dir("sample");
    // s1's two bad calls are masked (GQ < 20) rather than written missing, so
    // the mask stage does the work; s2/s3 are already missing in the file.
    let mut body = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=Chr01,length=100000>\n\
         ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
    );
    for i in 0..10_u32 {
        let s1 = if i < 2 { "0/0:5" } else { "0/0:99" };
        let s2 = if i == 0 { "./.:99" } else { "0/0:99" };
        let s3 = if i < 3 { "./.:99" } else { "0/0:99" };
        body.push_str(&format!(
            "Chr01\t{}\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/0:99\t{s1}\t{s2}\t{s3}\n",
            i + 1
        ));
    }
    let input = write_vcf(&dir, "samples.vcf", &body);
    let out = dir.join("out.vcf");
    let spec = FilterSpec {
        mind: Some(MissingRate::parse("0.1").expect("valid")),
        genotype: reed_filter::spec::GenotypeRules {
            mask_if: vec![reed_filter::spec::GenotypeMaskRule {
                field: "GQ".into(),
                when: "< 20".into(),
                on_missing: reed_filter::spec::OnMissing::Keep,
            }],
            ..reed_filter::spec::GenotypeRules::default()
        },
        ..FilterSpec::default()
    };
    let opts = FilterOptions {
        threads: 2,
        output: Some(out.clone()),
        write_index: reed::WriteIndex::None,
        compress_level: 6,
        strategy: None,
        command_line: "e2e_missing_boundary".to_string(),
        timing: None,
        report_json: None,
        force_samples: false,
        regions: Vec::new(),
        tmp_dir: None,
        output_compressed: None,
        output_bcf: None,
        mem_bytes: None,
        no_version: false,
        no_command: false,
    };
    run_filter(Some(&input), &spec, &opts).expect("invariant: run");
    let text = std::fs::read_to_string(&out).expect("invariant: read");
    let columns: Vec<&str> = text
        .lines()
        .find(|l| l.starts_with("#CHROM"))
        .expect("invariant: header")
        .split('\t')
        .collect();
    let kept: Vec<&str> = columns[9..].to_vec();
    assert_eq!(kept, vec!["s0", "s2"], "s2 sits exactly on the 0.1 boundary");

    if let Some(plink2) = tool("PLINK2", "biosoft/plink2") {
        let prefix = dir.join("plink2-mind");
        let out = Command::new(&plink2)
            .current_dir(&dir)
            .args(["--vcf"])
            .arg(&input)
            .args(["--vcf-min-gq", "20", "--mind", "0.1", "--export", "vcf", "--out"])
            .arg(&prefix)
            .output()
            .expect("invariant: spawn plink2");
        assert!(
            out.status.success(),
            "plink2 failed: {}",
            String::from_utf8_lossy(&out.stdout)
        );
        let text =
            std::fs::read_to_string(format!("{}.vcf", prefix.display())).expect("invariant: read");
        let columns: Vec<&str> = text
            .lines()
            .find(|l| l.starts_with("#CHROM"))
            .expect("invariant: header")
            .split('\t')
            .collect();
        assert_eq!(
            &columns[9..],
            &["s0", "s2"],
            "plink2 --vcf-min-gq 20 --mind 0.1 keeps the same samples"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

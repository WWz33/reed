//! Determinism gate (`docs/architecture-spec.md` §19 M3): the reports +
//! summary pass over the same BGZF input must produce byte-identical output
//! for every thread count (§7.6 rule 6).

use std::io::Write;
use std::path::PathBuf;

use sha2::Digest;
use reed::reports::{run_reports, ReportRequest};
use reed_engine::budget::{MemBudget, ThreadBudget};

/// Writes a synthetic BGZF VCF: mixed variant classes over 8 samples.
fn make_fixture(dir: &str) -> PathBuf {
    let mut text = String::from(
        "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100000>\n\
         ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for s in 0..8 {
        text.push_str(&format!("\ts{s}"));
    }
    text.push('\n');
    // 3000 records cycling through SNV / MNV / indel / symbolic / ref /
    // multiallelic, with a GT region that forces real sample parsing work.
    for i in 0..3_000 {
        let (reference, alt) = match i % 6 {
            0 => ("A", "G"),            // SNV
            1 => ("AG", "CT"),          // MNV
            2 => ("A", "ATT"),          // insertion
            3 => ("C", "<DEL>"),        // symbolic
            4 => ("T", "."),            // monomorphic
            _ => ("G", "GA,GC"),        // multiallelic
        };
        text.push_str(&format!("chr1\t{}\t.\t{reference}\t{alt}\t50\tPASS\tDP=9\tGT", 1000 + i * 7));
        for s in 0..8 {
            let gt = match (i + s) % 5 {
                0 => "0/0",
                1 => "0/1",
                2 => "1/1",
                3 => "./.",
                _ => "1|0",
            };
            text.push_str(&format!("\t{gt}"));
        }
        text.push('\n');
    }
    let path = PathBuf::from(format!("{dir}/reed-det-{}.vcf.gz", std::process::id()));
    let mut deflater = reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6");
    let file = std::fs::File::create(&path).expect("invariant: fixture file");
    let mut writer = std::io::BufWriter::new(file);
    let mut blocks = Vec::new();
    // 1.5 kB blocks force plenty of block boundaries (mid-line cuts are
    // fine — the chunk reader handles seams).
    for piece in text.as_bytes().chunks(1_500) {
        blocks.clear();
        deflater
            .deflate_block(piece, &mut blocks)
            .expect("invariant: fixture block fits");
        writer.write_all(&blocks).expect("invariant: fixture write");
    }
    writer
        .write_all(&reed_io::bgzf::block::EOF_BLOCK)
        .expect("invariant: fixture write");
    path
}

fn budget(threads: usize) -> ThreadBudget {
    ThreadBudget {
        threads,
        mem: MemBudget::from_bytes(u64::MAX),
    }
}

/// The report family plus the bcftools-stats summary must be byte-identical
/// for every thread count: the rendered summary text and every written
/// report file go into one digest per run.
#[test]
fn reports_and_summary_are_byte_identical_across_thread_counts() {
    let dir = std::env::temp_dir().join(format!("reed-det-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = make_fixture(dir.to_str().expect("invariant: temp path"));

    let mut digests = Vec::new();
    for threads in [1_usize, 2, 3, 7, 56] {
        let prefix = dir.join(format!("out-{threads}"));
        let request = ReportRequest {
            freq: true,
            hardy: true,
            het: true,
            missing: true,
            depth: true,
            summary: true,
            targets: Vec::new(),
            targets_exclude: Vec::new(),
            prefix: prefix.clone(),
        };
        let outcome = run_reports(Some(&input), &request, &budget(threads))
            .expect("invariant: reports run");
        assert_eq!(outcome.n_records, 3_000, "threads={threads}");
        let qc = outcome.qc.expect("invariant: summary requested");
        assert_eq!(qc.summary.n_records, 3_000, "threads={threads}");
        assert_eq!(qc.n_samples, 8, "threads={threads}");

        let mut hasher = sha2::Sha256::new();
        for path in &outcome.files {
            hasher.update(std::fs::read(path).expect("invariant: report file"));
        }
        hasher.update(reed::summary::render_bcftools_stats(&qc).as_bytes());
        digests.push(format!("{:x}", hasher.finalize()));
    }
    let unique: std::collections::HashSet<&String> = digests.iter().collect();
    assert_eq!(unique.len(), 1, "all thread counts must agree: {digests:?}");

    let _ = std::fs::remove_dir_all(&dir);
}

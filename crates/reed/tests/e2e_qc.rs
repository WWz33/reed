//! QC summary counts against a hand-computed fixture (§10.1), including
//! records whose primary ALT class differs from their extra alleles.

use std::io::Write;

use reed::qc::run_qc;
use reed_engine::budget::{MemBudget, ThreadBudget};

#[test]
fn summary_counts_match_hand_computation() {
    let dir = std::env::temp_dir().join(format!("reed-qc-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let text = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=100000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1
chr1\t200\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/0\t0/1
chr1\t300\t.\tC\t.\t50\tPASS\t.\tGT\t0/0\t0/0
chr1\t400\t.\tAC\tA\t50\tPASS\t.\tGT\t0/1\t./.
chr1\t500\t.\tAG\tCT\t50\tPASS\t.\tGT\t1|1\t0/0
";
    let path = dir.join("hand.vcf.gz");
    let mut deflater = reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6");
    let mut file = std::fs::File::create(&path).expect("invariant: fixture file");
    let mut compressed = Vec::new();
    deflater.deflate_block(text.as_bytes(), &mut compressed).expect("invariant: single block");
    file.write_all(&compressed).expect("invariant: write");
    file.write_all(&reed_io::bgzf::block::EOF_BLOCK).expect("invariant: write");
    drop(file);

    let report = run_qc(
        &path,
        &ThreadBudget { threads: 4, mem: MemBudget::from_bytes(u64::MAX) },
    )
    .expect("invariant: run_qc");

    assert_eq!(report.summary.n_records, 5);
    // The multiallelic row's primary ALT (G) is also an SNV.
    assert_eq!(report.summary.n_snv, 2);
    assert_eq!(report.summary.n_mnv, 1);
    assert_eq!(report.summary.n_indel, 1);
    assert_eq!(report.summary.n_ref, 1);
    assert_eq!(report.summary.n_multiallelic, 1);
    assert_eq!(report.summary.n_symbolic, 0);
    assert_eq!(report.n_samples, 2);
    assert_eq!(report.n_chunks, 1);

    let _ = std::fs::remove_dir_all(&dir);
}

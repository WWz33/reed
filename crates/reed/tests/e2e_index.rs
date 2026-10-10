//! `reed index` acceptance: the built CSI/TBI resolve region queries to the
//! same records a scan finds, and the guards (unsorted input, TBI
//! coordinate ceiling) fire.

use std::io::Write;
use std::path::{Path, PathBuf};

use reed::{WriteIndex, run_index};
use std::fmt::Write as _;

fn bgzf(text: &str) -> Vec<u8> {
    let mut deflater =
        reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(64) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block fits");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    out
}

fn write_fixture(name: &str, text: &str) -> (PathBuf, PathBuf) {
    let dir = std::env::temp_dir().join(format!("reed-index-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let path = dir.join("in.vcf.gz");
    let mut file = std::fs::File::create(&path).expect("invariant: fixture file");
    file.write_all(&bgzf(text)).expect("invariant: write");
    drop(file);
    (dir, path)
}

const HEADER: &str = "\
##fileformat=VCFv4.3
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0
";

fn rows() -> String {
    let mut text = String::from(HEADER);
    for pos in [10_u64, 100, 250, 900] {
        let _ = writeln!(text, "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1");
    }
    text.push_str("chr2\t500\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n");
    text
}

/// Reads every record line the index's chunks cover for `contig`.
fn query(index_path: &Path, vcf: &Path, contig: &str) -> Vec<String> {
    use noodles::csi::BinningIndex as _;
    let index = noodles::csi::fs::read(index_path).expect("invariant: index reads");
    let names: Vec<String> = index
        .header()
        .map(|h| {
            h.reference_sequence_names()
                .iter()
                .map(|n| n.to_string())
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    let tid = names
        .iter()
        .position(|n| n == contig)
        .expect("contig in index");
    let chunks: Vec<_> = index.reference_sequences()[tid]
        .bins()
        .values()
        .flat_map(|bin| bin.chunks())
        .copied()
        .collect();
    let mut reader =
        noodles::bgzf::io::Reader::new(std::fs::File::open(vcf).expect("invariant: vcf"));
    use std::io::BufRead as _;
    let query = noodles::csi::io::Query::new(&mut reader, chunks);
    query
        .split(b'\n')
        .map(|l| l.expect("invariant: query read"))
        .filter(|l| !l.is_empty() && l.starts_with(contig.as_bytes()))
        .map(|l| String::from_utf8_lossy(&l).into_owned())
        .collect()
}

#[test]
fn csi_index_resolves_queries() {
    let (dir, vcf) = write_fixture("csi", &rows());
    let outcome = run_index(&vcf, WriteIndex::Csi, None, 4).expect("invariant: index builds");
    assert_eq!(outcome.n_records, 5);
    assert_eq!(outcome.path, dir.join("in.vcf.gz.csi"));

    let got = query(&outcome.path, &vcf, "chr1");
    assert_eq!(got.len(), 4, "{got:?}");
    assert!(got.iter().all(|l| l.starts_with("chr1\t")));
    let got2 = query(&outcome.path, &vcf, "chr2");
    assert_eq!(got2.len(), 1, "{got2:?}");

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn tbi_is_written_and_refuses_wide_coordinates() {
    let (dir, vcf) = write_fixture("tbi", &rows());
    let outcome = run_index(&vcf, WriteIndex::Tbi, None, 4).expect("invariant: tbi builds");
    assert!(outcome.path.to_string_lossy().ends_with(".tbi"));
    // The TBI has its own magic: read it with the tabix reader.
    let index = noodles::tabix::fs::read(&outcome.path).expect("invariant: tbi reads");
    assert_eq!(index.reference_sequences().len(), 2);
    let _ = std::fs::remove_dir_all(&dir);

    // A position ≥ 2^29 cannot be represented in TBI.
    let mut wide = String::from(HEADER);
    wide.push_str("chr1\t600000000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n");
    let (dir2, vcf2) = write_fixture("tbi-wide", &wide);
    let err = run_index(&vcf2, WriteIndex::Tbi, None, 4).expect_err("invariant: TBI refuses 6e8");
    assert!(err.to_string().contains("TBI coordinates"), "{err}");
    let _ = std::fs::remove_dir_all(&dir2);
}

#[test]
fn unsorted_input_is_refused() {
    let mut text = String::from(HEADER);
    text.push_str("chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n");
    text.push_str("chr1\t50\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n");
    let (dir, vcf) = write_fixture("unsorted", &text);
    let err = run_index(&vcf, WriteIndex::Csi, None, 4).expect_err("invariant: unsorted refused");
    assert!(err.to_string().contains("ascending"), "{err}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// O4: the parallel index is byte-identical across thread counts (the chunk
/// framing and the seam chain make the record starts thread-independent).
#[test]
fn index_bytes_are_thread_count_independent() {
    let (dir, vcf) = write_fixture("indexthreads", &rows());
    let a = dir.join("a.csi");
    let b = dir.join("b.csi");
    run_index(&vcf, WriteIndex::Csi, Some(&a), 1).expect("invariant: one thread");
    run_index(&vcf, WriteIndex::Csi, Some(&b), 4).expect("invariant: four threads");
    assert_eq!(
        std::fs::read(&a).expect("invariant: a"),
        std::fs::read(&b).expect("invariant: b")
    );
    let _ = std::fs::remove_dir_all(&dir);
}

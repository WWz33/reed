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

/// Y6 batch G: the index `run_index` wrote is the same artifact
/// `open_regions` consumes for `-r`-style seeks — opening with a region must
/// yield exactly that region's records (the record predicate the filter and
/// stats drivers rely on).
#[test]
fn open_regions_consumes_its_own_index() {
    let (dir, vcf) = write_fixture("open-regions", &rows());
    let outcome = run_index(&vcf, WriteIndex::Csi, None, 4).expect("invariant: index builds");

    let mut opened = reed::open_regions(
        Some(&vcf),
        &[reed::chunk::region::Region {
            contig: "chr1".into(),
            start1: Some(50),
            end1: Some(300),
        }],
    )
    .expect("invariant: open with the built index");
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut text = String::new();
    while let Some(chunk) = opened.source.next_chunk().expect("invariant: chunk read") {
        let bytes = chunk.compressed.as_slice();
        let mut rest = bytes;
        while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
            let mut out = Vec::new();
            inflater
                .inflate_block(&rest[..total], &mut out)
                .expect("invariant: block inflates");
            text.push_str(&String::from_utf8_lossy(&out));
            rest = &rest[total..];
        }
    }
    // The seek is byte-range only (the record predicate is the job's target
    // filter, covered elsewhere); the contract here is that every record the
    // region names arrives through the index-sliced source. Bin granularity
    // may add out-of-region records, so assert on the region's own records.
    let records: Vec<(String, u64)> = text
        .lines()
        // The first chunk may start mid-header-line (the header tail rides
        // the region slice); record lines start with the contig name.
        .filter(|l| l.starts_with("chr1\t") || l.starts_with("chr2\t"))
        .map(|l| {
            let mut it = l.split('\t');
            (
                it.next().expect("invariant: contig").to_string(),
                it.next()
                    .expect("invariant: pos")
                    .parse::<u64>()
                    .expect("invariant: pos numeric"),
            )
        })
        .collect();
    let inside: Vec<u64> = records
        .iter()
        .filter(|(c, p)| c == "chr1" && (50..=300).contains(p))
        .map(|(_, p)| *p)
        .collect();
    assert_eq!(inside, [100, 250], "{records:?}");

    let _ = std::fs::remove_dir_all(&dir);
    let _ = outcome;
}

/// Regular gzip (no BGZF `BC` block) of the same text.
fn gzip(text: &str) -> Vec<u8> {
    use flate2::Compression;
    use flate2::write::GzEncoder;
    let mut enc = GzEncoder::new(Vec::new(), Compression::default());
    enc.write_all(text.as_bytes()).expect("invariant: gzip write");
    enc.finish().expect("invariant: gzip finish")
}

/// doc_v7 D7-3: tabix/CSI address BGZF blocks; a plain VCF or a regular gzip
/// stream has none, so indexing it must be refused — not silently written as
/// an all-zero-offset index.
#[test]
fn plain_and_gzip_inputs_are_refused() {
    let dir = std::env::temp_dir().join(format!("reed-index-plain-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");

    let plain = dir.join("in.vcf");
    std::fs::write(&plain, rows()).expect("invariant: write plain");
    let err = run_index(&plain, WriteIndex::Csi, None, 4).expect_err("plain VCF must be refused");
    assert!(err.to_string().contains("BGZF"), "{err}");

    let gz = dir.join("in.vcf.gz");
    std::fs::write(&gz, gzip(&rows())).expect("invariant: write gzip");
    let err = run_index(&gz, WriteIndex::Csi, None, 4).expect_err("regular gzip must be refused");
    assert!(err.to_string().contains("BGZF"), "{err}");

    // Control: the BGZF form still indexes.
    let (dir2, bgzf_path) = write_fixture("plain-control", &rows());
    run_index(&bgzf_path, WriteIndex::Csi, None, 4).expect("BGZF still indexes");

    let _ = std::fs::remove_dir_all(&dir);
    let _ = std::fs::remove_dir_all(&dir2);
}

/// doc_v7 D7-15: Tabix has no BCF flavour, so `reed index -t` on a BCF input
/// must be refused (it used to write a VCF-format .tbi for a BCF).
#[test]
fn bcf_input_refuses_tbi() {
    let bcf =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../reed-io/tests/fixtures/bcf/fix.bcf");
    let err = run_index(&bcf, WriteIndex::Tbi, None, 4).expect_err("BCF + TBI must be refused");
    assert!(err.to_string().contains("CSI"), "{err}");
    // CSI still works on the same BCF.
    let out = std::env::temp_dir().join(format!("reed-bcf-csi-{}.csi", std::process::id()));
    run_index(&bcf, WriteIndex::Csi, Some(&out), 4).expect("BCF + CSI must work");
    let _ = std::fs::remove_file(&out);
}

/// doc_v7 D7-16: a truncated BCF must produce an error, never an
/// out-of-bounds panic in `declared_total`.
#[test]
fn truncated_bcf_errors_without_panicking() {
    let src =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../reed-io/tests/fixtures/bcf/fix.bcf");
    let bytes = std::fs::read(&src).expect("invariant: fixture");
    let dir = std::env::temp_dir().join(format!("reed-bcf-trunc-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    // Several cut points exercise different crossing-record cases.
    for frac in [10_usize, 20, 50, 75, 90, 99] {
        let cut = bytes.len() * frac / 100;
        let out = dir.join(format!("in{frac}.bcf"));
        std::fs::write(&out, &bytes[..cut]).expect("invariant: write");
        let _ = run_index(&out, WriteIndex::Csi, None, 4)
            .expect_err("a truncated BCF must be reported, not panic");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

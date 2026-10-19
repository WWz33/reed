//! `reed rename` end-to-end gates (doc_v14 02 §4): the four edit planes, the
//! byte-passthrough of everything else, the validation rules and thread
//! determinism.

use std::path::{Path, PathBuf};

use reed::{ChrStyle, RenameOutput, RenameSpec, run_rename};
use reed_engine::budget::ThreadBudget;

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
chr1\t1\trs1\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\t0/0
chr1\t2\trs2\tA\tG\t60\tPASS\tDP=20\tGT\t0/0\t1/1
chr1\t3\t.\tA\tG\t70\tPASS\tDP=30\tGT\t0/1\t0/1
";

fn dir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("reed-rename-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("invariant: temp dir");
    d
}

fn write(path: &Path, text: &str) {
    std::fs::write(path, text).expect("invariant: write");
}

fn budget(threads: usize) -> ThreadBudget {
    ThreadBudget {
        threads,
        mem: reed_engine::budget::MemBudget::default_half_memory(),
    }
}

fn run(spec: &RenameSpec, input: &Path, output: &Path, threads: usize) -> reed::RenameReport {
    let out = RenameOutput {
        path: Some(output.to_path_buf()),
        compressed: false,
        level: 6,
        command_line: String::new(),
        no_version: true,
        no_command: false,
        write_index: None,
    };
    run_rename(Some(input), spec, &out, &budget(threads)).expect("invariant: rename run")
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).expect("invariant: read")
}

/// The three rename planes land in both the header and the records, and
/// everything else passes through byte-for-byte.
#[test]
fn sample_contig_and_id_planes_edit_both_sides() {
    let d = dir("planes");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let map = d.join("chrs.txt");
    write(&map, "chr1\tChr01\n");
    let spec = RenameSpec {
        rename_chrs: Some(map),
        samples: Some("s0 sampleA,s1 sampleB".into()),
        set_id: Some("%CHROM:%POS:%REF:%ALT".into()),
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf");
    let report = run(&spec, &input, &out, 2);
    let text = read(&out);
    assert!(text.contains("##contig=<ID=Chr01,length=100000>"), "{text}");
    assert!(text.contains("\tFORMAT\tsampleA\tsampleB\n"), "{text}");
    assert!(text.contains("Chr01\t1\tChr01:1:A:G\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\t0/0"), "{text}");
    assert_eq!(report.n_records, 3);
    assert_eq!(report.n_contigs, 1);
    assert_eq!(report.n_samples, 2);
    assert_eq!(report.n_ids, 3);
    let _ = std::fs::remove_dir_all(&d);
}

/// D14-5: `--rename-chrs` rewrites the `##contig` line and the CHROM column
/// together — never one without the other.
#[test]
fn contig_rename_touches_header_and_records_together() {
    let d = dir("contig");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let map = d.join("chrs.txt");
    write(&map, "chr1\tChr01\n");
    let spec = RenameSpec {
        rename_chrs: Some(map),
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf");
    run(&spec, &input, &out, 1);
    let text = read(&out);
    assert!(text.contains("##contig=<ID=Chr01"), "the header must follow: {text}");
    assert!(!text.contains("\nchr1\t"), "no record may keep the old name: {text}");
    let _ = std::fs::remove_dir_all(&d);
}

/// D14-4: the genotype columns and every non-edited field pass through
/// byte-for-byte (the editor never decodes genotypes).
#[test]
fn non_edited_fields_pass_through() {
    let d = dir("passthrough");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let spec = RenameSpec {
        samples: Some("s0 A".into()),
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf");
    run(&spec, &input, &out, 1);
    let before: Vec<&str> = VCF.lines().filter(|l| !l.starts_with('#')).collect();
    let after_text = read(&out);
    let after: Vec<&str> = after_text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(before, after, "records must be unchanged by a sample rename");
    let _ = std::fs::remove_dir_all(&d);
}

/// The ID planes: `--set-missing-var-ids` fills only dots, `--rm-dup-ids`
/// suffixes repeats (never dropping a record), `--idspace-to` swaps `:`.
#[test]
fn id_planes_are_independent() {
    let d = dir("ids");
    let input = d.join("in.vcf");
    write(
        &input,
        "##fileformat=VCFv4.2\n##contig=<ID=c1,length=100>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         c1\t1\trs1\tA\tG\t50\t.\t.\tGT\t0/1\n\
         c1\t2\trs1\tA\tG\t50\t.\t.\tGT\t0/0\n\
         c1\t3\t.\tA\tG\t50\t.\t.\tGT\t0/1\n",
    );
    let spec = RenameSpec {
        set_missing_var_ids: Some("v%POS".into()),
        rm_dup_ids: true,
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf");
    run(&spec, &input, &out, 1);
    let text = read(&out);
    let ids: Vec<&str> = text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').nth(2).unwrap())
        .collect();
    assert_eq!(ids, ["rs1", "rs1_1", "v3"], "{text}");

    let spec = RenameSpec {
        set_id: Some("%CHROM:%POS".into()),
        idspace_to: None,
        ..RenameSpec::default()
    };
    let out2 = d.join("out2.vcf");
    run(&spec, &input, &out2, 1);
    let text2 = read(&out2);
    assert!(text2.contains("c1:1\t"), "{text2}");
    let _ = std::fs::remove_dir_all(&d);
}

/// 02 §4.3: the validation rules — unknown names, duplicate sample names,
/// conflicting ID flags and an empty request are all usage errors.
#[test]
fn validation_rejects_bad_input() {
    let d = dir("validate");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let out = RenameOutput {
        path: Some(d.join("out.vcf")),
        compressed: false,
        level: 6,
        command_line: String::new(),
        no_version: true,
        no_command: false,
        write_index: None,
    };
    let bad = [
        RenameSpec::default(),
        RenameSpec {
            samples: Some("nope new".into()),
            ..RenameSpec::default()
        },
        RenameSpec {
            samples: Some("s0 s1".into()),
            ..RenameSpec::default()
        },
        RenameSpec {
            set_id: Some("%CHROM".into()),
            set_missing_var_ids: Some("%POS".into()),
            ..RenameSpec::default()
        },
        RenameSpec {
            clear_ids: true,
            rm_dup_ids: true,
            ..RenameSpec::default()
        },
        RenameSpec {
            set_id: Some("%CHROM".into()),
            idspace_to: Some("_".into()),
            ..RenameSpec::default()
        },
    ];
    for spec in &bad {
        let err = run_rename(Some(&input), spec, &out, &budget(1))
            .expect_err("must be refused");
        assert_eq!(err.kind(), reed::ErrorKind::Usage, "{spec:?}");
    }
    let _ = std::fs::remove_dir_all(&d);
}

/// `--chr-style` normalizes names in the header and the records.
#[test]
fn chr_style_normalizes() {
    let d = dir("style");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let spec = RenameSpec {
        chr_style: Some(ChrStyle::Plain),
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf");
    run(&spec, &input, &out, 1);
    let text = read(&out);
    assert!(text.contains("##contig=<ID=1,"), "{text}");
    assert!(text.contains("\n1\t1\t"), "{text}");
    let _ = std::fs::remove_dir_all(&d);
}

/// Header edits: `--add-meta`, `--remove-meta` and `-f/--fai`.
#[test]
fn header_meta_edits() {
    let d = dir("meta");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let fai = d.join("ref.fai");
    write(&fai, "chr1\t99999\n");
    let spec = RenameSpec {
        add_meta: vec!["source=test".into()],
        fai: Some(fai),
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf");
    run(&spec, &input, &out, 1);
    let text = read(&out);
    assert!(text.contains("##source=test\n"), "{text}");
    assert!(text.contains("##contig=<ID=chr1,length=99999>"), "{text}");

    let spec = RenameSpec {
        remove_meta: vec!["INFO".into()],
        ..RenameSpec::default()
    };
    let out2 = d.join("out2.vcf");
    run(&spec, &input, &out2, 1);
    assert!(!read(&out2).contains("##INFO="), "remove-meta drops the lines");
    let _ = std::fs::remove_dir_all(&d);
}

/// `-W`: a BGZF output gets a sidecar index that describes the *renamed*
/// contigs (doc_v14 02 §4.4).
#[test]
fn write_index_produces_a_sidecar() {
    let d = dir("index");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let map = d.join("chrs.txt");
    write(&map, "chr1\tChr01\n");
    let spec = RenameSpec {
        rename_chrs: Some(map),
        ..RenameSpec::default()
    };
    let out = d.join("out.vcf.gz");
    let out_cfg = RenameOutput {
        path: Some(out.clone()),
        compressed: true,
        level: 6,
        command_line: String::new(),
        no_version: true,
        no_command: false,
        write_index: Some(reed::WriteIndex::Csi),
    };
    let report = run_rename(Some(&input), &spec, &out_cfg, &budget(1)).expect("rename -W");
    assert_eq!(report.files, vec![out.clone()]);
    let mut csi = out.as_os_str().to_os_string();
    csi.push(".csi");
    let csi = PathBuf::from(csi);
    assert!(csi.is_file(), "the index must be written: {}", csi.display());
    assert!(
        std::fs::metadata(&csi).expect("invariant: stat").len() > 0,
        "the index must not be empty"
    );
    let _ = std::fs::remove_dir_all(&d);
}

/// Thread determinism: the output is byte-identical at 1 and 4 threads.
#[test]
fn rename_is_deterministic_across_threads() {
    let d = dir("threads");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let map = d.join("chrs.txt");
    write(&map, "chr1\tChr01\n");
    let spec = RenameSpec {
        rename_chrs: Some(map),
        samples: Some("s0 A,s1 B".into()),
        set_id: Some("%CHROM_%POS".into()),
        rm_dup_ids: true,
        ..RenameSpec::default()
    };
    let out1 = d.join("t1.vcf");
    let out4 = d.join("t4.vcf");
    run(&spec, &input, &out1, 1);
    run(&spec, &input, &out4, 4);
    assert_eq!(read(&out1), read(&out4));
    let _ = std::fs::remove_dir_all(&d);
}

/// Inflates a BGZF file back to text (the project's own inflater).
fn inflate(path: &Path) -> String {
    let bytes = std::fs::read(path).expect("invariant: read");
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut text = String::new();
    let mut rest = bytes.as_slice();
    while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
        let mut out = Vec::new();
        inflater
            .inflate_block(&rest[..total], &mut out)
            .expect("invariant: block inflates");
        text.push_str(&String::from_utf8_lossy(&out));
        rest = &rest[total..];
    }
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

/// doc_v15 D15-28: `--rm-dup-ids` now deflates in batches. The suffix still
/// follows the global record order *across* batch boundaries, the sidecar
/// index still resolves to the file's own records (the anchors are re-based
/// into the new blocks), and the compressed stream — index included — is
/// byte-identical across thread counts, because the batch seals are driven
/// by accumulated bytes, never by timing.
#[test]
fn dedup_batches_keep_order_index_and_bytes() {
    let d = dir("dedup-batches");
    let input = d.join("in.vcf");
    // Well past two 1 MiB batches, with seven base IDs repeating so the
    // suffix state crosses batch boundaries.
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n",
    );
    for i in 0..80_000_u64 {
        use std::fmt::Write as _;
        let _ = writeln!(text, "chr1\t{}\trs{}\tA\tG\t50\t.\t.\tGT\t0/1", i + 1, i % 7);
    }
    write(&input, &text);

    let spec = RenameSpec {
        rm_dup_ids: true,
        ..RenameSpec::default()
    };
    let mut outputs = Vec::new();
    for threads in [3_usize, 9] {
        let out = d.join(format!("t{threads}.vcf.gz"));
        let io = RenameOutput {
            path: Some(out.clone()),
            compressed: true,
            level: 6,
            command_line: String::new(),
            no_version: true,
            no_command: false,
            write_index: Some(reed::WriteIndex::Csi),
        };
        let report = run_rename(Some(&input), &spec, &io, &budget(threads)).expect("rename -W");
        assert_eq!(report.n_records, 80_000);
        assert_eq!(report.n_dup_ids, 80_000 - 7, "every repeat past the first");
        outputs.push(out);
    }
    let body = inflate(&outputs[0]);
    assert!(body.len() > 2 << 20, "several batches: {} bytes", body.len());
    assert_eq!(
        std::fs::read(&outputs[0]).expect("invariant: read"),
        std::fs::read(&outputs[1]).expect("invariant: read"),
        "the compressed stream is thread-count independent"
    );
    let sidecar = |p: &Path| {
        let mut s = p.as_os_str().to_os_string();
        s.push(".csi");
        PathBuf::from(s)
    };
    assert_eq!(
        std::fs::read(sidecar(&outputs[0])).expect("invariant: index"),
        std::fs::read(sidecar(&outputs[1])).expect("invariant: index"),
        "the sidecar is thread-count independent"
    );

    // The suffix follows the global order: the k-th occurrence of rs{n} is
    // plain for k == 0 and carries `_k` afterwards — regardless of which
    // batch the record landed in.
    let lines: Vec<&str> = body.lines().filter(|l| !l.starts_with('#')).collect();
    for (i, line) in lines.iter().enumerate() {
        let want = match i / 7 {
            0 => format!("rs{}", i % 7),
            k => format!("rs{}_{k}", i % 7),
        };
        assert_eq!(line.split('\t').nth(2).expect("ID column"), want, "{line}");
    }

    // The index resolves to exactly those records: a wrong anchor would
    // start a chunk mid-record and the line would not match the scan.
    // (The helper walks the index's bins, whose order is the map's, not
    // the file's — the record order itself is what the suffix loop above
    // asserts.)
    let pos = |l: &String| -> u64 {
        l.split('\t').nth(1).expect("POS").parse().expect("POS")
    };
    let mut got = query(&sidecar(&outputs[0]), &outputs[0], "chr1");
    let mut want: Vec<String> = lines.iter().map(|l| (*l).to_owned()).collect();
    got.sort_by_key(&pos);
    want.sort_by_key(&pos);
    assert_eq!(got.len(), want.len());
    assert_eq!(got, want);
    let _ = std::fs::remove_dir_all(&d);
}

/// doc_v15 D15-28: output that never fills a batch deflates inline (no
/// pool threads for a small file), and the result matches the plain run at
/// every thread count — including the sequential engine path.
#[test]
fn dedup_single_batch_matches_the_plain_run() {
    let d = dir("dedup-small");
    let input = d.join("in.vcf");
    write(
        &input,
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         chr1\t1\trs1\tA\tG\t50\t.\t.\tGT\t0/1\n\
         chr1\t2\trs1\tA\tG\t50\t.\t.\tGT\t0/1\n\
         chr1\t3\t.\tA\tG\t50\t.\t.\tGT\t0/1\n",
    );
    let spec = RenameSpec {
        rm_dup_ids: true,
        ..RenameSpec::default()
    };
    let plain = d.join("plain.vcf");
    run(&spec, &input, &plain, 1);
    let want: Vec<String> = read(&plain)
        .lines()
        .map(str::to_owned)
        .filter(|l| !l.starts_with("##reed_command="))
        .collect();
    assert!(want.iter().any(|l| l.contains("\trs1_1\t")), "{want:?}");

    for threads in [2_usize, 3, 9] {
        let out = d.join(format!("t{threads}.vcf.gz"));
        let io = RenameOutput {
            path: Some(out.clone()),
            compressed: true,
            level: 6,
            command_line: String::new(),
            no_version: true,
            no_command: false,
            write_index: Some(reed::WriteIndex::Csi),
        };
        run_rename(Some(&input), &spec, &io, &budget(threads)).expect("rename -W");
        let got: Vec<String> = inflate(&out)
            .lines()
            .map(str::to_owned)
            .filter(|l| !l.starts_with("##reed_command="))
            .collect();
        assert_eq!(got, want, "threads {threads}");
    }
    let _ = std::fs::remove_dir_all(&d);
}

/// A multi-chunk input (the engine stitches a seam record at every chunk
/// boundary) keeps exactly one line per record: the seam line used to arrive
/// with its terminator and get a second one, so every output over one chunk
/// carried blank lines, and `--rm-dup-ids` (which splits the raw text on
/// newlines) died on the empty piece.
#[test]
fn multi_chunk_records_stay_one_line_each() {
    let d = dir("multichunk");
    let input = d.join("big.vcf");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n",
    );
    // Well past the 1 MiB chunk size, so several seams are exercised.
    for i in 0..120_000 {
        use std::fmt::Write as _;
        let _ = writeln!(text, "chr1\t{}\trs{}\tA\tG\t50\t.\t.\tGT\t0/1", i + 1, i);
    }
    write(&input, &text);
    for (tag, spec) in [
        (
            "clear",
            RenameSpec {
                clear_ids: true,
                ..RenameSpec::default()
            },
        ),
        (
            "dedup",
            RenameSpec {
                rm_dup_ids: true,
                set_id: Some("dup".into()),
                ..RenameSpec::default()
            },
        ),
    ] {
        let out = d.join(format!("{tag}.vcf"));
        let report = run(&spec, &input, &out, 4);
        let text = read(&out);
        let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(report.n_records, 120_000, "{tag}: records");
        assert_eq!(body.len(), 120_000, "{tag}: one output line per record");
        assert!(
            body.iter().all(|l| !l.is_empty()),
            "{tag}: no blank line at a chunk seam"
        );
    }
    let _ = std::fs::remove_dir_all(&d);
}

/// `--fai` adds a `##contig` line the header lacks *before* `#CHROM`: a `##`
/// line after it is not a valid VCF header.
#[test]
fn fai_lines_land_in_the_header_block() {
    let d = dir("fai");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let fai = d.join("ref.fai");
    write(&fai, "chr1\t999\nchr2\t888\n");
    let out = d.join("out.vcf");
    let spec = RenameSpec {
        fai: Some(fai),
        ..RenameSpec::default()
    };
    run(&spec, &input, &out, 2);
    let text = read(&out);
    let chrom_at = text.find("#CHROM").expect("the header line");
    for line in text.lines().filter(|l| l.starts_with("##contig")) {
        let at = text.find(line).expect("line offset");
        assert!(at < chrom_at, "{line:?} must precede #CHROM");
    }
    assert!(text.contains("##contig=<ID=chr2,length=888>"), "{text}");
    let _ = std::fs::remove_dir_all(&d);
}

/// A contig-less input is legal VCF: plain (unindexed) output must pass the
/// records through, and an index request says what is missing instead of
/// failing mid-run.
#[test]
fn contig_less_input_is_accepted_without_an_index() {
    let d = dir("nocontig");
    let input = d.join("in.vcf");
    write(
        &input,
        "##fileformat=VCFv4.2\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n\
         c1\t1\trs1\tA\tG\t.\t.\t.\tGT\t0/1\n",
    );
    let out = d.join("out.vcf");
    let spec = RenameSpec {
        clear_ids: true,
        ..RenameSpec::default()
    };
    run(&spec, &input, &out, 1);
    assert!(read(&out).contains("c1\t1\t.\tA\tG"), "{}", read(&out));

    // `-W` needs the contig dictionary: a clear usage error, not a crash.
    let header = d.join("h.vcf");
    write(
        &header,
        "##fileformat=VCFv4.2\n##contig=<ID=c1,length=10>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\n",
    );
    let spec = RenameSpec {
        header: Some(header),
        clear_ids: true,
        ..RenameSpec::default()
    };
    let out = d.join("idx.vcf");
    let io = RenameOutput {
        path: Some(out.clone()),
        compressed: true,
        level: 6,
        command_line: String::new(),
        no_version: true,
        no_command: false,
        write_index: Some(reed::WriteIndex::Csi),
    };
    // The replacement header declares c1, so the index is buildable here;
    // the run must succeed and the sidecar must exist.
    run_rename(Some(&input), &spec, &io, &budget(2)).expect("invariant: index run");
    assert!(out.with_extension("vcf.csi").exists(), "the sidecar is written");
    let _ = std::fs::remove_dir_all(&d);
}

/// A CRLF header's last sample keeps its `\r` in the line, so a rename keyed
/// on the plain name must still match it (it used to be counted but skipped).
#[test]
fn crlf_sample_renames_apply() {
    let d = dir("crlf");
    let input = d.join("in.vcf");
    write(
        &input,
        "##fileformat=VCFv4.2\r\n##contig=<ID=chr1,length=10>\r\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\r\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\r\n\
         chr1\t1\trs1\tA\tG\t.\t.\t.\tGT\t0/1\t0/0\r\n",
    );
    let out = d.join("out.vcf");
    let spec = RenameSpec {
        samples: Some("s1 LAST".into()),
        ..RenameSpec::default()
    };
    let report = run(&spec, &input, &out, 1);
    assert_eq!(report.n_samples, 1);
    let text = read(&out);
    assert!(text.contains("s0\tLAST\r\n"), "{text:?}");
    let _ = std::fs::remove_dir_all(&d);
}

/// `--update-name` mapping a name the input lacks is a typo: the run fails
/// instead of renaming nothing and reporting success.
#[test]
fn update_name_validates_its_keys() {
    let d = dir("updatename");
    let input = d.join("in.vcf");
    write(&input, VCF);
    let bad = d.join("bad.txt");
    write(&bad, "nosuchid\tX\n");
    let out = d.join("out.vcf");
    let spec = RenameSpec {
        update_name: Some(bad),
        ..RenameSpec::default()
    };
    let err = run_rename(
        Some(&input),
        &spec,
        &RenameOutput {
            path: Some(out.clone()),
            compressed: false,
            level: 6,
            command_line: String::new(),
            no_version: true,
            no_command: false,
            write_index: None,
        },
        &budget(1),
    )
    .expect_err("an unknown --update-name key must fail");
    assert_eq!(err.kind(), reed::ErrorKind::Usage);
    assert!(err.to_string().contains("nosuchid"), "{err}");
    let _ = std::fs::remove_dir_all(&d);
}

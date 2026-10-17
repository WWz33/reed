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
    assert_eq!(ids, ["rs1", "rs1_2", "v3"], "{text}");

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

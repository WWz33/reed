//! `reed convert` acceptance (docs/rounds/v13): the expected bytes are the
//! *actual* outputs of the reference tools on the same fixture — vcf2phylip
//! 2.9 (PHYLIP/FASTA/NEXUS/binary-NEXUS/used_sites), vcftools 0.1.18
//! (`--012`) and vcf2msa (regions mapped to reed's inclusive convention) —
//! recorded in gold/popgen/run_convert_parity.sh.

use std::path::Path;

use reed::budget::{MemBudget, ThreadBudget};
use reed::{
    ConvertRequest, ConvertTarget, MsaRegion, PhyloFormats, PhyloOptions, Populations, SfsRequest,
    PopgenRequest, run_convert, run_popgen,
};

/// The phylogenetic fixture: 20 records over 8 samples, covering biallelic
/// and multiallelic SNPs, MNPs, indels, `*`/`<NON_REF>`/`.` ALT values,
/// phased and haploid genotypes, half calls, a fully missing sample and one
/// malformed (column-count) line.
const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=100000>
##contig=<ID=chr2,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\ts4\ts5\ts6\ts7
chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\t0/1\t1/1\t1/1\t./.
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/1\t0/1\t1/1\t1/1\t./.
chr1\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\t0/1\t1/1\t1/1\t./.
chr1\t4\t.\tC\tT\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\t1/1\t1/1\t./.
chr1\t5\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/0\t0/1\t2/2\t1/2\t0/2\t1/1\t2/2\t./.
chr1\t6\t.\tAT\tGC\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\t1/1\t0/0\t0/1\t./.
chr1\t7\t.\tA\tGT\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1\t1/1\t0/0\t0/1\t./.
chr1\t8\t.\tA\tG,*\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t2/2\t0/2\t1/1\t2/2\t./.
chr1\t9\t.\tA\t<NON_REF>\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\t0/1\t1/1\t1/1\t./.
chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t./.
chr1\t11\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\t1|0\t0|0\t1|1\t0/0\t0/1\t1/1\t./.
chr1\t12\t.\tA\tG\t50\tPASS\t.\tGT\t0\t0\t1\t1\t0\t1\t0\t1
chr1\t14\t.\tA\tG\t50\tPASS\t.\tGT\t0/.\t0/0\t0/0\t0/0\t0/1\t1/1\t1/1\t./.
chr1\t15\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t./.
chr1\t16\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1
chr2\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\t0/1\t1/1\t0/0\t./.
chr2\t2\t.\tC\tA\t50\tPASS\t.\tGT\t0/0\t0/0\t0/1\t1/1\t0/0\t0/1\t1/1\t./.
chr2\t3\t.\tT\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1\t1/1\t0/0\t0/1\t./.
chr2\t4\t.\tA\tC\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\t0/1\t1/1\t0/0\t./.
chr2\t5\t.\tG\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t0/1\t./.
";

/// vcf2phylip 2.9's `.min4.phy` for the fixture, `--outgroup s3`.
const REF_PHYLIP: &str = "\
8 17
s3   RRACK-AAGGAAAAKAK
s0   ARACAAAARAaAACKAK
s1   AARCRRAARAAARCGMK
s2   AAGCTGAAAGAAGMTCK
s4   RRRYWaAAAARARCGMK
s5   GGGTGGAARGGAGMTCK
s6   GGGTT-AAGAGAAAKAK
s7   NNNNNNNNNGNNNNNNN
";

/// vcf2phylip 2.9's `.bin.nexus` (NCHAR 14: the biallelic single-character
/// ALT stream), `--outgroup s3`.
const REF_BIN_NEXUS: &str = "\
#NEXUS

BEGIN DATA;
\tDIMENSIONS NTAX=8 NCHAR=14;
\tFORMAT DATATYPE=SNP MISSING=? GAP=- ;
MATRIX
s3   110002?0002101
s0   010001??000101
s1   001001?0010211
s2   002000?0021021
s4   111100?1010211
s5   222201?2021021
s6   222202?2002101
s7   ??????????????
;
END;
";

/// vcf2phylip 2.9's `.used_sites.tsv` (nucleotide-accepted sites only).
const REF_USED_SITES: &str = "\
#CHROM\tPOS\tNUM_SAMPLES
chr1\t1\t7
chr1\t2\t7
chr1\t3\t7
chr1\t4\t7
chr1\t5\t7
chr1\t8\t7
chr1\t9\t7
chr1\t10\t7
chr1\t11\t7
chr1\t12\t8
chr1\t14\t7
chr1\t15\t7
chr2\t1\t7
chr2\t2\t7
chr2\t3\t7
chr2\t4\t7
chr2\t5\t7
";

fn budget() -> ThreadBudget {
    ThreadBudget {
        threads: 2,
        mem: MemBudget::from_bytes(256 << 20),
    }
}

fn fixture_dir(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("reed-convert-{tag}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    dir
}

fn write_fixture(dir: &Path, name: &str, text: &str) -> std::path::PathBuf {
    let path = dir.join(name);
    std::fs::write(&path, text).expect("invariant: fixture write");
    path
}

fn read(path: &Path) -> String {
    std::fs::read_to_string(path).expect("invariant: output exists")
}

/// The full phylogenetic stream on the fixture reproduces vcf2phylip 2.9
/// byte for byte: PHYLIP, FASTA, NEXUS, binary NEXUS and used_sites
/// (doc_v13 §9.1/§9.2).
#[test]
fn phylo_outputs_match_vcf2phylip_reference() {
    let dir = fixture_dir("phylo");
    let vcf = write_fixture(&dir, "in.vcf", VCF);
    let request = ConvertRequest {
        targets: vec![
            ConvertTarget::Phylo {
                formats: PhyloFormats {
                    phylip: true,
                    fasta: true,
                    nexus: true,
                },
                options: PhyloOptions {
                    outgroup: Some("s3".into()),
                    resolve_iupac: false,
                    write_used_sites: true,
                },
            },
            ConvertTarget::NexusBin,
        ],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    let report = run_convert(Some(&vcf), &request, &budget()).expect("convert");
    assert_eq!(read(&dir.join("out.phy")), REF_PHYLIP);
    assert_eq!(read(&dir.join("out.bin.nexus")), REF_BIN_NEXUS);
    assert_eq!(read(&dir.join("out.used_sites.tsv")), REF_USED_SITES);

    // FASTA is the same model: header + unpadded sequence per sample.
    let fasta = read(&dir.join("out.fasta"));
    assert!(fasta.starts_with(">s3\nRRACK-AAGGAAAAKAK\n"), "{fasta}");
    // NEXUS is the padded model with the DNA header/tail.
    let nexus = read(&dir.join("out.nexus"));
    assert!(
        nexus.starts_with(
            "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=8 NCHAR=17;\n\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n"
        ),
        "{nexus}"
    );
    assert!(nexus.ends_with(";\nEND;\n"), "{nexus}");

    // The counters mirror vcf2phylip's stdout: 20 lines processed, 0
    // shallow, 2 MNPs, 17 accepted, 14 biallelic, 1 malformed (B13-2).
    let counters = report.counters.expect("row-target counters");
    assert_eq!(counters.total_records, 20);
    assert_eq!(counters.excess_missing, 0);
    assert_eq!(counters.mnps, 2);
    assert_eq!(counters.accepted, 17);
    assert_eq!(counters.biallelic_for_bin, 14);
    assert_eq!(counters.malformed, 1);
    let _ = std::fs::remove_dir_all(&dir);
}

/// B13-9: a record whose REF/ALT falls outside the AMBIG table, or whose GT
/// names an index beyond the declared alleles, is skipped and counted — the
/// reference script crashes on the first and skips the second; reed unifies
/// both onto the malformed path (B13-2's counter).
#[test]
fn ambig_miss_and_index_overflow_are_skipped_and_counted() {
    const VCF_BAD: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1
chr1\t2\t.\tR\tG\t.\tPASS\t.\tGT\t0/0\t0/1
chr1\t3\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t3/3
chr1\t4\t.\tA\tG\t.\tPASS\t.\tGT\t1/1\t0/1
";
    let dir = fixture_dir("b139");
    let vcf = write_fixture(&dir, "in.vcf", VCF_BAD);
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Phylo {
            formats: PhyloFormats {
                phylip: true,
                ..PhyloFormats::default()
            },
            options: PhyloOptions::default(),
        }],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    let report = run_convert(Some(&vcf), &request, &budget()).expect("convert");
    let phy = read(&dir.join("out.phy"));
    assert!(
        phy.starts_with("2 2\n"),
        "the two good sites survive: {phy}"
    );
    // The survivors are chr1:1 and chr1:4.
    assert!(phy.contains("s0   AG\n"), "{phy}");
    assert!(phy.contains("s1   RR\n"), "{phy}");
    let counters = report.counters.expect("counters");
    assert_eq!(counters.total_records, 4);
    assert_eq!(counters.malformed, 2, "REF=R and the index-3 GT");
    assert_eq!(counters.accepted, 2);
    let _ = std::fs::remove_dir_all(&dir);
}

/// B13-9: the same skip path on the hapmap→VCF side — a cell naming an
/// undeclared allele drops the whole row and counts it, rather than aborting
/// the run.
#[test]
fn hapmap_vcf_skips_rows_with_undeclared_alleles() {
    const HMP: &str = "\
rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	s0	s1
chr1:1	A/G	chr1	1	+	NA	NA	NA	NA	NA	NA	A	R
chr1:2	A/G	chr1	2	+	NA	NA	NA	NA	NA	NA	A	C
chr1:3	A/G	chr1	3	+	NA	NA	NA	NA	NA	NA	G	T
";
    let dir = fixture_dir("hmp_bad");
    let hmp = write_fixture(&dir, "in.hmp.txt", HMP);
    let out = dir.join("out.vcf");
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Vcf {
            compressed: false,
            level: 6,
            output: Some(out.clone()),
            command_line: String::new(),
            no_version: false,
        }],
        ..ConvertRequest::default()
    };
    let report = run_convert(Some(&hmp), &request, &budget()).expect("convert");
    let text = read(&out);
    // Only the good row is written; the `C` and `T` cells are not declared.
    assert_eq!(report.n_records, 1, "{text}");
    let counters = report.counters.expect("counters");
    assert_eq!(counters.malformed, 2);
    assert_eq!(
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .collect::<Vec<_>>(),
        vec!["chr1\t1\tchr1:1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// The `--min-samples-locus` gate (default 4) and its clamp to the sample
/// count behave like vcf2phylip: `-m 8` leaves 1 site (the all-called
/// haploid row), and `-m 100` clamps to 8 (same as `-m 8`).
#[test]
fn min_samples_locus_gates_and_clamps() {
    let dir = fixture_dir("gate");
    let vcf = write_fixture(&dir, "in.vcf", VCF);
    let run = |gate: u32, prefix: &str| {
        let request = ConvertRequest {
            targets: vec![ConvertTarget::Phylo {
                formats: PhyloFormats {
                    phylip: true,
                    ..PhyloFormats::default()
                },
                options: PhyloOptions::default(),
            }],
            prefix: dir.join(prefix),
            min_samples_locus: Some(gate),
            ..ConvertRequest::default()
        };
        run_convert(Some(&vcf), &request, &budget()).expect("convert");
        read(&dir.join(format!("{prefix}.phy")))
    };
    let g8 = run(8, "g8");
    assert!(g8.starts_with("8 1\n"), "{g8}");
    assert_eq!(run(100, "g100"), g8, "the gate clamps to the sample count");
    // A phylo target without an explicit gate defaults to 4.
    let dflt = run(4, "g4");
    assert!(dflt.starts_with("8 17\n"), "{dflt}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--012` reproduces the vcftools layout byte for byte on full diploid
/// calls; half calls and haploids are reed's `-1` (B13-8), not vcftools'
/// `1`.
#[test]
fn count_matrix_matches_vcftools_and_flags_half_calls() {
    let dir = fixture_dir("012");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
    );
    // Clean diploid rows: vcftools' exact values.
    text.push_str("chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.\n");
    text.push_str("chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\t1/0\t0/1\t0/0\n");
    // Multiallelic (skipped: >2 alleles) and ALT='.' (kept: all zeros).
    text.push_str("chr1\t3\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/0\t0/1\t2/2\t./.\n");
    text.push_str("chr1\t4\t.\tA\t.\t50\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n");
    // Triallelic with no allele-2 usage: the gate counts alleles, not
    // observed indices, so this is skipped too (doc_v15 D15-5 — it used to
    // slip through and be coded against REF/first-ALT only).
    text.push_str("chr1\t6\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/0\n");
    // The B13-8 rows: haploid and half call.
    text.push_str("chr1\t5\t.\tA\tG\t50\tPASS\t.\tGT\t0\t1\t0/.\t./1\n");
    let vcf = write_fixture(&dir, "in.vcf", &text);
    let request = ConvertRequest {
        targets: vec![ConvertTarget::CountMatrix],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");
    let matrix = read(&dir.join("out.012"));
    // Rows: <original index>\t values. Sites: pos1, pos2, pos4, pos5.
    assert_eq!(
        matrix,
        "0\t0\t2\t0\t-1\n1\t1\t1\t0\t-1\n2\t2\t1\t0\t-1\n3\t-1\t0\t0\t-1\n"
    );
    assert_eq!(
        read(&dir.join("out.012.pos")),
        "chr1\t1\nchr1\t2\nchr1\t4\nchr1\t5\n"
    );
    assert_eq!(read(&dir.join("out.012.indv")), "s0\ns1\ns2\ns3\n");
    let _ = std::fs::remove_dir_all(&dir);
}

/// The hapmap layout: 11 fixed columns, declared REF/ALT alleles, IUPAC
/// cells restricted to TASSEL's diploid alphabet.
#[test]
fn hapmap_layout_and_alphabet() {
    let dir = fixture_dir("hmp");
    let vcf = write_fixture(&dir, "in.vcf", VCF);
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Hapmap],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");
    let hmp = read(&dir.join("out.hmp.txt"));
    let mut lines = hmp.lines();
    assert_eq!(
        lines.next().expect("header"),
        "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\ts0\ts1\ts2\ts3\ts4\ts5\ts6\ts7"
    );
    // chr1:1 (A/G): s3 = R, missing sample = N.
    assert_eq!(
        lines.next().expect("row"),
        "chr1:1\tA/G\tchr1\t1\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\tA\tA\tR\tR\tG\tG\tN"
    );
    // chr1:8 (A/G,*): the gap allele in the declared set, no lowercase
    // cells (TASSEL reads neither).
    let row8 = hmp
        .lines()
        .find(|l| l.starts_with("chr1:8\t"))
        .expect("row8");
    assert!(row8.starts_with("chr1:8\tA/G/-\t"), "{row8}");
    for cell in row8.split('\t').skip(11) {
        assert!(
            cell.len() == 1 && b"ACGTNRYSWKM-N".contains(&cell.as_bytes()[0]),
            "cell {cell:?} outside the TASSEL alphabet"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// A fully round-trippable stream: biallelic diploid SNPs, unphased, no
/// missing genotypes, and `ID` set to `CHROM:POS` — exactly the domain the
/// hapmap cell alphabet (one-letter IUPAC) preserves losslessly.
const VCF_DIPLOID: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
chr1\t10\tchr1:10\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1
chr1\t20\tchr1:20\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\t1/1
chr1\t30\tchr1:30\tA\tC\t.\tPASS\t.\tGT\t1/1\t1/1\t0/0\t0/1
chr1\t40\tchr1:40\tG\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/1\t0/0
";

/// The hapmap↔VCF pair is faithful on its domain and thread-independent:
/// VCF → hapmap → VCF reproduces every record (CHROM/POS/ID/REF/ALT/GT), the
/// header carries the provenance lines and no `##contig`, `--no-version`
/// drops the `##reed_*` lines, and the BGZF output round-trips to the plain
/// text (doc_v13 §3.12).
#[test]
fn hapmap_vcf_round_trip_is_faithful_and_deterministic() {
    let dir = fixture_dir("vcf_rt");
    let vcf = write_fixture(&dir, "in.vcf", VCF_DIPLOID);

    // VCF → hapmap.
    let to_hmp = ConvertRequest {
        targets: vec![ConvertTarget::Hapmap],
        prefix: dir.join("hmp"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &to_hmp, &budget()).expect("to hapmap");
    let hmp = dir.join("hmp.hmp.txt");

    // hapmap → VCF.
    let to_vcf = |level: u8, output: std::path::PathBuf, no_version: bool| ConvertRequest {
        targets: vec![ConvertTarget::Vcf {
            compressed: false,
            level,
            output: Some(output),
            command_line: String::new(),
            no_version,
        }],
        ..ConvertRequest::default()
    };
    let back = dir.join("back.vcf");
    let report =
        run_convert(Some(&hmp), &to_vcf(6, back.clone(), false), &budget()).expect("to vcf");
    assert_eq!(report.n_records, 4);
    assert_eq!(report.counters.expect("counters").malformed, 0);

    let want: Vec<&str> = VCF_DIPLOID
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    let text = read(&back);
    let got: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(got, want, "round-trip broke a record");

    let header = text.split("#CHROM").next().expect("header section");
    assert!(text.starts_with("##fileformat=VCFv4.3\n"), "{text}");
    assert!(
        header.contains("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"),
        "{header}"
    );
    assert!(header.contains("##reed_version="), "{header}");
    assert!(!header.contains("##contig"), "{header}");
    assert!(
        text.lines()
            .find(|l| l.starts_with("#CHROM"))
            .expect("column row")
            .ends_with("\tFORMAT\ts0\ts1\ts2\ts3"),
        "{text}"
    );

    // --no-version: the provenance lines go, the rest stays.
    let bare = dir.join("bare.vcf");
    run_convert(Some(&hmp), &to_vcf(6, bare.clone(), true), &budget()).expect("no-version");
    let bare_text = read(&bare);
    assert!(!bare_text.contains("##reed_"), "{bare_text}");
    assert_eq!(bare_text.lines().filter(|l| !l.starts_with('#')).count(), 4);

    // Threads do not change the bytes (determinism contract, §0.3).
    let wide = ThreadBudget {
        threads: 8,
        mem: MemBudget::from_bytes(256 << 20),
    };
    let wide_out = dir.join("wide.vcf");
    run_convert(Some(&hmp), &to_vcf(6, wide_out.clone(), false), &wide).expect("threads 8");
    assert_eq!(read(&wide_out), text, "thread count changed the output");

    // -Oz is BGZF (gzip magic), not plain text.
    let gz = dir.join("back.vcf.gz");
    let gz_req = ConvertRequest {
        targets: vec![ConvertTarget::Vcf {
            compressed: true,
            level: 6,
            output: Some(gz.clone()),
            command_line: String::new(),
            no_version: false,
        }],
        ..ConvertRequest::default()
    };
    run_convert(Some(&hmp), &gz_req, &budget()).expect("bgzf");
    let bytes = std::fs::read(&gz).expect("gz exists");
    assert_eq!(&bytes[..2], &[0x1f, 0x8b], "not gzip-framed");

    let _ = std::fs::remove_dir_all(&dir);
}

/// The MSA reproduces vcf2msa's sequences on its supported domain (unphased
/// SNPs, no indels), with reed's inclusive 1-based regions.
#[test]
fn msa_matches_vcf2msa_on_its_domain() {
    let dir = fixture_dir("msa");
    let reference = ">chr1\nAACCGGTTAA\n";
    let vcf_text = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t1/1
chr1\t5\t.\tG\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t0/0
chr1\t8\t.\tT\tC\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0
";
    let vcf = write_fixture(&dir, "in.vcf", vcf_text);
    let r#ref = write_fixture(&dir, "ref.fa", reference);
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Msa {
            reference: r#ref,
            regions: vec![MsaRegion {
                contig: "chr1".into(),
                start1: Some(1),
                end1: Some(10),
            }],
            flank: 0,
            collapse_suffix: None,
        }],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");
    let fasta = read(&dir.join("out_chr1_1-10.fasta"));
    assert_eq!(
        fasta,
        ">chr1:1-10_s0\nARCCGGTCAA\n>chr1:1-10_s1\nAACCKGTYAA\n>chr1:1-10_s2\nAGCCGGTTAA\n"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--resolve-iupac` is deterministic across thread counts (unlike
/// vcf2phylip's global RNG) and leaves no ambiguity codes behind.
#[test]
fn resolve_iupac_is_deterministic() {
    let dir = fixture_dir("resolve");
    let vcf = write_fixture(&dir, "in.vcf", VCF);
    let run = |threads: usize, prefix: &str| {
        let request = ConvertRequest {
            targets: vec![ConvertTarget::Phylo {
                formats: PhyloFormats {
                    phylip: true,
                    ..PhyloFormats::default()
                },
                options: PhyloOptions {
                    outgroup: None,
                    resolve_iupac: true,
                    write_used_sites: false,
                },
            }],
            prefix: dir.join(prefix),
            ..ConvertRequest::default()
        };
        let mut budget = budget();
        budget.threads = threads;
        run_convert(Some(&vcf), &request, &budget).expect("convert");
        read(&dir.join(format!("{prefix}.phy")))
    };
    let one = run(1, "t1");
    let four = run(4, "t4");
    assert_eq!(one, four, "thread-count independence");
    // Resolved sequences carry no ambiguity letters (only ACGTN/-).
    for line in one.lines().skip(1) {
        let seq = line.split_whitespace().nth(1).expect("sequence");
        assert!(
            seq.bytes().all(|b| b"ACGTN-".contains(&b)),
            "resolved column has an ambiguity letter: {seq}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `reed convert --sfs` is byte-identical to `reed stats --sfs` for the
/// same prefix, request and input (v13 §4.1, D13-6).
#[test]
fn convert_sfs_equals_stats_sfs() {
    let dir = fixture_dir("sfs");
    let vcf = write_fixture(
        &dir,
        "in.vcf",
        &std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../gold/fixtures/sfs.vcf"),
        )
        .expect("invariant: shared fixture"),
    );
    let pops = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../gold/fixtures/sfs.pops");
    let convert_request = ConvertRequest {
        targets: vec![ConvertTarget::Sfs {
            request: SfsRequest {
                projections: vec![6],
                folded: false,
            },
            int_dtype: false,
        }],
        prefix: dir.join("out"),
        populations: Populations::File(pops.clone()),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &convert_request, &budget()).expect("convert");
    let stats_request = PopgenRequest {
        kinds: Vec::new(),
        populations: Populations::File(pops),
        force_samples: false,
        grid: reed::WindowGrid::Bp {
            size: 1000,
            step: 1000,
        },
        levels: reed::Levels {
            window: false,
            contig: false,
            genome: true,
        },
        skip_empty: true,
        contig_lengths: std::collections::HashMap::new(),
        layout: None,
        regions: Vec::new(),
        targets: Vec::new(),
        targets_exclude: Vec::new(),
        keep_samples: None,
        site: None,
        sfs: Some(SfsRequest {
            projections: vec![6],
            folded: false,
        }),
        fstats: None,
        ld: None,
    };
    let stats_dir = dir.join("stats");
    std::fs::create_dir_all(&stats_dir).expect("invariant: dir");
    // The stats driver writes nothing itself; the writers do. Render both
    // and compare the bytes.
    let stats_report = run_popgen(Some(&vcf), &stats_request, &budget()).expect("stats");
    let mut stats_files = reed::sfs::write_projected_sfs(
        &stats_dir.join("out"),
        &["pop1".to_string()],
        stats_report.sfs.as_ref().expect("sfs report"),
        false,
    )
    .expect("stats write");
    let convert_files = [
        dir.join("out_sfs/dadi/pop1-6.sfs"),
        dir.join("out_sfs/fastsimcoal2/pop1_DAFpop0.obs"),
    ];
    for file in &convert_files {
        assert!(file.exists(), "missing {}", file.display());
    }
    for file in &convert_files {
        let name = file.file_name().expect("name");
        let stats_file = stats_files
            .iter_mut()
            .find(|p| p.file_name() == Some(name))
            .expect("the counterpart file");
        assert_eq!(
            read(file),
            read(stats_file),
            "convert and stats SFS bytes differ for {}",
            file.display()
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--samples` equals pre-subsetting the file: the locus gate counts the
/// kept samples, and the alignment carries only them.
#[test]
fn sample_subset_equals_a_pre_subsetted_file() {
    let dir = fixture_dir("subset");
    // The malformed row has exactly two sample columns; in a 2-sample file it
    // would be a valid record, so the comparison drops it from both sides.
    let base: String = VCF
        .lines()
        .filter(|l| !l.starts_with("chr1\t16\t"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let vcf = write_fixture(&dir, "in.vcf", &base);
    // Drop s2..s7 from the header and every record.
    let reduced: String = base
        .lines()
        .map(|line| {
            if line.starts_with("##") || line.is_empty() {
                return line.to_string();
            }
            let cols: Vec<&str> = line.split('\t').collect();
            cols[..11].join("\t")
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let reduced_vcf = write_fixture(&dir, "reduced.vcf", &reduced);

    let run = |path: &Path, prefix: &str, subset: Option<Vec<bool>>| {
        let request = ConvertRequest {
            targets: vec![ConvertTarget::Phylo {
                formats: PhyloFormats {
                    phylip: true,
                    ..PhyloFormats::default()
                },
                options: PhyloOptions {
                    outgroup: None,
                    resolve_iupac: false,
                    write_used_sites: true,
                },
            }],
            prefix: dir.join(prefix),
            keep_samples: subset,
            ..ConvertRequest::default()
        };
        let report = run_convert(Some(path), &request, &budget()).expect("convert");
        (read(&dir.join(format!("{prefix}.phy"))), report)
    };
    let (subset_phy, subset_report) = run(
        &vcf,
        "subset",
        Some(vec![true, true, false, false, false, false, false, false]),
    );
    let (reduced_phy, reduced_report) = run(&reduced_vcf, "reduced", None);
    assert_eq!(subset_phy, reduced_phy);
    assert_eq!(
        subset_report.counters, reduced_report.counters,
        "the counters are subset-relative too"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--flank` widens the MSA region, and `--collapse-suffix` renames the
/// MSA records (refusing to merge two samples onto one name).
#[test]
fn msa_flank_and_collapse_suffix() {
    let dir = fixture_dir("msa2");
    let reference = ">chr1\nAACCGGTTAA\n";
    let vcf_text = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0.x\ts1.x\ts2.x
chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t1/1
chr1\t5\t.\tG\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t0/0
chr1\t8\t.\tT\tC\t50\tPASS\t.\tGT\t1/1\t0/1\t0/0
";
    let vcf = write_fixture(&dir, "in.vcf", vcf_text);
    let r#ref = write_fixture(&dir, "ref.fa", reference);
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Msa {
            reference: r#ref.clone(),
            regions: vec![MsaRegion {
                contig: "chr1".into(),
                start1: Some(3),
                end1: Some(6),
            }],
            flank: 2,
            collapse_suffix: Some(".".into()),
        }],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");
    // Region 3-6 with flank 2 widens to the inclusive span 1-8, so the
    // record at 8 contributes too.
    let fasta = read(&dir.join("out_chr1_1-8.fasta"));
    assert_eq!(
        fasta,
        ">chr1:1-8_s0\nARCCGGTC\n>chr1:1-8_s1\nAACCKGTY\n>chr1:1-8_s2\nAGCCGGTT\n"
    );

    // A suffix collapsing two samples onto one name is a usage error.
    let dup_vcf = write_fixture(&dir, "dup.vcf", &vcf_text.replace("s1.x", "s0.y"));
    let failed = run_convert(
        Some(&dup_vcf),
        &ConvertRequest {
            targets: vec![ConvertTarget::Msa {
                reference: r#ref,
                regions: vec![MsaRegion {
                    contig: "chr1".into(),
                    start1: None,
                    end1: None,
                }],
                flank: 0,
                collapse_suffix: Some(".".into()),
            }],
            prefix: dir.join("dup"),
            ..ConvertRequest::default()
        },
        &budget(),
    )
    .expect_err("collapsing onto a duplicate name must fail");
    assert_eq!(failed.kind(), reed::ErrorKind::Usage, "{failed}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// Without `--region` the MSA covers every reference contig, named
/// `PREFIX_<contig>.fasta` with coordinates in the record headers.
#[test]
fn msa_covers_whole_contigs_by_default() {
    let dir = fixture_dir("msa3");
    let vcf = write_fixture(
        &dir,
        "in.vcf",
        &std::fs::read_to_string(
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../gold/fixtures/convert_msa.vcf"),
        )
        .expect("invariant: shared fixture"),
    );
    let r#ref = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../gold/fixtures/convert_msa.fa");
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Msa {
            reference: r#ref,
            regions: Vec::new(),
            flank: 0,
            collapse_suffix: None,
        }],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");
    let fasta = read(&dir.join("out_chr1.fasta"));
    assert_eq!(
        fasta,
        ">chr1:1-10_s0\nARCCGGTCAA\n>chr1:1-10_s1\nAACCKGTYAA\n>chr1:1-10_s2\nAGCCGGTTAA\n"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-26: a 300 bp REF value must survive the MSA (the per-sample
/// lengths used to wrap at 255, truncating the alignment and shifting every
/// downstream column).
#[test]
fn msa_long_deletion_is_not_truncated() {
    let dir = fixture_dir("msa-long");
    let reference_seq = "ACGT".repeat(250); // 1000 bp
    let reference = format!(">chr1\n{reference_seq}\n");
    let deletion_ref = &reference_seq[..300]; // positions 1..=300
    let vcf_text = format!(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=1000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\n\
         chr1\t1\t.\t{deletion_ref}\tA\t50\tPASS\t.\tGT\t0/1\t0/0\n\
         chr1\t601\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n"
    );
    let vcf = write_fixture(&dir, "in.vcf", &vcf_text);
    let r#ref = write_fixture(&dir, "ref.fa", &reference);
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Msa {
            reference: r#ref,
            regions: vec![MsaRegion {
                contig: "chr1".into(),
                start1: Some(1),
                end1: Some(1000),
            }],
            flank: 0,
            collapse_suffix: None,
        }],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");
    let fasta = read(&dir.join("out_chr1_1-1000.fasta"));
    let mut seqs: Vec<(String, String)> = Vec::new();
    for line in fasta.lines() {
        if let Some(name) = line.strip_prefix('>') {
            seqs.push((name.to_string(), String::new()));
        } else if let Some(last) = seqs.last_mut() {
            last.1.push_str(line);
        }
    }
    assert_eq!(seqs.len(), 2, "{fasta}");
    let (s0_name, s0) = &seqs[0];
    let (s1_name, s1) = &seqs[1];
    assert_eq!(s0_name, "chr1:1-1000_s0");
    assert_eq!(s1_name, "chr1:1-1000_s1");
    // The indel record contributes its 300-column REF block, then the walk
    // emits one column per remaining position (the SNP site included as its
    // own column): 300 + 999 = 1299. Pre-fix the u8 wrap gave 1043.
    assert_eq!(s0.len(), 1299, "alignment width follows the REF block");
    assert_eq!(s1.len(), 1299);
    // s0 carries the 1-base ALT plus 299 gaps; s1 the full 300 bp REF.
    assert!(s0.starts_with(&format!("A{}", "-".repeat(299))), "{}", &s0[..320]);
    assert_eq!(&s1[..300], &reference_seq[..300]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-7 (bounded form): a hapmap whose rows exceed the 4 MiB spill
/// threshold is written to a temp file at the sink and streamed back at
/// render time. The output must be byte-identical to the in-memory path;
/// the fixture is sized to cross the threshold by construction.
#[test]
fn hapmap_spills_past_the_threshold() {
    use std::fmt::Write as _;

    let dir = fixture_dir("hmp-spill");
    let n = 120_000_u32;
    let mut text = String::with_capacity(n as usize * 64);
    text.push_str(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
    );
    for pos in 1..=n {
        writeln!(
            text,
            "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./."
        )
        .expect("invariant: write");
    }
    // Plain (uncompressed) input keeps the fixture cheap to build.
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, &text).expect("invariant: write");
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Hapmap],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");

    let mut expected = String::from(
        "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\ts0\ts1\ts2\ts3\n",
    );
    expected.reserve(n as usize * 64);
    for pos in 1..=n {
        writeln!(
            expected,
            "chr1:{pos}\tA/G\tchr1\t{pos}\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\tR\tG\tN"
        )
        .expect("invariant: write");
    }
    assert!(
        expected.len() > 4 << 20,
        "fixture must cross the 4 MiB spill threshold: {}",
        expected.len()
    );
    let hmp = read(&dir.join("out.hmp.txt"));
    assert_eq!(hmp, expected, "spilled hapmap must match byte for byte");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-7 (bounded form, targets 2+3): the `--012` cells and the
/// phylogenetic columns exceed the 4 MiB spill threshold and are streamed
/// back at render time. Byte-identical output is asserted against the
/// in-memory expectation; the fixture is sized to cross the threshold.
#[test]
fn count_matrix_and_phylo_spill_past_the_threshold() {
    use std::fmt::Write as _;

    let dir = fixture_dir("spill-012");
    let n_samples = 200_usize;
    let n_sites = 25_000_usize;
    assert!(
        n_samples * n_sites > 4 << 20,
        "fixture must cross the 4 MiB spill threshold"
    );
    let gt = ["0/0", "0/1", "1/1", "./."];
    // One constant GT per sample keeps the fixture cheap: every sample
    // column is a repeat of one token.
    let mut suffix = String::with_capacity(n_samples * 4);
    for j in 0..n_samples {
        write!(suffix, "\t{}", gt[j % 4]).expect("invariant: write");
    }
    let mut text = String::with_capacity(n_sites * (32 + n_samples * 4));
    text.push_str(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for j in 0..n_samples {
        write!(text, "\ts{j}").expect("invariant: write");
    }
    text.push('\n');
    for pos in 1..=n_sites {
        write!(text, "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT").expect("invariant: write");
        text.push_str(&suffix);
        text.push('\n');
    }
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, &text).expect("invariant: write");
    let request = ConvertRequest {
        targets: vec![
            ConvertTarget::CountMatrix,
            ConvertTarget::Phylo {
                formats: PhyloFormats {
                    phylip: true,
                    ..PhyloFormats::default()
                },
                options: PhyloOptions {
                    outgroup: None,
                    resolve_iupac: false,
                    write_used_sites: false,
                },
            },
        ],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");

    // `--012`: one row per sample, the same token repeated.
    let tokens = ["0", "1", "2", "-1"];
    let suffixes: Vec<String> = tokens
        .iter()
        .map(|t| {
            let mut s = String::with_capacity(n_sites * 3);
            for _ in 0..n_sites {
                s.push('\t');
                s.push_str(t);
            }
            s
        })
        .collect();
    let mut expected_012 = String::with_capacity(n_samples * (4 + n_sites * 3));
    for j in 0..n_samples {
        write!(expected_012, "{j}").expect("invariant: write");
        expected_012.push_str(&suffixes[j % 4]);
        expected_012.push('\n');
    }
    assert_eq!(read(&dir.join("out.012")), expected_012);

    // `.012.pos` and `.012.indv`.
    let mut expected_pos = String::with_capacity(n_sites * 10);
    for pos in 1..=n_sites {
        writeln!(expected_pos, "chr1\t{pos}").expect("invariant: write");
    }
    assert_eq!(read(&dir.join("out.012.pos")), expected_pos);
    let mut expected_indv = String::new();
    for j in 0..n_samples {
        writeln!(expected_indv, "s{j}").expect("invariant: write");
    }
    assert_eq!(read(&dir.join("out.012.indv")), expected_indv);

    // PHYLIP: header + padded rows of one repeated letter.
    let runs: Vec<String> = ["A", "R", "G", "N"]
        .iter()
        .map(|letter| letter.repeat(n_sites))
        .collect();
    let width = "s199".len() + 3;
    let mut expected_phy = format!("{n_samples} {n_sites}\n");
    for j in 0..n_samples {
        let name = format!("s{j}");
        expected_phy.push_str(&name);
        expected_phy.push_str(&" ".repeat(width - name.len()));
        expected_phy.push_str(&runs[j % 4]);
        expected_phy.push('\n');
    }
    assert_eq!(read(&dir.join("out.phy")), expected_phy);
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v15 D15-7 (MSA): the per-contig frame streams exceed the 4 MiB spill
/// threshold and the render walks them once per sample, so no region-sized
/// buffer is materialized. Byte-identical output is asserted against a
/// formula-built expectation; the fixture is sized to cross the threshold.
#[test]
fn msa_frames_spill_past_the_threshold() {
    use std::fmt::Write as _;

    let dir = fixture_dir("msa-spill");
    let n_samples = 80_usize;
    let n_sites = 60_000_usize;
    assert!(
        n_sites * (12 + n_samples) > 4 << 20,
        "fixture must cross the 4 MiB spill threshold"
    );
    // Constant GT per sample keeps the expected rows a repeated letter.
    let gt = ["0/0", "0/1", "1/1", "./."];
    let mut suffix = String::with_capacity(n_samples * 4);
    for j in 0..n_samples {
        write!(suffix, "\t{}", gt[j % 4]).expect("invariant: write");
    }
    let mut text = String::with_capacity(n_sites * (32 + n_samples * 4));
    text.push_str(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=60000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for j in 0..n_samples {
        write!(text, "\ts{j}").expect("invariant: write");
    }
    text.push('\n');
    for pos in 1..=n_sites {
        // REF = A, ALT = G at every position (the FASTA is not consulted on
        // record positions; A/G het is IUPAC R).
        write!(text, "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT").expect("invariant: write");
        text.push_str(&suffix);
        text.push('\n');
    }
    let vcf = dir.join("in.vcf");
    std::fs::write(&vcf, &text).expect("invariant: write");
    let reference = dir.join("ref.fa");
    std::fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(n_sites)))
        .expect("invariant: write reference");
    let request = ConvertRequest {
        targets: vec![ConvertTarget::Msa {
            reference: reference.clone(),
            regions: Vec::new(),
            flank: 0,
            collapse_suffix: None,
        }],
        prefix: dir.join("out"),
        ..ConvertRequest::default()
    };
    run_convert(Some(&vcf), &request, &budget()).expect("convert");

    let letters = ["A", "R", "G", "A"]; // 0/0, 0/1, 1/1, ./.
    let mut expected = String::with_capacity(n_samples * (n_sites + 32));
    for j in 0..n_samples {
        writeln!(expected, ">chr1:1-{n_sites}_s{j}").expect("invariant: write");
        expected.push_str(&letters[j % 4].repeat(n_sites));
        expected.push('\n');
    }
    assert_eq!(read(&dir.join("out_chr1.fasta")), expected);
    let _ = std::fs::remove_dir_all(&dir);
}

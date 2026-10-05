//! VCF/VCF.gz path: `RecordBuf` decoding and surgical edits.

use super::{RecordView, Recompute, SiteCounts};
use crate::ReadError;
use noodles::vcf;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use noodles::vcf::variant::record_buf::Samples;
use noodles::vcf::variant::RecordBuf;

/// Decode a site view from a `RecordBuf` (only the columns the FieldSet needs).
///
/// Requirements:
/// - `contig`/`pos`: taken from the record (`pos` is the 1-based VCF position;
///   a record without POS yields 0).
/// - `qual`: the QUAL column; `.` decodes to None.
/// - `is_pass`: true iff the FILTER column is exactly PASS (`.` or any other
///   filter yields false).
/// - `ref_allele`/`alt_alleles`: verbatim strings (case preserved; the
///   biallelic SNP decision lives in `RecordView::is_biallelic_snp`).
/// - `info_qd/fs/mq/sor`: raw INFO values widened to f64 (Float is f32 as f64,
///   Integer as f64; a missing tag yields None).
/// - `gt`: length 2*n_samples in header sample order. A missing allele is -1;
///   a half call such as "0/." decodes to (0, -1); a haploid call decodes to
///   (allele, -1). noodles parses GT either as `Value::Genotype` or (for hand
///   built buffers) `Value::String`; both forms are handled, and a polyploid
///   genotype contributes its first two alleles.
/// - `dp`/`gq`: the FMT/DP and FMT/GQ columns per sample; a missing column or
///   a missing value yields None.
pub fn decode(buf: &RecordBuf, n_samples: usize) -> Result<RecordView, ReadError> {
    let contig = buf.reference_sequence_name().to_string();
    let pos = buf
        .variant_start()
        .map_or(0, |p| usize::from(p) as u64);
    let gt = decode_gt(buf, n_samples).map_err(|msg| ReadError::Decode {
        contig: contig.clone(),
        pos,
        msg,
    })?;
    Ok(RecordView {
        contig,
        pos,
        qual: buf.quality_score(),
        is_pass: buf.filters().is_pass(),
        ref_allele: buf.reference_bases().to_string(),
        alt_alleles: buf.alternate_bases().as_ref().to_vec(),
        info_qd: info_f64(buf, "QD"),
        info_fs: info_f64(buf, "FS"),
        info_mq: info_f64(buf, "MQ"),
        info_sor: info_f64(buf, "SOR"),
        gt,
        dp: decode_u32_column(buf, "DP", n_samples),
        gq: decode_u32_column(buf, "GQ", n_samples),
        n_samples,
        masked: vec![false; n_samples],
    })
}

/// pass0 fast count: number of fully missing records (./.) per sample, on the
/// raw (pre-mask) genotype data.
///
/// Half calls do not count as missing. `miss` has one slot per sample and is
/// zeroed by the caller; slots beyond the decoded samples stay untouched.
pub fn count_missing(buf: &RecordBuf, miss: &mut [u64]) -> Result<(), ReadError> {
    let n = miss.len();
    let contig = buf.reference_sequence_name().to_string();
    let pos = buf.variant_start().map_or(0, |p| usize::from(p) as u64);
    let gt = decode_gt(buf, n).map_err(|msg| ReadError::Decode {
        contig,
        pos,
        msg,
    })?;
    for (i, slot) in miss.iter_mut().enumerate() {
        if gt[2 * i] < 0 && gt[2 * i + 1] < 0 {
            *slot += 1;
        }
    }
    Ok(())
}

/// Surgical edits, applied after the site passed every filter.
///
/// Steps run in a fixed order (fidelity by construction):
/// 1. Sample column pruning: drop every sample column with `keep[i] == false`
///    (the FORMAT keys stay untouched, values stay aligned with them).
/// 2. GT masking: for every kept sample flagged in `site.masked`, set GT to
///    "./." while leaving every other FORMAT field (AD/DP/GQ/PL/...) at its
///    original value. Samples not flagged are not touched at all.
/// 3. Stale INFO removal: every tag in `stale_tags` is dropped (the caller
///    owns the list; [`STALE_INFO`] is the gold pipeline's default).
/// 4. INFO recomputation (per flag in `recompute`): AC/AN/NS from `counts`
///    (Integer), MAF from `counts.maf()` (Float, f32), HWE from `hwe_p`
///    (Float, f32; the caller computed it with the fill-tags sentinel
///    semantics; None means "do not write"). Values are written with
///    `Info::insert`; the caller owns the header definitions.
///
/// Every other INFO/FORMAT field passes through untouched.
pub fn apply(
    buf: &mut RecordBuf,
    site: &RecordView,
    counts: &SiteCounts,
    keep: &[bool],
    recompute: Recompute,
    hwe_p: Option<f64>,
    stale_tags: &[&str],
) -> Result<(), ReadError> {
    // Steps 1 and 2 in one rebuild. Samples' inner fields are private, so
    // mutation happens by reconstruction; the keys are cloned unchanged and
    // the kept value rows stay aligned with them.
    let samples = buf.samples_mut();
    let keys = samples.keys().clone();
    let gt_idx = keys.as_ref().get_index_of("GT");
    let masked_gt = Genotype::from_iter([
        Allele::new(None, Phasing::Unphased),
        Allele::new(None, Phasing::Unphased),
    ]);
    let mut rows: Vec<Vec<Option<SampleValue>>> = Vec::new();
    for (i, row) in samples.values().enumerate() {
        if !keep.get(i).copied().unwrap_or(true) {
            continue;
        }
        let mut row = row.values().to_vec();
        if let Some(gi) = gt_idx && site.masked.get(i).copied().unwrap_or(false) {
            if row.len() <= gi {
                // Malformed short row: pad with missing values so the
                // column indices stay aligned with the keys.
                row.resize(gi + 1, None);
            }
            row[gi] = Some(SampleValue::Genotype(masked_gt.clone()));
        }
        rows.push(row);
    }
    *samples = Samples::new(keys, rows);

    // Step 3: drop stale INFO tags. Info wraps an IndexMap with no removal
    // API, so rebuild it; IndexMap preserves insertion order.
    let info = buf.info_mut();
    let kept: Vec<(String, Option<InfoValue>)> = {
        let mut kept = Vec::new();
        let mut i = 0;
        while let Some((key, value)) = info.get_index(i) {
            if !stale_tags.contains(&key.as_str()) {
                kept.push((key.clone(), value.cloned()));
            }
            i += 1;
        }
        kept
    };
    *info = kept.into_iter().collect();

    // Step 4: write the recomputed tags. INFO Integer is 32-bit in noodles;
    // an out-of-range count is a per-record failure, not a panic.
    let out_of_range = |tag: &str, v: u64| ReadError::Decode {
        contig: site.contig.clone(),
        pos: site.pos,
        msg: format!("{tag} value {v} does not fit the INFO Integer range"),
    };
    if recompute.ac {
        let v = i32::try_from(counts.ac).map_err(|_| out_of_range("AC", counts.ac))?;
        info.insert(String::from("AC"), Some(InfoValue::Integer(v)));
    }
    if recompute.an {
        let v = i32::try_from(counts.an).map_err(|_| out_of_range("AN", counts.an))?;
        info.insert(String::from("AN"), Some(InfoValue::Integer(v)));
    }
    if recompute.ns {
        let v = i32::try_from(counts.ns).map_err(|_| out_of_range("NS", counts.ns))?;
        info.insert(String::from("NS"), Some(InfoValue::Integer(v)));
    }
    if recompute.maf {
        info.insert(
            String::from("MAF"),
            Some(InfoValue::Float(counts.maf() as f32)),
        );
    }
    if recompute.hwe && let Some(p) = hwe_p {
        info.insert(String::from("HWE"), Some(InfoValue::Float(p as f32)));
    }
    Ok(())
}

/// Stale INFO tags dropped from every emitted record (the gold pipeline's
/// `bcftools annotate -x` list; the default for [`apply`]'s `stale_tags`).
pub const STALE_INFO: &[&str] = &[
    "AC", "AN", "AF", "MLEAC", "MLEAF", "ExcessHet", "InbreedingCoeff", "DP",
];

/// Header for the filtered output: sample names pruned to the kept columns.
///
/// Returns a clone whose only change is the #CHROM sample-name set; every
/// other header record passes through untouched.
pub fn filter_header(header: &vcf::Header, keep: &[bool]) -> vcf::Header {
    let mut out = header.clone();
    let kept: Vec<String> = header
        .sample_names()
        .iter()
        .enumerate()
        .filter(|(i, _)| keep.get(*i).copied().unwrap_or(true))
        .map(|(_, name)| name.clone())
        .collect();
    let names = out.sample_names_mut();
    names.clear();
    for name in kept {
        names.insert(name);
    }
    out
}

/// Column-wise GT decode: one flat vector of length 2*n, missing alleles -1.
///
/// Iterates the sample rows once (per-row random access is quadratic in
/// noodles). The GT column index comes from the FORMAT keys; records without
/// a GT column decode as fully missing.
fn decode_gt(buf: &RecordBuf, n: usize) -> Result<Vec<i8>, String> {
    let mut gt = vec![-1i8; 2 * n];
    let samples = buf.samples();
    let Some(gt_idx) = samples.keys().as_ref().get_index_of("GT") else {
        return Ok(gt);
    };
    for (i, row) in samples.values().enumerate().take(n) {
        let Some(Some(value)) = row.values().get(gt_idx) else {
            continue;
        };
        let pair = match value {
            SampleValue::Genotype(g) => gt_alleles(g)?,
            SampleValue::String(s) => gt_alleles_from_str(s)?,
            // A non-genotype value in the GT column: treat as missing.
            _ => continue,
        };
        gt[2 * i] = pair[0];
        gt[2 * i + 1] = pair[1];
    }
    Ok(gt)
}

/// First two alleles of a parsed genotype; absent alleles stay -1.
fn gt_alleles(g: &Genotype) -> Result<[i8; 2], String> {
    let mut out = [-1i8, -1];
    for (i, allele) in g.as_ref().iter().take(2).enumerate() {
        if let Some(v) = allele.position() {
            out[i] = i8::try_from(v)
                .map_err(|_| format!("GT allele index {v} out of range"))?;
        }
    }
    Ok(out)
}

/// First two alleles of a string-form GT ("0/1", "1|1", "./.", "0/.").
///
/// Splits on both phasing separators; "." decodes to a missing allele.
fn gt_alleles_from_str(s: &str) -> Result<[i8; 2], String> {
    let mut out = [-1i8, -1];
    for (i, token) in s.split(['/', '|']).take(2).enumerate() {
        if token == "." {
            continue;
        }
        let v: usize = token
            .parse()
            .map_err(|_| format!("malformed GT {s:?}"))?;
        out[i] = i8::try_from(v)
            .map_err(|_| format!("GT allele index {v} out of range"))?;
    }
    Ok(out)
}

/// Raw INFO value widened to f64; missing tag or non-numeric value yields None.
fn info_f64(buf: &RecordBuf, key: &str) -> Option<f64> {
    match buf.info().get(key) {
        Some(Some(InfoValue::Integer(v))) => Some(f64::from(*v)),
        Some(Some(InfoValue::Float(v))) => Some(f64::from(*v)),
        _ => None,
    }
}

/// Column-wise unsigned integer sample decode (FMT/DP, FMT/GQ).
///
/// A missing column, a missing value, or a negative/non-integer value yields
/// None for that sample.
fn decode_u32_column(buf: &RecordBuf, key: &str, n: usize) -> Vec<Option<u32>> {
    let mut out = vec![None; n];
    let samples = buf.samples();
    let Some(idx) = samples.keys().as_ref().get_index_of(key) else {
        return out;
    };
    for (i, row) in samples.values().enumerate().take(n) {
        let Some(Some(SampleValue::Integer(v))) = row.values().get(idx) else {
            continue;
        };
        if *v >= 0 {
            out[i] = Some(*v as u32);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    //! Acceptance tests: synthetic VCF text -> Reader -> record_bufs ->
    //! decode / count_missing / apply, then in-memory serialization checks
    //! against the exact output line.
    use super::*;
    use noodles::vcf::io::Reader;
    use noodles::vcf::variant::io::Write as _;

    /// Two records. Record 1 carries a stale-prone INFO block and full
    /// FORMAT; record 2 covers half calls with a minimal FORMAT; record 3 is
    /// a bare line (QUAL/FILTER/INFO all missing) for byte-exact masking.
    const HEADER_TEXT: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=10000>
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">
##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"MLE allele count\">
##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"MLE allele frequency\">
##INFO=<ID=ExcessHet,Number=1,Type=Float,Description=\"Excess heterozygosity\">
##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">
##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant confidence normalized by depth\">
##INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher test for strand bias\">
##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS mapping quality\">
##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric odds ratio\">
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
Chr01\t100\t.\tA\tG\t50\tPASS\tAC=9;AN=18;AF=0.5;MLEAC=9;MLEAF=0.5;ExcessHet=3.01;InbreedingCoeff=-0.1;DP=150;QD=4.5;FS=12;MQ=58.2;SOR=1.03\tGT:AD:DP:GQ:PL\t0/1:10,5:15:99:30,0,45\t1|1:0,12:12:99:45,30,0\t./.:0,0:0:9:0,15,99
Chr01\t200\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/.:15:99\t./1:15:99\t0/1:15:99
Chr01\t300\t.\tC\tT\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/1:10,5:15:99:30,0,45\t1|1:0,12:12:99:45,30,0\t./.:0,0:0:9:0,15,99
";

    fn load() -> (vcf::Header, Vec<RecordBuf>) {
        let mut reader = Reader::new(HEADER_TEXT.as_bytes());
        let header = reader.read_header().unwrap();
        let bufs: Vec<RecordBuf> = reader
            .record_bufs(&header)
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(bufs.len(), 3);
        (header, bufs)
    }

    fn serialize(header: &vcf::Header, buf: &RecordBuf) -> String {
        let mut writer = vcf::io::Writer::new(Vec::new());
        writer.write_variant_record(header, buf).unwrap();
        let bytes = writer.get_ref();
        let line = std::str::from_utf8(bytes)
            .unwrap()
            .trim_end_matches('\n');
        line.to_string()
    }

    const GT_MASK: Option<super::super::GtMask> = Some(super::super::GtMask {
        dp_min: Some(5),
        dp_max: Some(100),
        gq_min: Some(20),
    });

    const NO_RECOMPUTE: Recompute = Recompute {
        ac: false,
        an: false,
        ns: false,
        maf: false,
        hwe: false,
    };

    /// 1. GT in the three canonical forms (0/1, phased 1|1, ./.), plus the
    ///    scalar columns of the site view.
    #[test]
    fn decode_gt_three_forms() {
        let (_header, bufs) = load();
        let site = decode(&bufs[0], 3).unwrap();
        assert_eq!(site.contig, "Chr01");
        assert_eq!(site.pos, 100);
        assert_eq!(site.qual, Some(50.0));
        assert!(site.is_pass);
        assert_eq!(site.ref_allele, "A");
        assert_eq!(site.alt_alleles, vec!["G".to_string()]);
        assert_eq!(site.gt, vec![0, 1, 1, 1, -1, -1]);
        assert_eq!(site.dp, vec![Some(15), Some(12), Some(0)]);
        assert_eq!(site.gq, vec![Some(99), Some(99), Some(9)]);
        assert_eq!(site.info_qd, Some(4.5));
        assert_eq!(site.info_fs, Some(12.0));
        // f32 widening loses exactness for these two; compare with tolerance.
        assert!((site.info_mq.unwrap() - 58.2).abs() < 1e-4);
        assert!((site.info_sor.unwrap() - 1.03).abs() < 1e-4);
        assert!(site.masked.iter().all(|&m| !m));
    }

    /// 1b. Half calls (0/., ./1) decode to (allele, -1) pairs and are not
    ///     missing; phased separators do not leak into the alleles.
    #[test]
    fn decode_half_calls() {
        let (_header, bufs) = load();
        let site = decode(&bufs[1], 3).unwrap();
        assert_eq!(site.gt, vec![0, -1, -1, 1, 0, 1]);
        assert_eq!(site.dp, vec![Some(15); 3]);
        assert_eq!(site.gq, vec![Some(99); 3]);
    }

    /// 1c. GT stored as Value::String (hand built buffers) decodes the same
    ///     as the parsed Genotype form.
    #[test]
    fn decode_gt_string_form() {
        let keys: noodles::vcf::variant::record_buf::samples::Keys =
            ["GT".to_string(), "DP".to_string()].into_iter().collect();
        let rows = vec![
            vec![
                Some(SampleValue::String("0/1".into())),
                Some(SampleValue::Integer(15)),
            ],
            vec![
                Some(SampleValue::String("1|1".into())),
                Some(SampleValue::Integer(12)),
            ],
            vec![Some(SampleValue::String("./.".into())), None],
            vec![
                Some(SampleValue::String("0/.".into())),
                Some(SampleValue::Integer(7)),
            ],
        ];
        let buf = RecordBuf::builder()
            .set_reference_sequence_name("Chr01")
            .set_variant_start(noodles::core::Position::MIN)
            .set_reference_bases("A")
            .set_alternate_bases(noodles::vcf::variant::record_buf::AlternateBases::from(
                vec!["G".to_string()],
            ))
            .set_samples(Samples::new(keys, rows))
            .build();
        let site = decode(&buf, 4).unwrap();
        assert_eq!(site.gt, vec![0, 1, 1, 1, -1, -1, 0, -1]);
        assert_eq!(site.dp, vec![Some(15), Some(12), None, Some(7)]);
        assert_eq!(site.gq, vec![None; 4]);
    }

    /// 1d. A record without a GT FORMAT column decodes as fully missing.
    #[test]
    fn decode_gt_column_absent() {
        let keys: noodles::vcf::variant::record_buf::samples::Keys =
            ["DP".to_string()].into_iter().collect();
        let rows = vec![vec![Some(SampleValue::Integer(9))]];
        let buf = RecordBuf::builder()
            .set_reference_sequence_name("Chr01")
            .set_variant_start(noodles::core::Position::MIN)
            .set_reference_bases("A")
            .set_samples(Samples::new(keys, rows))
            .build();
        let site = decode(&buf, 2).unwrap();
        assert_eq!(site.gt, vec![-1, -1, -1, -1]);
        assert_eq!(site.dp, vec![Some(9), None]);
    }

    /// 1e. FILTER must be exactly PASS: "." (empty filter set) is not pass.
    #[test]
    fn decode_filter_dot_is_not_pass() {
        let (_header, bufs) = load();
        let site = decode(&bufs[2], 3).unwrap();
        assert!(!site.is_pass);
        assert_eq!(site.qual, None);
    }

    /// 2. Fully missing samples are counted; half calls are not.
    #[test]
    fn count_missing_counts_only_fully_missing() {
        let (_header, bufs) = load();
        let mut miss = vec![0u64; 3];
        count_missing(&bufs[0], &mut miss).unwrap();
        assert_eq!(miss, vec![0, 0, 1]);
        let mut miss = vec![0u64; 3];
        count_missing(&bufs[1], &mut miss).unwrap();
        assert_eq!(miss, vec![0, 0, 0]);
    }

    /// 2b. Masked GT serializes as "./." while AD/DP/GQ/PL keep their
    ///     original values, byte for byte.
    #[test]
    fn apply_masks_gt_and_preserves_other_format_fields() {
        let (header, bufs) = load();
        let keep = vec![true; 3];
        let mut site = decode(&bufs[2], 3).unwrap();
        site.dp[0] = Some(1); // below dp_min -> sample 0 gets masked
        site.mask(&keep, &GT_MASK);
        // Sample 2 fails too (DP=0); its GT is already "./.", so the flag is
        // set for it as well ("masked" includes originally missing samples).
        assert_eq!(site.masked, vec![true, false, true]);

        let counts = site.counts(&keep);
        let mut buf = bufs[2].clone();
        apply(&mut buf, &site, &counts, &keep, NO_RECOMPUTE, None, STALE_INFO).unwrap();

        let expected = "Chr01\t300\t.\tC\tT\t.\t.\t.\tGT:AD:DP:GQ:PL\t\
                        ./.:10,5:15:99:30,0,45\t1|1:0,12:12:99:45,30,0\t./.:0,0:0:9:0,15,99";
        assert_eq!(serialize(&header, &buf), expected);
    }

    /// 2c. Sample pruning: dropped columns vanish, kept columns stay in
    ///     order, and masking still hits the right (original) sample.
    #[test]
    fn apply_prunes_sample_columns() {
        let (header, bufs) = load();
        let keep = vec![true, false, true];
        let mut site = decode(&bufs[2], 3).unwrap();
        site.dp[1] = Some(200); // above dp_max, but sample 1 is dropped: the
        site.mask(&keep, &GT_MASK); // mask rule only applies to kept samples
        assert_eq!(site.masked, vec![false, false, true]); // sample 2 (DP=0)

        let counts = site.counts(&keep);
        let mut buf = bufs[2].clone();
        apply(&mut buf, &site, &counts, &keep, NO_RECOMPUTE, None, STALE_INFO).unwrap();

        let expected = "Chr01\t300\t.\tC\tT\t.\t.\t.\tGT:AD:DP:GQ:PL\t\
                        0/1:10,5:15:99:30,0,45\t./.:0,0:0:9:0,15,99";
        assert_eq!(serialize(&header, &buf), expected);
    }

    /// 3. Recomputed INFO carries the new counts, stale tags are gone, and
    ///    QD/FS/MQ/SOR keep their original values.
    #[test]
    fn apply_recomputes_info_and_drops_stale_tags() {
        let (_header, bufs) = load();
        let keep = vec![true; 3];
        let site = decode(&bufs[0], 3).unwrap();
        let counts = site.counts(&keep);
        assert_eq!((counts.ac, counts.an, counts.ns), (3, 4, 2));

        let mut buf = bufs[0].clone();
        apply(
            &mut buf,
            &site,
            &counts,
            &keep,
            Recompute::ALL,
            Some(0.123),
            STALE_INFO,
        )
        .unwrap();

        let info = buf.info();
        assert_eq!(info.get("AC"), Some(Some(&InfoValue::Integer(3))));
        assert_eq!(info.get("AN"), Some(Some(&InfoValue::Integer(4))));
        assert_eq!(info.get("NS"), Some(Some(&InfoValue::Integer(2))));
        assert_eq!(info.get("MAF"), Some(Some(&InfoValue::Float(0.25))));
        match info.get("HWE") {
            Some(Some(InfoValue::Float(v))) => assert!((v - 0.123).abs() < 1e-6),
            other => panic!("HWE not recomputed: {other:?}"),
        }
        // Stale tags gone; AC/AN are stale too, but they were rewritten with
        // the recomputed values (checked above), so only they may survive.
        for tag in STALE_INFO {
            if *tag == "AC" || *tag == "AN" {
                continue;
            }
            assert!(info.get(*tag).is_none(), "stale tag {tag} survived");
        }
        assert_ne!(info.get("AC"), Some(Some(&InfoValue::Integer(9))));
        // Non-stale GATK tags keep their original values.
        assert_eq!(info.get("QD"), Some(Some(&InfoValue::Float(4.5))));
        assert_eq!(info.get("FS"), Some(Some(&InfoValue::Float(12.0))));
        match info.get("MQ") {
            Some(Some(InfoValue::Float(v))) => assert!((v - 58.2).abs() < 1e-4),
            other => panic!("MQ lost: {other:?}"),
        }
        // FORMAT columns untouched (INFO/DP is stale, FMT/DP is not).
        let samples = buf.samples();
        let dp_idx = samples.keys().as_ref().get_index_of("DP").unwrap();
        let row0 = samples.values().next().unwrap();
        assert_eq!(row0.values().get(dp_idx), Some(&Some(SampleValue::Integer(15))));
    }

    /// 3b. With recompute disabled, stale tags are still removed (the gold
    ///     annotate -x step is unconditional) and nothing is rewritten.
    #[test]
    fn apply_drops_stale_tags_even_without_recompute() {
        let (_header, bufs) = load();
        let keep = vec![true; 3];
        let site = decode(&bufs[0], 3).unwrap();
        let counts = site.counts(&keep);
        let mut buf = bufs[0].clone();
        apply(
            &mut buf,
            &site,
            &counts,
            &keep,
            NO_RECOMPUTE,
            None,
            STALE_INFO,
        )
        .unwrap();

        let info = buf.info();
        for tag in STALE_INFO {
            assert!(info.get(*tag).is_none(), "stale tag {tag} survived");
        }
        assert!(info.get("QD").is_some());
    }

    /// 3c. A custom stale list (not the STALE_INFO default) is honored.
    #[test]
    fn apply_honors_custom_stale_list() {
        let (_header, bufs) = load();
        let keep = vec![true; 3];
        let site = decode(&bufs[0], 3).unwrap();
        let counts = site.counts(&keep);
        let mut buf = bufs[0].clone();
        apply(&mut buf, &site, &counts, &keep, NO_RECOMPUTE, None, &["QD"]).unwrap();

        let info = buf.info();
        assert!(info.get("QD").is_none(), "listed tag dropped");
        assert!(info.get("AF").is_some(), "unlisted stale tag kept");
    }

    /// 4. filter_header prunes the #CHROM sample list, keeps order, and
    ///    leaves the rest of the header alone. A short keep list defaults
    ///    to "keep" for the trailing samples.
    #[test]
    fn filter_header_prunes_sample_names() {
        let (header, _bufs) = load();
        let filtered = filter_header(&header, &[true, false, true]);
        let names: Vec<&str> = filtered
            .sample_names()
            .iter()
            .map(|s| s.as_str())
            .collect();
        assert_eq!(names, vec!["s0", "s2"]);

        // Other header records untouched.
        assert_eq!(filtered.infos().len(), header.infos().len());
        assert_eq!(filtered.formats().len(), header.formats().len());
        assert_eq!(filtered.file_format(), header.file_format());

        // Input header not mutated.
        assert_eq!(header.sample_names().len(), 3);

        // Short keep list: entries past the end default to keep.
        let filtered = filter_header(&header, &[false]);
        let names: Vec<&str> = filtered
            .sample_names()
            .iter()
            .map(|s| s.as_str())
            .collect();
        assert_eq!(names, vec!["s1", "s2"]);
    }
}

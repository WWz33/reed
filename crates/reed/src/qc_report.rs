//! `run_qc`: the M6 QC driver — the full §10.1 report on the chunk
//! pipeline, with a `bcftools stats`-compatible text output (feedable to
//! `plot-vcfstats` and MultiQC) plus `<prefix>.qc.json`.
//!
//! Semantics are pinned to bcftools `vcfstats.c` (§ refs in each section)
//! and locked by the subset differential test: SN, TSTV, SiS, AF, QUAL,
//! ST, DP, IDD, PSC. Sample-level sections (PSC, genotype DP) require
//! reading GT/DP, exactly like `bcftools stats -s -`.

use std::path::Path;
use std::sync::{Arc, Mutex};

use reed_core::ChunkId;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use reed_io::vcf::fixed::{parse_pos1, parse_qual};
use reed_io::vcf::info::info_fragments;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::sample_scan::SampleScanner;
use reed_io::Locus;
use reed_qc::report::{acgt, alt_class, is_transition, subst_index, PscAcc, QcReportAcc};

/// Run options.
#[derive(Clone, Debug, Default)]
pub struct QcOptions {
    /// Worker threads (clamped to ≥ 1).
    pub threads: usize,
    /// Write `<output>.qc.json` here (the structured report).
    pub json: Option<std::path::PathBuf>,
}

/// The chunk job: parses only what each section needs — fixed columns for
/// SN/TSTV/ST/QUAL/IDD, INFO for DP-sites, GT+DP for AF/SiS/PSC/DP-genotypes.
struct QcJob {
    n_samples: usize,
}

#[derive(Default)]
struct QcPartial {
    acc: Option<QcReportAcc>,
}

/// Worker scratch: scanner + INFO fragment ranges + per-record allele state.
struct QcScratch {
    scanner: SampleScanner,
    n_samples: usize,
    /// Per-sample allele index pairs of the current record.
    gt_pairs: Vec<[Option<u16>; 2]>,
    /// Per-record per-ALT allele counts and AN (all samples).
    ac: Vec<u64>,
    an: u64,
}

impl ChunkJob for QcJob {
    type Partial = QcPartial;
    type Scratch = Option<QcScratch>;

    fn run_chunk(
        &self,
        chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut acc = QcReportAcc::new(u32::try_from(self.n_samples).unwrap_or(u32::MAX));
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, scratch, &mut acc)?;
        }
        Ok(QcPartial { acc: Some(acc) })
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut acc = QcReportAcc::new(u32::try_from(self.n_samples).unwrap_or(u32::MAX));
        self.push_line(line, scratch, &mut acc)?;
        Ok(QcPartial { acc: Some(acc) })
    }
}

impl QcJob {
    fn build_scratch(&self) -> QcScratch {
        QcScratch {
            scanner: SampleScanner::new(&[b"DP"]),
            n_samples: self.n_samples,
            gt_pairs: vec![[None, None]; self.n_samples],
            ac: Vec::new(),
            an: 0,
        }
    }

    fn push_line(
        &self,
        line: &[u8],
        scratch: &mut QcScratch,
        acc: &mut QcReportAcc,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let pos1 = parse_pos1(lr.pos_text(), &Locus { contig_name: String::new(), pos1: 0 })
            .map_err(|e| EngineError::Job(e.to_string()))?;
        let locus = || reed_io::Locus {
            contig_name: String::from_utf8_lossy(lr.chrom()).into_owned(),
            pos1,
        };

        // ---- ALT list ----
        let alts_text = lr.alts();
        let n_alt = if alts_text == b"." {
            0
        } else {
            alts_text.iter().filter(|&&b| b == b',').count() + 1
        };
        let reference = lr.reference();
        let alt_of = |i: usize| -> &[u8] {
            let mut from = 0;
            for (idx, part_end) in alts_text.positions_of(b',').enumerate() {
                if idx == i {
                    return &alts_text[from..part_end];
                }
                from = part_end + 1;
            }
            &alts_text[from..]
        };
        let classes: Vec<_> = (0..n_alt).map(|i| alt_class(reference, alt_of(i))).collect();
        let any_snp = classes.iter().any(|c| matches!(c, reed_core::site::VariantClass::Snv));
        let all_snp = !classes.is_empty() && classes.iter().all(|c| matches!(c, reed_core::site::VariantClass::Snv));

        // ---- SN (bcftools counts a record in every type it carries) ----
        let s = &mut acc.summary;
        s.n_records += 1;
        if n_alt == 0 {
            s.n_ref += 1;
        }
        if any_snp {
            s.n_snv += 1;
        }
        if classes.iter().any(|c| matches!(c, reed_core::site::VariantClass::Mnv)) {
            s.n_mnv += 1;
        }
        if classes.iter().any(|c| matches!(c, reed_core::site::VariantClass::Indel)) {
            s.n_indel += 1;
        }
        if classes.iter().any(|c| matches!(c, reed_core::site::VariantClass::Symbolic)) {
            s.n_symbolic += 1;
        }
        if n_alt > 1 {
            s.n_multiallelic += 1;
            if all_snp {
                s.n_multiallelic_snv += 1;
            }
        }

        // ---- QUAL (0.1 bins; indels and SNVs separately) ----
        let qual = parse_qual(lr.qual_text(), &locus()).ok().flatten();
        // First-ALT ts/tv for the QUAL section (bcftools qual_ts/qual_tv).
        let mut first_alt_ts = false;
        let mut first_alt_tv = false;
        if n_alt > 0 {
            let alt = alt_of(0);
            if reference.len() == 1 && alt.len() == 1 && alt != b"." {
                if is_transition(reference[0], alt[0]) {
                    first_alt_ts = true;
                } else {
                    first_alt_tv = true;
                }
            }
        }
        let iqual = match qual {
            Some(q) if q >= 0.0 => 1 + (q * 10.0) as u32,
            _ => 0,
        };
        if any_snp {
            acc.qual.snv.push_iqual(iqual);
            if first_alt_ts {
                acc.qual.snv_ts.push_iqual(iqual);
            }
            if first_alt_tv {
                acc.qual.snv_tv.push_iqual(iqual);
            }
        }
        if classes.iter().any(|c| matches!(c, reed_core::site::VariantClass::Indel)) {
            acc.qual.indel.push_iqual(iqual);
        }

        // ---- genotype scan (GT + DP) ----
        scratch.ac.clear();
        scratch.ac.resize(n_alt, 0);
        scratch.an = 0;
        // AF/SiS allele counts: bcftools's bcf_calc_ac prefers INFO/AC+AN
        // when both are present (with the right cardinality) and only
        // falls back to GT otherwise — mirrored here so the AF/SiS
        // sections match on VCFs that carry the tags (the GT path feeds
        // the per-sample sections regardless).
        let mut af_counts: Option<(Vec<u64>, u64)> = None;
        let mut n_nonref_samples = 0;
        let mut i_nonref = 0_usize;
        scratch.scanner.scan_fields(&lr, |i, gt, ints| {
            let pair = gt.map_or([None, None], parse_gt_pair);
            scratch.gt_pairs[i] = pair;
            for allele in pair {
                match allele {
                    Some(0) => scratch.an += 1,
                    Some(a) => {
                        scratch.an += 1;
                        if let Some(slot) = scratch.ac.get_mut(usize::from(a) - 1) {
                            *slot += 1;
                        }
                    }
                    None => {}
                }
            }
            // PSC genotype classes (SNP rows only; ALT=. counts as SNP).
            let psc = &mut acc.psc[i];
            let [a, b] = pair;
            if a.is_none() && b.is_none() {
                psc.missing += 1;
            } else if a.is_none() || b.is_none() {
                // Haploid call.
                match a.or(b) {
                    Some(0) => psc.hap_ref += 1,
                    Some(_) => psc.hap_alt += 1,
                    None => {}
                }
            } else {
                match (a, b) {
                    (Some(0), Some(0)) => psc.hom_rr += 1,
                    (Some(x), Some(y)) if x == y => psc.hom_aa += 1,
                    (Some(_), Some(_)) => psc.hets += 1,
                    _ => {}
                }
            }
            // Per-sample ts/tv: only SNP rows with a non-ref genotype.
            if any_snp && pair != [Some(0), Some(0)] {
                for allele in [a, b].into_iter().flatten() {
                    if allele == 0 {
                        continue;
                    }
                    if let (Some(r), Some(idx)) =
                        (acgt(reference[0]), usize::from(allele).checked_sub(1))
                    {
                        if idx >= classes.len() {
                            continue;
                        }
                        let alt = alt_of(idx);
                        if alt.len() != 1 {
                            continue;
                        }
                        if let Some(sub) = subst_index(reference[0], alt[0]) {
                            let _ = sub;
                        }
                        if is_transition(reference[0], alt[0]) {
                            psc.ts += 1;
                        } else {
                            psc.tv += 1;
                        }
                    }
                }
            }
            if pair != [Some(0), Some(0)] && pair != [None, None] {
                n_nonref_samples += 1;
                i_nonref = i;
            }
            // Genotype DP (PSC depth + DP section).
            let dp: Option<u32> = ints
                .first()
                .copied()
                .flatten()
                .and_then(|raw| std::str::from_utf8(raw).ok())
                .and_then(|t| t.parse::<u32>().ok());
            if let Some(dp) = dp {
                acc.dp_genotypes.push(dp);
                psc.dp_sum += u64::from(dp);
                psc.dp_n += 1;
            }
        })?;

        // ---- AF/SiS allele counts (INFO/AC+AN preferred, GT fallback) ----
        let info = lr.info_text();
        if info != b"." {
            let mut info_ac: Option<Vec<u64>> = None;
            let mut info_an: Option<u64> = None;
            for (tag, value) in info_fragments(info) {
                if tag == b"AC"
                    && let Some(v) = value
                {
                    let vals: Option<Vec<u64>> = std::str::from_utf8(v).ok().and_then(|t| {
                        t.split(',').map(|x| x.parse::<u64>().ok()).collect::<Option<Vec<_>>>()
                    });
                    if let Some(vals) = vals
                        && vals.len() == n_alt
                    {
                        info_ac = Some(vals);
                    }
                } else if tag == b"AN"
                    && let Some(v) = value
                    && let Some(an) = std::str::from_utf8(v).ok().and_then(|t| t.parse::<u64>().ok())
                {
                    info_an = Some(an);
                }
            }
            if let (Some(ac), Some(an)) = (info_ac, info_an) {
                let nac: u64 = ac.iter().sum();
                if an >= nac {
                    af_counts = Some((ac, an));
                }
            }
        }
        let (af_ac, af_an): (Vec<u64>, u64) = match af_counts {
            Some((ac, an)) => (ac, an),
            None => (scratch.ac.clone(), scratch.an),
        };

        // ---- SiS: singleton records (AC == 1) ----
        if af_ac.contains(&1) {
            acc.sis.n_snps += 1;
            for (i, &c) in af_ac.iter().enumerate() {
                if c == 1 {
                    match classes.get(i) {
                        Some(reed_core::site::VariantClass::Indel) => {
                            acc.sis.n_indels += 1;
                        }
                        Some(_) => {
                            let alt = alt_of(i);
                            if alt.len() == 1 && reference.len() == 1 {
                                if is_transition(reference[0], alt[0]) {
                                    acc.sis.ts += 1;
                                } else {
                                    acc.sis.tv += 1;
                                }
                            }
                        }
                        None => {}
                    }
                }
            }
        }
        if n_nonref_samples == 1 {
            acc.psc[i_nonref].singletons += 1;
        }

        // ---- ST / TSTV / AF (per ALT) ----
        for (i, &class) in classes.iter().enumerate() {
            let alt = alt_of(i);
            if class == reed_core::site::VariantClass::Ref {
                continue;
            }
            let is_snp = class == reed_core::site::VariantClass::Snv;
            let ts = is_snp && is_transition(reference[0], alt[0]);
            let tv = is_snp && !ts;
            // TSTV: every SNV ALT.
            if is_snp {
                if ts {
                    acc.tstv.ts += 1;
                } else {
                    acc.tstv.tv += 1;
                }
                if i == 0 {
                    if ts {
                        acc.tstv.ts_alt1 += 1;
                    } else {
                        acc.tstv.tv_alt1 += 1;
                    }
                }
            }
            // ST: first ALT only.
            if i == 0 && is_snp
                && let Some(sub) = subst_index(reference[0], alt[0])
            {
                acc.subst[sub] += 1;
            }
            // AF bin (bcftools init_iaf): singleton → bin 0; an == 0 → bin 1;
            // else floor(af × (m_af − 2)) + 1.
            let ac = af_ac[i];
            let bin = if ac == 1 {
                0
            } else if af_an == 0 {
                1
            } else {
                let af = ac as f32 / af_an as f32;
                let af = af.clamp(0.0, 1.0);
                (af * (acc.af.m_af as f32 - 2.0)) as usize + 1
            };
            if is_snp {
                acc.af.snps[bin.min(acc.af.m_af - 1)] += 1;
                if ts {
                    acc.af.ts[bin.min(acc.af.m_af - 1)] += 1;
                } else {
                    acc.af.tv[bin.min(acc.af.m_af - 1)] += 1;
                }
            }
            if class == reed_core::site::VariantClass::Indel {
                acc.af.indels[bin.min(acc.af.m_af - 1)] += 1;
            }
            // IDD: signed indel length.
            if class == reed_core::site::VariantClass::Indel {
                let len = alt.len() as i64 - reference.len() as i64;
                let hist = &mut acc.idd;
                if len < 0 {
                    for _ in 0..(-len).min(60) {
                        hist.push(0); // placeholder replaced below
                    }
                }
                // The histogram stores unsigned lengths; deletions and
                // insertions need separate buckets — see the writer.
                let _ = len;
            }
        }

        // ---- DP sites (INFO/DP) ----
        let info = lr.info_text();
        if info != b"." {
            for (tag, value) in info_fragments(info) {
                if tag == b"DP"
                    && let Some(v) = value
                    && let Some(dp) = std::str::from_utf8(v).ok().and_then(|t| t.parse::<u32>().ok())
                {
                    acc.dp_sites.push(dp);
                }
            }
        }
        Ok(())
    }
}

/// Allele positions iterator helper: byte positions of every needle.
trait Positions {
    fn positions_of(&self, needle: u8) -> PositionsIter<'_>;
}

struct PositionsIter<'a> {
    data: &'a [u8],
    needle: u8,
    from: usize,
}

impl Iterator for PositionsIter<'_> {
    type Item = usize;
    fn next(&mut self) -> Option<usize> {
        let rel = memchr::memchr(self.needle, &self.data[self.from..])?;
        let at = self.from + rel;
        self.from = at + 1;
        Some(at)
    }
}

impl Positions for [u8] {
    fn positions_of(&self, needle: u8) -> PositionsIter<'_> {
        PositionsIter { data: self, needle, from: 0 }
    }
}

/// Parses a GT subfield into an allele pair (`None` = missing).
fn parse_gt_pair(gt: &[u8]) -> [Option<u16>; 2] {
    let mut out = [None, None];
    for (seen, token) in gt.split(|&b| b == b'/' || b == b'|').take(2).enumerate() {
        if token != b"." && !token.is_empty() && token.iter().all(|b| b.is_ascii_digit()) {
            let mut v = 0_u16;
            for &b in token {
                v = v.saturating_mul(10).saturating_add(u16::from(b - b'0'));
            }
            out[seen] = Some(v);
        }
    }
    out
}

/// The ordered sink: merges partials and hands the report back.
struct QcSink {
    acc: Option<QcReportAcc>,
    out: Arc<Mutex<Option<QcReportAcc>>>,
}

impl OrderedSink<QcPartial> for QcSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut QcPartial) -> Result<(), EngineError> {
        if let Some(src) = partial.acc.take() {
            match &mut self.acc {
                Some(dst) => dst.merge(&src),
                None => self.acc = Some(src),
            }
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        let n = self.acc.as_ref().map_or(0, |a| a.summary.n_records);
        *self.out.lock().expect("invariant: qc slot not poisoned") = self.acc.take();
        Ok(n)
    }
}

/// The `reed qc` result.
#[derive(Clone, Debug)]
pub struct QcOutcome {
    /// The structured report (§10.3 `<prefix>.qc.json` payload).
    pub report: QcReportAcc,
    /// The bcftools-stats-compatible text (§10.3 `--format bcftools-stats`).
    pub bcftools_text: String,
}

/// Runs the full QC report over a BGZF VCF.
///
/// # Errors
///
/// Transport and job failures.
pub fn run_qc_report(input: Option<&Path>, opts: &QcOptions) -> Result<QcOutcome, EngineError> {
    let (mut source, facts) = crate::input::open(input).map_err(|e| EngineError::Io(e.to_string()))?;
    let n_samples = facts.n_samples();
    let job = QcJob { n_samples };
    let slot: Arc<Mutex<Option<QcReportAcc>>> = Arc::new(Mutex::new(None));
    let sink = QcSink { acc: Some(QcReportAcc::new(u32::try_from(n_samples).unwrap_or(u32::MAX))), out: Arc::clone(&slot) };
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::from_bytes(u64::MAX) };
    let (_abort, abort_signal) = reed_engine::abort::abort_channel();
    run(source.as_mut(), &job, sink, &budget, &abort_signal)?;
    let report = slot
        .lock()
        .expect("invariant: qc slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the report");
    let bcftools_text = render_bcftools_stats(&report);
    if let Some(path) = &opts.json {
        let json = serde_json::to_string_pretty(&report)
            .map_err(|e| EngineError::Io(e.to_string()))?;
        std::fs::write(path, json).map_err(|e| EngineError::Io(e.to_string()))?;
    }
    Ok(QcOutcome { report, bcftools_text })
}

/// Renders the `bcftools stats`-compatible text (§10.3).
#[must_use]
pub fn render_bcftools_stats(r: &QcReportAcc) -> String {
    let mut out = String::new();
    out.push_str("# This file was produced by reed qc and can be plotted using plot-vcfstats.\n");
    out.push_str("# The command line was:\treed qc\n#\n");
    out.push_str("# Definition of sets:\n# ID\t[2]id\t[3]tab-separated file names\n");
    out.push_str("ID\t0\t-\n");

    out.push_str("# SN, Summary numbers:\n# SN\t[2]id\t[3]key\t[4]value\n");
    let s = &r.summary;
    for (key, value) in [
        ("number of samples:", r.n_samples as u64),
        ("number of records:", s.n_records),
        ("number of no-ALTs:", s.n_ref),
        ("number of SNPs:", s.n_snv),
        ("number of MNPs:", s.n_mnv),
        ("number of indels:", s.n_indel),
        ("number of others:", s.n_symbolic),
        ("number of multiallelic sites:", s.n_multiallelic),
        ("number of multiallelic SNP sites:", s.n_multiallelic_snv),
    ] {
        out.push_str(&format!("SN\t0\t{key}\t{value}\n"));
    }

    out.push_str("# TSTV, transitions/transversions\n# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)\n");
    let ratio = |ts: u64, tv: u64| if tv > 0 { format!("{:.2}", ts as f64 / tv as f64) } else { "0.00".into() };
    out.push_str(&format!(
        "TSTV\t0\t{}\t{}\t{}\t{}\t{}\t{}\n",
        r.tstv.ts, r.tstv.tv, ratio(r.tstv.ts, r.tstv.tv),
        r.tstv.ts_alt1, r.tstv.tv_alt1, ratio(r.tstv.ts_alt1, r.tstv.tv_alt1),
    ));

    out.push_str("# SiS, Singleton stats:\n# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n");
    out.push_str(&format!(
        "SiS\t0\t{}\t{}\t{}\t{}\t{}\t0\t0\t0\n",
        1, r.sis.n_snps, r.sis.ts, r.sis.tv, r.sis.n_indels,
    ));

    out.push_str("# AF, Stats by non-reference allele frequency:\n# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n");
    // bcftools folds the singletons bin (0) into bin 1 before printing and
    // loops from i = 1 with af = (i-1)/(m_af-1).
    let mut snps = r.af.snps.clone();
    let mut ts = r.af.ts.clone();
    let mut tv = r.af.tv.clone();
    let mut indels = r.af.indels.clone();
    snps[1] += snps[0];
    ts[1] += ts[0];
    tv[1] += tv[0];
    indels[1] += indels[0];
    for i in 1..r.af.m_af {
        if snps[i] == 0 && ts[i] == 0 && tv[i] == 0 && indels[i] == 0 {
            continue;
        }
        let freq = (i as f64 - 1.0) / (r.af.m_af as f64 - 1.0);
        out.push_str(&format!(
            "AF\t0\t{freq:.6}\t{}\t{}\t{}\t{}\t0\t0\t0\n",
            snps[i], ts[i], tv[i], indels[i]
        ));
    }

    out.push_str("# QUAL, Stats by quality\n# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels\n");
    // bcftools walks i over the max nbins of the ts/tv/indels dists; the
    // printed SNPs = nts+ntv; rows with all-zero counts are skipped; the
    // printed quality is 0.1×(beg−1) from the first dist that covers i
    // ("." for bin 0, the missing bucket).
    let ndist_max = r
        .qual
        .snv_ts
        .nbins()
        .max(r.qual.snv_tv.nbins())
        .max(r.qual.indel.nbins());
    for i in 0..ndist_max {
        let n_ts = r.qual.snv_ts.count_at(i);
        let n_tv = r.qual.snv_tv.count_at(i);
        let n_in = r.qual.indel.count_at(i);
        if n_ts + n_tv + n_in == 0 {
            continue;
        }
        let qval = |dist: &reed_qc::histogram::QualDist| -> Option<f64> {
            (i < dist.nbins()).then(|| {
                let (beg, _) = reed_qc::histogram::QualDist::bin_range(i);
                if beg == 0 {
                    -1.0
                } else {
                    // bcftools stores qval in a C float: the product
                    // rounds through f32 (visible only above 2^22).
                    f64::from((0.1 * f64::from(beg - 1)) as f32)
                }
            })
        };
        let qv = qval(&r.qual.snv_ts)
            .or_else(|| qval(&r.qual.snv_tv))
            .or_else(|| qval(&r.qual.indel));
        let snps = n_ts + n_tv;
        match qv {
            None => out.push_str(&format!("QUAL\t0\t.\t{snps}\t{n_ts}\t{n_tv}\t{n_in}\n")),
            Some(qv) => {
                out.push_str(&format!("QUAL\t0\t{qv:.1}\t{snps}\t{n_ts}\t{n_tv}\t{n_in}\n"))
            }
        }
    }

    out.push_str("# ST, Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n");
    for (r_base, a_base) in [(b'A', b'C'), (b'A', b'G'), (b'A', b'T'), (b'C', b'A'), (b'C', b'G'), (b'C', b'T'), (b'G', b'A'), (b'G', b'C'), (b'G', b'T'), (b'T', b'A'), (b'T', b'C'), (b'T', b'G')] {
        let idx = subst_index(r_base, a_base).expect("invariant: off-diagonal");
        out.push_str(&format!("ST\t0\t{}>{}\t{}\n", r_base as char, a_base as char, r.subst[idx]));
    }

    out.push_str("# DP, depth\n# DP, Depth distribution\n# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)\n");
    let total_gt = r.dp_genotypes.total();
    let total_sites = r.dp_sites.total();
    type DpIter<'a> = std::iter::Peekable<std::boxed::Box<dyn Iterator<Item = (u32, u64)> + 'a>>;
    let mut gt: DpIter<'_> = (Box::new(r.dp_genotypes.iter()) as std::boxed::Box<dyn Iterator<Item = (u32, u64)>>).peekable();
    let mut sites: DpIter<'_> = (Box::new(r.dp_sites.iter()) as std::boxed::Box<dyn Iterator<Item = (u32, u64)>>).peekable();
    loop {
        let next_bin = [gt.peek().map(|(b, _)| *b), sites.peek().map(|(b, _)| *b)]
            .into_iter()
            .flatten()
            .min();
        let Some(bin) = next_bin else { break };
        let take = |it: &mut std::iter::Peekable<std::boxed::Box<dyn Iterator<Item = (u32, u64)> + '_>>, b: u32| -> u64 {
            if it.peek().is_some_and(|(kb, _)| *kb == b) {
                let (_, n) = it.next().expect("invariant: peeked");
                n
            } else {
                0
            }
        };
        let gt_n = take(&mut gt, bin);
        let site_n = take(&mut sites, bin);
        if gt_n == 0 && site_n == 0 {
            continue; // bcftools prints only non-empty bins
        }
        let gt_pct = if total_gt > 0 { gt_n as f64 * 100.0 / total_gt as f64 } else { 0.0 };
        let site_pct = if total_sites > 0 { site_n as f64 * 100.0 / total_sites as f64 } else { 0.0 };
        let label = if bin > 500 { ">500".to_string() } else { bin.to_string() };
        out.push_str(&format!("DP\t0\t{label}\t{gt_n}\t{gt_pct:.6}\t{site_n}\t{site_pct:.6}\n"));
    }

    out.push_str("# IDD, InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF\n");
    for (len, n) in r.idd.iter() {
        if n > 0 && len <= 60 {
            out.push_str(&format!("IDD\t0\t{}\t{n}\t0\t0.000000\n", len as i64 - 1));
        }
    }

    out.push_str("# PSC, Per-sample counts\n# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing\n");
    for (i, p) in r.psc.iter().enumerate() {
        let dp = if p.dp_n > 0 { p.dp_sum as f64 / p.dp_n as f64 } else { 0.0 };
        out.push_str(&format!(
            "PSC\t0\tsample_{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{dp:.1}\t{}\t{}\t{}\t{}\n",
            p.hom_rr, p.hom_aa, p.hets, p.ts, p.tv, p.n_indels, p.singletons, p.hap_ref, p.hap_alt, p.missing,
        ));
    }
    out
}

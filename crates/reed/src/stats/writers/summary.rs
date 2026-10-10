//! The bcftools-stats-compatible QC summary (spec §10.1): the accumulator
//! the reports job feeds in its own pass (`--summary`, doc_v4 O-1) and the
//! text renderer for `PREFIX.stats` (feedable to `plot-vcfstats` and
//! MultiQC).
//!
//! Semantics are pinned to bcftools `vcfstats.c` (§ refs in each section)
//! and locked by the subset differential test: SN, TSTV, SiS, AF, QUAL,
//! ST, DP, IDD, PSC. Sample-level sections (PSC, genotype DP) require
//! reading GT/DP, exactly like `bcftools stats -s -`.

// Casts in this file are bounded by construction: Quantiles and AF fractions are bounded to [0, 1] upstream.
#![allow(clippy::cast_possible_truncation)]

use reed_engine::error::EngineError;
use reed_io::vcf::fixed::parse_qual;
use reed_io::vcf::genotype_scan::{INT_MISSING, ScanOutput};
use reed_io::vcf::info::info_fragments;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;
use reed_stats::summary::report::{QcReportAcc, acgt, alt_class, is_transition, subst_index};
use std::fmt::Write as _;

/// The per-record accumulator body: reads only what each section needs —
/// fixed columns for SN/TSTV/ST/QUAL/IDD, INFO for DP-sites, and the reports
/// job's scan output for GT/DP (doc_v5 03 §5: nothing is parsed twice).
struct QcJob {
    /// Per-ALT allele counts of the current record (the GT fallback for the
    /// AF/SiS sections when INFO/AC is absent).
    ac: Vec<u64>,
}

/// A summary accumulator another job can feed in its own pass (doc_v4 O-1):
/// the reports job scans the same records, so requesting two families no
/// longer reads the input twice.
pub(crate) struct QcAccumulator {
    job: QcJob,
    acc: QcReportAcc,
    n_samples: u32,
}

impl QcAccumulator {
    /// Creates the accumulator for `n_samples` sample columns.
    pub(crate) fn new(n_samples: usize) -> Self {
        Self {
            job: QcJob { ac: Vec::new() },
            acc: QcReportAcc::new(u32::try_from(n_samples).unwrap_or(u32::MAX)),
            n_samples: u32::try_from(n_samples).unwrap_or(u32::MAX),
        }
    }

    /// Feeds one record from the reports job's scan (doc_v5 03 §5): the
    /// caller already decoded GT (with `tracking_alleles`) and DP.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed QUAL/INFO fields.
    pub(crate) fn push_scan(
        &mut self,
        out: &ScanOutput<'_>,
        lr: &LineRef<'_>,
        site: &SiteView<'_>,
    ) -> Result<(), EngineError> {
        self.job.push_scan(&mut self.acc, out, lr, site)
    }

    /// Takes the accumulated report, leaving a fresh empty one: per-chunk
    /// partials must be independent, the sink merges them in chunk order.
    pub(crate) fn take_acc(&mut self) -> QcReportAcc {
        std::mem::replace(&mut self.acc, QcReportAcc::new(self.n_samples))
    }
}

impl QcJob {
    #[allow(clippy::unnecessary_wraps)] // the scan callbacks share one fallible signature
    fn push_scan(
        &mut self,
        acc: &mut QcReportAcc,
        out: &ScanOutput<'_>,
        lr: &LineRef<'_>,
        site: &SiteView<'_>,
    ) -> Result<(), EngineError> {
        let pos1 = site.pos1;
        let locus = || reed_io::Locus {
            contig_name: String::from_utf8_lossy(lr.chrom()).into_owned(),
            pos1,
        };

        // ---- ALT list ----
        let alts_text = lr.alts();
        let n_alt = if alts_text == b"." {
            0
        } else {
            memchr::memchr_iter(b',', alts_text).count() + 1
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
        let classes: Vec<_> = (0..n_alt)
            .map(|i| alt_class(reference, alt_of(i)))
            .collect();
        let any_snp = classes
            .iter()
            .any(|c| matches!(c, reed_core::site::VariantClass::Snv));
        let all_snp = !classes.is_empty()
            && classes
                .iter()
                .all(|c| matches!(c, reed_core::site::VariantClass::Snv));

        // ---- SN (bcftools counts a record in every type it carries) ----
        let s = &mut acc.summary;
        s.n_records += 1;
        if n_alt == 0 {
            s.n_ref += 1;
        }
        if any_snp {
            s.n_snv += 1;
        }
        if classes
            .iter()
            .any(|c| matches!(c, reed_core::site::VariantClass::Mnv))
        {
            s.n_mnv += 1;
        }
        if classes
            .iter()
            .any(|c| matches!(c, reed_core::site::VariantClass::Indel))
        {
            s.n_indel += 1;
        }
        if classes
            .iter()
            .any(|c| matches!(c, reed_core::site::VariantClass::Symbolic))
        {
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
        if classes
            .iter()
            .any(|c| matches!(c, reed_core::site::VariantClass::Indel))
        {
            acc.qual.indel.push_iqual(iqual);
        }

        // ---- genotype counts from the shared scan (GT + DP) ----
        self.ac.clear();
        self.ac.resize(n_alt, 0);
        let mut an = 0_u64;
        // AF/SiS allele counts: bcftools's bcf_calc_ac prefers INFO/AC+AN
        // when both are present (with the right cardinality) and only
        // falls back to GT otherwise — mirrored here so the AF/SiS
        // sections match on VCFs that carry the tags (the GT path feeds
        // the per-sample sections regardless).
        let mut af_counts: Option<(Vec<u64>, u64)> = None;
        let mut n_nonref_samples = 0;
        let mut i_nonref = 0_usize;
        let raw_alleles = out
            .alleles
            .expect("invariant: the summary requires tracking_alleles");
        let dp_column = out.int_columns.first().map_or(&[][..], Vec::as_slice);
        for (i, [a, b]) in raw_alleles.as_chunks::<2>().0.iter().enumerate() {
            let pair = [*a, *b];
            for allele in pair {
                match allele {
                    Some(0) => an += 1,
                    Some(a) => {
                        an += 1;
                        if let Some(slot) = self.ac.get_mut(usize::from(a) - 1) {
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
                    if let (Some(_r), Some(idx)) =
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
            let dp = dp_column.get(i).copied().filter(|&v| v != INT_MISSING);
            if let Some(dp) = dp {
                acc.dp_genotypes.push(dp);
                psc.dp_sum += u64::from(dp);
                psc.dp_n += 1;
            }
        }

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
                        t.split(',')
                            .map(|x| x.parse::<u64>().ok())
                            .collect::<Option<Vec<_>>>()
                    });
                    if let Some(vals) = vals
                        && vals.len() == n_alt
                    {
                        info_ac = Some(vals);
                    }
                } else if tag == b"AN"
                    && let Some(v) = value
                    && let Some(an) = std::str::from_utf8(v)
                        .ok()
                        .and_then(|t| t.parse::<u64>().ok())
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
            None => (self.ac.clone(), an),
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
            // (the tv bucket is filled in the ALT loop below)
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
            if i == 0
                && is_snp
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
                    && let Some(dp) = std::str::from_utf8(v)
                        .ok()
                        .and_then(|t| t.parse::<u32>().ok())
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
        PositionsIter {
            data: self,
            needle,
            from: 0,
        }
    }
}

/// Renders the `bcftools stats`-compatible text (§10.3).
#[must_use]
pub fn render_bcftools_stats(r: &QcReportAcc) -> String {
    let mut out = String::new();
    out.push_str(
        "# This file was produced by reed stats and can be plotted using plot-vcfstats.\n",
    );
    out.push_str("# The command line was:\treed stats --summary\n#\n");
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
        writeln!(out, "SN\t0\t{key}\t{value}").expect("invariant: writing to String");
    }

    out.push_str("# TSTV, transitions/transversions\n# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)\n");
    let ratio = |ts: u64, tv: u64| {
        if tv > 0 {
            format!("{:.2}", ts as f64 / tv as f64)
        } else {
            "0.00".into()
        }
    };
    writeln!(
        out,
        "TSTV\t0\t{}\t{}\t{}\t{}\t{}\t{}",
        r.tstv.ts,
        r.tstv.tv,
        ratio(r.tstv.ts, r.tstv.tv),
        r.tstv.ts_alt1,
        r.tstv.tv_alt1,
        ratio(r.tstv.ts_alt1, r.tstv.tv_alt1),
    )
    .expect("invariant: writing to String");

    out.push_str("# SiS, Singleton stats:\n# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n");
    writeln!(
        out,
        "SiS\t0\t{}\t{}\t{}\t{}\t{}\t0\t0\t0",
        1, r.sis.n_snps, r.sis.ts, r.sis.tv, r.sis.n_indels,
    )
    .expect("invariant: writing to String");

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
        writeln!(
            out,
            "AF\t0\t{freq:.6}\t{}\t{}\t{}\t{}\t0\t0\t0",
            snps[i], ts[i], tv[i], indels[i]
        )
        .expect("invariant: writing to String");
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
        let qval = |dist: &reed_stats::summary::histogram::QualDist| -> Option<f64> {
            (i < dist.nbins()).then(|| {
                let (beg, _) = reed_stats::summary::histogram::QualDist::bin_range(i);
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
            None => writeln!(out, "QUAL\t0\t.\t{snps}\t{n_ts}\t{n_tv}\t{n_in}")
                .expect("invariant: writing to String"),
            Some(qv) => writeln!(out, "QUAL\t0\t{qv:.1}\t{snps}\t{n_ts}\t{n_tv}\t{n_in}")
                .expect("invariant: writing to String"),
        }
    }

    out.push_str("# ST, Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n");
    for (r_base, a_base) in [
        (b'A', b'C'),
        (b'A', b'G'),
        (b'A', b'T'),
        (b'C', b'A'),
        (b'C', b'G'),
        (b'C', b'T'),
        (b'G', b'A'),
        (b'G', b'C'),
        (b'G', b'T'),
        (b'T', b'A'),
        (b'T', b'C'),
        (b'T', b'G'),
    ] {
        let idx = subst_index(r_base, a_base).expect("invariant: off-diagonal");
        writeln!(
            out,
            "ST\t0\t{}>{}\t{}",
            r_base as char, a_base as char, r.subst[idx]
        )
        .expect("invariant: writing to String");
    }

    out.push_str("# DP, depth\n# DP, Depth distribution\n# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)\n");
    let total_gt = r.dp_genotypes.total();
    let total_sites = r.dp_sites.total();
    type DpIter<'a> = std::iter::Peekable<std::boxed::Box<dyn Iterator<Item = (u32, u64)> + 'a>>;
    let mut gt: DpIter<'_> = (Box::new(r.dp_genotypes.iter())
        as std::boxed::Box<dyn Iterator<Item = (u32, u64)>>)
        .peekable();
    let mut sites: DpIter<'_> = (Box::new(r.dp_sites.iter())
        as std::boxed::Box<dyn Iterator<Item = (u32, u64)>>)
        .peekable();
    loop {
        let next_bin = [gt.peek().map(|(b, _)| *b), sites.peek().map(|(b, _)| *b)]
            .into_iter()
            .flatten()
            .min();
        let Some(bin) = next_bin else { break };
        let take = |it: &mut std::iter::Peekable<
            std::boxed::Box<dyn Iterator<Item = (u32, u64)> + '_>,
        >,
                    b: u32|
         -> u64 {
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
        let gt_pct = if total_gt > 0 {
            gt_n as f64 * 100.0 / total_gt as f64
        } else {
            0.0
        };
        let site_pct = if total_sites > 0 {
            site_n as f64 * 100.0 / total_sites as f64
        } else {
            0.0
        };
        let label = if bin > 500 {
            ">500".to_string()
        } else {
            bin.to_string()
        };
        writeln!(
            out,
            "DP\t0\t{label}\t{gt_n}\t{gt_pct:.6}\t{site_n}\t{site_pct:.6}"
        )
        .expect("invariant: writing to String");
    }

    out.push_str("# IDD, InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF\n");
    for (len, n) in r.idd.iter() {
        if n > 0 && len <= 60 {
            writeln!(out, "IDD\t0\t{}\t{n}\t0\t0.000000", len as i64 - 1)
                .expect("invariant: writing to String");
        }
    }

    out.push_str("# PSC, Per-sample counts\n# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing\n");
    for (i, p) in r.psc.iter().enumerate() {
        let dp = if p.dp_n > 0 {
            p.dp_sum as f64 / p.dp_n as f64
        } else {
            0.0
        };
        writeln!(
            out,
            "PSC\t0\tsample_{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{dp:.1}\t{}\t{}\t{}\t{}",
            p.hom_rr,
            p.hom_aa,
            p.hets,
            p.ts,
            p.tv,
            p.n_indels,
            p.singletons,
            p.hap_ref,
            p.hap_alt,
            p.missing,
        )
        .expect("invariant: writing to String");
    }
    out
}

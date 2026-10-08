//! The full QC report (§10.1): every section is an additive accumulator
//! with bcftools-stats semantics (SN / TSTV / SiS / AF / QUAL / ST / DP /
//! IDD / PSC), plus reed-only sections (per-TAG INFO distributions).
//!
//! bcftools semantics pinned by source (`vcfstats.c`) and fixture tests:
//! - AF bins: `m_af = max(101, n_samples + 1)`; bin =
//!   `floor(af × (m_af − 2))`; AC == 1 → the singletons bin (printed as
//!   0.000000); printed value is `bin / (m_af − 1)`.
//! - QUAL bins: 0.1 spaced; missing/NaN/negative → bin 0.
//! - Ts/Tv: `|ref_idx − alt_idx| == 2` on the A<C<G<T ring is a
//!   transition.
//! - SN: a record counts in every type it carries (a SNP+indel row
//!   increments both).
//! - DP (no `-s`): per-site INFO/DP; with `-s -`: per-genotype FMT/DP
//!   (AD fallback: the maximum non-ref AD).
//! - PSC: homRR/homAA/hets/ts/tv per sample over SNP rows; depth from
//!   FMT/DP only.

use reed_core::site::VariantClass;
use serde::Serialize;

use crate::summary::histogram::{IntHistogram, QualDist};

/// The whole-report accumulator (§10.1). All fields are additive.
#[derive(Clone, Debug, Serialize)]
pub struct QcReportAcc {
    /// Number of samples in the header.
    pub n_samples: u32,
    /// SN: record counts by class (§10.1 summary).
    pub summary: SummaryAcc,
    /// TSTV: total transitions / transversions (all ALTs) and 1st-ALT.
    pub tstv: TstvAcc,
    /// SiS: singleton records (AC == 1).
    pub sis: SiSAcc,
    /// AF: per-allele-frequency-bin site/ts/tv counts (bcftools AF).
    pub af: AfAcc,
    /// QUAL: the bcftools QUAL histogram per variant class.
    pub qual: QualAcc,
    /// ST: the 16-cell substitution matrix `ref << 2 | alt` over the
    /// A<C<G<T ring (the 12 off-diagonal cells are printed).
    pub subst: [u64; 16],
    /// DP: per-genotype depth histogram (FMT/DP; `-s -` semantics).
    pub dp_genotypes: IntHistogram,
    /// DP: per-site depth histogram (INFO/DP).
    pub dp_sites: IntHistogram,
    /// IDD: indel length distribution (negative = deletion).
    pub idd: IntHistogram,
    /// PSC: per-sample counters.
    pub psc: Vec<PscAcc>,
}

/// SN + Ts/Tv + SiS + AF accumulator parts.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct SummaryAcc {
    /// Total records seen.
    pub n_records: u64,
    /// Records whose primary ALT is a single-nucleotide substitution.
    pub n_snv: u64,
    /// Records whose primary ALT is a multi-nucleotide substitution.
    pub n_mnv: u64,
    /// Records whose primary ALT is an insertion or deletion.
    pub n_indel: u64,
    /// Records whose primary ALT is symbolic (`<...>`, breakends).
    pub n_symbolic: u64,
    /// Records with no ALT (`.` or REF-equal): monomorphic sites.
    pub n_ref: u64,
    /// Records carrying more than one ALT allele.
    pub n_multiallelic: u64,
    /// Multiallelic records whose every ALT is an SNV.
    pub n_multiallelic_snv: u64,
}

/// Transitions/transversions.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct TstvAcc {
    /// Transitions over every ALT of every record.
    pub ts: u64,
    /// Transversions over every ALT.
    pub tv: u64,
    /// Transitions of the first ALT.
    pub ts_alt1: u64,
    /// Transversions of the first ALT.
    pub tv_alt1: u64,
}

/// Singleton records (AC == 1 over the counted samples).
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct SiSAcc {
    /// Allele count of the singleton (always 1; kept for the table shape).
    pub allele_count: u64,
    /// Singleton records.
    pub n_snps: u64,
    /// Singleton transitions.
    pub ts: u64,
    /// Singleton transversions.
    pub tv: u64,
    /// Singleton indels.
    pub n_indels: u64,
}

/// Per-AF-bin site counts (bcftools AF section).
#[derive(Clone, Debug, Serialize)]
pub struct AfAcc {
    /// `m_af = max(101, n_samples + 1)`; bin = `floor(af × (m_af − 2))`;
    /// index 0 is the singletons bin.
    pub m_af: usize,
    /// Per bin: sites (SNPs), ts, tv, indels.
    pub snps: Vec<u64>,
    /// Per-bin transitions.
    pub ts: Vec<u64>,
    /// Per-bin transversions.
    pub tv: Vec<u64>,
    /// Per-bin indels.
    pub indels: Vec<u64>,
}

impl AfAcc {
    /// Bins for `n_samples`.
    #[must_use]
    pub fn new(n_samples: u32) -> Self {
        let m_af = 101_usize.max(n_samples as usize + 1);
        Self {
            m_af,
            snps: vec![0; m_af],
            ts: vec![0; m_af],
            tv: vec![0; m_af],
            indels: vec![0; m_af],
        }
    }
}

/// The QUAL histograms per variant class (bcftools QUAL section: per-bin
/// SNP counts with first-ALT ts/tv, plus per-bin indel counts).
#[derive(Clone, Debug, Serialize)]
pub struct QualAcc {
    /// SNV sites (0.1 bins).
    pub snv: QualDist,
    /// First-ALT transitions per SNV bin (0.1 bins).
    pub snv_ts: QualDist,
    /// First-ALT transversions per SNV bin (0.1 bins).
    pub snv_tv: QualDist,
    /// Indel sites (0.1 bins).
    pub indel: QualDist,
}

/// Per-sample counters (bcftools PSC).
#[derive(Clone, Debug, Default, Serialize)]
pub struct PscAcc {
    /// Homozygous reference.
    pub hom_rr: u64,
    /// Homozygous alt.
    pub hom_aa: u64,
    /// Heterozygous (same or different ALTs).
    pub hets: u64,
    /// Transitions in this sample's non-ref genotypes.
    pub ts: u64,
    /// Transversions.
    pub tv: u64,
    /// Indel genotypes.
    pub n_indels: u64,
    /// Sum of FMT/DP over genotypes with DP > 0.
    pub dp_sum: u64,
    /// Genotypes with DP > 0.
    pub dp_n: u64,
    /// Singletons: the sample is the only non-ref.
    pub singletons: u64,
    /// Haploid reference calls.
    pub hap_ref: u64,
    /// Haploid alt calls.
    pub hap_alt: u64,
    /// Missing genotypes.
    pub missing: u64,
}

impl QcReportAcc {
    /// An empty report for `n_samples`.
    #[must_use]
    pub fn new(n_samples: u32) -> Self {
        Self {
            n_samples,
            summary: SummaryAcc::default(),
            tstv: TstvAcc::default(),
            sis: SiSAcc::default(),
            af: AfAcc::new(n_samples),
            qual: QualAcc {
                snv: QualDist::new(),
                snv_ts: QualDist::new(),
                snv_tv: QualDist::new(),
                indel: QualDist::new(),
            },
            subst: [0; 16],
            dp_genotypes: IntHistogram::new(500),
            dp_sites: IntHistogram::new(500),
            idd: IntHistogram::new(60),
            psc: vec![PscAcc::default(); n_samples as usize],
        }
    }

    /// Component-wise merge (additive; §8.4 laws).
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.n_samples, other.n_samples, "invariant: same sample count");
        let SummaryAcc {
            n_records,
            n_snv,
            n_mnv,
            n_indel,
            n_symbolic,
            n_ref,
            n_multiallelic,
            n_multiallelic_snv,
        } = other.summary;
        self.summary.n_records += n_records;
        self.summary.n_snv += n_snv;
        self.summary.n_mnv += n_mnv;
        self.summary.n_indel += n_indel;
        self.summary.n_symbolic += n_symbolic;
        self.summary.n_ref += n_ref;
        self.summary.n_multiallelic += n_multiallelic;
        self.summary.n_multiallelic_snv += n_multiallelic_snv;

        self.tstv.ts += other.tstv.ts;
        self.tstv.tv += other.tstv.tv;
        self.tstv.ts_alt1 += other.tstv.ts_alt1;
        self.tstv.tv_alt1 += other.tstv.tv_alt1;

        self.sis.n_snps += other.sis.n_snps;
        self.sis.ts += other.sis.ts;
        self.sis.tv += other.sis.tv;
        self.sis.n_indels += other.sis.n_indels;

        for (dst, src) in self.af.snps.iter_mut().zip(&other.af.snps) {
            *dst += src;
        }
        for (dst, src) in self.af.ts.iter_mut().zip(&other.af.ts) {
            *dst += src;
        }
        for (dst, src) in self.af.tv.iter_mut().zip(&other.af.tv) {
            *dst += src;
        }
        for (dst, src) in self.af.indels.iter_mut().zip(&other.af.indels) {
            *dst += src;
        }

        self.qual.snv.merge(&other.qual.snv);
        self.qual.snv_ts.merge(&other.qual.snv_ts);
        self.qual.snv_tv.merge(&other.qual.snv_tv);
        self.qual.indel.merge(&other.qual.indel);
        for (dst, src) in self.subst.iter_mut().zip(&other.subst) {
            *dst += src;
        }
        self.dp_genotypes.merge(&other.dp_genotypes);
        self.dp_sites.merge(&other.dp_sites);
        self.idd.merge(&other.idd);
        for (dst, src) in self.psc.iter_mut().zip(&other.psc) {
            let PscAcc {
                hom_rr,
                hom_aa,
                hets,
                ts,
                tv,
                n_indels,
                dp_sum,
                dp_n,
                singletons,
                hap_ref,
                hap_alt,
                missing,
            } = src;
            dst.hom_rr += hom_rr;
            dst.hom_aa += hom_aa;
            dst.hets += hets;
            dst.ts += ts;
            dst.tv += tv;
            dst.n_indels += n_indels;
            dst.dp_sum += dp_sum;
            dst.dp_n += dp_n;
            dst.singletons += singletons;
            dst.hap_ref += hap_ref;
            dst.hap_alt += hap_alt;
            dst.missing += missing;
        }
    }
}

/// The A<C<G<T ring index (bcftools `bcf_acgt2int`).
#[must_use]
pub fn acgt(b: u8) -> Option<u8> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Transition test on the ring: `|ref − alt| == 2` (A↔G, C↔T).
#[must_use]
pub fn is_transition(ref_b: u8, alt_b: u8) -> bool {
    match (acgt(ref_b), acgt(alt_b)) {
        (Some(r), Some(a)) => r.abs_diff(a) == 2,
        _ => false,
    }
}

/// Substitution index: `ref << 2 | alt` over the A<C<G<T ring (0..16;
/// the diagonal is unreachable for distinct bases).
#[must_use]
pub fn subst_index(ref_b: u8, alt_b: u8) -> Option<usize> {
    match (acgt(ref_b), acgt(alt_b)) {
        (Some(r), Some(a)) if r != a => Some((usize::from(r) << 2) | usize::from(a)),
        _ => None,
    }
}

/// The variant class of one ALT allele (bcftools-style per-allele typing;
/// equal-length > 1 is MNP, unequal is indel, symbolic otherwise).
#[must_use]
pub fn alt_class(reference: &[u8], alt: &[u8]) -> VariantClass {
    if alt == b"." || alt == reference {
        return VariantClass::Ref;
    }
    if alt.starts_with(b"<")
        || alt.ends_with(b">")
        || alt.contains(&b'[')
        || alt.contains(&b']')
        || alt.contains(&b'*')
    {
        return VariantClass::Symbolic;
    }
    if reference.len() == alt.len() {
        if reference.len() == 1 {
            VariantClass::Snv
        } else {
            VariantClass::Mnv
        }
    } else {
        VariantClass::Indel
    }
}

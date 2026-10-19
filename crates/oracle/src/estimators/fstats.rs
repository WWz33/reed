//! f-statistics: f2, f3, f4, and the D-statistic (Patterson et al. 2012).
//!
//! All are computed over sites where every involved group has at least one
//! called allele; ALT-allele frequencies are taken from counts (per-allele
//! deletion). Block-jackknife standard errors are a v0.3 concern.
//!
//! Manifest gates (§4.2, judged per estimator): f2/f3/f4 are **site means**
//! that include monomorphic sites under all-sites input (they contribute 0 to
//! the sum but 1 to the count), so any variant-only input changes the
//! quantity — biased under every `VariantOnly` manifest. The D-statistic is a
//! **ratio of sums** (a monomorphic site contributes 0 to numerator and
//! denominator alike), so it is invariant to monomorphic-site presence:
//! `VariantOnly(Complete)` matches all-sites input, while `Filtered` — and
//! `Unknown`, which cannot be distinguished from it — loses pairs to removed
//! true variants → bias note.

use crate::population::PopulationView;
use crate::window::{ContigId, WindowBuffer};

use super::{samples_of, Unavailable};

fn alt_freq(window: &WindowBuffer, samples: &[u32], site: usize) -> Option<f64> {
    let counts = window.allele_counts(site, samples);
    let n = counts.ref_count + counts.alt_count;
    (n > 0).then(|| counts.alt_count as f64 / n as f64)
}

/// Bias note for the site-mean statistics (f2/f3/f4) under any
/// `VariantOnly` manifest.
const FSTATS_BIAS: &str = "variant-only input: monomorphic sites are absent, so this site mean is \
     inflated relative to all-sites semantics (site-density dependent); removed \
     true variants additionally bias the numerator";

/// Bias note for the D-statistic under `VariantOnly(Filtered/Unknown)`.
const DSTAT_BIAS: &str =
    "variant-only input: pairs involving removed true variants are missing; the \
     ratio of sums is otherwise invariant to monomorphic-site presence, and \
     Unknown cannot be distinguished from Filtered";

/// f2(A, B) = the bias-corrected Patterson sum
/// `Σ[(pA−pB)² − ĥA/anA − ĥB/anB]` over sites with anA, anB ≥ 2
/// (allel `patterson_f2`; doc_v15 D15-31).
#[derive(Clone, Copy, Debug)]
pub struct F2Result {
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group_a: u32,
    pub group_b: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub sum_sq: f64,
    pub bias: Option<&'static str>,
}

impl F2Result {
    /// Σ of the corrected per-site terms (allel `nansum`).
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then_some(self.sum_sq)
    }
}

/// f3(target; b, c) = the corrected sum
/// `Σ[(pT−pB)(pT−pC) − ĥT/anT]` over sites with anT ≥ 2 and anB, anC ≥ 1
/// (allel `patterson_f3`; doc_v15 D15-31).
#[derive(Clone, Copy, Debug)]
pub struct F3Result {
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub target: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub sum_prod: f64,
    pub bias: Option<&'static str>,
}

impl F3Result {
    /// Σ of the corrected per-site terms (allel `nansum`).
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then_some(self.sum_prod)
    }
}

/// f4(A, B; C, D) = `Σ(pA−pB)(pC−pD)` over sites with all four groups
/// called (allel `patterson_d` numerator; doc_v15 D15-31).
#[derive(Clone, Copy, Debug)]
pub struct F4Result {
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group_a: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub group_d: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub sum_prod: f64,
    pub bias: Option<&'static str>,
}

impl F4Result {
    /// Σ of the per-site numerator products (allel `nansum` of the f4 term).
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then_some(self.sum_prod)
    }
}

/// D-statistic = `Σ num / Σ den`, with
/// `num = (pA − pB)(pC − pD)` and
/// `den = (pA + pB − 2 pA pB)(pC + pD − 2 pC pD)`.
#[derive(Clone, Copy, Debug)]
pub struct DStatResult {
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group_a: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub group_d: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub num: f64,
    pub den: f64,
    pub bias: Option<&'static str>,
}

impl DStatResult {
    pub fn value(&self) -> Option<f64> {
        (self.den > 0.0).then(|| self.num / self.den)
    }
}

fn require_distinct(pairs: &[(u32, u32)]) -> Result<(), Unavailable> {
    for (i, &(a, b)) in pairs.iter().enumerate() {
        if a == b {
            return Err(Unavailable {
                reason: "f-statistics require distinct groups",
            });
        }
        for &(c, d) in &pairs[i + 1..] {
            if a == c && b == d {
                return Err(Unavailable {
                    reason: "f-statistics require distinct groups",
                });
            }
        }
    }
    Ok(())
}

/// f2 between two groups.
pub fn f2(
    window: &WindowBuffer,
    pops: &PopulationView,
    group_a: u32,
    group_b: u32,
) -> Result<F2Result, Unavailable> {
    require_distinct(&[(group_a, group_b)])?;
    let sa = samples_of(pops, window, group_a)?;
    let sb = samples_of(pops, window, group_b)?;
    let bias = super::variant_only_note(window.manifest, FSTATS_BIAS);
    let (mut n_sites, mut n_skipped, mut sum_sq) = (0u64, 0u64, 0.0f64);
    for site in 0..window.n_sites() {
        let a = window.allele_counts(site, sa);
        let b = window.allele_counts(site, sb);
        let an_a = a.ref_count + a.alt_count;
        let an_b = b.ref_count + b.alt_count;
        if an_a < 2 || an_b < 2 {
            n_skipped += 1;
            continue;
        }
        let pa = a.alt_count as f64 / an_a as f64;
        let pb = b.alt_count as f64 / an_b as f64;
        n_sites += 1;
        sum_sq += (pa - pb).powi(2)
            - (a.ref_count * a.alt_count) as f64 / (an_a * an_a * (an_a - 1)) as f64
            - (b.ref_count * b.alt_count) as f64 / (an_b * an_b * (an_b - 1)) as f64;
    }
    Ok(F2Result {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group_a,
        group_b,
        n_sites,
        n_skipped,
        sum_sq,
        bias,
    })
}

/// f3 with `target` as the outgroup-ish target (see formula in the type docs).
pub fn f3(
    window: &WindowBuffer,
    pops: &PopulationView,
    target: u32,
    group_b: u32,
    group_c: u32,
) -> Result<F3Result, Unavailable> {
    require_distinct(&[(target, group_b), (target, group_c), (group_b, group_c)])?;
    let st = samples_of(pops, window, target)?;
    let sb = samples_of(pops, window, group_b)?;
    let sc = samples_of(pops, window, group_c)?;
    let bias = super::variant_only_note(window.manifest, FSTATS_BIAS);
    let (mut n_sites, mut n_skipped, mut sum_prod) = (0u64, 0u64, 0.0f64);
    for site in 0..window.n_sites() {
        let t = window.allele_counts(site, st);
        let b = window.allele_counts(site, sb);
        let c = window.allele_counts(site, sc);
        let an_t = t.ref_count + t.alt_count;
        let an_b = b.ref_count + b.alt_count;
        let an_c = c.ref_count + c.alt_count;
        if an_t < 2 || an_b == 0 || an_c == 0 {
            n_skipped += 1;
            continue;
        }
        let pt = t.alt_count as f64 / an_t as f64;
        let pb = b.alt_count as f64 / an_b as f64;
        let pc = c.alt_count as f64 / an_c as f64;
        let h_hat = (t.ref_count * t.alt_count) as f64 / (an_t * (an_t - 1)) as f64;
        n_sites += 1;
        sum_prod += (pt - pb) * (pt - pc) - h_hat / an_t as f64;
    }
    Ok(F3Result {
        contig: window.contig,
        start: window.start,
        end: window.end,
        target,
        group_b,
        group_c,
        n_sites,
        n_skipped,
        sum_prod,
        bias,
    })
}

/// f4 between the pairs (A, B) and (C, D).
pub fn f4(
    window: &WindowBuffer,
    pops: &PopulationView,
    group_a: u32,
    group_b: u32,
    group_c: u32,
    group_d: u32,
) -> Result<F4Result, Unavailable> {
    require_distinct(&[
        (group_a, group_b),
        (group_a, group_c),
        (group_a, group_d),
        (group_b, group_c),
        (group_b, group_d),
        (group_c, group_d),
    ])?;
    let sa = samples_of(pops, window, group_a)?;
    let sb = samples_of(pops, window, group_b)?;
    let sc = samples_of(pops, window, group_c)?;
    let sd = samples_of(pops, window, group_d)?;
    let bias = super::variant_only_note(window.manifest, FSTATS_BIAS);
    let (mut n_sites, mut n_skipped, mut sum_prod) = (0u64, 0u64, 0.0f64);
    for site in 0..window.n_sites() {
        let a = window.allele_counts(site, sa);
        let b = window.allele_counts(site, sb);
        let c = window.allele_counts(site, sc);
        let d = window.allele_counts(site, sd);
        if a.ref_count + a.alt_count == 0
            || b.ref_count + b.alt_count == 0
            || c.ref_count + c.alt_count == 0
            || d.ref_count + d.alt_count == 0
        {
            n_skipped += 1;
            continue;
        }
        let pa = a.alt_count as f64 / (a.ref_count + a.alt_count) as f64;
        let pb = b.alt_count as f64 / (b.ref_count + b.alt_count) as f64;
        let pc = c.alt_count as f64 / (c.ref_count + c.alt_count) as f64;
        let pd = d.alt_count as f64 / (d.ref_count + d.alt_count) as f64;
        n_sites += 1;
        sum_prod += (pa - pb) * (pc - pd);
    }
    Ok(F4Result {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group_a,
        group_b,
        group_c,
        group_d,
        n_sites,
        n_skipped,
        sum_prod,
        bias,
    })
}

/// D-statistic (ABBA-BABA) between the pairs (A, B) and (C, D).
pub fn d_stat(
    window: &WindowBuffer,
    pops: &PopulationView,
    group_a: u32,
    group_b: u32,
    group_c: u32,
    group_d: u32,
) -> Result<DStatResult, Unavailable> {
    require_distinct(&[
        (group_a, group_b),
        (group_a, group_c),
        (group_a, group_d),
        (group_b, group_c),
        (group_b, group_d),
        (group_c, group_d),
    ])?;
    let sa = samples_of(pops, window, group_a)?;
    let sb = samples_of(pops, window, group_b)?;
    let sc = samples_of(pops, window, group_c)?;
    let sd = samples_of(pops, window, group_d)?;
    let bias = super::monomorphic_invariant_note(window.manifest, DSTAT_BIAS);
    let (mut n_sites, mut n_skipped, mut num, mut den) = (0u64, 0u64, 0.0f64, 0.0f64);
    for site in 0..window.n_sites() {
        match (
            alt_freq(window, sa, site),
            alt_freq(window, sb, site),
            alt_freq(window, sc, site),
            alt_freq(window, sd, site),
        ) {
            (Some(pa), Some(pb), Some(pc), Some(pd)) => {
                n_sites += 1;
                num += (pa - pb) * (pc - pd);
                den += (pa + pb - 2.0 * pa * pb) * (pc + pd - 2.0 * pc * pd);
            }
            _ => n_skipped += 1,
        }
    }
    Ok(DStatResult {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group_a,
        group_b,
        group_c,
        group_d,
        n_sites,
        n_skipped,
        num,
        den,
        bias,
    })
}

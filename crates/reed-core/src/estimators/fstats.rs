//! f-statistics: f2, f3, f4, and the D-statistic (Patterson et al. 2012).
//!
//! All are means over sites where every involved group has at least one called
//! allele; ALT-allele frequencies are taken from counts (per-allele deletion).
//! Block-jackknife standard errors are a v0.3 concern.

use crate::population::PopulationView;
use crate::window::WindowBuffer;

use super::Unavailable;

fn alt_freq(window: &WindowBuffer, samples: &[u32], site: usize) -> Option<f64> {
    let counts = window.allele_counts(site, samples);
    let n = counts.ref_count + counts.alt_count;
    (n > 0).then(|| counts.alt_count as f64 / n as f64)
}

/// f2(A, B) = mean over sites of `(pA − pB)²`.
#[derive(Clone, Copy, Debug)]
pub struct F2Result {
    pub group_a: u32,
    pub group_b: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub sum_sq: f64,
}

impl F2Result {
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then(|| self.sum_sq / self.n_sites as f64)
    }
}

/// f3(target; b, c) = mean over sites of `(pB − pTarget) * (pC − pTarget)`.
#[derive(Clone, Copy, Debug)]
pub struct F3Result {
    pub target: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub sum_prod: f64,
}

impl F3Result {
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then(|| self.sum_prod / self.n_sites as f64)
    }
}

/// f4(A, B; C, D) = mean over sites of `(pA − pB) * (pC − pD)`.
#[derive(Clone, Copy, Debug)]
pub struct F4Result {
    pub group_a: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub group_d: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub sum_prod: f64,
}

impl F4Result {
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then(|| self.sum_prod / self.n_sites as f64)
    }
}

/// D-statistic = `Σ num / Σ den`, with
/// `num = (pA − pB)(pC − pD)` and
/// `den = (pA + pB − 2 pA pB)(pC + pD − 2 pC pD)`.
#[derive(Clone, Copy, Debug)]
pub struct DStatResult {
    pub group_a: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub group_d: u32,
    pub n_sites: u64,
    pub n_skipped: u64,
    pub num: f64,
    pub den: f64,
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
    let sa = pops.samples(group_a);
    let sb = pops.samples(group_b);
    let (mut n_sites, mut n_skipped, mut sum_sq) = (0u64, 0u64, 0.0f64);
    for site in 0..window.n_sites() {
        match (alt_freq(window, sa, site), alt_freq(window, sb, site)) {
            (Some(pa), Some(pb)) => {
                n_sites += 1;
                sum_sq += (pa - pb).powi(2);
            }
            _ => n_skipped += 1,
        }
    }
    Ok(F2Result {
        group_a,
        group_b,
        n_sites,
        n_skipped,
        sum_sq,
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
    let st = pops.samples(target);
    let sb = pops.samples(group_b);
    let sc = pops.samples(group_c);
    let (mut n_sites, mut n_skipped, mut sum_prod) = (0u64, 0u64, 0.0f64);
    for site in 0..window.n_sites() {
        match (
            alt_freq(window, st, site),
            alt_freq(window, sb, site),
            alt_freq(window, sc, site),
        ) {
            (Some(pt), Some(pb), Some(pc)) => {
                n_sites += 1;
                sum_prod += (pb - pt) * (pc - pt);
            }
            _ => n_skipped += 1,
        }
    }
    Ok(F3Result {
        target,
        group_b,
        group_c,
        n_sites,
        n_skipped,
        sum_prod,
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
    let sa = pops.samples(group_a);
    let sb = pops.samples(group_b);
    let sc = pops.samples(group_c);
    let sd = pops.samples(group_d);
    let (mut n_sites, mut n_skipped, mut sum_prod) = (0u64, 0u64, 0.0f64);
    for site in 0..window.n_sites() {
        match (
            alt_freq(window, sa, site),
            alt_freq(window, sb, site),
            alt_freq(window, sc, site),
            alt_freq(window, sd, site),
        ) {
            (Some(pa), Some(pb), Some(pc), Some(pd)) => {
                n_sites += 1;
                sum_prod += (pa - pb) * (pc - pd);
            }
            _ => n_skipped += 1,
        }
    }
    Ok(F4Result {
        group_a,
        group_b,
        group_c,
        group_d,
        n_sites,
        n_skipped,
        sum_prod,
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
    let sa = pops.samples(group_a);
    let sb = pops.samples(group_b);
    let sc = pops.samples(group_c);
    let sd = pops.samples(group_d);
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
        group_a,
        group_b,
        group_c,
        group_d,
        n_sites,
        n_skipped,
        num,
        den,
    })
}

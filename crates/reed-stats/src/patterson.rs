//! Patterson f-statistics (`docs/spec/architecture.md` §9.5): the
//! Patterson (2012) unbiased estimators as implemented by scikit-allel's
//! `admixture.py` (MIT), in `Corrected` and `Naive` sealed definitions.
//! All accumulators carry per-site (numerator, denominator) sums so the
//! §8.6 block jackknife can reuse them directly.
//!
//! - f2(A,B) = (pA−pB)² − ĥA/anA − ĥB/anB, where ĥ = ac0·ac1/(an(an−1));
//!   sites with anA < 2 or anB < 2 are skipped (ĥ undefined).
//! - f3(C;A,B): T = (c−a)(c−b) − ĥC/anC, B = 2·ĥC (the normalized f3* is
//!   ΣT/ΣB); sites with anC < 2 skipped.
//! - f4(A,B;C,D) = (pA−pB)(pC−pD) — unbiased, no correction; D is the
//!   normalized form Σnum/Σden with den = (a+b−2ab)(c+d−2cd); sites with
//!   any an == 0 skipped (p undefined).
//! - Naive variants drop the ĥ correction terms (f2/f3 only; f4/D are
//!   identical under both).

use reed_core::count::SiteGroupCounts;
use reed_core::math::sum::NeumaierSum;

/// Which Patterson definition a value comes from (§9.5).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PattersonKind {
    /// The unbiased estimator (default).
    Corrected,
    /// The raw moment without the ĥ correction terms.
    Naive,
}

impl PattersonKind {
    /// Stable identifier.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::Corrected => "corrected",
            Self::Naive => "naive",
        }
    }
}

/// Per-site f2(A, B) (Patterson 2012; allel `patterson_f2`): the unbiased
/// `(pA−pB)² − ĥA/anA − ĥB/anB`, or the raw moment for `Naive`. `None`
/// when the site is skipped (anA or anB zero; corrected also needs both ≥ 2).
#[must_use]
pub fn f2_site(
    row: &[SiteGroupCounts],
    a: usize,
    b: usize,
    kind: PattersonKind,
) -> Option<f64> {
    let (g1, g2) = (&row[a], &row[b]);
    let an1 = u64::from(g1.n_ref_alleles + g1.n_alt_alleles);
    let an2 = u64::from(g2.n_ref_alleles + g2.n_alt_alleles);
    let corrected = kind == PattersonKind::Corrected;
    if corrected && (an1 < 2 || an2 < 2) {
        return None;
    }
    if an1 == 0 || an2 == 0 {
        return None;
    }
    let p1 = f64::from(g1.n_alt_alleles) / an1 as f64;
    let p2 = f64::from(g2.n_alt_alleles) / an2 as f64;
    let mut v = (p1 - p2) * (p1 - p2);
    if corrected {
        v -= f64::from(g1.n_ref_alleles) * f64::from(g1.n_alt_alleles)
            / (an1 as f64 * an1 as f64 * (an1 - 1) as f64);
        v -= f64::from(g2.n_ref_alleles) * f64::from(g2.n_alt_alleles)
            / (an2 as f64 * an2 as f64 * (an2 - 1) as f64);
    }
    Some(v)
}

/// Per-site f3(C; A, B): `(T, B)` — the un-normalized estimate and the
/// heterozygosity term the normalized f3* divides by (allel
/// `patterson_f3`). `None` when the site is skipped.
#[must_use]
pub fn f3_site(
    row: &[SiteGroupCounts],
    c: usize,
    a: usize,
    b: usize,
    kind: PattersonKind,
) -> Option<(f64, f64)> {
    let (gc, ga, gb) = (&row[c], &row[a], &row[b]);
    let anc = u64::from(gc.n_ref_alleles + gc.n_alt_alleles);
    let ana = u64::from(ga.n_ref_alleles + ga.n_alt_alleles);
    let anb = u64::from(gb.n_ref_alleles + gb.n_alt_alleles);
    let corrected = kind == PattersonKind::Corrected;
    if corrected && anc < 2 {
        return None;
    }
    if anc == 0 || ana == 0 || anb == 0 {
        return None;
    }
    let pc = f64::from(gc.n_alt_alleles) / anc as f64;
    let pa = f64::from(ga.n_alt_alleles) / ana as f64;
    let pb = f64::from(gb.n_alt_alleles) / anb as f64;
    let mut t = (pc - pa) * (pc - pb);
    let h_hat = f64::from(gc.n_ref_alleles) * f64::from(gc.n_alt_alleles)
        / (anc as f64 * (anc as f64 - 1.0));
    if corrected {
        t -= h_hat / anc as f64;
    }
    Some((t, 2.0 * h_hat))
}

/// Per-site f4(A, B; C, D): `(num, den)` — `num = (pA−pB)(pC−pD)` and the
/// D-denominator `(pA+pB−2pApB)(pC+pD−2pCpD)` (allel `patterson_d`). `None`
/// when any of the four populations has no called allele.
#[must_use]
pub fn f4_site(
    row: &[SiteGroupCounts],
    a: usize,
    b: usize,
    c: usize,
    d: usize,
) -> Option<(f64, f64)> {
    let g = [&row[a], &row[b], &row[c], &row[d]];
    let mut an = [0_u64; 4];
    for (i, gi) in g.iter().enumerate() {
        an[i] = u64::from(gi.n_ref_alleles + gi.n_alt_alleles);
        if an[i] == 0 {
            return None;
        }
    }
    let p = [
        f64::from(g[0].n_alt_alleles) / an[0] as f64,
        f64::from(g[1].n_alt_alleles) / an[1] as f64,
        f64::from(g[2].n_alt_alleles) / an[2] as f64,
        f64::from(g[3].n_alt_alleles) / an[3] as f64,
    ];
    let num = (p[0] - p[1]) * (p[2] - p[3]);
    let den = (p[0] + p[1] - 2.0 * p[0] * p[1]) * (p[2] + p[3] - 2.0 * p[2] * p[3]);
    Some((num, den))
}

/// f2(A, B) between two groups.
#[derive(Clone, Debug)]
pub struct F2Acc {
    a: usize,
    b: usize,
    corrected: bool,
    sum: NeumaierSum,
    n_sites: u64,
    n_skipped: u64,
}

impl F2Acc {
    /// Accumulator over the `a`/`b` group pair.
    #[must_use]
    pub fn new(a: usize, b: usize, kind: PattersonKind) -> Self {
        Self {
            a,
            b,
            corrected: kind == PattersonKind::Corrected,
            sum: NeumaierSum::new(),
            n_sites: 0,
            n_skipped: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let kind = if self.corrected { PattersonKind::Corrected } else { PattersonKind::Naive };
        match f2_site(row, self.a, self.b, kind) {
            Some(v) => {
                self.sum.add(v);
                self.n_sites += 1;
            }
            None => self.n_skipped += 1,
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.sum.merge(&other.sum);
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// The estimate Σf2 over contributing sites; `None` when empty.
    #[must_use]
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then_some(self.sum.value())
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }
}

/// f3(C; A, B) — the three-population test. Carries ΣT and ΣB so the
/// normalized f3* = ΣT/ΣB is available alongside the raw estimate.
#[derive(Clone, Debug)]
pub struct F3Acc {
    c: usize,
    a: usize,
    b: usize,
    corrected: bool,
    t: NeumaierSum,
    b_sum: NeumaierSum,
    n_sites: u64,
    n_skipped: u64,
}

impl F3Acc {
    /// Accumulator over the `(c; a, b)` triple.
    #[must_use]
    pub fn new(c: usize, a: usize, b: usize, kind: PattersonKind) -> Self {
        Self {
            c,
            a,
            b,
            corrected: kind == PattersonKind::Corrected,
            t: NeumaierSum::new(),
            b_sum: NeumaierSum::new(),
            n_sites: 0,
            n_skipped: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let kind = if self.corrected { PattersonKind::Corrected } else { PattersonKind::Naive };
        match f3_site(row, self.c, self.a, self.b, kind) {
            Some((t, b_term)) => {
                self.t.add(t);
                self.b_sum.add(b_term);
                self.n_sites += 1;
            }
            None => self.n_skipped += 1,
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.t.merge(&other.t);
        self.b_sum.merge(&other.b_sum);
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// The raw f3 estimate ΣT.
    #[must_use]
    pub fn value(&self) -> Option<f64> {
        (self.n_sites > 0).then_some(self.t.value())
    }

    /// The normalized f3* = ΣT / ΣB (Patterson 2012: B = 2·ĥC).
    #[must_use]
    pub fn value_normalized(&self) -> Option<f64> {
        let b = self.b_sum.value();
        (self.n_sites > 0 && b != 0.0).then(|| self.t.value() / b)
    }

    /// (ΣT, ΣB) for the block jackknife.
    #[must_use]
    pub const fn sums(&self) -> (&NeumaierSum, &NeumaierSum) {
        (&self.t, &self.b_sum)
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }
}

/// f4(A, B; C, D) and D(A, B; C, D) — the same per-site numerator and
/// denominator; f4 is Σnum, D is Σnum/Σden. Unbiased (no correction).
#[derive(Clone, Debug)]
pub struct F4Acc {
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    num: NeumaierSum,
    den: NeumaierSum,
    n_sites: u64,
    n_skipped: u64,
}

impl F4Acc {
    /// Accumulator over the `(a, b; c, d)` quadruple.
    #[must_use]
    pub fn new(a: usize, b: usize, c: usize, d: usize) -> Self {
        Self {
            a,
            b,
            c,
            d,
            num: NeumaierSum::new(),
            den: NeumaierSum::new(),
            n_sites: 0,
            n_skipped: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        match f4_site(row, self.a, self.b, self.c, self.d) {
            Some((num, den)) => {
                self.num.add(num);
                self.den.add(den);
                self.n_sites += 1;
            }
            None => self.n_skipped += 1,
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.num.merge(&other.num);
        self.den.merge(&other.den);
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// f4 = Σnum.
    #[must_use]
    pub fn value_f4(&self) -> Option<f64> {
        (self.n_sites > 0).then_some(self.num.value())
    }

    /// D = Σnum / Σden.
    #[must_use]
    pub fn value_d(&self) -> Option<f64> {
        let den = self.den.value();
        (self.n_sites > 0 && den != 0.0).then(|| self.num.value() / den)
    }

    /// (Σnum, Σden) for the block jackknife and the f4-ratio.
    #[must_use]
    pub const fn sums(&self) -> (&NeumaierSum, &NeumaierSum) {
        (&self.num, &self.den)
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }
}

/// f4-ratio: Σ(num quadruples) / Σ(den quadruples) over two sets of f4
/// statistics (§9.5; snputils semantics — the ratio of sums, with the
/// block jackknife providing the SE).
#[must_use]
pub fn f4_ratio(num: &F4Acc, den: &F4Acc) -> Option<f64> {
    let d = den.den.value();
    (num.n_sites > 0 && den.n_sites > 0 && d != 0.0).then(|| num.num.value() / d)
}

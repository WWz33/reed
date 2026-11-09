//! Between-group Fst: Hudson (1992) as Bhatia et al. 2013's ratio of sums
//! and Weir & Cockerham (1984) variance components
//! (`docs/rounds/v2/04-algorithms.md` §3). The components follow scikit-allel
//! and are verified against the reference implementation in `gold/popgen`
//! (the parity script names the tool). Nei 1986 is not implemented (04 §1).

use reed_core::count::SiteGroupCounts;
use reed_core::math::sum::NeumaierSum;

use crate::pi::PairAcc;

/// Which sealed Fst definition a value comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FstKind {
    /// Hudson 1992 / Bhatia 2013 (allel `hudson_fst`).
    Hudson,
    /// Weir & Cockerham 1984 variance components (allel `fst`).
    WeirCockerham84,
}

impl FstKind {
    /// Stable identifier.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::Hudson => "hudson",
            Self::WeirCockerham84 => "weir_cockerham_84",
        }
    }
}

/// WC84 variance components summed over sites (all exact float adds go
/// through `NeumaierSum`).
#[derive(Clone, Debug, Default)]
pub struct WcAcc {
    /// Between-population component.
    pub a: NeumaierSum,
    /// Within-population component.
    pub b: NeumaierSum,
    /// Inbreeding component.
    pub c: NeumaierSum,
    /// Contributing sites.
    pub n_sites: u64,
    /// Sites skipped by the guards.
    pub n_skipped: u64,
}

impl WcAcc {
    fn merge(&mut self, other: &Self) {
        self.a.merge(&other.a);
        self.b.merge(&other.b);
        self.c.merge(&other.c);
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }
}

/// One site's WC84 outcome at the **multi-population scale** (the scale
/// `WcMultiAcc::components` reports). The pair path is the `r = 2` instance,
/// so there is one formula and one set of tests (doc_v24 V24-05).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Wc84Site {
    /// A guard rejected the site (`n̄ ≤ 1`, `n_C ≤ 0`, or no called group).
    Skip,
    /// A group has no called diploid individual (`an = 0`, or a lone half
    /// call with `n = an // 2 = 0` when `n_C` collapses): a/b are undefined
    /// and only `c = h̄/2` is kept, matching allel's `nansum` (doc_v7 D7-20).
    KeepC { c: f64 },
    /// Full contribution; `Fst = a / (a + b + c)`.
    Full { a: f64, b: f64, c: f64 },
}

/// One population's WC84 inputs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Wc84Group {
    /// Called alleles (`an`).
    pub an: u64,
    /// Reference-allele count.
    pub refs: u64,
    /// Fully-called heterozygotes (`n_het`).
    pub het: u64,
}

/// The single WC84 formula. Pair path: two groups; multi path: r groups.
/// Allocation-free; callers own the group slice.
pub(crate) fn wc84_site(groups: &[Wc84Group]) -> Wc84Site {
    let r = groups.len();
    if r < 2 {
        return Wc84Site::Skip;
    }
    let (mut an_total, mut refs_total, mut h_sum) = (0_u64, 0_u64, 0_u64);
    let (mut n_total, mut n_sq_sum) = (0.0_f64, 0.0_f64);
    for group in groups {
        an_total += group.an;
        refs_total += group.refs;
        h_sum += group.het;
        let ni = (group.an / 2) as f64;
        n_total += ni;
        n_sq_sum += ni * ni;
    }
    if n_total == 0.0 {
        return Wc84Site::Skip;
    }
    let r_f = r as f64;
    let n_bar = n_total / r_f;
    let n_c = (n_total - n_sq_sum / n_total) / (r_f - 1.0);
    let any_empty = groups.iter().any(|g| g.an == 0);
    let any_half_only = groups.iter().any(|g| g.an < 2);
    if any_empty || (any_half_only && n_c <= 0.0) {
        return Wc84Site::KeepC {
            c: (h_sum as f64 / n_total) / 2.0,
        };
    }
    if n_bar <= 1.0 || n_c <= 0.0 {
        return Wc84Site::Skip;
    }
    let p_bar = refs_total as f64 / an_total as f64;
    let mut s2 = 0.0_f64;
    for group in groups {
        let pi = group.refs as f64 / group.an as f64;
        let ni = (group.an / 2) as f64;
        s2 += ni * (pi - p_bar).powi(2);
    }
    s2 /= n_bar * (r_f - 1.0);
    let h_bar = h_sum as f64 / n_total;
    let a = (n_bar / n_c)
        * (s2
            - (1.0 / (n_bar - 1.0))
                * (p_bar * (1.0 - p_bar) - (r_f - 1.0) * s2 / r_f - h_bar / 4.0));
    let b = (n_bar / (n_bar - 1.0))
        * (p_bar * (1.0 - p_bar)
            - (r_f - 1.0) * s2 / r_f
            - (2.0 * n_bar - 1.0) * h_bar / (4.0 * n_bar));
    Wc84Site::Full {
        a,
        b,
        c: h_bar / 2.0,
    }
}

/// Fst accumulator for one group pair. One pass fills every estimator's
/// components (the integer `PairAcc` audit trail is estimator-independent,
/// exactly as in the oracle).
#[derive(Clone, Debug)]
pub struct FstAcc {
    group_a: usize,
    group_b: usize,
    within_a: PairAcc,
    within_b: PairAcc,
    between: PairAcc,
    hudson_num: NeumaierSum,
    hudson_den: NeumaierSum,
    hudson_any: bool,
    wc: WcAcc,
    /// Variant sites pushed (the reference `no_snps`; invariant records do not count).
    n_sites: u64,
}

impl FstAcc {
    /// Accumulator for the `group_a`/`group_b` pair (indices into the row).
    #[must_use]
    pub fn new(group_a: usize, group_b: usize) -> Self {
        Self {
            group_a,
            group_b,
            within_a: PairAcc::default(),
            within_b: PairAcc::default(),
            between: PairAcc::default(),
            hudson_num: NeumaierSum::new(),
            hudson_den: NeumaierSum::new(),
            hudson_any: false,
            wc: WcAcc::default(),
            n_sites: 0,
        }
    }

    /// Absorbs one site.
    /// Absorbs one site for callers that do not track site-wide
    /// variant/invariant status (the legacy unit-test and matrix APIs).
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        self.push_with_variant(row, true);
    }

    /// Absorbs one site and counts it in `no_snps` only when it is variant
    /// (the reference convention; doc_v27 V27-27).
    pub fn push_with_variant(&mut self, row: &[SiteGroupCounts], all_variant: bool) {
        let (a, b) = (&row[self.group_a], &row[self.group_b]);
        let n_a = u64::from(a.n_ref_alleles + a.n_alt_alleles);
        let n_b = u64::from(b.n_ref_alleles + b.n_alt_alleles);

        // Audit trail (every site, skipped or not) — oracle fst.rs-aligned.
        let comps_a = n_a * n_a.saturating_sub(1) / 2;
        let comps_b = n_b * n_b.saturating_sub(1) / 2;
        self.within_a.add_partial(
            u64::from(a.n_ref_alleles) * u64::from(a.n_alt_alleles),
            comps_a,
        );
        self.within_b.add_partial(
            u64::from(b.n_ref_alleles) * u64::from(b.n_alt_alleles),
            comps_b,
        );
        let between_comps = n_a * n_b;
        self.between.add_partial(
            u64::from(a.n_ref_alleles) * u64::from(b.n_alt_alleles)
                + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles),
            between_comps,
        );

        // Hudson (reference rule): both groups need at least two called alleles for
        // BOTH the numerator and denominator. scikit-allel admits n = 1
        // into the denominator alone; that convention is deliberately not
        // used here (doc_v27 V27-03).
        let mpd_a = (n_a >= 2).then(|| {
            u64::from(a.n_ref_alleles) as f64 * u64::from(a.n_alt_alleles) as f64
                / (n_a * (n_a - 1) / 2) as f64
        });
        let mpd_b = (n_b >= 2).then(|| {
            u64::from(b.n_ref_alleles) as f64 * u64::from(b.n_alt_alleles) as f64
                / (n_b * (n_b - 1) / 2) as f64
        });
        let mpd_ab = (n_a >= 2 && n_b >= 2).then(|| {
            (u64::from(a.n_ref_alleles) * u64::from(b.n_alt_alleles)
                + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles)) as f64
                / (n_a * n_b) as f64
        });
        if let Some(pab) = mpd_ab {
            self.hudson_den.add(pab);
            self.hudson_any = true;
            if let (Some(pa), Some(pb)) = (mpd_a, mpd_b) {
                self.hudson_num.add(pab - f64::midpoint(pa, pb));
            }
        }

        // WC84: the r = 2 instance of the one core formula (doc_v24 V24-05).
        // The pair path used to carry its own copy at 2× the multi scale;
        // Fst is a ratio, so the estimate is unchanged and `components()`
        // now has one documented scale.
        let groups = [
            Wc84Group {
                an: n_a,
                refs: u64::from(a.n_ref_alleles),
                het: u64::from(a.n_het),
            },
            Wc84Group {
                an: n_b,
                refs: u64::from(b.n_ref_alleles),
                het: u64::from(b.n_het),
            },
        ];
        match wc84_site(&groups) {
            Wc84Site::Skip => self.wc.n_skipped += 1,
            Wc84Site::KeepC { c } => {
                self.wc.c.add(c);
                self.wc.n_sites += 1;
            }
            Wc84Site::Full { a, b, c } => {
                self.wc.a.add(a);
                self.wc.b.add(b);
                self.wc.c.add(c);
                self.wc.n_sites += 1;
            }
        }

        if all_variant {
            self.n_sites += 1;
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.within_a.merge(&other.within_a);
        self.within_b.merge(&other.within_b);
        self.between.merge(&other.between);
        self.hudson_num.merge(&other.hudson_num);
        self.hudson_den.merge(&other.hudson_den);
        self.hudson_any |= other.hudson_any;
        self.wc.merge(&other.wc);
        self.n_sites += other.n_sites;
    }

    /// Integer audit trail (within_a, within_b, between).
    #[must_use]
    pub const fn audit(&self) -> (&PairAcc, &PairAcc, &PairAcc) {
        (&self.within_a, &self.within_b, &self.between)
    }

    /// Hudson estimate.
    #[must_use]
    pub fn value_hudson(&self) -> Option<f64> {
        if !self.hudson_any {
            return None;
        }
        let den = self.hudson_den.value();
        (den != 0.0).then(|| self.hudson_num.value() / den)
    }

    /// WC84 estimate.
    #[must_use]
    pub fn value_wc84(&self) -> Option<f64> {
        if self.wc.n_sites == 0 {
            return None;
        }
        let (a, b, c) = (self.wc.a.value(), self.wc.b.value(), self.wc.c.value());
        let total = a + b + c;
        (total != 0.0).then(|| a / total)
    }

    /// The estimate for a kind.
    #[must_use]
    pub fn value(&self, kind: FstKind) -> Option<f64> {
        match kind {
            FstKind::Hudson => self.value_hudson(),
            FstKind::WeirCockerham84 => self.value_wc84(),
        }
    }

    /// Sites pushed (the reference `no_snps`).
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }

    /// WC84 site tally (contributing, skipped) — the pair-path counterpart
    /// of [`WcMultiAcc::site_tally`].
    #[must_use]
    pub const fn wc_tally(&self) -> (u64, u64) {
        (self.wc.n_sites, self.wc.n_skipped)
    }

    /// WC84 variance components `(a, b, c)` at the **multi-population
    /// scale** shared with [`WcMultiAcc::components`] (doc_v24 V24-05).
    /// Before v24 the pair path carried its own 2×-scaled copy; `Fst` is a
    /// ratio and is unchanged.
    #[must_use]
    pub fn wc_components(&self) -> (f64, f64, f64) {
        (self.wc.a.value(), self.wc.b.value(), self.wc.c.value())
    }
}

// ---------------------------------------------------- multi-population WC84

/// Weir-Cockerham (1984) variance components over r populations (§9.2:
/// generalizes the pair-only path above; formulas per allel's
/// `weir_cockerham_fst` — the r-population ANOVA with n_i = an_i // 2
/// floor semantics and h̄ from the count table's n_het).
///
/// Sums (a, b, c) are per-site NeumaierSums over the contributing sites.
/// A site where some (not all) groups have a called-diploid count of 0 —
/// `an_i = 0` (allel `nansum` semantics, doc_v7 D7-20) or a lone half-call
/// whose `n_i = an_i // 2` floor leaves `n_C = 0` (doc_v20 S46) —
/// contributes only the finite c component; at `r ≥ 3` such a site keeps
/// its full a/b/c contribution while `n_C > 0` (allel parity). The
/// `n̄ ≤ 1` / `n_C ≤ 0` guards still skip otherwise.
#[derive(Clone, Debug)]
pub struct WcMultiAcc {
    groups: Vec<usize>,
    /// Per-site input scratch, reused across `push` (doc_v24 V24-05/07).
    site: Vec<Wc84Group>,
    a: NeumaierSum,
    b: NeumaierSum,
    c: NeumaierSum,
    n_sites: u64,
    n_skipped: u64,
}

impl WcMultiAcc {
    /// Accumulator over the given group indices.
    #[must_use]
    pub fn new(groups: Vec<usize>) -> Self {
        Self {
            groups,
            site: Vec::new(),
            a: NeumaierSum::new(),
            b: NeumaierSum::new(),
            c: NeumaierSum::new(),
            n_sites: 0,
            n_skipped: 0,
        }
    }

    /// Absorbs one site through the same core the pair path uses.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        if self.groups.len() < 2 {
            self.n_skipped += 1;
            return;
        }
        self.site.clear();
        for &g in &self.groups {
            let c = &row[g];
            self.site.push(Wc84Group {
                an: u64::from(c.n_ref_alleles + c.n_alt_alleles),
                refs: u64::from(c.n_ref_alleles),
                het: u64::from(c.n_het),
            });
        }
        match wc84_site(&self.site) {
            Wc84Site::Skip => self.n_skipped += 1,
            Wc84Site::KeepC { c } => {
                self.c.add(c);
                self.n_sites += 1;
            }
            Wc84Site::Full { a, b, c } => {
                self.a.add(a);
                self.b.add(b);
                self.c.add(c);
                self.n_sites += 1;
            }
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.a.merge(&other.a);
        self.b.merge(&other.b);
        self.c.merge(&other.c);
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// Fst = Σa / (Σa + Σb + Σc).
    #[must_use]
    pub fn value(&self) -> Option<f64> {
        if self.n_sites == 0 {
            return None;
        }
        let (a, b, c) = (self.a.value(), self.b.value(), self.c.value());
        let total = a + b + c;
        (total != 0.0).then(|| a / total)
    }

    /// Variance-component sums (a, b, c).
    #[must_use]
    pub const fn components(&self) -> (&NeumaierSum, &NeumaierSum, &NeumaierSum) {
        (&self.a, &self.b, &self.c)
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }
}

// ------------------------------------------------------------- Fst matrix

/// One matrix entry: the pair (i, j) and its estimate.
#[derive(Clone, Debug)]
pub struct FstMatrixEntry {
    /// Row group (index into the request's group list).
    pub i: usize,
    /// Column group.
    pub j: usize,
    /// The estimate for the pair (None = no contributing site).
    pub value: Option<f64>,
}

/// The all-pairs Fst matrix (§9.2 `Fst::matrix`): one pass fills every
/// group pair of the given kind, with each population's within-pair terms
/// computed once per site (the accumulator-union shape — the O(G²) pair
/// terms are scalar work over per-group quantities already computed).
#[derive(Clone, Debug)]
pub struct FstMatrixAcc {
    n_groups: usize,
    kind: FstKind,
    /// Pair (i, j), i < j, flattened as i·n + j − i(i+1)/2 − 1.
    pairs: Vec<FstAcc>,
}

impl FstMatrixAcc {
    /// A matrix accumulator over `n_groups` groups.
    #[must_use]
    pub fn new(n_groups: usize, kind: FstKind) -> Self {
        let pairs = (0..n_groups)
            .flat_map(|i| (i + 1..n_groups).map(move |j| FstAcc::new(i, j)))
            .collect();
        Self {
            n_groups,
            kind,
            pairs,
        }
    }

    /// The estimator kind.
    #[must_use]
    pub const fn kind(&self) -> FstKind {
        self.kind
    }

    /// Absorbs one site into every pair.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        for pair in &mut self.pairs {
            pair.push(row);
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        crate::merge_zip_eq(
            "FstMatrixEntry pairs",
            &mut self.pairs,
            &other.pairs,
            |dst, src| {
                dst.merge(src);
            },
        );
    }

    /// The matrix entries in (i, j) row-major order, i < j.
    #[must_use]
    pub fn entries(&self) -> Vec<FstMatrixEntry> {
        (0..self.n_groups)
            .flat_map(|i| (i + 1..self.n_groups).map(move |j| (i, j)))
            .zip(&self.pairs)
            .map(|((i, j), acc)| FstMatrixEntry {
                i,
                j,
                value: acc.value(self.kind),
            })
            .collect()
    }
}

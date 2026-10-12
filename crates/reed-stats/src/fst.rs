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
    hudson_num_any: bool,
    wc: WcAcc,
    /// Every site pushed (the reference `no_snps` column).
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
            hudson_num_any: false,
            wc: WcAcc::default(),
            n_sites: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let (a, b) = (&row[self.group_a], &row[self.group_b]);
        let n_a = u64::from(a.n_ref_alleles + a.n_alt_alleles);
        let n_b = u64::from(b.n_ref_alleles + b.n_alt_alleles);

        // Audit trail (every site, skipped or not) — oracle fst.rs-aligned.
        let comps_a = n_a * n_a.saturating_sub(1) / 2;
        let comps_b = n_b * n_b.saturating_sub(1) / 2;
        self.within_a.add(
            u64::from(a.n_ref_alleles) * u64::from(a.n_alt_alleles),
            comps_a,
            comps_a,
        );
        self.within_b.add(
            u64::from(b.n_ref_alleles) * u64::from(b.n_alt_alleles),
            comps_b,
            comps_b,
        );
        let between_comps = n_a * n_b;
        self.between.add(
            u64::from(a.n_ref_alleles) * u64::from(b.n_alt_alleles)
                + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles),
            between_comps,
            between_comps,
        );

        // Hudson (allel inclusion rule): den whenever both groups have ≥1
        // called allele; num also requires n ≥ 2 within each group.
        let mpd_a = (n_a >= 2).then(|| {
            u64::from(a.n_ref_alleles) as f64 * u64::from(a.n_alt_alleles) as f64
                / (n_a * (n_a - 1) / 2) as f64
        });
        let mpd_b = (n_b >= 2).then(|| {
            u64::from(b.n_ref_alleles) as f64 * u64::from(b.n_alt_alleles) as f64
                / (n_b * (n_b - 1) / 2) as f64
        });
        let mpd_ab = (n_a >= 1 && n_b >= 1).then(|| {
            (u64::from(a.n_ref_alleles) * u64::from(b.n_alt_alleles)
                + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles)) as f64
                / (n_a * n_b) as f64
        });
        if let Some(pab) = mpd_ab {
            self.hudson_den.add(pab);
            self.hudson_any = true;
            if let (Some(pa), Some(pb)) = (mpd_a, mpd_b) {
                self.hudson_num.add(pab - f64::midpoint(pa, pb));
                self.hudson_num_any = true;
            }
        }

        // WC84 (r = 2 populations; n = called // 2 floor, allel fst.py:163).
        // Skipping a group with no called alleles matches vcftools
        // (--weir-fst-pop prints `-nan` for the site and drops it from the
        // weighted mean); pixy / scikit-allel instead keep the finite `c`
        // component (a paired demo with an all-missing group: reed/vcftools
        // 0.150 vs pixy/allel 0.125). Documented divergence, doc_v7 D7-20.
        if n_a == 0 || n_b == 0 {
            self.wc.n_skipped += 1;
        } else {
            let (n1, n2) = ((n_a / 2) as f64, (n_b / 2) as f64);
            let n_total = n1 + n2;
            let n_bar = n_total / 2.0;
            let n_c = (n_total - (n1 * n1 + n2 * n2) / n_total) / (2.0 - 1.0);
            if n1 == 0.0 || n2 == 0.0 || n_bar <= 1.0 || n_c <= 0.0 {
                self.wc.n_skipped += 1;
            } else {
                let p_ref_1 = u64::from(a.n_ref_alleles) as f64 / n_a as f64;
                let p_ref_2 = u64::from(b.n_ref_alleles) as f64 / n_b as f64;
                let p_bar =
                    u64::from(a.n_ref_alleles + b.n_ref_alleles) as f64 / (n_a + n_b) as f64;
                let s2 = (n1 * (p_ref_1 - p_bar).powi(2) + n2 * (p_ref_2 - p_bar).powi(2))
                    / (n_bar * (2.0 - 1.0));
                let h_bar = (u64::from(a.n_het + b.n_het) as f64) / n_total;
                let term = p_bar * (1.0 - p_bar) - (2.0 - 1.0) * s2 / 2.0 - h_bar / 4.0;
                let a_allele = (n_bar / n_c) * (s2 - term / (n_bar - 1.0));
                let b_allele = (n_bar / (n_bar - 1.0))
                    * (p_bar * (1.0 - p_bar)
                        - (2.0 - 1.0) * s2 / 2.0
                        - (2.0 * n_bar - 1.0) * h_bar / (4.0 * n_bar));
                let c_allele = h_bar / 2.0;
                self.wc.a.add(2.0 * a_allele);
                self.wc.b.add(2.0 * b_allele);
                self.wc.c.add(2.0 * c_allele);
                self.wc.n_sites += 1;
            }
        }

        self.n_sites += 1;
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.within_a.merge(&other.within_a);
        self.within_b.merge(&other.within_b);
        self.between.merge(&other.between);
        self.hudson_num.merge(&other.hudson_num);
        self.hudson_den.merge(&other.hudson_den);
        self.hudson_any |= other.hudson_any;
        self.hudson_num_any |= other.hudson_num_any;
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
        if !self.hudson_any || !self.hudson_num_any {
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
}

// ---------------------------------------------------- multi-population WC84

/// Weir-Cockerham (1984) variance components over r populations (§9.2:
/// generalizes the pair-only path above; formulas per allel's
/// `weir_cockerham_fst` — the r-population ANOVA with n_i = an_i // 2
/// floor semantics and h̄ from the count table's n_het).
///
/// Sums (a, b, c) are per-site NeumaierSums over the contributing sites
/// (all populations must have n_i ≥ 1 and the totals must define n_C > 0).
#[derive(Clone, Debug)]
pub struct WcMultiAcc {
    groups: Vec<usize>,
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
            a: NeumaierSum::new(),
            b: NeumaierSum::new(),
            c: NeumaierSum::new(),
            n_sites: 0,
            n_skipped: 0,
        }
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let r = self.groups.len();
        if r < 2 {
            self.n_skipped += 1;
            return;
        }
        let counts: Vec<&SiteGroupCounts> = self.groups.iter().map(|&g| &row[g]).collect();
        // n_i = an_i // 2 (allel floor semantics), an_i = called alleles.
        let an: Vec<u64> = counts
            .iter()
            .map(|c| u64::from(c.n_ref_alleles + c.n_alt_alleles))
            .collect();
        if an.contains(&0) {
            self.n_skipped += 1;
            return;
        }
        // No per-group `n_i == 0` skip here (only `an_i == 0` above): matches
        // vcftools, which drops a site whose WEIR_AND_COCKERHAM_FST is `-nan`.
        // NOTE (doc_v7 D7-20): this differs from pixy / scikit-allel, which
        // keep allel's finite `c` component for an all-missing group
        // (pairwise demo, one group all-missing: reed/vcftools 0.15 vs
        // pixy/allel 0.125; r=3 likewise 0.2195 vs 0.1698). reed follows the
        // vcftools/hierfstat reference the docs name for the multi-group
        // WC84, not pixy (which has no r-population merged form).
        let n: Vec<f64> = an.iter().map(|&a| (a / 2) as f64).collect();
        let n_total: f64 = n.iter().sum();
        let n_bar = n_total / r as f64;
        let n_c = (n_total - n.iter().map(|x| x * x).sum::<f64>() / n_total) / (r as f64 - 1.0);
        if n_bar <= 1.0 || n_c <= 0.0 {
            self.n_skipped += 1;
            return;
        }
        // Allele frequencies and the total (reference-allele convention,
        // like allel: p_i of the reference allele).
        let p: Vec<f64> = an
            .iter()
            .zip(&counts)
            .map(|(&a, c)| f64::from(c.n_ref_alleles) / a as f64)
            .collect();
        let ac_total: u64 = counts.iter().map(|c| u64::from(c.n_ref_alleles)).sum();
        let an_total: u64 = an.iter().sum();
        let p_bar = ac_total as f64 / an_total as f64;
        let s2 = n
            .iter()
            .zip(&p)
            .map(|(&ni, &pi)| ni * (pi - p_bar) * (pi - p_bar))
            .sum::<f64>()
            / (n_bar * (r as f64 - 1.0));
        // h̄ = Σ n_het_i / n_total (from the count table). `n_het` counts
        // genotypes with both alleles called and different (a half call is
        // not a full heterozygote) — the same rule as scikit-allel's
        // `count_het`, so h̄ matches allel (doc_v7 D7-27).
        let h_bar = counts.iter().map(|c| u64::from(c.n_het)).sum::<u64>() as f64 / n_total;
        let a = (n_bar / n_c)
            * (s2
                - (1.0 / (n_bar - 1.0))
                    * (p_bar * (1.0 - p_bar) - (r as f64 - 1.0) * s2 / r as f64 - h_bar / 4.0));
        let b = (n_bar / (n_bar - 1.0))
            * (p_bar * (1.0 - p_bar)
                - (r as f64 - 1.0) * s2 / r as f64
                - (2.0 * n_bar - 1.0) * h_bar / (4.0 * n_bar));
        let c = h_bar / 2.0;
        self.a.add(a);
        self.b.add(b);
        self.c.add(c);
        self.n_sites += 1;
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
        for (dst, src) in self.pairs.iter_mut().zip(&other.pairs) {
            dst.merge(src);
        }
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

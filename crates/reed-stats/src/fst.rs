//! Between-group Fst: Hudson / Weir-Cockerham-84 / Nei-86
//! (`docs/architecture-spec.md` §9.2), components identical to the oracle's
//! `estimators::fst` (scikit-allel semantics), including the allel inclusion
//! rule for Hudson and the `n = called // 2` floor of WC84.

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
    /// Nei 1986 GST.
    Nei86,
}

impl FstKind {
    /// Stable identifier.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::Hudson => "hudson",
            Self::WeirCockerham84 => "weir_cockerham_84",
            Self::Nei86 => "nei_86",
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

/// Nei (1986) components summed over sites.
#[derive(Clone, Debug, Default)]
pub struct NeiAcc {
    /// Within-group gene diversity summed over sites.
    pub hs: NeumaierSum,
    /// Total gene diversity summed over sites.
    pub ht: NeumaierSum,
    /// Contributing sites.
    pub n_sites: u64,
    /// Sites skipped.
    pub n_skipped: u64,
}

impl NeiAcc {
    fn merge(&mut self, other: &Self) {
        self.hs.merge(&other.hs);
        self.ht.merge(&other.ht);
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
    nei: NeiAcc,
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
            nei: NeiAcc::default(),
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
                + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles))
                as f64
                / (n_a * n_b) as f64
        });
        if let Some(pab) = mpd_ab {
            self.hudson_den.add(pab);
            self.hudson_any = true;
            if let (Some(pa), Some(pb)) = (mpd_a, mpd_b) {
                self.hudson_num.add(pab - (pa + pb) / 2.0);
                self.hudson_num_any = true;
            }
        }

        // WC84 (r = 2 populations; n = called // 2 floor, allel fst.py:163).
        // Its guards are WC84-only: Nei86 below counts the same sites the
        // oracle's independent per-estimator loop counts.
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
                let p_bar = u64::from(a.n_ref_alleles + b.n_ref_alleles) as f64 / (n_a + n_b) as f64;
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

        // Nei 86 (independent of the WC84 guards).
        if n_a == 0 || n_b == 0 {
            self.nei.n_skipped += 1;
        } else {
            let p1 = u64::from(a.n_alt_alleles) as f64 / n_a as f64;
            let p2 = u64::from(b.n_alt_alleles) as f64 / n_b as f64;
            let hs = (2.0 * p1 * (1.0 - p1) + 2.0 * p2 * (1.0 - p2)) / 2.0;
            let p_bar_alt = (p1 + p2) / 2.0;
            let ht = 2.0 * p_bar_alt * (1.0 - p_bar_alt);
            self.nei.hs.add(hs);
            self.nei.ht.add(ht);
            self.nei.n_sites += 1;
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
        self.hudson_num_any |= other.hudson_num_any;
        self.wc.merge(&other.wc);
        self.nei.merge(&other.nei);
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

    /// Nei-86 estimate.
    #[must_use]
    pub fn value_nei86(&self) -> Option<f64> {
        if self.nei.n_sites == 0 {
            return None;
        }
        let ht = self.nei.ht.value();
        (ht != 0.0).then(|| 1.0 - self.nei.hs.value() / ht)
    }

    /// The estimate for a kind.
    #[must_use]
    pub fn value(&self, kind: FstKind) -> Option<f64> {
        match kind {
            FstKind::Hudson => self.value_hudson(),
            FstKind::WeirCockerham84 => self.value_wc84(),
            FstKind::Nei86 => self.value_nei86(),
        }
    }
}

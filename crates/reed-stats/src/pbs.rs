//! Population branch statistic (PBS) among three groups
//! (`docs/architecture-spec.md` §9.3): three Hudson pairs, window Fst
//! clipped to allel's `[0, 0.99999]` before the log transform.

use reed_core::count::SiteGroupCounts;
use reed_core::math::sum::NeumaierSum;

use crate::pi::PairAcc;

/// allel's window-Fst clip (selection.py:1322): caps T at ≈ 11.5.
const FST_CLIP_MAX: f64 = 0.99999;

/// Which sealed PBS definition a value comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PbsKind {
    /// Yi et al. 2010.
    Yi2010,
    /// Normalised PBSn1 (Malaspinas 2016; scikit-allel default).
    PbsN1,
}

impl PbsKind {
    /// Stable identifier.
    #[must_use]
    pub const fn id(self) -> &'static str {
        match self {
            Self::Yi2010 => "yi_2010",
            Self::PbsN1 => "pbs_n1",
        }
    }
}

#[derive(Clone, Debug, Default)]
struct PairSums {
    num: NeumaierSum,
    den: NeumaierSum,
    den_any: bool,
    num_any: bool,
}

impl PairSums {
    fn merge(&mut self, other: &Self) {
        self.num.merge(&other.num);
        self.den.merge(&other.den);
        self.den_any |= other.den_any;
        self.num_any |= other.num_any;
    }

    /// The pair's Hudson Fst, unclipped (the same value `--fst hudson`
    /// reports; the PBS transform clips it separately).
    fn fst(&self) -> Option<f64> {
        if !self.den_any || !self.num_any {
            return None;
        }
        let den = self.den.value();
        (den != 0.0).then(|| self.num.value() / den)
    }

    /// The Fst the PBS transform uses: clipped to `[0, 0.99999]` like
    /// scikit-allel's `pbs` (04-algorithms.md §3).
    fn fst_clipped(&self) -> Option<f64> {
        self.fst().map(|f| f.clamp(0.0, FST_CLIP_MAX))
    }
}

/// PBS accumulator for the triple (a, b, c); `a` is the outgroup position in
/// `t_ab + t_ac − t_bc` (Yi 2010 convention).
#[derive(Clone, Debug)]
pub struct PbsAcc {
    group_a: usize,
    group_b: usize,
    group_c: usize,
    within_a: PairAcc,
    within_b: PairAcc,
    within_c: PairAcc,
    between_ab: PairAcc,
    between_ac: PairAcc,
    between_bc: PairAcc,
    pair_ab: PairSums,
    pair_ac: PairSums,
    pair_bc: PairSums,
    n_sites: u64,
    n_skipped: u64,
}

impl PbsAcc {
    /// Accumulator for the triple.
    #[must_use]
    pub fn new(group_a: usize, group_b: usize, group_c: usize) -> Self {
        Self {
            group_a,
            group_b,
            group_c,
            within_a: PairAcc::default(),
            within_b: PairAcc::default(),
            within_c: PairAcc::default(),
            between_ab: PairAcc::default(),
            between_ac: PairAcc::default(),
            between_bc: PairAcc::default(),
            pair_ab: PairSums::default(),
            pair_ac: PairSums::default(),
            pair_bc: PairSums::default(),
            n_sites: 0,
            n_skipped: 0,
        }
    }

    fn push_pair(pair: &mut PairSums, x: &SiteGroupCounts, y: &SiteGroupCounts) {
        let n_x = u64::from(x.n_ref_alleles + x.n_alt_alleles);
        let n_y = u64::from(y.n_ref_alleles + y.n_alt_alleles);
        let within_x = (n_x >= 2).then(|| {
            u64::from(x.n_ref_alleles) as f64 * u64::from(x.n_alt_alleles) as f64
                / (n_x * (n_x - 1) / 2) as f64
        });
        let within_y = (n_y >= 2).then(|| {
            u64::from(y.n_ref_alleles) as f64 * u64::from(y.n_alt_alleles) as f64
                / (n_y * (n_y - 1) / 2) as f64
        });
        let between = (n_x >= 1 && n_y >= 1).then(|| {
            (u64::from(x.n_ref_alleles) * u64::from(y.n_alt_alleles)
                + u64::from(x.n_alt_alleles) * u64::from(y.n_ref_alleles))
                as f64
                / (n_x * n_y) as f64
        });
        if let Some(bet) = between {
            pair.den.add(bet);
            pair.den_any = true;
            if let (Some(wx), Some(wy)) = (within_x, within_y) {
                pair.num.add(bet - (wx + wy) / 2.0);
                pair.num_any = true;
            }
        }
    }

    /// Absorbs one site: each of the three Hudson pairs decides on its own
    /// coverage (04-algorithms.md §1 retires the v0.2 three-way
    /// intersection; the pairs are independent, as in Yi et al. 2010).
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let (a, b, c) = (&row[self.group_a], &row[self.group_b], &row[self.group_c]);
        let n_a = u64::from(a.n_ref_alleles + a.n_alt_alleles);
        let n_b = u64::from(b.n_ref_alleles + b.n_alt_alleles);
        let n_c = u64::from(c.n_ref_alleles + c.n_alt_alleles);
        if n_a > 0 && n_b > 0 {
            Self::push_pair(&mut self.pair_ab, a, b);
        }
        if n_a > 0 && n_c > 0 {
            Self::push_pair(&mut self.pair_ac, a, c);
        }
        if n_b > 0 && n_c > 0 {
            Self::push_pair(&mut self.pair_bc, b, c);
        }
        self.n_sites += 1;
        if n_a == 0 || n_b == 0 || n_c == 0 {
            self.n_skipped += 1;
        }
        self.push_audit(row);
    }

    /// The pooled integer audit trail over every site (both site sets).
    fn push_audit(&mut self, row: &[SiteGroupCounts]) {
        let (a, b, c) = (&row[self.group_a], &row[self.group_b], &row[self.group_c]);
        let n_a = u64::from(a.n_ref_alleles + a.n_alt_alleles);
        let n_b = u64::from(b.n_ref_alleles + b.n_alt_alleles);
        let n_c = u64::from(c.n_ref_alleles + c.n_alt_alleles);

        // Audit trail (every site, skipped or not).
        let comps_a = n_a * n_a.saturating_sub(1) / 2;
        let comps_b = n_b * n_b.saturating_sub(1) / 2;
        let comps_c = n_c * n_c.saturating_sub(1) / 2;
        self.within_a.add(u64::from(a.n_ref_alleles) * u64::from(a.n_alt_alleles), comps_a, comps_a);
        self.within_b.add(u64::from(b.n_ref_alleles) * u64::from(b.n_alt_alleles), comps_b, comps_b);
        self.within_c.add(u64::from(c.n_ref_alleles) * u64::from(c.n_alt_alleles), comps_c, comps_c);
        let ab = n_a * n_b;
        let ac = n_a * n_c;
        let bc = n_b * n_c;
        self.between_ab.add(
            u64::from(a.n_ref_alleles) * u64::from(b.n_alt_alleles)
                + u64::from(a.n_alt_alleles) * u64::from(b.n_ref_alleles),
            ab,
            ab,
        );
        self.between_ac.add(
            u64::from(a.n_ref_alleles) * u64::from(c.n_alt_alleles)
                + u64::from(a.n_alt_alleles) * u64::from(c.n_ref_alleles),
            ac,
            ac,
        );
        self.between_bc.add(
            u64::from(b.n_ref_alleles) * u64::from(c.n_alt_alleles)
                + u64::from(b.n_alt_alleles) * u64::from(c.n_ref_alleles),
            bc,
            bc,
        );

    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.within_a.merge(&other.within_a);
        self.within_b.merge(&other.within_b);
        self.within_c.merge(&other.within_c);
        self.between_ab.merge(&other.between_ab);
        self.between_ac.merge(&other.between_ac);
        self.between_bc.merge(&other.between_bc);
        self.pair_ab.merge(&other.pair_ab);
        self.pair_ac.merge(&other.pair_ac);
        self.pair_bc.merge(&other.pair_bc);
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// The estimate for a kind; `None` when any pair denominator is empty.
    #[must_use]
    pub fn value(&self, kind: PbsKind) -> Option<f64> {
        let (Some(f_ab), Some(f_ac), Some(f_bc)) = (
            self.pair_ab.fst_clipped(),
            self.pair_ac.fst_clipped(),
            self.pair_bc.fst_clipped(),
        ) else {
            return None;
        };
        let t = |f: f64| -(1.0 - f).ln();
        let (t_ab, t_ac, t_bc) = (t(f_ab), t(f_ac), t(f_bc));
        let raw = (t_ab + t_ac - t_bc) / 2.0;
        Some(match kind {
            PbsKind::Yi2010 => raw,
            PbsKind::PbsN1 => raw / (1.0 + (t_ab + t_ac + t_bc) / 2.0),
        })
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }

    /// The three pairwise Hudson Fst values `(f_ab, f_ac, f_bc)`, unclipped
    /// (the PBS transform applies the clip itself).
    #[must_use]
    pub fn pair_fsts(&self) -> (Option<f64>, Option<f64>, Option<f64>) {
        (self.pair_ab.fst(), self.pair_ac.fst(), self.pair_bc.fst())
    }

    /// Sites absorbed (`no_snps`).
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }
}

//! Statistics requests (`docs/spec/architecture.md` §8.1–§8.2): the sealed
//! set of `(statistic, estimator, groups)` triples and the matching sealed
//! accumulator union — static dispatch, no `dyn`. [`build_accs`] validates
//! every request up front (group bounds, distinctness, chromosome counts),
//! so a rejected combination fails at build time, never mid-stream.

use reed_core::CoreError;
use reed_core::count::SiteGroupCounts;

use crate::dxy::DxyAcc;
use crate::fst::{FstAcc, FstKind, WcMultiAcc};
use crate::pbs::{PbsAcc, PbsKind};
use crate::pi::PiAcc;
use crate::sfs::SfsAcc;
use crate::tajima::TajimaAcc;
use crate::theta::ThetaAcc;

/// One requested statistic over group indices (a group index is a position
/// in the count row handed to every `push`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StatKind {
    /// π for one group; `per_bp` selects the `AccessibleBp` convention
    /// (divides by window bp) over the `PerPairComparable` ratio.
    Pi {
        /// Group index.
        group: usize,
        /// Divide by window bp instead of comparable pairs.
        per_bp: bool,
    },
    /// Absolute divergence dxy for one group pair; `per_bp` as in
    /// [`StatKind::Pi`].
    Dxy {
        /// First group index.
        a: usize,
        /// Second group index.
        b: usize,
        /// Divide by window bp instead of comparable pairs.
        per_bp: bool,
    },
    /// Fst for one group pair; one accumulator pass fills all three
    /// estimators' components, `kind` picks the reported one.
    Fst {
        /// First group index.
        a: usize,
        /// Second group index.
        b: usize,
        /// Which estimator's value to report.
        kind: FstKind,
    },
    /// Global (r-population) WC84 Fst over every group (docs/rounds/v5/05
    /// §3.2): one value, `Σa / (Σa + Σb + Σc)`. Only WC84 has a global
    /// form; Hudson is a two-population quantity.
    FstGlobal,
    /// Population branch statistic over one group triple (`a` is the branch
    /// of interest: `t_ab + t_ac − t_bc`).
    Pbs {
        /// Branch-of-interest group index.
        a: usize,
        /// Second group index.
        b: usize,
        /// Third group index.
        c: usize,
        /// Which flavour to report.
        kind: PbsKind,
    },
    /// Site frequency spectrum for one group (complete-case sites).
    Sfs {
        /// Group index.
        group: usize,
        /// Fold to the minor allele count.
        folded: bool,
    },
    /// Watterson's θ for one group; `per_bp` divides by window bp.
    ThetaW {
        /// Group index.
        group: usize,
        /// Divide by window bp.
        per_bp: bool,
    },
    /// Tajima's D for one group with a segregating-site guard.
    TajimaD {
        /// Group index.
        group: usize,
        /// Minimum segregating sites (below → `None`).
        min_sites: u64,
    },
}

/// Stable output identifier of one request, e.g. `fst.hudson.0.1`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StatKey(pub String);

impl StatKey {
    /// The identifier as a string slice.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::fmt::Display for StatKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl StatKind {
    /// The stable key this request reports under.
    #[must_use]
    pub fn key(&self) -> StatKey {
        match *self {
            Self::Pi { group, per_bp } => StatKey(format!(
                "pi.{}.{}",
                if per_bp { "bp" } else { "pair" },
                group
            )),
            Self::Dxy { a, b, per_bp } => StatKey(format!(
                "dxy.{}.{}.{}",
                if per_bp { "bp" } else { "pair" },
                a,
                b
            )),
            Self::Fst { a, b, kind } => StatKey(format!("fst.{}.{a}.{b}", kind.id())),
            Self::FstGlobal => StatKey("fst.wc84.global".to_string()),
            Self::Pbs { a, b, c, kind } => StatKey(format!("pbs.{}.{a}.{b}.{c}", kind.id())),
            Self::Sfs { group, folded } => StatKey(format!(
                "sfs.{}.{group}",
                if folded { "folded" } else { "unfolded" }
            )),
            Self::ThetaW { group, per_bp } => StatKey(format!(
                "theta_w.{}.{}",
                if per_bp { "bp" } else { "total" },
                group
            )),
            Self::TajimaD { group, .. } => StatKey(format!("tajima_d.{group}")),
        }
    }

    /// Validates the request against the group sizes and builds its
    /// accumulator.
    ///
    /// # Errors
    ///
    /// [`CoreError::LengthMismatch`] when a group index is out of bounds,
    /// [`CoreError::InvalidWindow`] when an estimator's chromosome
    /// requirement is not met (e.g. SFS/θ need `n_chrom ≥ 2`) or groups are
    /// not distinct where required.
    pub fn build(&self, n_group: &[u32]) -> Result<StatAcc, CoreError> {
        let chroms = |g: usize| -> Result<u64, CoreError> {
            let n = n_group.get(g).copied().ok_or(CoreError::LengthMismatch {
                expected: g + 1,
                got: n_group.len(),
            })?;
            Ok(u64::from(n) * 2)
        };
        let distinct = |gs: &[usize]| -> Result<(), CoreError> {
            for (i, &x) in gs.iter().enumerate() {
                for &y in &gs[i + 1..] {
                    if x == y {
                        return Err(CoreError::InvalidWindow);
                    }
                }
            }
            Ok(())
        };
        match *self {
            Self::Pi { group, .. } => {
                let n_group = u32::try_from(chroms(group)? / 2).unwrap_or(u32::MAX);
                Ok(StatAcc::Pi(PiAcc::new(group, n_group)))
            }
            Self::Dxy { a, b, .. } => {
                distinct(&[a, b])?;
                let (na, nb) = (chroms(a)?, chroms(b)?);
                let (na, nb) = (
                    u32::try_from(na / 2).unwrap_or(u32::MAX),
                    u32::try_from(nb / 2).unwrap_or(u32::MAX),
                );
                Ok(StatAcc::Dxy(DxyAcc::new(a, b, na, nb)))
            }
            Self::Fst { a, b, .. } => {
                distinct(&[a, b])?;
                Ok(StatAcc::Fst(FstAcc::new(a, b)))
            }
            Self::FstGlobal => {
                if n_group.len() < 2 {
                    return Err(CoreError::LengthMismatch {
                        expected: 2,
                        got: n_group.len(),
                    });
                }
                Ok(StatAcc::FstGlobal(WcMultiAcc::new(
                    (0..n_group.len()).collect(),
                )))
            }
            Self::Pbs { a, b, c, .. } => {
                distinct(&[a, b, c])?;
                Ok(StatAcc::Pbs(PbsAcc::new(a, b, c)))
            }
            Self::Sfs { group, folded } => {
                let n = chroms(group)?;
                let n_group = u32::try_from(n / 2).unwrap_or(u32::MAX);
                Ok(StatAcc::Sfs(SfsAcc::new(group, n_group, folded)?))
            }
            Self::ThetaW { group, .. } => {
                let n = chroms(group)?;
                let n_group = u32::try_from(n / 2).unwrap_or(u32::MAX);
                let acc = ThetaAcc::new(group, n_group).ok_or(CoreError::InvalidWindow)?;
                Ok(StatAcc::Theta(acc))
            }
            Self::TajimaD { group, .. } => {
                let n = chroms(group)?;
                let n_group = u32::try_from(n / 2).unwrap_or(u32::MAX);
                let acc = TajimaAcc::new(group, n_group).ok_or(CoreError::InvalidWindow)?;
                Ok(StatAcc::Tajima(acc))
            }
        }
    }
}

/// Builds one accumulator per request, in request order.
///
/// # Errors
///
/// The first failing request (see [`StatKind::build`]).
pub fn build_accs(kinds: &[StatKind], n_group: &[u32]) -> Result<Vec<StatAcc>, CoreError> {
    kinds.iter().map(|kind| kind.build(n_group)).collect()
}

/// The sealed accumulator union — one variant per [`StatKind`] family.
#[derive(Clone, Debug)]
pub enum StatAcc {
    /// π accumulator.
    Pi(PiAcc),
    /// dxy accumulator.
    Dxy(DxyAcc),
    /// Fst accumulator (all three estimators in one pass).
    Fst(FstAcc),
    /// PBS accumulator (both flavours in one pass).
    Pbs(PbsAcc),
    /// SFS accumulator.
    Sfs(SfsAcc),
    /// Global (r-population) WC84 accumulator.
    FstGlobal(WcMultiAcc),
    /// Watterson θ accumulator.
    Theta(ThetaAcc),
    /// Tajima's D accumulator.
    Tajima(TajimaAcc),
}

/// Reference `_pbs.txt` columns: `(value, no_snps, f_ab, f_ac, f_bc)`.
pub type PbsColumns = (Option<f64>, u64, (Option<f64>, Option<f64>, Option<f64>));

/// Reference `_tajima_d.txt` components: `(raw_pi, raw_theta, stdev, no_sites)`.
pub type TajimaColumns = (Option<f64>, Option<f64>, Option<f64>, u64);

impl StatAcc {
    /// Absorbs one site's count row; `all_variant` is the site-wide
    /// "any ALT called" flag (variant/invariant site tallies).
    pub fn push(&mut self, row: &[SiteGroupCounts], all_variant: bool) {
        match self {
            Self::Pi(acc) => acc.push(row, all_variant),
            Self::Dxy(acc) => acc.push(row, all_variant),
            Self::Fst(acc) => acc.push(row),
            Self::FstGlobal(acc) => acc.push(row),
            Self::Pbs(acc) => acc.push(row),
            Self::Sfs(acc) => acc.push(row),
            Self::Theta(acc) => acc.push(row),
            Self::Tajima(acc) => acc.push(row),
        }
    }

    /// Registers a present-but-uncountable site (a skipped multiallelic
    /// record): only the variant-only pi/dxy denominator reads it, every
    /// other statistic ignores it (doc_v7 D7-11).
    pub fn note_skipped(&mut self) {
        match self {
            Self::Pi(acc) => acc.note_skipped_site(),
            Self::Dxy(acc) => acc.note_skipped_site(),
            _ => {}
        }
    }

    /// Component-wise merge (associative; exact integer components).
    ///
    /// # Panics
    ///
    /// Debug builds panic when the two variants differ (a caller bug: parts
    /// of one request are never merged across requests).
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(
            std::mem::discriminant(self),
            std::mem::discriminant(other),
            "invariant: merge combines identical accumulator kinds"
        );
        match (self, other) {
            (Self::Pi(a), Self::Pi(b)) => a.merge(b),
            (Self::Dxy(a), Self::Dxy(b)) => a.merge(b),
            (Self::Fst(a), Self::Fst(b)) => a.merge(b),
            (Self::FstGlobal(a), Self::FstGlobal(b)) => a.merge(b),
            (Self::Pbs(a), Self::Pbs(b)) => a.merge(b),
            (Self::Sfs(a), Self::Sfs(b)) => a.merge(b),
            (Self::Theta(a), Self::Theta(b)) => a.merge(b),
            (Self::Tajima(a), Self::Tajima(b)) => a.merge(b),
            _ => {}
        }
    }

    /// The scalar estimate this accumulator holds for its request kind;
    /// `None` for an empty/skipped window. The SFS has no scalar value.
    #[must_use]
    pub fn value(&self, kind: &StatKind, window_bp: u64) -> Option<f64> {
        debug_assert_eq!(self.kind_id(), kind_key_family(kind));
        match (self, kind) {
            (Self::Pi(acc), StatKind::Pi { per_bp: false, .. }) => acc.value_per_pair(),
            (Self::Pi(acc), StatKind::Pi { per_bp: true, .. }) => {
                acc.value_accessible_bp(window_bp)
            }
            (Self::Dxy(acc), StatKind::Dxy { per_bp: false, .. }) => acc.value_per_pair(),
            (Self::Dxy(acc), StatKind::Dxy { per_bp: true, .. }) => {
                acc.value_accessible_bp(window_bp)
            }
            (Self::Fst(acc), StatKind::Fst { kind, .. }) => acc.value(*kind),
            (Self::FstGlobal(acc), StatKind::FstGlobal) => acc.value(),
            (Self::Pbs(acc), StatKind::Pbs { kind, .. }) => acc.value(*kind),
            (Self::Theta(acc), StatKind::ThetaW { per_bp: false, .. }) => acc.value_avg(),
            (Self::Theta(acc), StatKind::ThetaW { per_bp: true, .. }) => {
                acc.value_per_bp(window_bp)
            }
            (Self::Tajima(acc), StatKind::TajimaD { min_sites, .. }) => acc.value(*min_sites),
            _ => None,
        }
    }

    /// π/dxy/θ with the `--variant-only` denominator completion
    /// (04-algorithms.md §3); ratio statistics (Fst, PBS, Tajima's D) are
    /// unaffected.
    ///
    /// `per_bp` is intentionally ignored here: the design fixes one
    /// variant-only estimator per statistic (the `absent`-completion
    /// denominator changes, the point estimate does not), and the CLI only
    /// builds `per_bp: false` kinds for this route, so the `pi.bp` /
    /// `theta_w.bp` columns never carry a per-pair quantity in practice
    /// (doc_v7 D7-23).
    #[must_use]
    pub fn value_variant_only(&self, kind: &StatKind, window_bp: u64) -> Option<f64> {
        match (self, kind) {
            (Self::Pi(acc), StatKind::Pi { .. }) => acc.value_variant_only(window_bp),
            (Self::Dxy(acc), StatKind::Dxy { .. }) => acc.value_variant_only(window_bp),
            (Self::Theta(acc), StatKind::ThetaW { .. }) => acc.value_variant_only(window_bp),
            _ => self.value(kind, window_bp),
        }
    }

    /// Reference `_watterson_theta.txt` columns:
    /// `(raw_theta, no_sites, no_var_sites, weighted_no_sites)`.
    #[must_use]
    pub fn theta_columns(&self) -> Option<(f64, u64, u64, f64)> {
        match self {
            Self::Theta(acc) => Some((
                acc.raw_theta(),
                acc.n_sites(),
                acc.n_var_sites(),
                acc.weighted_no_sites(),
            )),
            _ => None,
        }
    }

    /// `PREFIX_pbs.txt` columns: `(value, no_snps, f_ab, f_ac, f_bc)`.
    #[must_use]
    pub fn pbs_columns(&self, kind: &StatKind) -> Option<PbsColumns> {
        match (self, kind) {
            (Self::Pbs(acc), StatKind::Pbs { kind, .. }) => {
                Some((acc.value(*kind), acc.n_sites(), acc.pair_fsts()))
            }
            _ => None,
        }
    }

    /// Reference `_fst.txt` columns: `(value, no_snps)`.
    #[must_use]
    pub fn fst_columns(&self, kind: &StatKind) -> Option<(Option<f64>, u64)> {
        match (self, kind) {
            (Self::Fst(acc), StatKind::Fst { kind, .. }) => Some((acc.value(*kind), acc.n_sites())),
            _ => None,
        }
    }

    /// Global-Fst columns: `(value, no_snps)` (`PREFIX_fst_global.txt`).
    #[must_use]
    pub fn fst_global_columns(&self, kind: &StatKind) -> Option<(Option<f64>, u64)> {
        match (self, kind) {
            (Self::FstGlobal(acc), StatKind::FstGlobal) => Some((acc.value(), acc.site_tally().0)),
            _ => None,
        }
    }

    /// Reference `_tajima_d.txt` components: `(raw_pi, raw_theta, stdev,
    /// no_sites)`. The values are `None` where the denominator is undefined
    /// (the reference prints `NA` there); the site count is always present.
    #[must_use]
    pub fn tajima_columns(&self) -> Option<TajimaColumns> {
        match self {
            Self::Tajima(acc) => Some((
                Some(acc.raw_pi()),
                Some(acc.raw_theta()),
                acc.d_stdev(),
                acc.site_tally().0,
            )),
            _ => None,
        }
    }

    /// Sites present in the window for π/dxy (`|S|`); `None` for others.
    #[must_use]
    pub const fn n_present(&self) -> Option<u64> {
        match self {
            Self::Pi(acc) => Some(acc.n_present()),
            Self::Dxy(acc) => Some(acc.n_present()),
            _ => None,
        }
    }

    /// Integer audit components `(diffs, comparable, missing)` for π/dxy.
    #[must_use]
    pub const fn pair_counts(&self) -> Option<(u64, u64, u64)> {
        match self {
            Self::Pi(acc) => {
                let c = acc.counts();
                Some((c.diffs, c.comparable, c.missing))
            }
            Self::Dxy(acc) => {
                let c = acc.counts();
                Some((c.diffs, c.comparable, c.missing))
            }
            _ => None,
        }
    }

    /// The SFS bins when this is the SFS accumulator.
    #[must_use]
    pub fn sfs_bins(&self) -> Option<&[u64]> {
        match self {
            Self::Sfs(acc) => Some(acc.bins()),
            _ => None,
        }
    }

    fn kind_id(&self) -> &'static str {
        match self {
            Self::Pi(_) => "pi",
            Self::Dxy(_) => "dxy",
            Self::Fst(_) | Self::FstGlobal(_) => "fst",
            Self::Pbs(_) => "pbs",
            Self::Sfs(_) => "sfs",
            Self::Theta(_) => "theta",
            Self::Tajima(_) => "tajima",
        }
    }
}

fn kind_key_family(kind: &StatKind) -> &'static str {
    match kind {
        StatKind::Pi { .. } => "pi",
        StatKind::Dxy { .. } => "dxy",
        StatKind::Fst { .. } | StatKind::FstGlobal => "fst",
        StatKind::Pbs { .. } => "pbs",
        StatKind::Sfs { .. } => "sfs",
        StatKind::ThetaW { .. } => "theta",
        StatKind::TajimaD { .. } => "tajima",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keys_are_stable_and_distinct() {
        let kinds = [
            StatKind::Pi {
                group: 0,
                per_bp: false,
            },
            StatKind::Pi {
                group: 0,
                per_bp: true,
            },
            StatKind::Dxy {
                a: 0,
                b: 1,
                per_bp: false,
            },
            StatKind::Fst {
                a: 0,
                b: 1,
                kind: FstKind::Hudson,
            },
            StatKind::Fst {
                a: 0,
                b: 1,
                kind: FstKind::WeirCockerham84,
            },
            StatKind::Pbs {
                a: 0,
                b: 1,
                c: 2,
                kind: PbsKind::PbsN1,
            },
            StatKind::Sfs {
                group: 1,
                folded: true,
            },
            StatKind::ThetaW {
                group: 0,
                per_bp: false,
            },
            StatKind::TajimaD {
                group: 0,
                min_sites: 3,
            },
        ];
        let keys: Vec<_> = kinds.iter().map(StatKind::key).collect();
        for (i, a) in keys.iter().enumerate() {
            for b in &keys[i + 1..] {
                assert_ne!(a, b, "keys must be distinct");
            }
        }
        assert_eq!(keys[3].as_str(), "fst.hudson.0.1");
    }

    #[test]
    fn build_validates_bounds_and_distinctness() {
        let ok = [StatKind::Fst {
            a: 0,
            b: 1,
            kind: FstKind::Hudson,
        }];
        assert!(build_accs(&ok, &[3, 2]).is_ok());
        let oob = [StatKind::Pi {
            group: 2,
            per_bp: false,
        }];
        assert!(build_accs(&oob, &[3, 2]).is_err());
        let same = [StatKind::Fst {
            a: 1,
            b: 1,
            kind: FstKind::Hudson,
        }];
        assert!(build_accs(&same, &[3, 2]).is_err());
        let pbs = [StatKind::Pbs {
            a: 0,
            b: 1,
            c: 1,
            kind: PbsKind::Yi2010,
        }];
        assert!(build_accs(&pbs, &[3, 2, 2]).is_err());
        // SFS with zero chromosomes has no spectrum.
        let sfs = [StatKind::Sfs {
            group: 0,
            folded: false,
        }];
        assert!(build_accs(&sfs, &[0]).is_err());
    }
}

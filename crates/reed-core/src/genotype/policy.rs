//! Explicit input policies (`docs/architecture-spec.md` §5.4).
//!
//! Every policy exists so that reed never silently drops, keeps or rewrites
//! input: an unsupported genotype either errors (the conservative default)
//! or is transformed exactly as the declared policy prescribes, with every
//! trigger counted into a [`PolicyCounters`] that lands in the run report.

/// How non-diploid genotypes are handled during decoding.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PloidyPolicy {
    /// Reject any haploid (or higher-ploid) genotype. Default.
    #[default]
    RequireDiploid,
    /// A haploid genotype becomes a half call: one called slot, one missing
    /// slot (then [`HalfCallPolicy`] applies).
    HaploidAsHalfCall,
    /// A haploid genotype is duplicated into both slots.
    HaploidAsHomozygous,
}

/// How allele indices above 1 (multi-allelic sites) are handled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MultiallelicPolicy {
    /// Reject the record. Default.
    #[default]
    Reject,
    /// Any ALT index counts as the (single) ALT allele.
    CollapseAlt,
    /// Keep per-ALT planes (one plane set per ALT allele; not yet implemented,
    /// see `reed-stats` multi-allelic support).
    SplitPerAlt,
}

/// How half calls (`0/.`, `./1`) are handled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HalfCallPolicy {
    /// Keep the called slot as called and the missing slot as missing
    /// (per-allele deletion semantics). Default.
    #[default]
    AlleleLevel,
    /// Turn the whole genotype into a missing call.
    TreatAsMissing,
}

/// The bundle of decoding policies a pass runs under.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GenotypePolicies {
    /// Non-diploid handling.
    pub ploidy: PloidyPolicy,
    /// Multi-allelic handling.
    pub multiallelic: MultiallelicPolicy,
    /// Half-call handling.
    pub half_call: HalfCallPolicy,
}

/// Counts of every policy trigger, accumulated across the run.
///
/// All fields are exact integers and the struct satisfies the accumulator
/// laws (merge is addition), so per-chunk partials sum exactly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PolicyCounters {
    /// Haploid genotypes converted to half calls.
    pub n_haploid_as_half_call: u64,
    /// Haploid genotypes duplicated into homozygous calls.
    pub n_haploid_as_homozygous: u64,
    /// Haploid genotypes zeroed to missing.
    pub n_haploid_as_missing: u64,
    /// Half calls (including haploid-as-half-call) zeroed to missing.
    pub n_half_call_as_missing: u64,
    /// ALT alleles with index > 1 collapsed onto the single ALT bit.
    pub n_multiallelic_collapsed: u64,
}

impl PolicyCounters {
    /// Merges another counter block into this one (component-wise addition).
    pub fn merge(&mut self, other: &Self) {
        self.n_haploid_as_half_call += other.n_haploid_as_half_call;
        self.n_haploid_as_homozygous += other.n_haploid_as_homozygous;
        self.n_haploid_as_missing += other.n_haploid_as_missing;
        self.n_half_call_as_missing += other.n_half_call_as_missing;
        self.n_multiallelic_collapsed += other.n_multiallelic_collapsed;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_conservative() {
        let p = GenotypePolicies::default();
        assert_eq!(p.ploidy, PloidyPolicy::RequireDiploid);
        assert_eq!(p.multiallelic, MultiallelicPolicy::Reject);
        assert_eq!(p.half_call, HalfCallPolicy::AlleleLevel);
    }

    #[test]
    fn counters_merge_additively() {
        let mut a = PolicyCounters {
            n_haploid_as_half_call: 2,
            ..PolicyCounters::default()
        };
        let b = PolicyCounters {
            n_haploid_as_half_call: 3,
            n_multiallelic_collapsed: 7,
            ..PolicyCounters::default()
        };
        a.merge(&b);
        assert_eq!(a.n_haploid_as_half_call, 5);
        assert_eq!(a.n_multiallelic_collapsed, 7);
    }
}

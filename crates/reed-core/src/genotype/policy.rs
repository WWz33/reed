//! Explicit input policies (`docs/spec/architecture.md` §5.4).
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
///
/// The field list is the canonical trigger table (doc_v25 V25-07c): the BCF
/// fast table's trigger bits and the filter's report both derive from
/// [`Self::FIELDS`], so adding a trigger touches this file and the enum that
/// produces it — not three crates. `haploid_as_missing` was removed: no
/// `PloidyPolicy` variant produces it, so the column was zero forever.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub struct PolicyCounters {
    /// Haploid genotypes converted to half calls.
    #[serde(rename = "haploid_as_half_call")]
    pub n_haploid_as_half_call: u64,
    /// Haploid genotypes duplicated into homozygous calls.
    #[serde(rename = "haploid_as_homozygous")]
    pub n_haploid_as_homozygous: u64,
    /// Half calls (including haploid-as-half-call) zeroed to missing.
    #[serde(rename = "half_call_as_missing")]
    pub n_half_call_as_missing: u64,
    /// ALT alleles with index > 1 collapsed onto the single ALT bit.
    #[serde(rename = "multiallelic_collapsed")]
    pub n_multiallelic_collapsed: u64,
}

/// One row of the canonical policy-trigger table.
pub struct PolicyCounterField {
    /// Report name (also the trigger's identity).
    pub name: &'static str,
    /// Read accessor.
    pub get: fn(&PolicyCounters) -> u64,
    /// Write accessor (the BCF fast table adds one trigger per entry).
    pub get_mut: fn(&mut PolicyCounters) -> &mut u64,
}

impl PolicyCounters {
    /// The canonical table, in trigger-bit order: bit `i` of the BCF fast
    /// table's trigger byte means `FIELDS[i] > 0`.
    pub const FIELDS: [PolicyCounterField; 4] = [
        PolicyCounterField {
            name: "haploid_as_half_call",
            get: |c| c.n_haploid_as_half_call,
            get_mut: |c| &mut c.n_haploid_as_half_call,
        },
        PolicyCounterField {
            name: "haploid_as_homozygous",
            get: |c| c.n_haploid_as_homozygous,
            get_mut: |c| &mut c.n_haploid_as_homozygous,
        },
        PolicyCounterField {
            name: "half_call_as_missing",
            get: |c| c.n_half_call_as_missing,
            get_mut: |c| &mut c.n_half_call_as_missing,
        },
        PolicyCounterField {
            name: "multiallelic_collapsed",
            get: |c| c.n_multiallelic_collapsed,
            get_mut: |c| &mut c.n_multiallelic_collapsed,
        },
    ];

    /// Merges another counter block into this one (component-wise addition).
    pub fn merge(&mut self, other: &Self) {
        for field in &Self::FIELDS {
            *(field.get_mut)(self) += (field.get)(other);
        }
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

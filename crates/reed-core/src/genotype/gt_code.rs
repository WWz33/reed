//! [`GtCode`]: one byte per sample, the single genotype representation every
//! consumer decodes to (`docs/spec/architecture.md` §3.3 of docs/rounds/v2/03).
//!
//! Layout: the low four bits are `slot0_state | slot1_state << 2` with
//! 0 = REF, 1 = ALT, 2 = missing (9 of the 16 states occur); bit 4 is the
//! phasing bit and bit 5 records that the input text was haploid. Consumers
//! that need bit planes convert once with [`PlaneBlock::push_codes`].
//!
//! The genotype policies (ploidy, multi-allelic, half calls) are executed in
//! exactly one place — [`GtCode::from_alleles`] — with every trigger counted
//! into a [`PolicyCounters`].

use crate::error::CoreError;
use crate::genotype::{
    GenotypePolicies, HalfCallPolicy, MultiallelicPolicy, Ploidy, PloidyPolicy, PolicyCounters,
};

/// Slot state: REF allele.
pub const STATE_REF: u8 = 0;
/// Slot state: ALT allele (any index above 0 under `CollapseAlt`).
pub const STATE_ALT: u8 = 1;
/// Slot state: missing allele (`.`).
pub const STATE_MISSING: u8 = 2;

/// Bit 4: the genotype was phased (`|`).
const PHASED: u8 = 0x10;
/// Bit 5: the genotype text held one allele (haploid).
const HAPLOID: u8 = 0x20;

/// One sample's decoded genotype.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(transparent)]
pub struct GtCode(u8);

impl GtCode {
    /// `./.`: both slots missing.
    pub const MISSING: Self = Self(STATE_MISSING | STATE_MISSING << 2);

    /// Builds a code from raw bits (decoders and tests).
    #[must_use]
    pub const fn from_raw(raw: u8) -> Self {
        Self(raw)
    }

    /// The raw byte.
    #[must_use]
    pub const fn raw(self) -> u8 {
        self.0
    }

    /// Builds a code from two slot states.
    #[must_use]
    pub const fn from_states(slot0: u8, slot1: u8, phased: bool, haploid: bool) -> Self {
        Self(slot0 | slot1 << 2 | if phased { PHASED } else { 0 } | if haploid { HAPLOID } else { 0 })
    }

    /// The low four bits: `slot0_state | slot1_state << 2`.
    #[must_use]
    pub const fn state(self) -> u8 {
        self.0 & 0x0F
    }

    /// Whether the input genotype was phased.
    #[must_use]
    pub const fn phased(self) -> bool {
        self.0 & PHASED != 0
    }

    /// Whether the input genotype text held a single allele.
    #[must_use]
    pub const fn is_haploid(self) -> bool {
        self.0 & HAPLOID != 0
    }

    /// The state of one slot (0/1/2).
    #[must_use]
    pub const fn slot_state(self, slot: usize) -> u8 {
        (self.0 >> (2 * slot)) & 0b11
    }

    /// Whether one slot carries the ALT allele.
    #[must_use]
    pub const fn slot_alt(self, slot: usize) -> bool {
        self.slot_state(slot) == STATE_ALT
    }

    /// Whether one slot is missing.
    #[must_use]
    pub const fn slot_missing(self, slot: usize) -> bool {
        self.slot_state(slot) == STATE_MISSING
    }

    /// Both slots missing (`./.`).
    #[must_use]
    pub const fn is_fully_missing(self) -> bool {
        self.state() == Self::MISSING.0
    }

    /// At least one slot called.
    #[must_use]
    pub const fn is_called_any(self) -> bool {
        !self.is_fully_missing()
    }

    /// The masked form: both slots missing, the haploid flag preserved
    /// (a masked haploid writes `.`, not `./.`).
    #[must_use]
    pub const fn masked(self) -> Self {
        Self(Self::MISSING.0 | self.0 & HAPLOID)
    }

    /// Decodes one genotype against the declared policies.
    ///
    /// `slots` holds the allele index per slot (`None` = `.`), `ploidy` what
    /// the input text showed, `phased` whether it used `|`. Multi-allelic
    /// allele indices above 1 go through [`MultiallelicPolicy`]; haploid
    /// genotypes through [`PloidyPolicy`]; one-called-slot genotypes through
    /// [`HalfCallPolicy`]. Every transformation is counted.
    ///
    /// # Errors
    ///
    /// [`CoreError::PolicyRejected`] when the declared policy forbids the
    /// genotype, [`CoreError::Unimplemented`] for `SplitPerAlt`.
    pub fn from_alleles(
        slots: [Option<u16>; 2],
        ploidy: Ploidy,
        phased: bool,
        policies: &GenotypePolicies,
        counters: &mut PolicyCounters,
    ) -> Result<Self, CoreError> {
        match ploidy {
            Ploidy::Haploid => {
                let Some(allele) = slots[0] else {
                    // A haploid `.` is a fully missing call; no policy involved.
                    return Ok(Self::from_states(
                        STATE_MISSING,
                        STATE_MISSING,
                        phased,
                        true,
                    ));
                };
                match policies.ploidy {
                    PloidyPolicy::RequireDiploid => Err(CoreError::PolicyRejected(
                        "haploid genotype under PloidyPolicy::RequireDiploid",
                    )),
                    PloidyPolicy::HaploidAsHomozygous => {
                        counters.n_haploid_as_homozygous += 1;
                        let state = allele_state(allele, policies, counters)?;
                        Ok(Self::from_states(state, state, phased, true))
                    }
                    PloidyPolicy::HaploidAsHalfCall => {
                        counters.n_haploid_as_half_call += 1;
                        half_call(0, allele, phased, true, policies, counters)
                    }
                }
            }
            Ploidy::Diploid => match (slots[0], slots[1]) {
                (None, None) => Ok(Self::from_states(STATE_MISSING, STATE_MISSING, phased, false)),
                (Some(allele), None) => half_call(0, allele, phased, false, policies, counters),
                (None, Some(allele)) => half_call(1, allele, phased, false, policies, counters),
                (Some(a), Some(b)) => {
                    let slot0 = allele_state(a, policies, counters)?;
                    let slot1 = allele_state(b, policies, counters)?;
                    Ok(Self::from_states(slot0, slot1, phased, false))
                }
            },
        }
    }
}

/// One genotype with exactly one called slot (haploid-as-half-call included).
fn half_call(
    called_slot: usize,
    allele: u16,
    phased: bool,
    haploid: bool,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
) -> Result<GtCode, CoreError> {
    match policies.half_call {
        HalfCallPolicy::TreatAsMissing => {
            counters.n_half_call_as_missing += 1;
            Ok(GtCode::from_states(
                STATE_MISSING,
                STATE_MISSING,
                phased,
                haploid,
            ))
        }
        HalfCallPolicy::AlleleLevel => {
            let state = allele_state(allele, policies, counters)?;
            let (slot0, slot1) = if called_slot == 0 {
                (state, STATE_MISSING)
            } else {
                (STATE_MISSING, state)
            };
            Ok(GtCode::from_states(slot0, slot1, phased, haploid))
        }
    }
}

/// Resolves one called allele index to its slot state (`0` = REF, `1` = ALT,
/// `>= 2` through the multi-allelic policy).
fn allele_state(
    allele: u16,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
) -> Result<u8, CoreError> {
    match allele {
        0 => Ok(STATE_REF),
        1 => Ok(STATE_ALT),
        _ => match policies.multiallelic {
            MultiallelicPolicy::Reject => Err(CoreError::PolicyRejected(
                "allele index above 1 under MultiallelicPolicy::Reject",
            )),
            MultiallelicPolicy::CollapseAlt => {
                counters.n_multiallelic_collapsed += 1;
                Ok(STATE_ALT)
            }
            MultiallelicPolicy::SplitPerAlt => Err(CoreError::Unimplemented(
                "MultiallelicPolicy::SplitPerAlt plane sets",
            )),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn states_and_bits_round_trip() {
        let code = GtCode::from_states(STATE_ALT, STATE_MISSING, true, false);
        assert_eq!(code.state(), STATE_ALT | STATE_MISSING << 2);
        assert!(code.slot_alt(0) && code.slot_missing(1));
        assert!(code.phased() && !code.is_haploid());
        assert!(!code.is_fully_missing() && code.is_called_any());
        assert_eq!(GtCode::MISSING.state(), 0b1010);
        assert!(GtCode::MISSING.is_fully_missing());
    }

    #[test]
    fn masking_keeps_only_the_haploid_bit() {
        let call = GtCode::from_states(STATE_ALT, STATE_REF, true, true);
        let masked = call.masked();
        assert!(masked.is_fully_missing());
        assert!(masked.is_haploid());
        assert!(!masked.phased());
        assert_eq!(masked.state(), GtCode::MISSING.state());
    }

    fn policies() -> GenotypePolicies {
        GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        }
    }

    #[test]
    fn diploid_full_calls_decode() {
        let mut c = PolicyCounters::default();
        let p = policies();
        let code = GtCode::from_alleles([Some(0), Some(1)], Ploidy::Diploid, true, &p, &mut c).unwrap();
        assert_eq!(code.state(), STATE_REF | STATE_ALT << 2);
        assert!(code.phased());
        assert_eq!(c, PolicyCounters::default());
    }

    #[test]
    fn half_call_and_haploid_policies() {
        let mut c = PolicyCounters::default();
        let p = policies();
        let half = GtCode::from_alleles([None, Some(0)], Ploidy::Diploid, false, &p, &mut c).unwrap();
        assert_eq!(half.state(), STATE_MISSING | STATE_REF << 2);
        let hap = GtCode::from_alleles([Some(1), None], Ploidy::Haploid, false, &p, &mut c).unwrap();
        assert_eq!(hap.state(), STATE_ALT | STATE_MISSING << 2);
        assert!(hap.is_haploid());
        assert_eq!(c.n_haploid_as_half_call, 1);

        // Default rejects a called haploid.
        let mut c2 = PolicyCounters::default();
        let err = GtCode::from_alleles(
            [Some(1), None],
            Ploidy::Haploid,
            false,
            &GenotypePolicies::default(),
            &mut c2,
        )
        .unwrap_err();
        assert!(matches!(err, CoreError::PolicyRejected(_)));

        // A haploid `.` is plain missing, whatever the policy says.
        let mut c3 = PolicyCounters::default();
        let dot = GtCode::from_alleles(
            [None, None],
            Ploidy::Haploid,
            false,
            &GenotypePolicies::default(),
            &mut c3,
        )
        .unwrap();
        assert!(dot.is_fully_missing() && dot.is_haploid());
        assert_eq!(c3, PolicyCounters::default());
    }

    #[test]
    fn multiallelic_policies_count_and_reject() {
        let mut collapse = PolicyCounters::default();
        let collapsed =
            GtCode::from_alleles([Some(2), Some(3)], Ploidy::Diploid, false, &policies(), &mut collapse).unwrap();
        assert_eq!(collapsed.state(), STATE_ALT | STATE_ALT << 2);
        assert_eq!(collapse.n_multiallelic_collapsed, 2);

        let mut reject = PolicyCounters::default();
        let err = GtCode::from_alleles(
            [Some(0), Some(2)],
            Ploidy::Diploid,
            false,
            &GenotypePolicies::default(),
            &mut reject,
        )
        .unwrap_err();
        assert!(matches!(err, CoreError::PolicyRejected(_)));
    }

    #[test]
    fn half_call_as_missing_counts() {
        let policies = GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: MultiallelicPolicy::CollapseAlt,
            half_call: HalfCallPolicy::TreatAsMissing,
        };
        let mut c = PolicyCounters::default();
        let code = GtCode::from_alleles([Some(1), None], Ploidy::Diploid, false, &policies, &mut c).unwrap();
        assert!(code.is_fully_missing());
        assert_eq!(c.n_half_call_as_missing, 1);
    }
}

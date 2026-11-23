//! INFO/AC + INFO/AN: the one validation of the cohort allele counts.
//!
//! htslib's `bcf_calc_ac` (vcfutils.c:56-64) is the contract every consumer
//! shares: AC carries one value per ALT allele, AN is present and positive,
//! and ΣAC must not exceed AN. The REF count is `AN − ΣAC`. Before
//! doc_v30 R02 the filter's reader and `stats --summary` each implemented
//! this rule, with different duplicate-tag and overflow behavior: the same
//! record could yield two different AC values in one command.

/// Per-allele AC counts with REF filled in, plus AN.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InfoCounts {
    /// REF first, then one count per ALT: `n_alt + 1` entries.
    pub per_allele: Vec<u64>,
    /// INFO/AN (always > 0; a zero AN is [`InfoCountsOutcome::Absent`]).
    pub an: u64,
}

/// What [`InfoCounts::from_info`] found, keeping the malformed cases apart.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum InfoCountsOutcome {
    /// AC and AN are present and consistent.
    Usable(InfoCounts),
    /// AC or AN is absent, unparsable, or AN is zero: the caller falls back
    /// to the genotype counts, which is what bcftools does when no INFO
    /// value carries the cohort.
    Absent,
    /// AC is present with the wrong number of values.
    AcLength {
        /// Values AC carried (a prefix count when the list is oversized).
        got: usize,
        /// ALT alleles the record declares.
        want: u32,
    },
    /// ΣAC > AN. htslib prints `[E::bcf_calc_ac] Incorrect AN/AC counts`
    /// and aborts the run; callers decide whether that is fatal.
    AcExceedsAn {
        /// Sum of the AC values.
        alt_sum: u64,
        /// INFO/AN.
        an: u64,
    },
    /// The AC values sum past `u64::MAX` (doc_v26 V26-02): a malformed INFO
    /// value, never usable counts. The old `sum()` panicked in debug and
    /// wrapped in release.
    AcOverflow {
        /// Values AC carried.
        got: usize,
    },
}

impl InfoCounts {
    /// Validates the raw `INFO/AC` and `INFO/AN` **values** (not the whole
    /// INFO column: tag lookup belongs to the caller's index, so this module
    /// does not grow a second INFO parser).
    ///
    /// `n_alt` is the record's ALT count; `AC=.` or a non-integer element
    /// yields [`InfoCountsOutcome::Absent`], matching bcftools' missing-value
    /// path rather than a length error.
    #[must_use]
    pub fn from_info(ac: Option<&[u8]>, an: Option<&[u8]>, n_alt: u32) -> InfoCountsOutcome {
        let (Some(ac), Some(an)) = (ac, an) else {
            return InfoCountsOutcome::Absent;
        };
        let Ok(ac) = std::str::from_utf8(ac) else {
            return InfoCountsOutcome::Absent;
        };
        let Some(an) = std::str::from_utf8(an)
            .ok()
            .and_then(|text| text.parse::<u64>().ok())
            .filter(|&an| an > 0)
        else {
            return InfoCountsOutcome::Absent;
        };

        // One pass (doc_v30 §4.7): parse, check the length, and sum as we go.
        // The old shape counted first (to avoid a proportional allocation,
        // V26-02) and then re-parsed, which left a `return Absent` in the
        // second loop that could never fire. The bound is enforced here
        // before each push: once `got` exceeds `n_alt` the list is already
        // too long, so no further element is examined or stored.
        let mut per_allele = Vec::with_capacity((n_alt as usize).saturating_add(1));
        let mut got = 0_usize;
        let mut alt_sum = 0_u64;
        for part in ac.split(',') {
            // One-based count of the element now being handled: `AcLength`
            // reports it as "seen", and `AcOverflow` reports how many were
            // parsed, matching the old two-pass shape's values.
            got += 1;
            if got > n_alt as usize {
                return InfoCountsOutcome::AcLength { got, want: n_alt };
            }
            let Ok(count) = part.parse::<u64>() else {
                return InfoCountsOutcome::Absent;
            };
            let Some(next) = alt_sum.checked_add(count) else {
                return InfoCountsOutcome::AcOverflow { got };
            };
            alt_sum = next;
            per_allele.push(count);
        }
        if got != n_alt as usize {
            return InfoCountsOutcome::AcLength { got, want: n_alt };
        }
        if alt_sum > an {
            return InfoCountsOutcome::AcExceedsAn { alt_sum, an };
        }
        per_allele.insert(0, an - alt_sum);
        InfoCountsOutcome::Usable(InfoCounts { per_allele, an })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn usable_counts_fill_ref_from_an() {
        let InfoCountsOutcome::Usable(counts) = InfoCounts::from_info(Some(b"1,2"), Some(b"10"), 2)
        else {
            panic!("consistent AC/AN must be usable");
        };
        assert_eq!(counts.per_allele, [7, 1, 2]);
        assert_eq!(counts.an, 10);
    }

    #[test]
    fn absent_cases_fall_back() {
        for (ac, an) in [
            (None, Some(&b"10"[..])),
            (Some(&b"1"[..]), None),
            (Some(&b"1"[..]), Some(&b"0"[..])),
            (Some(&b"."[..]), Some(&b"10"[..])),
            (Some(&b"x"[..]), Some(&b"10"[..])),
            (Some(&b"1"[..]), Some(&b"x"[..])),
        ] {
            assert_eq!(
                InfoCounts::from_info(ac, an, 1),
                InfoCountsOutcome::Absent,
                "{ac:?}/{an:?}"
            );
        }
    }

    #[test]
    fn ac_length_mismatch_is_kept_apart() {
        assert_eq!(
            InfoCounts::from_info(Some(b"1,2"), Some(b"3"), 1),
            InfoCountsOutcome::AcLength { got: 2, want: 1 }
        );
        // The oversized list is cut off before it is materialised.
        assert_eq!(
            InfoCounts::from_info(Some(b"1,2,3,4,5"), Some(b"0"), 1),
            InfoCountsOutcome::Absent,
            "AN=0 wins over the length check"
        );
        assert_eq!(
            InfoCounts::from_info(Some(b"1,2,3,4,5"), Some(b"100"), 1),
            InfoCountsOutcome::AcLength { got: 2, want: 1 }
        );
    }

    /// doc_v30 §4.7: the single pass stops at the first excess element, so a
    /// pathologically long AC field is never materialised. The old shape
    /// parsed every element twice; this pins the "first excess element wins"
    /// boundary and that the allocation stays bounded by `n_alt + 1`.
    #[test]
    fn a_pathological_ac_field_is_cut_off_at_the_first_excess_element() {
        // 10,000 elements against n_alt = 1: the loop must return after the
        // second element rather than scanning the whole field.
        let long = vec!["1"; 10_000].join(",");
        assert_eq!(
            InfoCounts::from_info(Some(long.as_bytes()), Some(b"100000"), 1),
            InfoCountsOutcome::AcLength { got: 2, want: 1 }
        );
        // The valid case stays exact for a multi-allelic record.
        assert_eq!(
            InfoCounts::from_info(Some(b"1,2"), Some(b"10"), 2),
            InfoCountsOutcome::Usable(InfoCounts {
                per_allele: vec![7, 1, 2],
                an: 10
            })
        );
    }

    #[test]
    fn sum_above_an_is_an_error() {
        assert_eq!(
            InfoCounts::from_info(Some(b"200"), Some(b"100"), 1),
            InfoCountsOutcome::AcExceedsAn {
                alt_sum: 200,
                an: 100
            }
        );
    }

    #[test]
    fn sum_overflow_is_reported_not_wrapped() {
        assert_eq!(
            InfoCounts::from_info(
                Some(b"18446744073709551615,1"),
                Some(b"18446744073709551615"),
                2
            ),
            InfoCountsOutcome::AcOverflow { got: 2 }
        );
    }
}

//! Genotype (GT) text scanning for the chunk fast path: per sample column,
//! alleles and genotype classes are counted straight from the GT bytes —
//! no noodles record is built (`docs/architecture-spec.md` §6.4). The
//! tallies map 1:1 onto `reed_core::count::SiteGroupCounts` fields.

/// Per-sample tallies of one GT subfield (all fields saturate into the
/// matching `SiteGroupCounts` field at the group level).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GtTally {
    /// Called REF alleles.
    pub n_ref: u32,
    /// Called ALT alleles (any nonzero index).
    pub n_alt: u32,
    /// Full calls with two REF alleles.
    pub n_hom_ref: u32,
    /// Full calls with two different alleles.
    pub n_het: u32,
    /// Full calls with two identical ALT alleles.
    pub n_hom_alt: u32,
    /// Calls with exactly one called allele (the called allele is counted
    /// in `n_ref`/`n_alt`).
    pub n_half: u32,
}

/// Tallies one sample column's GT subfield (text up to the first `:`;
/// trailing `DP:GQ…` subfields are ignored).
///
/// Returns `None` for polyploid GT (more than two alleles) — the count
/// kernel is diploid-only, and the caller must surface that as an error
/// rather than silently miscount.
#[must_use]
pub fn tally_gt(sample: &[u8]) -> Option<GtTally> {
    let gt_end = memchr::memchr(b':', sample).unwrap_or(sample.len());
    let gt = &sample[..gt_end];
    let mut alleles = [None::<u32>; 2];
    let mut n_alleles = 0_usize;
    for field in gt.split(|&b| b == b'/' || b == b'|') {
        if n_alleles >= 2 {
            return None;
        }
        if field == b"." || field.is_empty() {
            n_alleles += 1;
            continue;
        }
        let mut value: u32 = 0;
        for &b in field {
            let digit = (b as char).to_digit(10)?;
            value = value.checked_mul(10)?.checked_add(digit)?;
        }
        alleles[n_alleles] = Some(value);
        n_alleles += 1;
    }
    let mut t = GtTally::default();
    match (alleles[0], alleles[1]) {
        (None, None) => {}
        (Some(x), None) | (None, Some(x)) => {
            t.n_half += 1;
            if x == 0 {
                t.n_ref += 1;
            } else {
                t.n_alt += 1;
            }
        }
        (Some(x), Some(y)) => {
            t.n_ref += u32::from(x == 0) + u32::from(y == 0);
            t.n_alt += u32::from(x != 0) + u32::from(y != 0);
            if x == 0 && y == 0 {
                t.n_hom_ref += 1;
            } else if x != y {
                t.n_het += 1;
            } else {
                t.n_hom_alt += 1;
            }
        }
    }
    Some(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tally(gt: &[u8]) -> (u32, u32, u32, u32, u32, u32) {
        let t = tally_gt(gt).expect("invariant: diploid GT");
        (t.n_ref, t.n_alt, t.n_hom_ref, t.n_het, t.n_hom_alt, t.n_half)
    }

    #[test]
    fn classifies_diploid_genotypes() {
        assert_eq!(tally(b"0/0"), (2, 0, 1, 0, 0, 0));
        assert_eq!(tally(b"0/1"), (1, 1, 0, 1, 0, 0));
        assert_eq!(tally(b"1|0"), (1, 1, 0, 1, 0, 0));
        assert_eq!(tally(b"1/1"), (0, 2, 0, 0, 1, 0));
        assert_eq!(tally(b"./."), (0, 0, 0, 0, 0, 0));
        assert_eq!(tally(b"0/."), (1, 0, 0, 0, 0, 1));
        assert_eq!(tally(b"./1"), (0, 1, 0, 0, 0, 1));
    }

    #[test]
    fn reads_first_subfield_only() {
        assert_eq!(tally(b"0/1:99:60"), (1, 1, 0, 1, 0, 0));
        assert_eq!(tally(b"1/1"), (0, 2, 0, 0, 1, 0));
    }

    #[test]
    fn haploid_counts_one_allele() {
        assert_eq!(tally(b"1"), (0, 1, 0, 0, 0, 1));
        assert_eq!(tally(b"0"), (1, 0, 0, 0, 0, 1));
        assert_eq!(tally(b"."), (0, 0, 0, 0, 0, 0));
    }

    #[test]
    fn polyploid_and_garbage_are_none() {
        assert_eq!(tally_gt(b"0/0/1"), None);
        assert_eq!(tally_gt(b"x/1"), None);
    }
}

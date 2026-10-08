//! The SketchThenEmit sketch (§11.4): pass0 stores, per head-surviving
//! record, the raw 4-plane genotype bitmap over samples plus one firing
//! bitmap per mask rule — `(4 + R) × W` words per site (per-rule bitmaps
//! because the cohort keep mask is only known after pass0: the mid stage
//! intersects them with `keep` to attribute masks exactly like a
//! SinglePass recount). The mid stage derives every post-mask count from
//! `keep ∧ planes` with word operations only, so the decision equals the
//! inline recount bit for bit while pass1 re-parses survivors only.

// Casts in this file are bounded by construction: Record counts fit `usize` on every supported target.
#![allow(clippy::cast_possible_truncation)]

/// Bit `r` of a head bitmap: record `r` of the fragment passed all head
/// rules.
#[derive(Clone, Debug, Default)]
pub struct SketchFragment {
    /// Records covered by this fragment (head-passing or not).
    pub n_records: u64,
    /// Bit per record (see above).
    pub head_bits: Vec<u64>,
    /// Words per sample bitmap, `ceil(n_samples / 64)`.
    pub w: usize,
    /// Head-surviving records (sites); `data` holds `n_sites × (4 + R) ×
    /// w` words, stride `[p0_alt, p0_miss, p1_alt, p1_miss, rule_0, …,
    /// rule_{R-1}]`.
    pub n_sites: usize,
    /// Number of mask rules (`R` in the stride).
    pub n_rules: usize,
    pub data: Vec<u64>,
}

impl SketchFragment {
    /// An empty fragment for `w` words per bitmap and `n_rules` mask rules.
    #[must_use]
    pub fn new(w: usize, n_rules: usize) -> Self {
        Self { n_records: 0, head_bits: Vec::new(), w, n_sites: 0, n_rules, data: Vec::new() }
    }

    /// Marks record `r`'s head outcome and advances the record counter.
    pub fn push_record(&mut self, passed_head: bool) {
        let r = self.n_records as usize;
        let word = r / 64;
        if word >= self.head_bits.len() {
            self.head_bits.resize(word + 1, 0);
        }
        if passed_head {
            self.head_bits[word] |= 1 << (r % 64);
        }
        self.n_records += 1;
    }

    /// Whether record `r` passed the head rules.
    #[must_use]
    pub fn head_passed(&self, r: usize) -> bool {
        self.head_bits.get(r / 64).is_some_and(|word| word >> (r % 64) & 1 == 1)
    }

    /// Appends another fragment's records and sites (bit offsets shift by
    /// `self.n_records`; the sink uses this to close a chunk with its seam
    /// record).
    pub fn append(&mut self, other: &Self) {
        debug_assert_eq!(self.w, other.w, "invariant: same bitmap width");
        debug_assert_eq!(self.n_rules, other.n_rules, "invariant: same rule count");
        for r in 0..other.n_records {
            let passed = other.head_passed(r as usize);
            self.push_record(passed);
        }
        self.data.extend_from_slice(&other.data);
        self.n_sites += other.n_sites;
    }

    /// Appends one head-surviving site's bitmaps.
    ///
    /// `planes` is `[p0_alt, p0_miss, p1_alt, p1_miss]`, each `w` words;
    /// `rule_masks` holds one firing bitmap per mask rule (all samples —
    /// the mid stage intersects with the keep mask).
    pub fn push_site(&mut self, planes: &[&[u64]], rule_masks: &[&[u64]]) {
        debug_assert_eq!(planes.len(), 4, "invariant: 4 genotype planes");
        debug_assert_eq!(rule_masks.len(), self.n_rules, "invariant: one bitmap per rule");
        self.data.reserve((4 + self.n_rules) * self.w);
        for plane in planes {
            self.data.extend_from_slice(plane);
        }
        for rule in rule_masks {
            self.data.extend_from_slice(rule);
        }
        self.n_sites += 1;
    }

    /// One genotype plane of site `s` (`plane` ∈ 0..4).
    #[must_use]
    pub fn plane_of(&self, s: usize, plane: usize) -> &[u64] {
        let base = s * (4 + self.n_rules) * self.w + self.w * plane;
        &self.data[base..base + self.w]
    }

    /// The firing bitmap of mask rule `r` at site `s`.
    #[must_use]
    pub fn rule_mask_of(&self, s: usize, r: usize) -> &[u64] {
        let base = s * (4 + self.n_rules) * self.w + self.w * (4 + r);
        &self.data[base..base + self.w]
    }

    /// The effective mask of site `s`: the OR of every rule's firing
    /// bitmap (the mid stage intersects with `keep`).
    #[must_use]
    pub fn mask_of(&self, s: usize) -> Vec<u64> {
        let mut out = vec![0_u64; self.w];
        for r in 0..self.n_rules {
            for (dst, src) in out.iter_mut().zip(self.rule_mask_of(s, r)) {
                *dst |= src;
            }
        }
        out
    }

    /// Estimated in-memory size in bytes (the §11.4 budget check).
    #[must_use]
    pub fn size_bytes(&self) -> u64 {
        u64::try_from(
            self.data.len() * 8
                + self.head_bits.len() * 8
                + std::mem::size_of_val(&self.n_records)
                + std::mem::size_of_val(&self.n_sites),
        )
        .unwrap_or(u64::MAX)
    }

}

/// Post-mask per-site counts derived from the sketch planes
/// (`docs/architecture-spec.md` §11.4 mid stage). Identical by
/// construction to a SinglePass recount over `keep ∧ ¬mask` samples.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SketchCounts {
    /// ALT allele count.
    pub ac: u64,
    /// Called allele number.
    pub an: u64,
    /// Samples with data (half-calls included).
    pub ns: u64,
    /// Kept samples without data (masked or raw-missing).
    pub n_missing: u64,
    /// Full calls, both alleles REF.
    pub hom_ref: u64,
    /// Full calls, exactly one ALT.
    pub het: u64,
    /// Full calls, both ALT.
    pub hom_alt: u64,
    /// Called-to-missing transitions this site's mask performed.
    pub masked_gt: u64,
}

/// Counts one sketch site over the kept samples (`keep`: one word per
/// sample bitmap word; a sample is counted when its keep bit is set and
/// its mask bit is clear — masked samples are missing by definition).
#[must_use]
pub fn count_site_sketch(fragment: &SketchFragment, s: usize, keep: &[u64]) -> SketchCounts {
    let w = fragment.w;
    let mask = &fragment.mask_of(s);
    let p0a = fragment.plane_of(s, 0);
    let p0m = fragment.plane_of(s, 1);
    let p1a = fragment.plane_of(s, 2);
    let p1m = fragment.plane_of(s, 3);
    let mut c = SketchCounts::default();

    let mut kept_raw_total = 0_u64;
    for i in 0..w {
        let keep_word = keep.get(i).copied().unwrap_or(0);
        let mask_word = mask.get(i).copied().unwrap_or(0);
        let p0a = p0a.get(i).copied().unwrap_or(0);
        let p0m = p0m.get(i).copied().unwrap_or(0);
        let p1a = p1a.get(i).copied().unwrap_or(0);
        let p1m = p1m.get(i).copied().unwrap_or(0);

        let kept_raw = keep_word & !mask_word;
        kept_raw_total += kept_raw.count_ones() as u64;

        // ALT alleles over called slots of kept, unmasked samples.
        c.ac += (kept_raw & p0a).count_ones() as u64;
        c.ac += (kept_raw & p1a).count_ones() as u64;
        // Called slots: everything but the missing bit.
        c.an += 2 * (kept_raw.count_ones() as u64)
            - (kept_raw & p0m).count_ones() as u64
            - (kept_raw & p1m).count_ones() as u64;
        // Samples with data: at least one called slot.
        c.ns += (kept_raw & !(p0m & p1m)).count_ones() as u64;

        // Full calls: both slots called.
        let full = kept_raw & !p0m & !p1m;
        c.het += (full & (p0a ^ p1a)).count_ones() as u64;
        c.hom_alt += (full & p0a & p1a).count_ones() as u64;
        c.hom_ref += (full & !p0a & !p1a).count_ones() as u64;

        // Masked samples that were called before the mask: the called→
        // missing transitions.
        c.masked_gt += (keep_word & mask_word & !(p0m & p1m)).count_ones() as u64;
    }
    let kept_total: u32 = keep.iter().map(|word| word.count_ones()).sum();
    c.n_missing = u64::from(kept_total) - c.ns;
    c
}

/// Per-rule mask attribution for site `s`: `hits[r]` = kept samples whose
/// first firing rule is `r`… no — the raw kept firing count of rule `r`
/// (the caller folds them in funnel order; a sample may fire several
/// rules but SinglePass attributes the FIRST one — mirrored by folding
/// `keep ∧ rule_r ∧ ¬(rule_0 ∨ … ∨ rule_{r-1})`).
#[must_use]
pub fn rule_hits_of(fragment: &SketchFragment, s: usize, keep: &[u64]) -> Vec<u64> {
    let w = fragment.w;
    let mut or_so_far = vec![0_u64; w];
    let mut hits = vec![0_u64; fragment.n_rules];
    for (r, hit) in hits.iter_mut().enumerate() {
        let rule = fragment.rule_mask_of(s, r);
        for (i, or_word) in or_so_far.iter_mut().enumerate().take(w) {
            let keep_word = keep.get(i).copied().unwrap_or(0);
            let first_hit = keep_word & rule[i] & !*or_word;
            *hit += first_hit.count_ones() as u64;
            *or_word |= rule[i];
        }
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 4 samples, w = 1, one mask rule. s0 = 0/1, s1 = 1/1, s2 = ./., s3 = 0/.
    fn fragment() -> SketchFragment {
        let mut f = SketchFragment::new(1, 1);
        f.push_record(true);
        // rule 0 fires on sample 2 (mask).
        let rule0 = [0b0100_u64];
        let p0a = [0b0010_u64]; // s0 slot0=0, s1 slot0=1, s2 ., s3 slot0=0
        let p0m = [0b0100_u64];
        let p1a = [0b0011_u64]; // s0 slot1=1, s1 slot1=1, s2 ., s3 slot1 missing
        let p1m = [0b1100_u64];
        f.push_site(&[&p0a, &p0m, &p1a, &p1m], &[&rule0]);
        f
    }

    #[test]
    fn counts_match_hand_computation() {
        let f = fragment();
        // keep all 4 samples.
        let keep = [0b1111_u64];
        let c = count_site_sketch(&f, 0, &keep);
        // s0 0/1: ac 1, an 2, het; s1 1/1: ac 2, an 2, hom_alt; s2 masked →
        // missing; s3 0/. half call: ac 0, an 1, ns.
        assert_eq!(c.ac, 3);
        assert_eq!(c.an, 5);
        assert_eq!(c.ns, 3);
        assert_eq!(c.n_missing, 1);
        assert_eq!(c.het, 1);
        assert_eq!(c.hom_alt, 1);
        assert_eq!(c.hom_ref, 0);
        // masked_gt: only s2 was masked and s2 was raw-missing (./.) → the
        // called→missing transition count is 0.
        assert_eq!(c.masked_gt, 0);
    }

    #[test]
    fn keep_restricts_the_denominator() {
        let f = fragment();
        // keep s1 and s3 only.
        let keep = [0b1010_u64];
        let c = count_site_sketch(&f, 0, &keep);
        // s1 1/1 (hom_alt), s3 0/. (half call): ac 2, an 3, ns 2, missing 0.
        assert_eq!(c.ac, 2);
        assert_eq!(c.an, 3);
        assert_eq!(c.ns, 2);
        assert_eq!(c.n_missing, 0);
        assert_eq!(c.hom_alt, 1);
    }

    #[test]
    fn masked_called_sample_counts_as_transition() {
        let mut f = SketchFragment::new(1, 1);
        f.push_record(true);
        // mask sample 0 (0/1 → called): called→missing transition.
        // Samples 1–3 are raw-missing (./.).
        let rule0 = [0b0001_u64];
        let p0a = [0b0001_u64];
        let p0m = [0b1110_u64];
        let p1a = [0b0001_u64];
        let p1m = [0b1110_u64];
        f.push_site(&[&p0a, &p0m, &p1a, &p1m], &[&rule0]);
        let keep = [0b1111_u64];
        let c = count_site_sketch(&f, 0, &keep);
        assert_eq!(c.masked_gt, 1);
        assert_eq!(c.n_missing, 4);
        assert_eq!(c.ns, 0);
    }

    #[test]
    fn head_bits_round_trip() {
        let mut f = SketchFragment::new(1, 0);
        for passed in [true, false, true, true, false] {
            f.push_record(passed);
        }
        assert_eq!(f.n_records, 5);
        assert!(f.head_passed(0));
        assert!(!f.head_passed(1));
        assert!(f.head_passed(2));
        assert!(f.head_passed(3));
        assert!(!f.head_passed(4));
    }
}

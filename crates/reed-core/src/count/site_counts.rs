//! Per-site exact allele counting for biallelic sites.
//!
//! This is the canonical population-genetics counting primitive shared by the
//! filter and stats layers: one definition of AC/AN/NS, the genotype classes
//! and missingness, with floats derived from exact integer counts only
//! afterwards (the reed rule). The semantics are extracted verbatim from
//! reed-filter's validated implementation (subset parity PASS against the
//! three-tool gold standard), so both layers count identically by
//! construction.
//!
//! These are **site-local** counts: every returned value depends only on the
//! genotype vector handed in, so there is no cross-site denominator and no
//! manifest gate (§4.2, site-local statistics exemption) — sites absent from
//! the input are simply not counted.

use crate::math::hwe::hwe_exact_p;

/// Exact per-site allele and genotype counts over a sample subset.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SiteCounts {
    /// ALT allele count (a half-call counts its present allele).
    pub ac: u64,
    /// Total called alleles (2 for a full call, 1 for a half-call; samples
    /// outside the retain mask excluded).
    pub an: u64,
    /// Samples with data (half-calls included).
    pub ns: u64,
    /// Fully missing (`./.`) samples among retained samples.
    pub n_missing: u64,
    /// Samples homozygous REF (full calls only).
    pub hom_ref: u64,
    /// Samples heterozygous (full calls only).
    pub het: u64,
    /// Samples homozygous ALT (full calls only).
    pub hom_alt: u64,
    /// Half-calls (`0/.`, `./1`; haploid-as-half-call included): one called
    /// allele, already counted in [`SiteCounts::ns`]/`an`/`ac`, never in a
    /// genotype class. bcftools' `bcf_gt_type` reports these as UNKN, so
    /// the `-g/--genotype` rule reads them as missing (`miss` sees
    /// `n_missing + n_half_called`).
    pub n_half_called: u64,
    /// The record carried a GT field at all. `false` when the FORMAT lacks
    /// the GT key (or the record has no FORMAT): every sample then decodes
    /// as missing, but that is **not** the same as a GT whose values are
    /// all `./.` — bcftools skips `-g/--genotype` entirely on such records
    /// (`bcf_get_fmt` finds no GT), and so does the filter's tail stage.
    pub has_gt: bool,
}

impl SiteCounts {
    /// Minor allele frequency, `min(AC, AN-AC)/AN` (bcftools fill-tags MAF).
    ///
    /// 0 when AN is 0 or when every called allele sits on one side:
    /// fill-tags gives MAF 0 when either the nref or the nalt side is empty.
    //
    // AC/AN are per-site called-allele counts (≪ 2^53); the ratio cast is exact.
    #[allow(clippy::cast_precision_loss)]
    pub fn maf(&self) -> f64 {
        if self.an == 0 {
            return 0.0;
        }
        let ac = self.ac.min(self.an - self.ac);
        ac as f64 / self.an as f64
    }

    /// HWE input sentinel: fill-tags reports p = 1 ("no evidence against
    /// HWE") when either allele side is absent, so detect that before the
    /// exact test. True when AN is 0, or all called alleles are REF, or all
    /// are ALT.
    pub fn hwe_is_degenerate(&self) -> bool {
        let nref = self.an.saturating_sub(self.ac);
        self.an == 0 || self.ac == 0 || nref == 0
    }

    /// The exact-test inputs `(n_called, n_minor, n_het)` behind
    /// [`Self::hwe_p`], `None` for the degenerate sites that report the
    /// fill-tags sentinel `p = 1`. Cached consumers (the `--fill-tags HWE`
    /// fill, the `--hwe` tail rule, the `--hardy` report) feed these to
    /// [`crate::math::hwe::HweCache`] instead of paying the full
    /// distribution per site (A6-16).
    #[must_use]
    pub fn hwe_test_input(&self) -> Option<(u32, u32, u32)> {
        if self.hwe_is_degenerate() {
            return None;
        }
        let n_called = self.hom_ref + self.het + self.hom_alt;
        // Minor allele count = 2·(rarer homozygote) + hets.
        let minor_hom = self.hom_ref.min(self.hom_alt);
        let mac = 2 * minor_hom + self.het;
        let cast = |v: u64| u32::try_from(v).unwrap_or(u32::MAX);
        Some((cast(n_called), cast(mac), cast(self.het)))
    }

    /// Per-site exact HWE p-value, paired with
    /// [`crate::math::hwe::hwe_exact_p`]
    /// (Wigginton, Cutler & Abecasis 2005) over the genotype classes.
    ///
    /// Degenerate sites (see [`SiteCounts::hwe_is_degenerate`]) return 1.0 —
    /// the fill-tags sentinel — which also covers the empty site, where the
    /// raw test would return NaN. Half-calls contribute to AC/AN but not to
    /// the genotype classes, so they never enter the test.
    pub fn hwe_p(&self) -> f64 {
        match self.hwe_test_input() {
            None => 1.0,
            Some((n_called, mac, het)) => hwe_exact_p(n_called, mac, het),
        }
    }
}

/// Counts one biallelic site from a flat genotype vector: two alleles per
/// sample in sample order, `-1` = missing allele, any allele value `> 0`
/// counts as ALT.
///
/// `retain` is the per-sample keep mask (`None` = keep all); a sample past
/// the end of a too-short mask counts as retained, matching the filter
/// layer's keep-table default. Dropped samples are never counted: they add
/// nothing to any denominator (their output column is removed too).
///
/// Half-calls (`0/.`, `./1`, ...) count the single present allele toward
/// AN/AC and the sample toward NS, never toward a genotype class and never
/// toward [`SiteCounts::n_missing`] — they land in
/// [`SiteCounts::n_half_called`]. A trailing odd element of `gt`
/// (no partner allele) is ignored.
pub fn count_site(gt: &[i8], retain: Option<&[bool]>) -> SiteCounts {
    let mut c = SiteCounts::default();
    for i in 0..gt.len() / 2 {
        let retained = retain.and_then(|keep| keep.get(i).copied()).unwrap_or(true);
        if !retained {
            continue; // dropped sample: never counted, never in the denominator
        }
        let a1 = gt[2 * i];
        let a2 = gt[2 * i + 1];
        match (a1 >= 0, a2 >= 0) {
            (true, true) => {
                c.ns += 1;
                c.an += 2;
                let n_alt = u64::from(a1 > 0) + u64::from(a2 > 0);
                c.ac += n_alt;
                match n_alt {
                    0 => c.hom_ref += 1,
                    1 => c.het += 1,
                    _ => c.hom_alt += 1,
                }
            }
            (true, false) | (false, true) => {
                // Half-call: fill-tags counts the one present allele (a
                // hemizygous call), so it feeds AN/AC/NS plus the
                // half-call tally, but no genotype class.
                c.ns += 1;
                c.an += 1;
                c.ac += u64::from(a1.max(a2) > 0);
                c.n_half_called += 1;
            }
            (false, false) => c.n_missing += 1,
        }
    }
    // Counting a GT vector is counting a record that has one.
    c.has_gt = true;
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Mirrors reed-filter's `counts_full` reference example (the validated
    /// three-tool parity implementation): s0 = 0/1, s1 = 1/1, s2 = ./.
    #[test]
    fn mixed_genotypes_with_missing() {
        let gt = [0, 1, 1, 1, -1, -1];
        let c = count_site(&gt, None);
        assert_eq!((c.ac, c.an, c.ns), (3, 4, 2));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 1, 1));
        assert_eq!(c.n_missing, 1);
        // MAF = min(3, 1)/4 = 0.25.
        assert!((c.maf() - 0.25).abs() < 1e-12);
    }

    /// Mirrors reed-filter's `counts_with_sample_drop`: dropping a sample
    /// removes it from every denominator, including the missing count.
    #[test]
    fn retained_subset_excludes_dropped_samples() {
        let gt = [0, 1, 1, 1, -1, -1];
        let c = count_site(&gt, Some(&[true, true, false]));
        assert_eq!((c.ac, c.an, c.ns), (3, 4, 2));
        assert_eq!(c.n_missing, 0);

        // Dropping a called sample instead: only s1 = 1/1 remains.
        let c = count_site(&gt, Some(&[false, true, false]));
        assert_eq!((c.ac, c.an, c.ns), (2, 2, 1));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 0, 1));
        assert_eq!(c.n_missing, 0);
    }

    /// A mask shorter than the sample count keeps the unlisted samples
    /// (the filter layer's keep-table default).
    #[test]
    fn short_retain_mask_defaults_to_keep() {
        let gt = [0, 1, 1, 1, -1, -1];
        let c = count_site(&gt, Some(&[false]));
        assert_eq!((c.ac, c.an, c.ns), (2, 2, 1));
        assert_eq!(c.n_missing, 1); // s2 is kept and still missing
    }

    /// End-state of reed-filter's `mask_by_dp` reference example: after a
    /// DP rule masks s0, the site is 1/1 plus two missing samples.
    #[test]
    fn masked_site_end_state() {
        let gt = [-1, -1, 1, 1, -1, -1];
        let c = count_site(&gt, None);
        assert_eq!((c.ac, c.an, c.ns), (2, 2, 1));
        assert_eq!(c.n_missing, 2);
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 0, 1));
    }

    /// Half-calls: one present allele feeds AN/AC/NS, never a genotype
    /// class, never `n_missing` — the half-call tally holds them (the
    /// `-g/--genotype` rule reads it as missing, bcftools `bcf_gt_type`
    /// UNKN). Mirrors reed-filter's hemizygous handling.
    #[test]
    fn half_calls_count_one_allele() {
        // 0/. contributes one REF allele; ./1 contributes one ALT allele.
        let c = count_site(&[0, -1, -1, 1], None);
        assert_eq!((c.ac, c.an, c.ns), (1, 2, 2));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 0, 0));
        assert_eq!(c.n_missing, 0);
        assert_eq!(c.n_half_called, 2);
        assert!((c.maf() - 0.5).abs() < 1e-12);

        // The allele side is irrelevant to the AC contribution: 1/. and ./1.
        let c = count_site(&[1, -1, -1, 1], None);
        assert_eq!((c.ac, c.an), (2, 2));
        assert_eq!(c.n_half_called, 2);

        // ./0: present REF allele, no ALT contribution.
        let c = count_site(&[-1, 0], None);
        assert_eq!((c.ac, c.an, c.ns), (0, 1, 1));
        assert_eq!(c.n_half_called, 1);
        assert_eq!(c.maf(), 0.0);
    }

    /// MAF boundaries: 0 on an empty/fully-missing site and on fixed sites
    /// (either side empty), 0.5 at parity.
    #[test]
    fn maf_boundaries() {
        // AN = 0: no called alleles at all.
        assert_eq!(count_site(&[], None).maf(), 0.0);
        assert_eq!(count_site(&[-1, -1, -1, -1], None).maf(), 0.0);
        // Fixed REF (AC = 0) and fixed ALT (AC = AN): one side is empty.
        assert_eq!(count_site(&[0, 0, 0, 0], None).maf(), 0.0);
        assert_eq!(count_site(&[1, 1, 1, 1], None).maf(), 0.0);
        // Parity: min(2, 2)/4.
        let c = count_site(&[0, 1, 1, 0], None);
        assert!((c.maf() - 0.5).abs() < 1e-12);
        // The filter layer's an = 10, ac = 5 boundary case.
        let c = SiteCounts {
            ac: 5,
            an: 10,
            ..SiteCounts::default()
        };
        assert!((c.maf() - 0.5).abs() < 1e-12);
    }

    /// HWE degenerate sentinels, mirroring the filter layer: degenerate
    /// exactly when AN = 0 or an allele side is empty.
    #[test]
    fn hwe_degenerate_sentinels() {
        let mut c = SiteCounts {
            ac: 0,
            an: 10,
            ..SiteCounts::default()
        };
        assert!(c.hwe_is_degenerate());
        c.ac = 10;
        assert!(c.hwe_is_degenerate());
        c.ac = 5;
        assert!(!c.hwe_is_degenerate());
        assert!(!count_site(&[0, 1, 0, 1], None).hwe_is_degenerate());
        // Empty site: AN = 0.
        assert!(count_site(&[], None).hwe_is_degenerate());
    }

    /// Degenerate sites report the fill-tags sentinel p = 1 (good), never
    /// NaN — including the empty site, where the raw test would give NaN.
    #[test]
    fn hwe_p_degenerate_is_one() {
        assert_eq!(count_site(&[0, 0, 0, 0], None).hwe_p(), 1.0); // no ALT side
        assert_eq!(count_site(&[1, 1], None).hwe_p(), 1.0); // no REF side
        assert_eq!(count_site(&[], None).hwe_p(), 1.0); // empty
        assert_eq!(count_site(&[-1, -1], None).hwe_p(), 1.0); // all missing
    }

    /// Non-degenerate sites pass straight through to the exact test with
    /// (het, hom_ref, hom_alt) in that order.
    #[test]
    fn hwe_p_pairs_with_exact_test() {
        // One of each genotype class: the observed configuration is the
        // mode for mac = 3, so the two-sided p is 1.
        let c = count_site(&[0, 0, 0, 1, 1, 1], None);
        assert!((c.hwe_p() - 1.0).abs() < 1e-12);

        // 2 hom_ref + 2 hom_alt, 0 het: the observed configuration is the
        // least likely one (hets = 0 out of {0, 2, 4}), p = 3/35.
        let c = count_site(&[0, 0, 0, 0, 1, 1, 1, 1], None);
        assert!((c.hwe_p() - 3.0 / 35.0).abs() < 1e-12);

        // Wigginton two-sided exact test on a known site: n = 4 calls,
        // mac = 3 → valid het counts {1, 3} with P(1) = 3/7, P(3) = 4/7;
        // the observed configuration is the less likely one, so p = 3/7.
        let gt = [0, 0, 0, 1, 1, 1, 1, 1]; // hom_ref = 1, het = 1, hom_alt = 2
        let c = count_site(&gt, None);
        assert!((c.hwe_p() - 3.0 / 7.0).abs() < 1e-12, "{}", c.hwe_p());
    }

    /// A site whose only calls are half-calls has zero genotype
    /// observations: nothing is testable, so the p-value is the fill-tags
    /// convention for degenerate sites (HWE = 1, the site is kept).
    #[test]
    fn half_call_only_site_hwe_is_one() {
        let c = count_site(&[0, -1, -1, 1], None);
        assert!(!c.hwe_is_degenerate());
        assert!((c.hwe_p() - 1.0).abs() < 1e-12);
    }

    /// A trailing unpaired element carries no sample and is ignored.
    #[test]
    fn odd_gt_tail_is_ignored() {
        let c = count_site(&[0, 1, 0], None);
        assert_eq!((c.ac, c.an, c.ns), (1, 2, 1));
        assert_eq!(c.het, 1);
    }
}

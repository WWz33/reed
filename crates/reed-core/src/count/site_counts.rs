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

use crate::bitset::Bitset;
use crate::genotype::{GtClass, GtCode};
use crate::math::hwe::hwe_exact_p;
use crate::sample::SampleMask;

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

impl SiteCounts {
    /// Counts one site from decoded genotypes (doc_v25 V25-07a): the one
    /// counting semantic for the filter, stats and report paths.
    ///
    /// `keep` is the static sample mask (its length is a type guarantee) and
    /// `masked` the genotype mask: a masked call counts as missing, exactly
    /// like a raw `./.`. `has_gt` is a record property (whether the FORMAT
    /// carried GT), not a genotype one — the caller sets it.
    #[must_use]
    pub fn from_codes(codes: &[GtCode], keep: &SampleMask, masked: Option<&Bitset>) -> Self {
        assert_eq!(
            codes.len(),
            keep.n_samples(),
            "invariant: codes and keep-mask lengths differ"
        );
        Self::count(codes, |i| keep.get(i), masked)
    }

    /// [`Self::from_codes`] for callers that hold the full-length `&[bool]`
    /// pipeline mask. The length is asserted, never defaulted: the old
    /// `unwrap_or(true)` short-mask default was S68's root API.
    #[must_use]
    pub fn from_codes_kept(codes: &[GtCode], keep: &[bool], masked: Option<&Bitset>) -> Self {
        assert_eq!(
            codes.len(),
            keep.len(),
            "invariant: codes and keep-mask lengths differ"
        );
        Self::count(codes, |i| keep[i], masked)
    }

    fn count(codes: &[GtCode], keep: impl Fn(usize) -> bool, masked: Option<&Bitset>) -> Self {
        if let Some(masked) = masked {
            assert_eq!(
                masked.len(),
                codes.len(),
                "invariant: codes and genotype-mask lengths differ"
            );
        }
        let mut counts = Self::default();
        for (i, &code) in codes.iter().enumerate() {
            if !keep(i) {
                continue;
            }
            let code = if masked.is_some_and(|m| m.get(i)) {
                GtCode::MISSING
            } else {
                code
            };
            let class = code.class();
            counts.ac += u64::from(code.n_alt());
            counts.an += u64::from(code.n_ref() + code.n_alt());
            if class != GtClass::Missing {
                counts.ns += 1;
            }
            match class {
                GtClass::HomRef => counts.hom_ref += 1,
                GtClass::Het => counts.het += 1,
                GtClass::HomAlt => counts.hom_alt += 1,
                GtClass::HalfCall => counts.n_half_called += 1,
                GtClass::Missing => counts.n_missing += 1,
            }
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::{STATE_ALT, STATE_MISSING, STATE_REF};

    fn code(slot0: u8, slot1: u8) -> GtCode {
        GtCode::from_states(slot0, slot1, false, false)
    }

    fn all(n: usize) -> SampleMask {
        SampleMask::all(n)
    }

    /// s0 = 0/1, s1 = 1/1, s2 = ./., s3 = 0/0.
    #[test]
    fn mixed_genotypes_with_missing() {
        let codes = [
            code(STATE_REF, STATE_ALT),
            code(STATE_ALT, STATE_ALT),
            GtCode::MISSING,
            code(STATE_REF, STATE_REF),
        ];
        let c = SiteCounts::from_codes(&codes, &all(4), None);
        assert_eq!((c.ac, c.an, c.ns), (3, 6, 3));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (1, 1, 1));
        assert_eq!(c.n_missing, 1);
        assert_eq!(c.n_half_called, 0);
    }

    /// Dropping a sample removes it from every denominator, missing count
    /// included.
    #[test]
    fn retained_subset_excludes_dropped_samples() {
        let codes = [
            code(STATE_REF, STATE_ALT),
            code(STATE_ALT, STATE_ALT),
            GtCode::MISSING,
        ];
        let mask = SampleMask::from_ids(3, [crate::SampleId(0), crate::SampleId(1)]).expect("ids");
        let c = SiteCounts::from_codes(&codes, &mask, None);
        assert_eq!((c.ac, c.an, c.ns), (3, 4, 2));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 1, 1));
        assert_eq!(c.n_missing, 0);
    }

    /// Half-calls feed AC/AN/NS but no genotype class; a masked called
    /// sample becomes fully missing (the filter's count attribution).
    #[test]
    fn half_calls_and_masks_follow_the_shared_class() {
        let codes = [code(STATE_ALT, STATE_MISSING), code(STATE_REF, STATE_ALT)];
        let c = SiteCounts::from_codes(&codes, &all(2), None);
        assert_eq!((c.ac, c.an, c.ns), (2, 3, 2));
        assert_eq!(c.n_half_called, 1);
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 1, 0));

        let mut masked = Bitset::new(2);
        masked.set(1, true);
        let c = SiteCounts::from_codes(&codes, &all(2), Some(&masked));
        assert_eq!((c.ac, c.an, c.ns), (1, 1, 1));
        assert_eq!(c.n_missing, 1);
        assert_eq!(c.n_half_called, 1);
    }

    /// One site from slot-state pairs (keeps the old `count_site` fixtures
    /// readable after its deletion).
    fn site(pairs: &[(u8, u8)]) -> SiteCounts {
        let codes: Vec<GtCode> = pairs.iter().map(|&(a, b)| code(a, b)).collect();
        SiteCounts::from_codes(&codes, &SampleMask::all(codes.len()), None)
    }

    /// After a DP mask, s0 is masked (missing) and s1 = 1/1 remains; the
    /// transition accounting itself lives in the filter, the counts here.
    #[test]
    fn masked_site_end_state() {
        let c = site(&[
            (STATE_MISSING, STATE_MISSING),
            (STATE_ALT, STATE_ALT),
            (STATE_MISSING, STATE_MISSING),
        ]);
        assert_eq!((c.ac, c.an, c.ns), (2, 2, 1));
        assert_eq!(c.n_missing, 2);
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 0, 1));
    }

    /// Half-calls feed AC/AN/NS, never a genotype class and never
    /// `n_missing`; the allele side is irrelevant to the AC contribution.
    #[test]
    fn half_calls_count_one_allele() {
        let c = site(&[(STATE_REF, STATE_MISSING), (STATE_MISSING, STATE_ALT)]);
        assert_eq!((c.ac, c.an, c.ns), (1, 2, 2));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 0, 0));
        assert_eq!(c.n_missing, 0);
        assert_eq!(c.n_half_called, 2);
        assert!((c.maf() - 0.5).abs() < 1e-12);

        let c = site(&[(STATE_ALT, STATE_MISSING), (STATE_MISSING, STATE_ALT)]);
        assert_eq!((c.ac, c.an), (2, 2));
        assert_eq!(c.n_half_called, 2);

        let c = site(&[(STATE_MISSING, STATE_REF)]);
        assert_eq!((c.ac, c.an, c.ns), (0, 1, 1));
        assert_eq!(c.n_half_called, 1);
        assert_eq!(c.maf(), 0.0);
    }

    /// MAF boundaries: 0 on an empty/fully-missing site and on fixed sites
    /// (either side empty), 0.5 at parity.
    #[test]
    fn maf_boundaries() {
        assert_eq!(site(&[]).maf(), 0.0);
        assert_eq!(
            site(&[
                (STATE_MISSING, STATE_MISSING),
                (STATE_MISSING, STATE_MISSING)
            ])
            .maf(),
            0.0
        );
        assert_eq!(
            site(&[(STATE_REF, STATE_REF), (STATE_REF, STATE_REF)]).maf(),
            0.0
        );
        assert_eq!(
            site(&[(STATE_ALT, STATE_ALT), (STATE_ALT, STATE_ALT)]).maf(),
            0.0
        );
        let c = site(&[(STATE_REF, STATE_ALT), (STATE_ALT, STATE_REF)]);
        assert!((c.maf() - 0.5).abs() < 1e-12);
        let c = SiteCounts {
            ac: 5,
            an: 10,
            ..SiteCounts::default()
        };
        assert!((c.maf() - 0.5).abs() < 1e-12);
    }

    /// HWE degenerate sentinels: AN = 0 or an empty allele side.
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
        assert!(!site(&[(STATE_REF, STATE_ALT), (STATE_REF, STATE_ALT)]).hwe_is_degenerate());
        assert!(site(&[]).hwe_is_degenerate());
    }

    /// Degenerate sites report the fill-tags sentinel p = 1 (good), never
    /// NaN — including the empty site, where the raw test would give NaN.
    #[test]
    fn hwe_p_degenerate_is_one() {
        assert_eq!(
            site(&[(STATE_REF, STATE_REF), (STATE_REF, STATE_REF)]).hwe_p(),
            1.0
        );
        assert_eq!(site(&[(STATE_ALT, STATE_ALT)]).hwe_p(), 1.0);
        assert_eq!(site(&[]).hwe_p(), 1.0);
        assert_eq!(
            site(&[
                (STATE_MISSING, STATE_MISSING),
                (STATE_MISSING, STATE_MISSING)
            ])
            .hwe_p(),
            1.0
        );
    }

    /// Non-degenerate sites pass through to the exact test with
    /// (n_called, n_minor, n_het).
    #[test]
    fn hwe_p_pairs_with_exact_test() {
        let c = site(&[
            (STATE_REF, STATE_REF),
            (STATE_REF, STATE_ALT),
            (STATE_ALT, STATE_ALT),
        ]);
        assert!((c.hwe_p() - 1.0).abs() < 1e-12);

        let c = site(&[
            (STATE_REF, STATE_REF),
            (STATE_REF, STATE_REF),
            (STATE_ALT, STATE_ALT),
            (STATE_ALT, STATE_ALT),
        ]);
        assert!((c.hwe_p() - 3.0 / 35.0).abs() < 1e-12);

        // n = 4 calls, mac = 3: hom_ref = 1, het = 1, hom_alt = 2.
        let c = site(&[
            (STATE_REF, STATE_REF),
            (STATE_REF, STATE_ALT),
            (STATE_ALT, STATE_ALT),
            (STATE_ALT, STATE_ALT),
        ]);
        assert!((c.hwe_p() - 3.0 / 7.0).abs() < 1e-12, "{}", c.hwe_p());
    }

    /// A site whose only calls are half-calls has zero genotype observations:
    /// the p-value is the degenerate sentinel (HWE = 1, the site is kept).
    #[test]
    fn half_call_only_site_hwe_is_one() {
        let c = site(&[(STATE_REF, STATE_MISSING), (STATE_MISSING, STATE_ALT)]);
        assert!(!c.hwe_is_degenerate());
        assert!((c.hwe_p() - 1.0).abs() < 1e-12);
    }

    /// A short keep mask is a caller invariant violation, not an implicit
    /// "keep the rest" default (S68's root API).
    #[test]
    #[should_panic(expected = "keep-mask lengths differ")]
    fn short_keep_mask_panics() {
        let codes = [code(STATE_REF, STATE_ALT)];
        let _ = SiteCounts::from_codes_kept(&codes, &[], None);
    }
}

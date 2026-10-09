//! Popcount counting kernel (`docs/spec/architecture.md` §5.6).
//!
//! Computes all genotype/allele counts of one (site, group) pair from the bit
//! planes with AND/XOR/NOT plus population counts — no per-sample loop.
//!
//! Two implementations with byte-identical integer results:
//! 1. the scalar canonical path (always correct, parity oracle for the SIMD
//!    path);
//! 2. an AVX-512 path using `VPOPCNTDQ`, chosen at runtime when the CPU
//!    supports it. All arithmetic is exact integer popcounting, so dispatch
//!    never affects determinism.
//!
//! This module is the sanctioned `unsafe` exception in reed-core (spec
//! §14.4); every block carries a `// SAFETY:` comment.

#![allow(unsafe_code)]

use crate::genotype::SitePlanes;

use super::SiteGroupCounts;

/// Raw per-word popcount tallies before the `SiteGroupCounts` derivation.
#[derive(Debug, Default, Clone, Copy)]
struct RawTallies {
    alt0: u64,
    alt1: u64,
    miss0: u64,
    miss1: u64,
    full: u64,
    het: u64,
    hom_alt: u64,
    both_miss: u64,
}

/// Counts one (site, group) pair into `out`.
///
/// `group_words` is the group's sample mask (`SampleMask::words()`), of the
/// same word length as the site planes; `n_group` is its selected-sample
/// count (the caller already knows it, so the kernel does not recount).
pub fn count_site_group(
    site: &SitePlanes<'_>,
    group_words: &[u64],
    n_group: u32,
    out: &mut SiteGroupCounts,
) {
    debug_assert_eq!(site.n_words(), group_words.len(), "invariant: mask/plane word length");
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512f")
            && std::arch::is_x86_feature_detected!("avx512vpopcntdq")
        {
            let mut tallies = RawTallies::default();
            // SAFETY: both required CPU features were just detected at
            // runtime; the callee is annotated with `#[target_feature]`
            // enabling exactly those features and reads only in-bounds
            // `u64`-aligned memory.
            unsafe { count_words_avx512(site, group_words, &mut tallies) };
            finish(n_group, &tallies, out);
            return;
        }
    }
    count_site_group_scalar(site, group_words, n_group, out);
}

/// Scalar canonical implementation (parity oracle for the SIMD path).
pub fn count_site_group_scalar(
    site: &SitePlanes<'_>,
    group_words: &[u64],
    n_group: u32,
    out: &mut SiteGroupCounts,
) {
    let mut tallies = RawTallies::default();
    count_words_scalar(site, group_words, &mut tallies);
    finish(n_group, &tallies, out);
}

/// Accumulates the popcount tallies of every word with scalar operations.
fn count_words_scalar(site: &SitePlanes<'_>, group_words: &[u64], t: &mut RawTallies) {
    for (w, &g) in group_words.iter().enumerate().take(site.n_words()) {
        let a0 = site.slot0_alt[w] & g;
        let m0 = site.slot0_miss[w] & g;
        let a1 = site.slot1_alt[w] & g;
        let m1 = site.slot1_miss[w] & g;
        t.alt0 += u64::from(a0.count_ones());
        t.alt1 += u64::from(a1.count_ones());
        t.miss0 += u64::from(m0.count_ones());
        t.miss1 += u64::from(m1.count_ones());
        let full = !m0 & !m1 & g;
        t.full += u64::from(full.count_ones());
        t.het += u64::from(((a0 ^ a1) & full).count_ones());
        t.hom_alt += u64::from((a0 & a1).count_ones());
        t.both_miss += u64::from((m0 & m1).count_ones());
    }
}

/// Vector prefix (8 words per iteration) of the AVX-512 path; the caller
/// handles any trailing words with [`count_words_scalar`].
///
/// # Safety
///
/// `avx512f` and `avx512vpopcntdq` must be available on the current CPU
/// (checked with `is_x86_feature_detected!` by the dispatcher). Every slice
/// must be at least `group_words.len()` words long, which the dispatcher
/// asserts.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512vpopcntdq")]
// `_mm512_loadu_si512` is the unaligned load; the pointer cast never
// dereferences with alignment requirements.
#[allow(clippy::cast_ptr_alignment)]
unsafe fn count_words_avx512(
    site: &SitePlanes<'_>,
    group_words: &[u64],
    t: &mut RawTallies,
) {
    use std::arch::x86_64::{
        _mm512_add_epi64, _mm512_and_si512, _mm512_andnot_si512, _mm512_loadu_si512,
        _mm512_popcnt_epi64, _mm512_setzero_si512, _mm512_storeu_si512, _mm512_xor_si512,
        __m512i,
    };

    let n_words = site.n_words();
    let vector_words = n_words - (n_words % 8);
    let (mut acc_alt0, mut acc_alt1) = (_mm512_setzero_si512(), _mm512_setzero_si512());
    let (mut acc_miss0, mut acc_miss1) = (_mm512_setzero_si512(), _mm512_setzero_si512());
    let (mut acc_full, mut acc_het) = (_mm512_setzero_si512(), _mm512_setzero_si512());
    let (mut acc_hom_alt, mut acc_both_miss) = (_mm512_setzero_si512(), _mm512_setzero_si512());

    let mut w = 0;
    while w < vector_words {
        // SAFETY: `w + 8 <= n_words` and every slice is `n_words` long, so
        // each 64-byte load stays in bounds.
        unsafe {
            let g = _mm512_loadu_si512(group_words.as_ptr().add(w).cast::<__m512i>());
            let a0 = _mm512_and_si512(
                _mm512_loadu_si512(site.slot0_alt.as_ptr().add(w).cast::<__m512i>()),
                g,
            );
            let m0 = _mm512_and_si512(
                _mm512_loadu_si512(site.slot0_miss.as_ptr().add(w).cast::<__m512i>()),
                g,
            );
            let a1 = _mm512_and_si512(
                _mm512_loadu_si512(site.slot1_alt.as_ptr().add(w).cast::<__m512i>()),
                g,
            );
            let m1 = _mm512_and_si512(
                _mm512_loadu_si512(site.slot1_miss.as_ptr().add(w).cast::<__m512i>()),
                g,
            );
            // full = !m0 & !m1 & g (`andnot(a, b) = !a & b`).
            let full = _mm512_andnot_si512(m0, _mm512_andnot_si512(m1, g));
            acc_alt0 = _mm512_add_epi64(acc_alt0, _mm512_popcnt_epi64(a0));
            acc_alt1 = _mm512_add_epi64(acc_alt1, _mm512_popcnt_epi64(a1));
            acc_miss0 = _mm512_add_epi64(acc_miss0, _mm512_popcnt_epi64(m0));
            acc_miss1 = _mm512_add_epi64(acc_miss1, _mm512_popcnt_epi64(m1));
            acc_full = _mm512_add_epi64(acc_full, _mm512_popcnt_epi64(full));
            acc_het = _mm512_add_epi64(
                acc_het,
                _mm512_popcnt_epi64(_mm512_and_si512(_mm512_xor_si512(a0, a1), full)),
            );
            acc_hom_alt = _mm512_add_epi64(
                acc_hom_alt,
                _mm512_popcnt_epi64(_mm512_and_si512(a0, a1)),
            );
            acc_both_miss = _mm512_add_epi64(
                acc_both_miss,
                _mm512_popcnt_epi64(_mm512_and_si512(m0, m1)),
            );
        }
        w += 8;
    }

    let mut lanes = [0_u64; 8];
    let mut horizontal = |acc: __m512i| {
        // SAFETY: `lanes` is a local 64-byte buffer; storeu needs no alignment.
        unsafe { _mm512_storeu_si512(lanes.as_mut_ptr().cast::<__m512i>(), acc) };
        lanes.iter().sum::<u64>()
    };
    t.alt0 += horizontal(acc_alt0);
    t.alt1 += horizontal(acc_alt1);
    t.miss0 += horizontal(acc_miss0);
    t.miss1 += horizontal(acc_miss1);
    t.full += horizontal(acc_full);
    t.het += horizontal(acc_het);
    t.hom_alt += horizontal(acc_hom_alt);
    t.both_miss += horizontal(acc_both_miss);

    // Trailing words (n_words % 8) go through the scalar path.
    if vector_words < n_words {
        let tail_site = SitePlanes {
            slot0_alt: &site.slot0_alt[vector_words..],
            slot0_miss: &site.slot0_miss[vector_words..],
            slot1_alt: &site.slot1_alt[vector_words..],
            slot1_miss: &site.slot1_miss[vector_words..],
            phased: None,
        };
        count_words_scalar(&tail_site, &group_words[vector_words..], t);
    }
}

/// Derives [`SiteGroupCounts`] from the raw tallies (shared by both paths so
/// the §5.6 derivation rules live in exactly one place).
///
/// Cast safety: every tally is bounded by `2 * n_group <= 2 * u32::MAX`, and
/// the per-group sample count fits `u32` by construction of `SampleMask`.
#[allow(clippy::cast_possible_truncation)]
fn finish(n_group: u32, t: &RawTallies, out: &mut SiteGroupCounts) {
    let n_alt = t.alt0 + t.alt1;
    let n_missing = t.miss0 + t.miss1;
    let n_group64 = u64::from(n_group);
    let n_called_alleles = 2 * n_group64 - n_missing;
    let n_full_called = t.full;
    let n_hom_ref = n_full_called - t.het - t.hom_alt;
    let n_half_called = n_group64 - n_full_called - t.both_miss;
    out.n_ref_alleles = (n_called_alleles - n_alt) as u32;
    out.n_alt_alleles = n_alt as u32;
    out.n_hom_ref = n_hom_ref as u32;
    out.n_het = t.het as u32;
    out.n_hom_alt = t.hom_alt as u32;
    out.n_half_called = n_half_called as u32;
}

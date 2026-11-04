//! GT text → BCF code helpers (doc_v25 V25-02).
//!
//! The *syntax* of a text genotype — `(allele+1)<<1 | phased`, htslib's
//! `vcf.c:3496-3558` — lived inside `encode_native.rs` next to the typed
//! value encoders, while the domain rule (V25-01) lives in
//! `gt.rs`/`genotype_scan.rs`. This module is the single home for the
//! syntax side, shared by both encoder entry points; `bcf/gt.rs` remains the
//! table-based *decoder* and `vcf/genotype_scan.rs` the text scanner.
//!
//! Pure move: the functions are byte-for-byte the ones from
//! `encode_native.rs` (doc_v25 08 §3 ①).

/// `(ploidy, largest allele index)` of one text genotype.
pub(crate) fn gt_shape(sub: &[u8]) -> (usize, i64) {
    let mut ploidy = 0_usize;
    let mut max_allele = 0_i64;
    for allele in strip_gt_prefix(sub).split(|&b| b == b'/' || b == b'|') {
        ploidy += 1;
        if allele != b"."
            && !allele.is_empty()
            && let Some(v) = parse_i64(allele)
        {
            max_allele = max_allele.max(v);
        }
    }
    (ploidy, max_allele)
}

/// Drops a VCF 4.4 phasing prefix (htslib consumes it only for ≥ 4.4
/// headers; accepted unconditionally here, as in the scanner's decoder).
fn strip_gt_prefix(sub: &[u8]) -> &[u8] {
    match sub.first() {
        Some(b'|' | b'/') => &sub[1..],
        _ => sub,
    }
}

/// One text genotype's BCF codes in `out` (`(allele+1)<<1 | phased`),
/// exactly htslib's parser (vcf.c:3496-3558): allele `i > 0` carries the bit
/// of the separator before it; allele 0's bit is the explicit prefix when
/// present, else inferred (a haploid is phased, else phased iff every later
/// allele is); a missing allele stores the bare phasing bit.
pub(crate) fn gt_codes_into(sub: &[u8], out: &mut Vec<i32>) {
    out.clear();
    let (prefix, mut rest) = match sub.first() {
        Some(b'|') => (Some(true), &sub[1..]),
        Some(b'/') => (Some(false), &sub[1..]),
        _ => (None, sub),
    };
    let mut ploidy = 0_usize;
    let mut anyunphased = false;
    let mut bit = prefix.unwrap_or(false);
    loop {
        match rest.iter().position(|&b| b == b'/' || b == b'|') {
            Some(i) => {
                push_gt_code(out, &rest[..i], bit);
                ploidy += 1;
                bit = rest[i] == b'|';
                anyunphased |= !bit;
                rest = &rest[i + 1..];
            }
            None => {
                push_gt_code(out, rest, bit);
                ploidy += 1;
                break;
            }
        }
    }
    if prefix.is_none() {
        // No explicit phasing for allele 0: infer it the v44 way (a missing
        // haploid allele — code 0 — stays unphased, htslib's `unknown1`).
        let set = if ploidy == 1 {
            out.first().is_some_and(|&c| c != 0)
        } else {
            !anyunphased
        };
        if set && let Some(first) = out.first_mut() {
            *first |= 1;
        }
    }
}

/// One allele token's code: the index plus the phasing bit, or the bare bit
/// for a missing allele (htslib stores `.` as `is_phased`, vcf.c:3519).
fn push_gt_code(out: &mut Vec<i32>, allele: &[u8], phased: bool) {
    let bit = i32::from(phased);
    let code = match parse_i64(allele) {
        // An allele index beyond i32 cannot be represented in any BCF GT
        // encoding; treat it as missing rather than truncating. The bound is
        // the last allele whose `(a+1)<<1 | bit` still fits i32
        // (a = 1073741822; 1073741823 would wrap, doc_v10 D10-1).
        Some(a) if (0..=(i64::from(i32::MAX >> 1) - 1)).contains(&a) => {
            (((a + 1) << 1) as i32) | bit
        }
        _ => bit,
    };
    out.push(code);
}

/// Strict integer parse (optional leading '-', digits only).
/// The one INFO integer grammar (doc_v12 D12-14); see
/// [`crate::vcf::info::parse_info_int`].
pub(crate) fn parse_i64(text: &[u8]) -> Option<i64> {
    crate::vcf::info::parse_info_int(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_codes_follow_the_htslib_parser() {
        let cases: &[(&[u8], &[i32])] = &[
            (b"0/1", &[2, 4]),
            (b"0|1", &[3, 5]),
            (b"|0/1", &[3, 4]), // 4.4 prefix: allele 0 explicitly phased
            (b"/0|1", &[2, 5]), // explicit unphased first allele
            (b"1", &[5]),       // haploid: implicitly phased
            (b".", &[0]),       // missing haploid stays unphased
            (b"|.", &[1]),      // prefixed missing: the bare phasing bit
            (b"0|.", &[3, 1]),  // all phased: allele 0 inferred; missing carries the bit
            (b"0/.", &[2, 0]),
            (b"./.", &[0, 0]),
            (b".|.", &[1, 1]), // all phased -> first allele inferred phased
        ];
        let mut codes = Vec::new();
        for (text, want) in cases {
            gt_codes_into(text, &mut codes);
            assert_eq!(codes, *want, "GT {}", String::from_utf8_lossy(text));
        }
    }

    #[test]
    fn gt_shape_counts_alleles_after_the_prefix() {
        assert_eq!(gt_shape(b"|0/1"), (2, 1));
        assert_eq!(gt_shape(b"0/1"), (2, 1));
        assert_eq!(gt_shape(b"|1"), (1, 1));
        assert_eq!(gt_shape(b"."), (1, 0));
    }
}

//! Randomized robustness test for the unified decoder
//! (`docs/rounds/v2/03-architecture.md` §3.4: "fuzz 目标改为它").
//!
//! The scanner must never panic on adversarial sample regions: every input
//! either decodes or returns an error. Valid-but-unusual inputs (multi-digit
//! alleles, half calls, haploid runs, absent GT, short/long FORMAT strings)
//! must agree with an independent text-level reference. Property cases are
//! generated deterministically by proptest, so a failure reproduces from the
//! seed.

use proptest::prelude::*;

use reed_core::genotype::{
    GenotypePolicies, GtCode, MultiallelicPolicy, Ploidy, PloidyPolicy, PolicyCounters,
    STATE_ALT, STATE_MISSING, STATE_REF,
};
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;

fn policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

/// Text-level reference for one GT subfield: `None` = the text is not a
/// decodable GT (the scanner must error).
fn reference_gt(gt: &str) -> Option<(GtCode, [Option<u16>; 2])> {
    if gt.is_empty() {
        return None;
    }
    if gt == "." {
        return Some((
            GtCode::from_states(STATE_MISSING, STATE_MISSING, false, true),
            [None, None],
        ));
    }
    let mut slots = [None::<u16>; 2];
    let mut n = 0;
    let mut phased = false;
    let mut from = 0;
    let bytes = gt.as_bytes();
    loop {
        let end = bytes[from..]
            .iter()
            .position(|&b| b == b'/' || b == b'|')
            .map_or(bytes.len(), |rel| from + rel);
        if n >= 2 {
            return None; // polyploid: rejected
        }
        let token = &gt[from..end];
        slots[n] = if token == "." {
            None
        } else {
            let value: u16 = token.parse().ok()?;
            Some(value)
        };
        n += 1;
        if end == bytes.len() {
            break;
        }
        phased |= bytes[end] == b'|';
        from = end + 1;
    }
    let ploidy = if n == 1 { Ploidy::Haploid } else { Ploidy::Diploid };
    let mut scratch = PolicyCounters::default();
    let code = GtCode::from_alleles(slots, ploidy, phased, &policies(), &mut scratch).ok()?;
    Some((code, slots))
}

/// A sample column: a GT value plus trailing subfields.
fn arb_sample() -> impl Strategy<Value = String> {
    let allele = prop_oneof![Just(".".to_string()), (0_u16..5).prop_map(|i| i.to_string())];
    let gt = prop_oneof![
        (allele.clone(), allele.clone()).prop_map(|(a, b)| format!("{a}/{b}")),
        (allele.clone(), allele.clone()).prop_map(|(a, b)| format!("{a}|{b}")),
        allele.clone(),
        Just(".".to_string()),
        Just(String::new()),
    ];
    (gt, prop::option::of(0_u32..200), prop::option::of(0_u32..100))
        .prop_map(|(gt, dp, gq)| {
            let dp = dp.map_or(".".to_string(), |v| v.to_string());
            let gq = gq.map_or(".".to_string(), |v| v.to_string());
            format!("{gt}:{dp}:{gq}")
        })
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 2_000,
        max_shrink_iters: 1_000,
        ..ProptestConfig::default()
    })]

    /// Decoding is total: for any generated record the scanner either
    /// produces one code per sample or reports an error — it never panics,
    /// and where the reference decodes, the codes agree.
    #[test]
    fn scanner_is_total_and_matches_the_text_reference(
        samples in prop::collection::vec(arb_sample(), 1..8),
        format in prop_oneof![
            Just("GT:DP:GQ".to_string()),
            Just("GT".to_string()),
            Just("DP:GQ".to_string()),
            Just("GQ:GT".to_string()),
            Just("GT::DP".to_string()),
        ],
    ) {
        let line_text = format!("c\t1\t.\tA\tG\t.\tPASS\t.\t{format}\t{}\n", samples.join("\t"));
        let line = LineRef::parse(line_text.as_bytes()).expect("invariant: 10 columns");
        let site = SiteView::parse(&line).expect("invariant: site parses");
        let n_samples = samples.len();
        let mut scanner = GenotypeScanner::new(n_samples, &[b"DP", b"GQ"], &[], policies()).tracking_alleles();
        let mut counters = PolicyCounters::default();
        let out = scanner.scan(&line, &site, &mut counters);

        let gt_key_at: Option<usize> = format.split(':').position(|k| k == "GT");
        // The reference decides whether the record is decodable at all.
        let mut expected_ok = true;
        for sample in &samples {
            let subfields: Vec<&str> = sample.split(':').collect();
            let gt = gt_key_at.and_then(|i| subfields.get(i)).copied().unwrap_or(".");
            if reference_gt(gt).is_none() {
                expected_ok = false;
            }
        }
        match out {
            Ok(out) => {
                prop_assert!(expected_ok, "scanner decoded a record the reference rejects");
                prop_assert_eq!(out.codes.len(), n_samples);
                let alleles = out.alleles.expect("invariant: tracking is on");
                for (i, sample) in samples.iter().enumerate() {
                    let subfields: Vec<&str> = sample.split(':').collect();
                    let Some(gt) = gt_key_at.and_then(|k| subfields.get(k)).copied() else {
                        // GT absent from the FORMAT string: plain missing,
                        // no ploidy claim (there is no GT text to show one).
                        prop_assert!(out.codes[i].is_fully_missing());
                        prop_assert!(!out.codes[i].is_haploid());
                        prop_assert_eq!(&alleles[2 * i..2 * i + 2], &[None, None]);
                        continue;
                    };
                    let (code, slots) = reference_gt(gt).expect("reference decodes");
                    prop_assert_eq!(out.codes[i].state(), code.state(), "sample {} state", i);
                    prop_assert_eq!(out.codes[i].is_haploid(), code.is_haploid(), "sample {} ploidy", i);
                    prop_assert_eq!(&alleles[2 * i..2 * i + 2], &slots[..], "sample {} alleles", i);
                }
            }
            Err(_) => {
                prop_assert!(!expected_ok, "scanner rejected a record the reference decodes");
            }
        }
    }

    /// Mutated record text (inserted separators, garbage) must not panic.
    #[test]
    fn mutated_records_never_panic(
        samples in prop::collection::vec(arb_sample(), 1..5),
        mutation in prop::collection::vec(prop_oneof![
            Just(b'0'), Just(b'1'), Just(b'.'), Just(b'/'), Just(b'|'), Just(b':'),
            Just(b'x'), Just(b'9'),
        ], 0..6),
        at in 0_usize..64,
    ) {
        let text = format!("c\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t{}\n", samples.join("\t"));
        let mut bytes = text.into_bytes();
        for (k, &m) in mutation.iter().enumerate() {
            let idx = (at + k) % bytes.len();
            bytes[idx] = m;
        }
        let Some(line) = LineRef::parse(&bytes) else { return Ok(()); };
        let Ok(site) = SiteView::parse(&line) else { return Ok(()); };
        let mut scanner = GenotypeScanner::new(samples.len(), &[b"DP"], &[], policies());
        let mut counters = PolicyCounters::default();
        let _ = scanner.scan(&line, &site, &mut counters);
    }
}

#[test]
fn reference_and_scanner_agree_on_the_documented_edge_cases() {
    let cases = [
        ("0/0", 0b0000_u8),
        ("0/1", 0b0100),
        ("1/0", 0b0001),
        ("1/1", 0b0101),
        ("./.", 0b1010),
        ("0/.", 0b1000),
        ("./1", 0b0110),
        ("1", 0b1001),
        ("0", 0b1000),
        (".", 0b1010),
        ("10/11", 0b0101),
        ("1|0", 0b0001),
    ];
    let line_text = format!("c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t{}\n", cases.iter().map(|(g, _)| *g).collect::<Vec<_>>().join("\t"));
    let line = LineRef::parse(line_text.as_bytes()).expect("invariant: 10 columns");
    let site = SiteView::parse(&line).expect("invariant: site parses");
    let mut scanner = GenotypeScanner::new(cases.len(), &[], &[], policies());
    let mut counters = PolicyCounters::default();
    let out = scanner.scan(&line, &site, &mut counters).expect("invariant: decodes");
    for (i, (gt, expected)) in cases.iter().enumerate() {
        assert_eq!(out.codes[i].state(), *expected, "GT {gt}");
        // The text reference agrees with the pinned expectation.
        let (code, _) = reference_gt(gt).expect("invariant: reference decodes");
        assert_eq!(code.state(), *expected, "reference GT {gt}");
    }
    assert_eq!(counters.n_haploid_as_half_call, 2, "`1` and `0` are half calls");
    assert_eq!(counters.n_multiallelic_collapsed, 2, "`10` and `11` collapse");
    let _ = STATE_REF;
    let _ = STATE_ALT;
}

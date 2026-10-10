//! O6 gate (doc_v3 03-performance.md): time the genotype scanner over one
//! realistic 348-sample record, for the GT-only and the GT+DP+GQ plans.
// The criterion macro expands to a `pub fn main()` with nothing to document.
#![allow(missing_docs)]

use criterion::{Criterion, criterion_group, criterion_main};

use reed_core::genotype::{GenotypePolicies, PolicyCounters};
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;

/// A 348-sample record shaped like the clean.vcf.gz data.
fn record() -> Vec<u8> {
    let mut line = String::from("Chr01\t36\t.\tG\tA\t4596.04\tPASS\t.\tGT:DP:GQ:PL");
    for i in 0..348 {
        let (gt, dp, gq) = match i % 50 {
            0 => ("./.", 0, 0),
            7 => ("0/1", 3, 5),
            11 => ("1/1", 150, 99),
            _ => ("0/1", 30, 99),
        };
        use std::fmt::Write as _;
        write!(line, "\t{gt}:{dp}:{gq}:0,10,100").expect("invariant: writing to String");
    }
    line.push('\n');
    line.into_bytes()
}

fn bench(c: &mut Criterion) {
    let line = record();
    let lr = LineRef::parse(&line).expect("invariant: record parses");
    let site = SiteView::parse(&lr).expect("invariant: site parses");
    let policies = GenotypePolicies::default();

    let mut gt_only = GenotypeScanner::new(348, &[], &[], policies);
    let mut counters = PolicyCounters::default();
    c.bench_function("scan_gt_only_348_samples", |b| {
        b.iter(|| {
            let out = gt_only
                .scan(&lr, &site, &mut counters)
                .expect("invariant: scan");
            std::hint::black_box(out.codes.len());
        });
    });

    let mut gt_dp_gq = GenotypeScanner::new(348, &[b"DP", b"GQ"], &[], policies);
    c.bench_function("scan_gt_dp_gq_348_samples", |b| {
        b.iter(|| {
            let out = gt_dp_gq
                .scan(&lr, &site, &mut counters)
                .expect("invariant: scan");
            std::hint::black_box((out.codes.len(), out.int_columns.len()));
        });
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);

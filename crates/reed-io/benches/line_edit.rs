//! O7 gate (doc_v3 03-performance.md): time the sample-column edit of one
//! realistic 348-sample record, with the gold mask rate applied.
// The criterion macro expands to a `pub fn main()` with nothing to document.
#![allow(missing_docs)]

use criterion::{Criterion, criterion_group, criterion_main};

use reed_core::Bitset;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{EditPlan, InfoValueOut, edit_vcf_line};

/// A 348-sample record whose GT+DP+GQ profile mirrors the clean.vcf.gz data.
fn record() -> Vec<u8> {
    let mut line = String::from(
        "Chr01\t36\tChr01_36\tG\tA\t4596.04\tPASS\tAC=55;AN=696\tGT:DP:GQ:PL",
    );
    for i in 0..348 {
        let (gt, dp, gq) = match i % 50 {
            0 => ("./.", 0, 0),
            7 => ("0/1", 3, 5),
            11 => ("1/1", 150, 99),
            _ => ("0/1", 30, 99),
        };
        line.push_str(&format!("\t{gt}:{dp}:{gq}:0,10,100"));
    }
    line.push('\n');
    line.into_bytes()
}

fn bench(c: &mut Criterion) {
    let line = record();
    let lr = LineRef::parse(&line).expect("invariant: record parses");
    // The data's masking rate: DP<5 or DP>100 or GQ<20.
    let mut mask = Bitset::new(348);
    for i in (0..348).filter(|i| matches!(i % 50, 7 | 11)) {
        mask.set(i, true);
    }
    let stale: [&[u8]; 2] = [b"AC", b"AN"];
    let info: Vec<(&[u8], InfoValueOut)> = vec![
        (b"AC", InfoValueOut::Int(51)),
        (b"AN", InfoValueOut::Int(688)),
    ];
    let plan = EditPlan {
        keep_samples: None,
        mask_gt: Some(&mask),
        info_remove: &stale,
        info_set: &info,
    };
    let mut out = Vec::with_capacity(line.len());
    c.bench_function("edit_vcf_line_348_samples", |b| {
        b.iter(|| {
            out.clear();
            edit_vcf_line(&lr, &plan, Some(0), &mut out).expect("invariant: edit");
        });
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);

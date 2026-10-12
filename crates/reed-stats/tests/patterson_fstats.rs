//! Patterson f-statistics canonical tests: hand-computed per-site values
//! against allel's `admixture.py` formulas (read directly), estimator
//! identities, and the block jackknife against a brute-force delete-one
//! recomputation.

use reed_core::count::SiteGroupCounts;
use reed_core::math::sum::NeumaierSum;
use reed_stats::jackknife::{BlockJackknife, BlockPartial};
use reed_stats::patterson::{F2Acc, F3Acc, F4Acc, PattersonKind, f4_ratio};

fn row(an: [u32; 4], alt: [u32; 4]) -> Vec<SiteGroupCounts> {
    an.iter()
        .zip(alt)
        .map(|(an, alt)| SiteGroupCounts {
            n_ref_alleles: an - alt,
            n_alt_alleles: alt,
            ..SiteGroupCounts::default()
        })
        .collect()
}

/// allel admixture.py reference formulas (f64, per site).
fn ref_f2(an: [u32; 4], alt: [u32; 4]) -> f64 {
    let (p1, p2) = (
        f64::from(alt[0]) / f64::from(an[0]),
        f64::from(alt[1]) / f64::from(an[1]),
    );
    let h = |a: u32, alt: u32| -> f64 {
        let na = f64::from(a);
        (f64::from(a - alt) * f64::from(alt)) / (na * (na - 1.0))
    };
    (p1 - p2) * (p1 - p2)
        - h(an[0], alt[0]) / f64::from(an[0])
        - h(an[1], alt[1]) / f64::from(an[1])
}

fn ref_f3(an: [u32; 4], alt: [u32; 4]) -> (f64, f64) {
    let (pa, pb, pc) = (
        f64::from(alt[0]) / f64::from(an[0]),
        f64::from(alt[1]) / f64::from(an[1]),
        f64::from(alt[2]) / f64::from(an[2]),
    );
    let hc =
        (f64::from(an[2] - alt[2]) * f64::from(alt[2])) / (f64::from(an[2]) * f64::from(an[2] - 1));
    let t = (pc - pa) * (pc - pb) - hc / f64::from(an[2]);
    (t, 2.0 * hc)
}

fn ref_d(an: [u32; 4], alt: [u32; 4]) -> (f64, f64) {
    let p: Vec<f64> = (0..4)
        .map(|i| f64::from(alt[i]) / f64::from(an[i]))
        .collect();
    let num = (p[0] - p[1]) * (p[2] - p[3]);
    let den = (p[0] + p[1] - 2.0 * p[0] * p[1]) * (p[2] + p[3] - 2.0 * p[2] * p[3]);
    (num, den)
}

#[test]
fn f2_corrected_matches_allel_formula() {
    let mut acc = F2Acc::new(0, 1, PattersonKind::Corrected);
    let sites = [
        ([10_u32, 10, 10, 10], [2_u32, 7, 1, 5]),
        ([8, 6, 12, 4], [4, 3, 6, 2]),
        ([6, 12, 8, 10], [3, 9, 2, 8]),
    ];
    let mut total = 0.0;
    for (an, alt) in sites {
        acc.push(&row(an, alt));
        total += ref_f2(an, alt);
    }
    assert_eq!(acc.site_tally(), (3, 0));
    assert!((acc.value().expect("invariant: sites") - total).abs() < 1e-12);
}

#[test]
fn f2_naive_drops_the_correction() {
    let mut corrected = F2Acc::new(0, 1, PattersonKind::Corrected);
    let mut naive = F2Acc::new(0, 1, PattersonKind::Naive);
    let (an, alt) = ([10_u32, 10, 10, 10], [2_u32, 7, 1, 5]);
    corrected.push(&row(an, alt));
    naive.push(&row(an, alt));
    let (p1, p2) = (0.2, 0.7);
    assert!((naive.value().expect("invariant: sites") - (p1 - p2) * (p1 - p2)).abs() < 1e-12);
    assert!(
        corrected.value().expect("invariant: sites") < naive.value().expect("invariant: sites")
    );
}

#[test]
fn f2_corrected_skips_single_chromosome_sites() {
    let mut acc = F2Acc::new(0, 1, PattersonKind::Corrected);
    acc.push(&row([1_u32, 10, 10, 10], [1_u32, 5, 3, 2])); // anA = 1
    assert_eq!(acc.site_tally(), (0, 1));
    assert!(acc.value().is_none());
}

#[test]
fn f3_carries_the_normalized_variant() {
    let mut acc = F3Acc::new(2, 0, 1, PattersonKind::Corrected);
    let sites = [
        ([10_u32, 8, 12, 10], [3_u32, 4, 6, 5]),
        ([8, 10, 6, 12], [2, 6, 1, 7]),
    ];
    let (mut t, mut b) = (0.0, 0.0);
    for (an, alt) in sites {
        acc.push(&row(an, alt));
        let (t_i, b_i) = ref_f3(an, alt);
        t += t_i;
        b += b_i;
    }
    assert!((acc.value().expect("invariant: sites") - t).abs() < 1e-12);
    let (t_sum, b_sum) = acc.sums();
    assert_eq!(t_sum.value(), t);
    assert_eq!(b_sum.value(), b);
    assert!((acc.value_normalized().expect("invariant: B ≠ 0") - t / b).abs() < 1e-12);
}

#[test]
fn f4_and_d_share_the_sums() {
    let mut acc = F4Acc::new(0, 1, 2, 3);
    let sites = [
        ([10_u32, 12, 8, 14], [4_u32, 6, 2, 8]),
        ([6, 10, 12, 8], [1, 5, 7, 3]),
        ([12, 8, 10, 10], [6, 2, 4, 6]),
    ];
    let (mut num, mut den) = (0.0, 0.0);
    for (an, alt) in sites {
        acc.push(&row(an, alt));
        let (n_i, d_i) = ref_d(an, alt);
        num += n_i;
        den += d_i;
    }
    assert!((acc.value_f4().expect("invariant: sites") - num).abs() < 1e-12);
    assert!((acc.value_d().expect("invariant: den ≠ 0") - num / den).abs() < 1e-12);
    let (n_sum, d_sum) = acc.sums();
    assert_eq!(n_sum.value(), num);
    assert_eq!(d_sum.value(), den);
}

#[test]
fn f4_ratio_is_the_ratio_of_sums() {
    let mut num_acc = F4Acc::new(0, 1, 2, 3);
    let mut den_acc = F4Acc::new(0, 3, 2, 1);
    num_acc.push(&row([10_u32, 12, 8, 14], [4, 6, 2, 8]));
    den_acc.push(&row([10_u32, 12, 8, 14], [4, 6, 2, 8]));
    den_acc.push(&row([6, 10, 12, 8], [1, 5, 7, 3]));
    let ratio = f4_ratio(&num_acc, &den_acc).expect("invariant: den ≠ 0");
    // Both operands are f4 statistics, so the denominator is the second
    // f4's *numerator* sum, not its D-statistic denominator (doc_v7 D7-10).
    let (n, d) = (num_acc.sums().0.value(), den_acc.sums().0.value());
    assert!((ratio - n / d).abs() < 1e-12);
}

#[test]
fn jackknife_se_matches_brute_force_delete_one() {
    // Deterministic pseudo-data: (num, den) pairs.
    let pairs: Vec<(f64, f64)> = (0..200)
        .map(|i| {
            let i = f64::from(i);
            ((i * 7.0 % 13.0) - 3.0, 5.0 + (i * 11.0 % 17.0))
        })
        .collect();
    let block_size = 20_u64;
    let mut jk = BlockJackknife::new(block_size);
    for (num, den) in &pairs {
        jk.push(*num, *den);
    }
    let (est, se, vj) = jk.jackknife().expect("invariant: ≥ 2 blocks");

    // Brute force: recompute the sums per block and delete each.
    let n_blocks = usize::try_from((pairs.len() as u64).div_ceil(block_size)).unwrap();
    let mut block_sums = vec![(0.0, 0.0); n_blocks];
    for (i, (num, den)) in pairs.iter().enumerate() {
        let b = u64::try_from(i).unwrap() / block_size;
        block_sums[usize::try_from(b).unwrap()].0 += num;
        block_sums[usize::try_from(b).unwrap()].1 += den;
    }
    let total_num: f64 = block_sums.iter().map(|(n, _)| n).sum();
    let total_den: f64 = block_sums.iter().map(|(_, d)| d).sum();
    assert!((est - total_num / total_den).abs() < 1e-9);

    let mut vj_brute = Vec::new();
    for (bn, bd) in &block_sums {
        let num = total_num - bn;
        let den = total_den - bd;
        if den != 0.0 {
            vj_brute.push(num / den);
        }
    }
    let n = vj_brute.len() as f64;
    let mean = vj_brute.iter().sum::<f64>() / n;
    let se_brute = (((n - 1.0) / n)
        * vj_brute
            .iter()
            .map(|v| (v - mean) * (v - mean))
            .sum::<f64>())
    .sqrt();
    assert_eq!(vj.len(), vj_brute.len());
    for (a, b) in vj.iter().zip(&vj_brute) {
        assert!((a - b).abs() < 1e-9);
    }
    assert!((se - se_brute).abs() < 1e-9, "{se} vs {se_brute}");
}

#[test]
fn jackknife_needs_two_blocks() {
    // block_size = 10: two sites are one open block.
    let mut jk = BlockJackknife::new(10);
    jk.push(1.0, 2.0);
    jk.push(3.0, 4.0);
    assert_eq!(jk.n_blocks(), 1);
    assert!(jk.jackknife().is_none());
    // block_size = 1: two sealed singleton blocks.
    let mut jk = BlockJackknife::new(1);
    jk.push(1.0, 2.0);
    jk.push(3.0, 4.0);
    assert_eq!(jk.n_blocks(), 2);
    assert!(jk.jackknife().is_some());
}

#[test]
fn chunked_accept_reproduces_the_whole() {
    // The block grid is the global ordinal grid: feeding two chunk
    // partials through one accumulator (any split point) reproduces the
    // single-partial run exactly.
    let pairs: Vec<(f64, f64)> = (0..100)
        .map(|i| (f64::from(i), 10.0 + f64::from(i)))
        .collect();
    let mut whole = BlockJackknife::new(10);
    whole.accept(&BlockPartial {
        values: pairs.clone(),
    });

    let mut split = BlockJackknife::new(10);
    split.accept(&BlockPartial {
        values: pairs[..55].to_vec(),
    });
    split.accept(&BlockPartial {
        values: pairs[55..].to_vec(),
    });

    assert_eq!(split.n_blocks(), whole.n_blocks());
    let (est_w, se_w, vj_w) = whole.jackknife().expect("invariant: blocks");
    let (est_s, se_s, vj_s) = split.jackknife().expect("invariant: blocks");
    assert!((est_w - est_s).abs() < 1e-12);
    assert!((se_w - se_s).abs() < 1e-12);
    for (x, y) in vj_w.iter().zip(&vj_s) {
        assert!((x - y).abs() < 1e-12);
    }
}

#[test]
fn neumaier_is_what_the_sums_hold() {
    // A visible guard: the sums are compensated, so the block totals are
    // order-stable (§8.4).
    let mut a = NeumaierSum::new();
    for v in [1.0e16, 1.0, 1.0, -1.0e16] {
        a.add(v);
    }
    let mut b = NeumaierSum::new();
    for v in [1.0, -1.0e16, 1.0e16, 1.0] {
        b.add(v);
    }
    assert_eq!(a.value(), 2.0);
    assert_eq!(b.value(), 2.0);
}

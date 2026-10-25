//! Does using the GRM geometry improve accuracy?
//!
//! Controlled experiment: simulate genotypes from a known three-population
//! model, then run each reduction method twice -- once on the raw genotype
//! geometry (Euclidean distance between dosage rows) and once on GRM geometry
//! (`d^2 = a_ii + a_jj - 2 a_ij`) -- and score both against the known labels.
//!
//! Metrics:
//!   * `silhouette`: (between - within) / max(between, within) on the known
//!     three groups; higher is better, 1.0 means perfectly separated.
//!   * `knn_acc`: fraction of samples whose k nearest neighbours in the
//!     embedding share their true label.
//!
//! Usage: `cargo run --release --example grm_accuracy -- [n_per] [n_snp] [reps]`

#![allow(clippy::needless_range_loop)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_borrow)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]
#![allow(clippy::explicit_counter_loop)]
#![allow(clippy::useless_vec)]
#![allow(clippy::redundant_closure)]

use reed_reduce::extra::{lle_from_sq_distances, tsne_from_sq_distances, TsneConfig};
use reed_reduce::graph::{kernel_pca, laplacian_eigenmaps, pairwise_distances, pcoa, Kernel};
use reed_reduce::pca::{pca_dense, standardize, Grm, Stand};
use reed_reduce::umap::{umap_from_distances, UmapConfig};
use reed_reduce::{grm_distances, Mat};
use std::time::Instant;

/// Allele-frequency spectrum to draw the ancestral frequencies from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Spectrum {
    /// Skewed toward rare variants (`0.01 + 0.98 u^2`).
    Rare,
    /// Uniform over [0.05, 0.95].
    Uniform,
    /// Common variants only ([0.2, 0.8]).
    Common,
}

/// Simulate three populations with Fst-like divergence.
fn simulate(n_per: usize, n_snp: usize, seed: u64, spec: Spectrum) -> (Mat, Vec<usize>) {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    // Ancestral frequencies, skewed toward rare variants as in real data.
    let base: Vec<f64> = (0..n_snp)
        .map(|_| {
            let u = next();
            match spec {
                Spectrum::Rare => 0.01 + 0.98 * u * u,
                Spectrum::Uniform => 0.05 + 0.90 * u,
                Spectrum::Common => 0.20 + 0.60 * u,
            }
        })
        .collect();
    // Per-population drift.
    let drift = [0.0, 0.16, -0.13];
    let mut labels = Vec::with_capacity(3 * n_per);
    for g in 0..3 {
        labels.extend(std::iter::repeat_n(g, n_per));
    }
    let x = Mat::from_fn(3 * n_per, n_snp, move |i, j| {
        let g = i / n_per;
        let p = (base[j] + drift[g]).clamp(0.001, 0.999);
        let u = next();
        if u < p * p {
            2.0
        } else if u < p * p + 2.0 * p * (1.0 - p) {
            1.0
        } else {
            0.0
        }
    });
    (x, labels)
}

/// Silhouette on three known groups, using the embedding's Euclidean metric.
fn silhouette(coords: &Mat, labels: &[usize]) -> f64 {
    let n = coords.rows;
    let dim = coords.cols;
    if dim == 0 || n < 4 {
        return 0.0;
    }
    let mut total = 0.0;
    let mut counted = 0usize;
    for i in 0..n {
        let li = labels[i];
        let mut same = 0usize;
        let mut sum_same = 0.0;
        let mut sum_other = vec![0.0f64; 3];
        let mut cnt_other = vec![0usize; 3];
        for j in 0..n {
            if i == j {
                continue;
            }
            let mut s = 0.0;
            for d in 0..dim {
                s += (coords.get(i, d) - coords.get(j, d)).powi(2);
            }
            let dist = s.sqrt();
            if labels[j] == li {
                sum_same += dist;
                same += 1;
            } else {
                sum_other[labels[j]] += dist;
                cnt_other[labels[j]] += 1;
            }
        }
        if same == 0 {
            continue;
        }
        let a = sum_same / same as f64;
        let mut b = f64::INFINITY;
        for g in 0..3 {
            if g != li && cnt_other[g] > 0 {
                b = b.min(sum_other[g] / cnt_other[g] as f64);
            }
        }
        if !b.is_finite() {
            continue;
        }
        total += (b - a) / a.max(b);
        counted += 1;
    }
    if counted == 0 {
        0.0
    } else {
        total / counted as f64
    }
}

/// Fraction of the 10 nearest neighbours sharing the sample's label.
fn knn_acc(coords: &Mat, labels: &[usize], k: usize) -> f64 {
    let n = coords.rows;
    let dim = coords.cols;
    let mut hits = 0usize;
    let mut total = 0usize;
    for i in 0..n {
        let mut order: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (coords.get(i, d) - coords.get(j, d)).powi(2);
                }
                (s, j)
            })
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for &(_, j) in order.iter().take(k) {
            if labels[j] == labels[i] {
                hits += 1;
            }
            total += 1;
        }
    }
    hits as f64 / total.max(1) as f64
}

struct Score {
    sil: f64,
    acc: f64,
}

fn score(m: &Mat, labels: &[usize]) -> Score {
    Score { sil: silhouette(m, labels), acc: knn_acc(m, labels, 10) }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n_per: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(40);
    let n_snp: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(500);
    let reps: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(5);

    println!(
        "experiment: {reps} replicates, {} samples (3 x {n_per}) x {n_snp} SNPs\n",
        3 * n_per
    );
    println!(
        "{:<22} {:>18} {:>18} {:>10}",
        "method", "raw geometry", "GRM geometry", "delta(sil)"
    );
    println!("{:<22} {:>18} {:>18} {:>10}", "", "sil / knn_acc", "sil / knn_acc", "");

    for spec in [Spectrum::Rare, Spectrum::Uniform, Spectrum::Common] {
    // Accumulators.
    type Acc = (f64, f64, f64, f64);
    let mut acc: std::collections::BTreeMap<&str, Acc> = std::collections::BTreeMap::new();

    for rep in 0..reps {
        let seed = 0x9E3779B97F4A7C15u64.wrapping_mul(rep as u64 + 1);
        let (x, labels) = simulate(n_per, n_snp, seed, spec);
        let st = standardize(&x, Stand::Binom2);
        let grm = Grm::mean_impute(&st.w);

        let d_raw = pairwise_distances(&x);
        let d_grm = grm_distances(&grm);
        // Squared distances for t-SNE and LLE.
        let sq = |d: &Mat| {
            let mut out = Mat::zeros(d.rows, d.cols);
            for i in 0..d.rows {
                for j in 0..d.cols {
                    let v = d.get(i, j);
                    out.set(i, j, v * v);
                }
            }
            out
        };
        let d_raw2 = sq(&d_raw);
        let d_grm2 = sq(&d_grm);

        let w_raw = reed_reduce::graph::knn_graph(&x, 10, 0.0);
        let w_grm = reed_reduce::grm_knn_weights(&grm, 10);

        let tsne_cfg = TsneConfig { perplexity: 10.0, n_iter: 400, seed: 3, ..Default::default() };
        let umap_cfg = UmapConfig { n_neighbors: 10, n_epochs: 250, seed: 7, ..Default::default() };

        let pairs: Vec<(&str, Mat, Mat)> = vec![
            ("pca", pca_dense(&x, 2).scores(), pca_dense(&x, 2).scores()),
            ("pcoa", pcoa(&d_raw, 2).coords, pcoa(&d_grm, 2).coords),
            (
                "kernel_pca(rbf)",
                kernel_pca(&x, 2, Kernel::Rbf, 0.0).coords,
                reed_reduce::grm_kernel(&grm, 2, Kernel::Rbf).coords,
            ),
            (
                "laplacian",
                laplacian_eigenmaps(&w_raw, 2).coords,
                laplacian_eigenmaps(&w_grm, 2).coords,
            ),
            (
                "diffusion",
                reed_reduce::graph::diffusion_map(&w_raw, 2, 1.0).coords,
                reed_reduce::graph::diffusion_map(&w_grm, 2, 1.0).coords,
            ),
            (
                "isomap",
                reed_reduce::graph::isomap(&x, 10, 2).coords,
                reed_reduce::grm_isomap(&grm, 2, 10).coords,
            ),
            (
                "tsne",
                tsne_from_sq_distances(&d_raw2, &tsne_cfg).expect("tsne (raw)"),
                tsne_from_sq_distances(&d_grm2, &tsne_cfg).expect("tsne (grm)"),
            ),
            (
                "umap",
                umap_from_distances(&d_raw, &umap_cfg),
                umap_from_distances(&d_grm, &umap_cfg),
            ),
            (
                "lle",
                lle_from_sq_distances(&x, &d_raw2, 8, 2),
                lle_from_sq_distances(&st.w, &d_grm2, 8, 2),
            ),
        ];

        for (name, raw, grm_m) in pairs {
            let t = Instant::now();
            let s_raw = score(&raw, &labels);
            let s_grm = score(&grm_m, &labels);
            let e = acc.entry(name).or_insert((0.0, 0.0, 0.0, 0.0));
            e.0 += s_raw.sil;
            e.1 += s_raw.acc;
            e.2 += s_grm.sil;
            e.3 += s_grm.acc;
            let _ = t;
        }
    }

    println!("=== spectrum: {:?} ===", spec);
    for (name, (rs, ra, gs, ga)) in &acc {
        let n = reps as f64;
        let (rs, ra, gs, ga) = (rs / n, ra / n, gs / n, ga / n);
        println!(
            "{:<22} {:8.3} / {:6.3}  {:8.3} / {:6.3}  {:>+9.3}",
            name,
            rs,
            ra,
            gs,
            ga,
            gs - rs
        );
    }
    println!();
    }
    println!("pca uses the same input on both sides (it is already GRM-based).");
    println!("delta > 0 means GRM geometry separated the known populations better.");
}

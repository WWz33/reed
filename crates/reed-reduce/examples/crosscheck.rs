//! Cross-check example: emits deterministic embeddings for the shared "blobs"
//! fixture so the numbers can be compared with scikit-learn / umap-learn.
//!
//! Usage: `cargo run --release --example crosscheck`

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

use reed_reduce::extra::{lle, nmf, ppca, tsne, TsneConfig};
use reed_reduce::graph::{
    diffusion_map, isomap, kernel_pca, knn_graph, laplacian_eigenmaps, pcoa, Kernel,
};
use reed_reduce::pca::pca_dense;
use reed_reduce::umap::{fuzzy_simplicial_set, find_ab_params, umap, UmapConfig};
use reed_reduce::Mat;

/// Genotype-like fixture: 3 populations with different allele frequencies on
/// 400 SNPs. This is the shape the PCA path is actually designed for.
fn geno_blobs(n_per: usize, n_snp: usize) -> Mat {
    let mut state = 0x243F6A8885A308D3u64;
    let mut rnd = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let base: Vec<f64> = (0..n_snp).map(|_| 0.05 + 0.9 * rnd()).collect();
    Mat::from_fn(3 * n_per, n_snp, move |i, j| {
        let group = i / n_per;
        let shift = [0.0, 0.18, -0.15][group];
        let p = (base[j] + shift).clamp(0.001, 0.999);
        let u = rnd();
        if u < p * p {
            2.0
        } else if u < p * p + 2.0 * p * (1.0 - p) {
            1.0
        } else {
            0.0
        }
    })
}

fn centroids(m: &Mat, n_per: usize) -> Vec<Vec<f64>> {
    (0..3)
        .map(|g| {
            (0..m.cols)
                .map(|d| {
                    (g * n_per..(g + 1) * n_per).map(|i| m.get(i, d)).sum::<f64>()
                        / n_per as f64
                })
                .collect()
        })
        .collect()
}

fn within_spread(m: &Mat, n_per: usize) -> f64 {
    let cent = centroids(m, n_per);
    let mut s = 0.0;
    for g in 0..3 {
        for i in g * n_per..(g + 1) * n_per {
            let mut acc = 0.0;
            for d in 0..m.cols {
                acc += (m.get(i, d) - cent[g][d]).powi(2);
            }
            s += acc;
        }
    }
    (s / (3 * n_per) as f64).sqrt()
}

fn min_separation(m: &Mat, n_per: usize) -> f64 {
    let cent = centroids(m, n_per);
    let mut best = f64::MAX;
    for a in 0..3 {
        for b in (a + 1)..3 {
            let mut s = 0.0;
            for d in 0..m.cols {
                s += (cent[a][d] - cent[b][d]).powi(2);
            }
            best = best.min(s.sqrt());
        }
    }
    best
}

fn report(name: &str, m: &Mat, n_per: usize) {
    let w = within_spread(m, n_per);
    let s = min_separation(m, n_per);
    println!(
        "{name:<22} dims={} within={:.4} min_sep={:.4} ratio={:.2}",
        m.cols,
        w,
        s,
        s / w.max(1e-12)
    );
}

fn main() {
    let n_per = 20;
    let x = geno_blobs(n_per, 400);
    println!("fixture: {} samples, {} SNPs, 3 groups of {n_per}", x.rows, x.cols);

    // --- PCA family -------------------------------------------------------
    let p = pca_dense(&x, 2);
    let scores = p.scores();
    report("pca (binom2)", &scores, n_per);

    let d = reed_reduce::graph::pairwise_distances(&x);
    report("pcoa", &pcoa(&d, 2).coords, n_per);
    report("kernel_pca(rbf)", &kernel_pca(&x, 2, Kernel::Rbf, 0.05).coords, n_per);
    report("kernel_pca(ibs)", &kernel_pca(&x, 2, Kernel::Ibs, 1.0).coords, n_per);

    // --- graph family -----------------------------------------------------
    let g = knn_graph(&x, 10, 0.05);
    report("laplacian", &laplacian_eigenmaps(&g, 2).coords, n_per);
    report("diffusion", &diffusion_map(&g, 2, 1.0).coords, n_per);
    report("isomap", &isomap(&x, 10, 2).coords, n_per);

    // --- UMAP -------------------------------------------------------------
    let (a, b) = find_ab_params(1.0, 0.1);
    println!("umap find_ab_params(1.0, 0.1) = ({a:.9}, {b:.9})");
    let cfg = UmapConfig { n_neighbors: 10, n_epochs: 200, seed: 7, ..Default::default() };
    let fg = fuzzy_simplicial_set(&x, &cfg);
    let nnz = fg.weights.data.iter().filter(|v| **v > 0.0).count();
    println!("umap fuzzy graph: {} non-zeros", nnz);
    report("umap", &umap(&x, &cfg), n_per);

    // --- others -----------------------------------------------------------
    let tc = TsneConfig { perplexity: 10.0, n_iter: 300, seed: 3, ..Default::default() };
    report("tsne", &tsne(&x, &tc).expect("tsne"), n_per);
    report("lle", &lle(&x, 8, 2), n_per);
    let (w, h) = nmf(&x, 3, 200, 1);
    let recon = w.matmul(&h);
    let mut err = 0.0;
    let mut base = 0.0;
    for i in 0..x.rows {
        for j in 0..x.cols {
            err += (x.get(i, j) - recon.get(i, j)).powi(2);
            base += x.get(i, j).powi(2);
        }
    }
    println!("nmf rank=3: relative reconstruction error {:.4}", err / base);
    report("ppca", &ppca(&x, 2, 60, 5), n_per);
}

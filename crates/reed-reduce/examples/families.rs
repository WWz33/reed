//! Uniform conformance and timing harness for every reduction family.
//!
//! Each method runs on the same fixture and reports:
//!   * a rotation-invariant separation ratio (between-centroid / within spread)
//!   * wall time
//!
//! The ratios can be compared against scikit-learn / umap-learn on the same
//! fixture (see `docs/rounds/v16/STATUS.md` for the reference numbers).
//!
//! Usage: `cargo run --release --example families -- [small|medium]`

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
use reed_reduce::pca::{pca_dense, pca_dense_with, standardize, Grm, Stand};
use reed_reduce::umap::{umap, UmapConfig};
use reed_reduce::{lanczos_top_k, Mat, StreamConfig};
use std::time::Instant;

/// Genotype-like fixture: three populations with shifted allele frequencies.
fn geno(n_per: usize, n_snp: usize, seed: u64) -> Mat {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let base: Vec<f64> = (0..n_snp).map(|_| 0.05 + 0.9 * next()).collect();
    Mat::from_fn(3 * n_per, n_snp, move |i, j| {
        let group = i / n_per;
        let shift = [0.0, 0.18, -0.15][group];
        let p = (base[j] + shift).clamp(0.001, 0.999);
        let u = next();
        if u < p * p {
            2.0
        } else if u < p * p + 2.0 * p * (1.0 - p) {
            1.0
        } else {
            0.0
        }
    })
}

fn separation(m: &Mat, n_per: usize) -> f64 {
    let dim = m.cols;
    if dim == 0 {
        return 0.0;
    }
    let mut cent = vec![vec![0.0f64; dim]; 3];
    for g in 0..3 {
        for i in g * n_per..(g + 1) * n_per {
            for d in 0..dim {
                cent[g][d] += m.get(i, d) / n_per as f64;
            }
        }
    }
    let mut within = 0.0f64;
    for g in 0..3 {
        for i in g * n_per..(g + 1) * n_per {
            let mut s = 0.0;
            for d in 0..dim {
                s += (m.get(i, d) - cent[g][d]).powi(2);
            }
            within += s;
        }
    }
    let within = (within / (3 * n_per) as f64).sqrt();
    let mut min_sep = f64::MAX;
    for a in 0..3 {
        for b in (a + 1)..3 {
            let mut s = 0.0;
            for d in 0..dim {
                s += (cent[a][d] - cent[b][d]).powi(2);
            }
            min_sep = min_sep.min(s.sqrt());
        }
    }
    min_sep / within.max(1e-12)
}

fn report(name: &str, m: &Mat, n_per: usize, secs: f64) {
    println!("{name:<26} {:8.3}s  ratio={:10.2}", secs, separation(m, n_per));
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let scale = args.get(1).map(|s| s.as_str()).unwrap_or("small");
    let (n_per, n_snp) = match scale {
        "medium" => (100, 2000),
        _ => (20, 400),
    };
    let x = geno(n_per, n_snp, 0x243F6A8885A308D3);
    println!(
        "fixture: {} samples x {} SNPs, 3 groups of {n_per}",
        x.rows, x.cols
    );
    println!("backend: {}", reed_reduce::COMPILED_BACKEND.name());
    println!();

    // --- PCA family -------------------------------------------------------
    let t = Instant::now();
    let p = pca_dense(&x, 2);
    report("pca_dense (binom2,mi)", &p.scores(), n_per, t.elapsed().as_secs_f64());

    let t = Instant::now();
    let pb = pca_dense_with(&x, 2, Stand::Binom);
    report("pca_dense (binom,mi)", &pb.scores(), n_per, t.elapsed().as_secs_f64());

    let st = standardize(&x, Stand::Binom2);
    let cfg = StreamConfig { k: 2, tol: 1e-10, max_cycles: 200, ..Default::default() };
    let t = Instant::now();
    let s = lanczos_top_k(&st.w, &cfg);
    let scores = Mat::from_fn(s.vectors.rows, s.vectors.cols, |i, c| {
        s.vectors.get(i, c) * s.values[c].max(0.0).sqrt()
    });
    report("stream lanczos (binom2)", &scores, n_per, t.elapsed().as_secs_f64());

    let t = Instant::now();
    let g = Grm::mean_impute(&st.w);
    let gp = reed_reduce::grm_pcoa(&g, 2);
    report("grm_pcoa", &gp.coords, n_per, t.elapsed().as_secs_f64());

    // --- distance / kernel family ----------------------------------------
    let d = reed_reduce::graph::pairwise_distances(&x);
    let t = Instant::now();
    report("pcoa (euclid)", &pcoa(&d, 2).coords, n_per, t.elapsed().as_secs_f64());

    let t = Instant::now();
    report(
        "kernel_pca (rbf)",
        &kernel_pca(&x, 2, Kernel::Rbf, 0.0).coords,
        n_per,
        t.elapsed().as_secs_f64(),
    );
    let t = Instant::now();
    report(
        "kernel_pca (ibs)",
        &kernel_pca(&x, 2, Kernel::Ibs, 1.0).coords,
        n_per,
        t.elapsed().as_secs_f64(),
    );
    let t = Instant::now();
    report(
        "grm_kernel_pca",
        &reed_reduce::grm_kernel_pca(&g, 2).coords,
        n_per,
        t.elapsed().as_secs_f64(),
    );

    // --- graph family -----------------------------------------------------
    let w = knn_graph(&x, 10, 0.0);
    let t = Instant::now();
    report("laplacian", &laplacian_eigenmaps(&w, 2).coords, n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    report("diffusion", &diffusion_map(&w, 2, 1.0).coords, n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    report("isomap", &isomap(&x, 10, 2).coords, n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let gl = reed_reduce::grm_laplacian(&g, 2, 10);
    report("grm_laplacian", &gl.coords, n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    report("grm_diffusion", &reed_reduce::grm_diffusion(&g, 2, 10, 1.0).coords, n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let grm_tsne_cfg = TsneConfig { perplexity: 10.0, n_iter: 300, seed: 3, ..Default::default() };
    report("grm_tsne", &reed_reduce::grm_tsne(&g, &grm_tsne_cfg), n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let grm_umap_cfg = UmapConfig { n_neighbors: 10, n_epochs: 200, seed: 7, ..Default::default() };
    report("grm_umap", &reed_reduce::grm_umap(&g, &grm_umap_cfg), n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    report("grm_lle", &reed_reduce::grm_lle(&g, &st.w, 8, 2), n_per, t.elapsed().as_secs_f64());

    // --- manifold learning ------------------------------------------------
    let t = Instant::now();
    let ll = lle(&x, 8, 2);
    report("lle", &ll, n_per, t.elapsed().as_secs_f64());

    // --- probabilistic / matrix factorisation -----------------------------
    let t = Instant::now();
    report("ppca (em)", &ppca(&x, 2, 60, 5), n_per, t.elapsed().as_secs_f64());
    let t = Instant::now();
    let (wm, hm) = nmf(&x, 3, 200, 1);
    let recon = wm.matmul(&hm);
    let mut err = 0.0;
    let mut base = 0.0;
    for i in 0..x.rows {
        for j in 0..x.cols {
            err += (x.get(i, j) - recon.get(i, j)).powi(2);
            base += x.get(i, j).powi(2);
        }
    }
    println!(
        "{:<26} {:8.3}s  rel-recon-error={:.4}",
        "nmf (rank 3)",
        t.elapsed().as_secs_f64(),
        err / base
    );

    // --- nonlinear embeddings ---------------------------------------------
    let cfg = UmapConfig { n_neighbors: 10, n_epochs: 200, seed: 7, ..Default::default() };
    let t = Instant::now();
    report("umap", &umap(&x, &cfg), n_per, t.elapsed().as_secs_f64());
    let perplexity = (x.rows as f64 / 8.0).clamp(5.0, 50.0);
    let tc = TsneConfig { perplexity, n_iter: 500, seed: 3, ..Default::default() };
    let t = Instant::now();
    report("tsne", &tsne(&x, &tc), n_per, t.elapsed().as_secs_f64());
}

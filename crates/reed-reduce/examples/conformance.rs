//! Conformance example: prints the eigenvalue vector for each semantic
//! combination so it can be diffed against plink2 / flashpca output.
//!
//! Usage: `cargo run --release --example conformance -- <matrix> <k>`
//!
//! The four combinations are `--stand {binom2,binom} x --missing
//! {pairwise,meanimpute}`. Reference commands:
//!
//! ```text
//! plink2 --bfile X --pca k                 # pairwise, binom2
//! plink2 --bfile X --pca k meanimpute      # meanimpute, binom2
//! flashpca --bfile X --ndim k --standx binom2 --div p
//! flashpca --bfile X --ndim k --standx binom  --div p   # = 2x the binom2 values
//! ```

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

use reed_reduce::pca::{
    pca_dense_pairwise_with, pca_dense_with, standardize, Grm, Stand,
};
use reed_reduce::{lanczos_top_k, Mat, StreamConfig};
use std::time::Instant;

fn read_matrix(path: &str) -> Mat {
    let txt = std::fs::read_to_string(path).expect("read matrix");
    let mut it = txt.split_whitespace();
    let n: usize = it.next().unwrap().parse().unwrap();
    let m: usize = it.next().unwrap().parse().unwrap();
    let mut mat = Mat::zeros(n, m);
    for i in 0..n {
        for j in 0..m {
            let tok = it.next().unwrap();
            mat.set(i, j, if tok == "NA" { f64::NAN } else { tok.parse().unwrap() });
        }
    }
    mat
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let k: usize = args[2].parse().unwrap();
    let x = read_matrix(path);
    println!("matrix {} x {}, k={}", x.rows, x.cols, k);
    println!("backend: {}", reed_reduce::COMPILED_BACKEND.name());

    for stand in [Stand::Binom2, Stand::Binom] {
        // mean imputation
        let t = Instant::now();
        let p = pca_dense_with(&x, k, stand);
        println!(
            "\nstand={:<7} missing=meanimpute  {:.3}s\n  lam = {:?}",
            stand.name(),
            t.elapsed().as_secs_f64(),
            p.values
        );
        // pairwise
        let t = Instant::now();
        let pp = pca_dense_pairwise_with(&x, k, stand);
        println!(
            "stand={:<7} missing=pairwise     {:.3}s\n  lam = {:?}",
            stand.name(),
            t.elapsed().as_secs_f64(),
            pp.values
        );
    }

    // Cross-check: binom should be exactly binom2 / 2 (same eigenvectors).
    let a = pca_dense_with(&x, k, Stand::Binom2);
    let b = pca_dense_with(&x, k, Stand::Binom);
    let ratio: Vec<f64> = a.values.iter().zip(&b.values).map(|(u, v)| u / v).collect();
    println!("\nlambda(binom2)/lambda(binom) per PC = {:?}", ratio);

    // Stream path on the default convention.
    let st = standardize(&x, Stand::Binom2);
    let cfg = StreamConfig { k, tol: 1e-10, max_cycles: 200, ..Default::default() };
    let t = Instant::now();
    let s = lanczos_top_k(&st.w, &cfg);
    println!(
        "\nstream (binom2, meanimpute): {:.3}s, {} cycles, residual {:.2e}\n  lam = {:?}",
        t.elapsed().as_secs_f64(),
        s.cycles,
        s.max_residual,
        s.values
    );
    let worst = a
        .values
        .iter()
        .zip(&s.values)
        .map(|(u, v)| (u - v).abs() / u.abs().max(1e-12))
        .fold(0.0f64, f64::max);
    println!("stream vs dense max rel diff: {:.3e}", worst);

    // GRM pair counts (the pairwise-complete denominator).
    let g = Grm::pairwise(&x);
    if let Some(counts) = g.n_pairs {
        let min = counts.data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = counts.data.iter().cloned().fold(0.0f64, f64::max);
        println!("pairwise counts: min {min}, max {max}, m = {}", g.m);
    }
}

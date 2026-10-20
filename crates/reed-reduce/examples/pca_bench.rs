//! Benchmark: PCA (dense and Lanczos) on a matrix file, with timings and
//! residual diagnostics. Input format: `n m` then n rows of m values, `NA`
//! for missing.
//!
//! Usage: `cargo run --release --example pca_bench -- <matrix> <k> [dense|lanczos]`

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

use reed_reduce::lanczos::{lanczos_top_k, LanczosConfig};
use reed_reduce::pca::{pca_dense, pca_dense_pairwise, pca_top_k, standardize_binom2};
use reed_reduce::Mat;
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
            let v = if tok == "NA" { f64::NAN } else { tok.parse().unwrap() };
            mat.set(i, j, v);
        }
    }
    mat
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: pca_bench <matrix> <k> [dense|lanczos|pairwise]");
        return;
    }
    let path = &args[1];
    let k: usize = args[2].parse().unwrap();
    let mode = args.get(3).map(|s| s.as_str()).unwrap_or("both");

    let t = Instant::now();
    let x = read_matrix(path);
    println!("read: {:.3} s  ({} x {})", t.elapsed().as_secs_f64(), x.rows, x.cols);

    let t = Instant::now();
    let st = standardize_binom2(&x);
    println!("standardise (binom2): {:.3} s", t.elapsed().as_secs_f64());

    if mode == "dense" || mode == "both" {
        let t = Instant::now();
        let p = pca_dense(&x, k);
        let secs = t.elapsed().as_secs_f64();
        println!("dense PCA: {:.3} s", secs);
        println!("  lam: {:?}", &p.values[..k.min(p.values.len())]);
        println!("EXACT {:?}", &p.values[..k.min(p.values.len())]);
        println!("  pve: {:?}", &p.pve[..k.min(p.pve.len())]);
    }

    if mode == "topk" || mode == "both" {
        let t = Instant::now();
        let p = pca_top_k(&x, k);
        let secs = t.elapsed().as_secs_f64();
        println!("top-k PCA (bisection+invit): {:.3} s", secs);
        println!("  lam: {:?}", &p.values[..k.min(p.values.len())]);
    }

    if mode == "pairwise" {
        let t = Instant::now();
        let p = pca_dense_pairwise(&x, k);
        println!("pairwise PCA: {:.3} s", t.elapsed().as_secs_f64());
        println!("  lam: {:?}", &p.values[..k.min(p.values.len())]);
        println!("EXACT {:?}", &p.values[..k.min(p.values.len())]);
    }

    if mode == "lanczos" || mode == "both" {
        let cfg = LanczosConfig { k, block: 8, ..Default::default() };
        let t = Instant::now();
        let part = lanczos_top_k(&st.w, &cfg);
        let secs = t.elapsed().as_secs_f64();
        println!("stream PCA: {:.3} s  ({} restart cycles)", secs, part.cycles);
        println!("  lam: {:?}", part.values);
        println!("  max rel residual: {:.3e}", part.max_residual);
    }
}

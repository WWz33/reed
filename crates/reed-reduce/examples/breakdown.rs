//! Timing breakdown of the PCA pipeline, to attribute cost to each stage.
//!
//! Usage: `cargo run --release --example breakdown -- <matrix> <k>`

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

use reed_reduce::backend::{sym_eig, sym_eig_top_k};
use reed_reduce::pca::{standardize_binom2, Grm};
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
            mat.set(i, j, if tok == "NA" { f64::NAN } else { tok.parse().unwrap() });
        }
    }
    mat
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let k: usize = args[2].parse().unwrap();

    let t = Instant::now();
    let x = read_matrix(path);
    println!("read:              {:8.3} s  ({} x {})", t.elapsed().as_secs_f64(), x.rows, x.cols);

    let t = Instant::now();
    let st = standardize_binom2(&x);
    println!("standardise:       {:8.3} s", t.elapsed().as_secs_f64());

    let t = Instant::now();
    let grm = Grm::mean_impute(&st.w);
    let grm_secs = t.elapsed().as_secs_f64();
    println!("GRM = W W'/m:      {:8.3} s", grm_secs);

    let t = Instant::now();
    let (vals, _vecs, _be) = sym_eig_top_k(&grm.a, k);
    let partial_secs = t.elapsed().as_secs_f64();
    println!("partial eig:       {:8.3} s   lam1 = {:.9}", partial_secs, vals[0]);

    // Dense full spectrum on the same matrix, for the comparison.
    let t = Instant::now();
    let full = sym_eig(&grm.a);
    let dense_secs = t.elapsed().as_secs_f64();
    println!("dense full eig:    {:8.3} s   lam1 = {:.9}", dense_secs, full.values[full.values.len() - 1]);

    println!();
    println!("GRM share of pipeline: {:.1}%", 100.0 * grm_secs / (grm_secs + partial_secs));
    println!("partial vs dense solver: {:.2}x", dense_secs / partial_secs);
}

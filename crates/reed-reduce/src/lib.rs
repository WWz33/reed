//! Dimension reduction for population-genetics work.
//!
//! PCA with the binom2 standardization (the convention shared by GCTA, plink2
//! and flashpca), plus PCoA, kernel PCA, Laplacian eigenmaps, diffusion maps,
//! Isomap, LLE, NMF, PPCA-EM, UMAP and t-SNE.
//!
//! The dense linear algebra goes through [`backend`], which selects between
//! `faer` (default) and a dependency-free fallback. Both produce identical
//! results; the fallback exists so every number can be audited against a
//! from-scratch implementation, and so the crate still builds where no
//! dependency is allowed.
//!
//! # Lint policy
//!
//! The algorithms here are array mathematics, so the loops are written as
//! `for i in 0..n { a[i] ... }` over several parallel arrays. Clippy's
//! `needless_range_loop` would rewrite those into zipped iterators, which is
//! strictly harder to read against the published formulas they implement
//! (tred2/tql2, the UMAP gradient, the Laplacian normalization). The lint is
//! therefore allowed for the crate, together with the arithmetic-cast lints
//! that fire on index-to-float conversions. `needless_return` and
//! `useless_vec` are allowed for the same reason: the `#[cfg]`-selected
//! backend branches need explicit returns to type-check identically in both
//! configurations.
#![allow(clippy::needless_range_loop)]
#![allow(clippy::needless_return)]
#![allow(clippy::needless_borrow)]
#![allow(clippy::useless_vec)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::similar_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::suboptimal_flops)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::excessive_precision)]
#![allow(clippy::if_same_then_else)]

pub mod backend;
pub mod eig;
pub mod extra;
pub mod graph;
pub mod grm_reduce;

pub mod lanczos;
pub mod lda;
pub mod pca;
pub mod umap;
pub mod xval;

pub use backend::{gram_divide, sym_eig, sym_eig_top_k, Backend, SymEigResult, COMPILED_BACKEND};
pub use eig::SymEig;
pub use extra::{lle, lle_from_sq_distances, nmf, ppca, tsne, tsne_from_sq_distances, TsneConfig};
pub use grm_reduce::{grm_diffusion, grm_distances, grm_isomap, grm_kernel, grm_kernel_pca, grm_knn_weights, grm_laplacian, grm_lle, grm_pcoa, grm_similarity, grm_tsne, grm_umap};
pub use graph::{diffusion_map, isomap, kernel_pca, knn_graph, laplacian_eigenmaps, pcoa, pairwise_distances, Embedding, Kernel};
pub use reed_core::linalg::Mat;
pub use lanczos::{lanczos_top_k, LanczosConfig, PartialEig, StreamConfig, StreamEig};
pub use lda::{fit as fit_lda, Lda, LdaError, LdaProjection};
pub use xval::{candidate_counts, XvalConfig, XvalError, XvalReport, XvalResult};
pub use umap::{fuzzy_simplicial_set, fuzzy_simplicial_set_from_distances, umap, umap_from_distances, UmapConfig};
pub use pca::{pca_dense, pca_dense_pairwise, pca_dense_pairwise_with, pca_dense_with, pca_top_k, pca_top_k_with, standardize, standardize_binom2, Grm, Pca, Stand, Standardized};

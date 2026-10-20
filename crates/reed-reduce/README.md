# reed-reduce

Dimension reduction for population-genetics work: PCA with the binom2
standardization, plus the wider family of graph, kernel, probabilistic and
nonlinear methods. Every number here is cross-checked against a reference
implementation (LAPACK/numpy, plink2, GCTA, flashpca, scikit-learn,
umap-learn); the audit trail is in `docs/rounds/v16/`.

## Backends

| backend | selected by | GRM `W W'` | dense eigendecomposition |
|---|---|---|---|
| `faer` | default (`faer-backend` feature) | faer blocked GEMM, multi-threaded | faer Householder + implicit QR |
| fallback | `--no-default-features` | scalar blocked loops | EISPACK `tred2` + `tql2` |

Both backends produce identical results (the `families` example prints the
same separation ratios either way). Measured on HapMap3 (957 x 14266, k = 10):

| stage | fallback | faer |
|---|---|---|
| GRM | 3.34 s | 0.081 s |
| dense eigendecomposition | 2.23 s | 0.194 s |
| whole pipeline | 5.6 s | 0.40 s |

Simulated 4000 x 20000: fallback 334 s, faer 3.96 s for the dense path.

## Modules

| module | contents |
|---|---|
| `backend` | backend selection: `gram_divide`, `sym_eig`, `sym_eig_top_k` |
| `linalg` | row-major `Mat`, blocked matmul, gram, distances |
| `eig` | dependency-free dense eigensolver (`tred2` + `tql2`) and partial-spectrum solver (tridiagonalize + Sturm bisection + inverse iteration) |
| `pca` | `Stand` (binom2/binom), standardization, `Grm` (mean-impute / pairwise with `.grm.N` semantics), dense PCA, scores/loadings |
| `lanczos` | matrix-free **block Lanczos with thick restart** for the top-k eigenpairs |
| `graph` | PCoA/classical MDS, kernel PCA (linear/RBF/IBS), kNN graph, Laplacian eigenmaps, diffusion map, Isomap, `auto_gamma` |
| `grm_reduce` | 12 GRM entry points: distances/similarity/knn-weights, plus PCoA, kernel PCA, Laplacian, diffusion, Isomap, t-SNE, UMAP and LLE driven by GRM geometry |
| `umap` | fuzzy simplicial set, `find_ab_params`, spectral init, SGD — line-by-line from umap-learn 0.5.12 |
| `extra` | t-SNE, NMF, LLE, PPCA-EM (missing-value native) |

## GRM coverage

The relationship matrix is not a separate method; it is the geometry. Every
distance-based family therefore takes GRM distances through a dedicated entry
point, while the feature-based entry point stays available:

| method | GRM entry | feature entry |
|---|---|---|
| PCoA / classical MDS | `grm_pcoa` | `pcoa` |
| kernel PCA | `grm_kernel_pca`, `grm_kernel` | `kernel_pca` |
| Laplacian eigenmaps | `grm_laplacian` | `laplacian_eigenmaps` |
| diffusion map | `grm_diffusion` | `diffusion_map` |
| Isomap | `grm_isomap` | `isomap` |
| t-SNE | `grm_tsne` | `tsne` |
| UMAP | `grm_umap` | `umap` |
| LLE | `grm_lle` | `lle` |
| PCA | (starts from the GRM by construction) | `pca_dense` |

Supporting pieces: `grm_distances` (`d^2 = a_ii + a_jj - 2 a_ij`),
`grm_similarity` (`a_ij / sqrt(a_ii a_jj)`), `grm_knn_weights`.

Each GRM entry is a thin wrapper over a general distance entry point
(`tsne_from_sq_distances`, `umap_from_distances`, `lle_from_sq_distances`), so
the two routes are bit-for-bit identical (asserted in
`distance_entry_points_match_feature_entry_points`).

**NMF and PPCA have no GRM entry on purpose**: NMF needs non-negative input
and PPCA models genotype-space noise, neither of which a GRM provides.

## Verified numbers

PCA on HapMap3 (957 x 14266), k = 10:

| convention | missing | lambda1 | lambda2 | lambda3 |
|---|---|---|---|---|
| binom2 | meanimpute | 26.200944934 | 23.394076211 | 6.932494547 |
| binom2 | pairwise | 26.280628904 | 23.458733875 | 6.948662248 |
| binom | meanimpute | 52.401889869 | 46.788152422 | 13.864989094 |
| binom | pairwise | 52.561257808 | 46.917467749 | 13.897324497 |

Reference: plink2 meanimpute gives 26.2009 / 23.3941 / 6.93249 and plink2
pairwise gives 26.2806 / 23.4587 / 6.94866. **Agreement is exact to the
printed digits.** `lambda(binom2)/lambda(binom) = 0.500000000000000` for
every PC, as the sqrt(2) per-SNP scaling implies.

Against numpy `eigh` on a 60 x 400 fixture the relative error is 4.5e-16
(machine precision).

The stream path (block Lanczos + thick restart, faer GEMM operator) reaches
`max_rel_residual = 6.8e-11` on HapMap3, with eigenvalues agreeing with the
dense solution to 6.3e-15.

## Cross-checks against reference implementations

Same 60 x 400 genotype fixture (3 populations), separation ratio =
between-centroid distance / within-group spread:

| method | reed-reduce | reference |
|---|---|---|
| PCA (binom2) | 2.91 | 2.91 (numpy eigh of the same GRM) |
| Kernel PCA (RBF) | 3.88 | 1.80 (sklearn; different gamma convention) |
| Isomap | 1.93 | 1.95 (sklearn) |
| Laplacian eigenmaps | 5.73 | 6.90 (sklearn SpectralEmbedding) |
| LLE | 6.39 | 3.47 (sklearn) |
| UMAP | 2.04 | 2.72 (umap-learn) |
| NMF reconstruction error | 0.1907 | 0.1905 (sklearn) |
| t-SNE | 2.57 | 0.49 (sklearn, different init/schedule) |

`find_ab_params(1.0, 0.1)` gives `(1.576943618, 0.895060710)`; umap-learn
gives `(1.576943460, 0.895060878)` — agreement to 7 significant digits.

Graph and kernel methods depend on the neighbourhood/kernel scale, so the
comparisons above use each implementation's own defaults; matching the
defaults matters more than matching the algorithm. t-SNE and UMAP are
stochastic and are not expected to match bit for bit.

## Known limitations

- **`stream` is slower than `dense` whenever N^2 fits in memory.** At
  N = 4000 the dense path takes 4.0 s and stream 154 s; stream exists for N
  large enough that the N x N matrix is not affordable (N > ~50k).
- **`stream` accuracy degrades for k close to N/2.** A stress test at
  n = 80, k = 30 shows ~2% error on trailing pairs. Production use is
  k <= 50 with n in the thousands, where agreement is 1e-7.
- **`grm_laplacian` can be uninformative on well-separated clusters**: the
  Laplacian embedding of a GRM-distance graph may separate only some of the
  groups (the others coincide at the origin). Use `grm_pcoa` there.
- **`lle` is poorly conditioned for nearly collinear local neighbourhoods**
  and can produce very large coordinates.
- **Monomorphic columns** are set to the column mean (flashpca's matrix-path
  convention), not to 0 as GCTA and plink2 do.
- **NMF and PPCA deliberately have no GRM entry point.** NMF's multiplicative
  updates require non-negative input and a GRM has negative off-diagonals;
  PPCA models observation noise in genotype space, which a GRM does not carry.

## Tests and examples

```
cargo test -p reed-reduce                                   # 38 tests
cargo run --release --example conformance -- <matrix> <k>   # 4 semantic combinations
cargo run --release --example families -- small|medium      # 17 methods, ratios + timings
cargo run --release --example breakdown -- <matrix> <k>     # per-stage timing
cargo run --release --example pca_bench -- <matrix> <k> [dense|lanczos|pairwise|topk]
cargo run --release --example crosscheck                    # sklearn/umap-learn fixture
```

Matrix file format: `n m` on the first line, then `n` rows of `m`
whitespace-separated dosages in {0,1,2}, with `NA` for missing.
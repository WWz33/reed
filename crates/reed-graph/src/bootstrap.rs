//! Bootstrap support for tree branches.
//!
//! Resamples **loci** (columns), rebuilds the tree from the resampled distance
//! matrix, and counts how often each branch of the original tree reappears.
//! This is `ape::boot.phylo`'s procedure, which `poppr::aboot` delegates to.
//!
//! # Why loci, not samples
//!
//! Bootstrapping individuals would answer "how stable is the tree if I had
//! sampled different individuals", which is not the question a phylogeny asks.
//! Resampling loci answers "how stable is the tree given the markers I chose".
//! The two give very different support values on the same data, and only the
//! second is what `aboot` reports.

use std::collections::HashMap;

use reed_core::linalg::Mat;
use reed_core::rng::Rng;

use crate::error::GraphError;
use crate::tree::{nj, upgma, Tree};

/// Bootstrap configuration.
#[derive(Clone, Debug)]
pub struct BootstrapConfig {
    /// Replicates.
    pub replicates: usize,
    /// Seed for the per-replicate sub-streams.
    pub seed: u64,
    /// Below this percentage the support is reported as `None`.
    pub cutoff: f64,
}

impl Default for BootstrapConfig {
    fn default() -> Self {
        Self {
            replicates: 100,
            seed: 42,
            cutoff: 0.0,
        }
    }
}

/// Bootstrap outcome for one tree.
#[derive(Clone, Debug)]
pub struct BootstrapReport {
    /// Support percentage per internal node, in the same order as
    /// [`Tree::node_labels`] would be filled: ascending node id.
    pub support: Vec<f64>,
    /// Replicates actually built.
    pub replicates: usize,
}

/// Which tree the support is computed for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TreeMethod {
    /// UPGMA (rooted, so support is by clade).
    Upgma,
    /// Neighbor-joining (unrooted, so support is by bipartition).
    Nj,
}

impl TreeMethod {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// A message naming the accepted values.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "upgma" => Ok(Self::Upgma),
            "nj" => Ok(Self::Nj),
            other => Err(format!("unknown tree method {other:?} (upgma|nj)")),
        }
    }

    /// Builds a tree with this method.
    fn build(self, d: &Mat) -> Result<Tree, GraphError> {
        match self {
            Self::Upgma => upgma(d),
            Self::Nj => nj(d),
        }
    }

    /// Whether the resulting tree is rooted.
    #[must_use]
    pub fn rooted(self) -> bool {
        matches!(self, Self::Upgma)
    }
}

/// Computes bootstrap support for `tree`, rebuilding from `genotypes`.
///
/// `genotypes` is `n_samples x n_sites`; `tip_labels[i]` names row `i` of the
/// tree, and the rows of `genotypes` must be in the same order. `distance_fn`
/// turns a resampled genotype matrix into a distance matrix — the caller
/// supplies it so this module does not depend on `reed-distance`.
///
/// # Errors
///
/// `TooFewTaxa` when there are fewer than four samples (a tree over fewer has
/// no internal branch to support), and whatever the tree builder reports.
pub fn support(
    tree: &Tree,
    genotypes: &Mat,
    method: TreeMethod,
    distance_fn: impl Fn(&Mat) -> Mat,
    config: &BootstrapConfig,
) -> Result<BootstrapReport, GraphError> {
    let n = genotypes.rows;
    if n < 4 {
        return Err(GraphError::TooFewTaxa { got: n, need: 4 });
    }
    if genotypes.cols == 0 {
        return Err(GraphError::TooFewTaxa { got: 0, need: 1 });
    }

    // The branches of the original tree, keyed for matching.
    let target = branch_keys(tree, method.rooted());

    // Count how often each target branch appears among the replicates.
    let mut counts: HashMap<Vec<u64>, usize> = HashMap::new();
    let base = Rng::from_seed(config.seed);
    let mut built = 0usize;

    for rep in 0..config.replicates {
        let mut rng = base.derive(rep as u64);
        // Resample columns with replacement: the same locus can appear twice,
        // and another can be absent.
        let picks = rng.resample_n(genotypes.cols, genotypes.cols);
        let mut resampled = Mat::zeros(n, picks.len());
        for (dst, &src) in picks.iter().enumerate() {
            for r in 0..n {
                resampled.set(r, dst, genotypes.get(r, src));
            }
        }

        let dist = distance_fn(&resampled);
        // A replicate whose distance matrix is degenerate (every resampled
        // locus monomorphic) has no tree; skip it rather than counting a
        // spurious zero. `ape` would error here; skipping keeps the run going
        // and is reported through `replicates`.
        let Ok(replicate) = method.build(&dist) else {
            continue;
        };
        built += 1;
        for key in branch_keys(&replicate, method.rooted()).into_keys() {
            if target.contains_key(&key) {
                *counts.entry(key).or_insert(0) += 1;
            }
        }
    }

    // Percentage per internal node, in ascending node-id order.
    let mut support = Vec::with_capacity(target.len());
    let mut ids: Vec<u32> = target.values().copied().collect();
    ids.sort_unstable();
    let by_id: HashMap<u32, &Vec<u64>> = target.iter().map(|(k, &v)| (v, k)).collect();
    for id in ids {
        let key = by_id.get(&id).expect("invariant: id came from the map");
        let hit = counts.get(*key).copied().unwrap_or(0);
        let pct = if built == 0 {
            0.0
        } else {
            100.0 * hit as f64 / built as f64
        };
        support.push(if pct >= config.cutoff { pct } else { f64::NAN });
    }

    Ok(BootstrapReport { support, replicates: built })
}

/// Maps each internal node of `tree` to its branch key.
///
/// The returned map is key -> node id, so the caller can order the support
/// values by node.
fn branch_keys(tree: &Tree, rooted: bool) -> HashMap<Vec<u64>, u32> {
    let n = tree.n_tips();
    let words = n.div_ceil(64);

    let max_id = tree.max_node_id() as usize;
    let mut children: Vec<Vec<u32>> = vec![Vec::new(); max_id + 1];
    for &(p, c) in &tree.edges {
        children[p as usize].push(c);
    }

    // Descendant leaves of every node, computed bottom-up by node id. Both
    // constructors create children before parents, so an ascending scan works.
    let mut leaves: Vec<Vec<u64>> = vec![vec![0u64; words]; max_id + 1];
    for id in 1..=max_id {
        if id <= n {
            let idx = id - 1;
            leaves[id][idx / 64] |= 1u64 << (idx % 64);
        } else {
            for &c in &children[id] {
                // Copy the child's words out first: `leaves` is borrowed
                // immutably by the read and mutably by the write.
                let child = leaves[c as usize].clone();
                for (dst, src) in leaves[id].iter_mut().zip(child.iter()) {
                    *dst |= src;
                }
            }
        }
    }

    let mut out = HashMap::new();
    for id in 1..=max_id {
        if id <= n {
            continue;
        }
        let mut key = leaves[id].clone();
        if !rooted {
            // Normalise the bipartition: record the side containing tip 0.
            // Exactly one side of a split does, so the key no longer depends
            // on which side the traversal happened to collect.
            if key[0] & 1 == 0 {
                for (w, word) in key.iter_mut().enumerate() {
                    let remaining = n.saturating_sub(w * 64).min(64);
                    let mask = if remaining == 64 {
                        u64::MAX
                    } else {
                        (1u64 << remaining) - 1
                    };
                    *word = !*word & mask;
                }
            }
        }
        out.insert(key, id as u32);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn levels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("s{i}")).collect()
    }

    /// A genotype matrix with a strong two-block structure: sites 0..k favour
    /// block A, the rest favour block B, so the tree is stable under
    /// resampling.
    fn structured(n: usize, m: usize) -> Mat {
        let mut x = Mat::zeros(n, m);
        for r in 0..n {
            let block = if r < n / 2 { 0.0 } else { 1.0 };
            for c in 0..m {
                // The first half of the sites is informative, the rest noise.
                let v = if c < m / 2 {
                    block * 2.0
                } else {
                    ((r * 7 + c * 13) % 3) as f64
                };
                x.set(r, c, v);
            }
        }
        x
    }

    /// Euclidean distance between rows of a dosage matrix.
    fn euclid(x: &Mat) -> Mat {
        let n = x.rows;
        let mut d = Mat::zeros(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let mut acc = 0.0;
                for c in 0..x.cols {
                    let dv = x.get(i, c) - x.get(j, c);
                    acc += dv * dv;
                }
                let v = acc.sqrt();
                d.set(i, j, v);
                d.set(j, i, v);
            }
        }
        d
    }

    fn tree_of(x: &Mat, method: TreeMethod) -> Tree {
        let mut t = method.build(&euclid(x)).expect("tree");
        t.tip_labels = levels(x.rows);
        t
    }

    #[test]
    fn strong_branches_get_high_support() {
        let x = structured(12, 40);
        let tree = tree_of(&x, TreeMethod::Upgma);
        let cfg = BootstrapConfig {
            replicates: 50,
            ..Default::default()
        };
        let report = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
        assert_eq!(report.support.len(), tree.n_internal());
        // The two-block split is the root's child branch; it must be near 100.
        let best = report
            .support
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(best > 90.0, "no branch reached high support: {:?}", report.support);
    }

    #[test]
    fn support_is_a_percentage() {
        let x = structured(8, 30);
        let tree = tree_of(&x, TreeMethod::Upgma);
        let cfg = BootstrapConfig {
            replicates: 20,
            ..Default::default()
        };
        let report = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
        for &s in &report.support {
            assert!(
                s.is_nan() || (0.0..=100.0).contains(&s),
                "support {s} out of range"
            );
        }
        assert_eq!(report.replicates, 20);
    }

    #[test]
    fn nj_uses_bipartitions_not_clades() {
        // On an unrooted tree the same split must be counted regardless of
        // which side the traversal collects. Two runs with different seeds
        // must both find the strong split.
        let x = structured(10, 40);
        let tree = tree_of(&x, TreeMethod::Nj);
        let cfg = BootstrapConfig {
            replicates: 30,
            seed: 7,
            ..Default::default()
        };
        let a = support(&tree, &x, TreeMethod::Nj, euclid, &cfg).unwrap();
        let cfg2 = BootstrapConfig { seed: 99, ..cfg };
        let b = support(&tree, &x, TreeMethod::Nj, euclid, &cfg2).unwrap();
        assert_eq!(a.support.len(), tree.n_internal());
        assert_eq!(b.support.len(), tree.n_internal());
        let best_a = a.support.iter().copied().fold(0.0f64, f64::max);
        let best_b = b.support.iter().copied().fold(0.0f64, f64::max);
        assert!(best_a > 80.0 && best_b > 80.0, "{best_a} {best_b}");
    }

    #[test]
    fn cutoff_turns_low_support_into_nan() {
        let x = structured(10, 40);
        let tree = tree_of(&x, TreeMethod::Upgma);
        let cfg = BootstrapConfig {
            replicates: 20,
            cutoff: 99.0,
            ..Default::default()
        };
        let report = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
        // With a 99% cutoff, branches below it are NaN; at least one should be.
        assert!(
            report.support.iter().any(|s| s.is_nan()),
            "cutoff 99 left everything above it: {:?}",
            report.support
        );
    }

    #[test]
    fn is_reproducible_for_a_fixed_seed() {
        let x = structured(10, 30);
        let tree = tree_of(&x, TreeMethod::Upgma);
        let cfg = BootstrapConfig {
            replicates: 15,
            seed: 5,
            ..Default::default()
        };
        let a = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
        let b = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
        assert_eq!(a.support, b.support);
    }

    #[test]
    fn different_seeds_move_the_support() {
        let x = structured(10, 30);
        let tree = tree_of(&x, TreeMethod::Upgma);
        let cfg = BootstrapConfig {
            replicates: 15,
            seed: 1,
            ..Default::default()
        };
        let a = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
        let cfg2 = BootstrapConfig { seed: 2, ..cfg };
        let b = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg2).unwrap();
        // Different resamples; at least one branch should differ (exact
        // equality across all of them would mean the seed is ignored).
        assert_ne!(a.support, b.support);
    }

    /// Deterministic parity against `ape::prop.clades`.
    ///
    /// The reference values come from R on a fixed 6x5 fixture with three
    /// hand-written resamples (`/tmp/bootdet.R`):
    ///
    /// ```text
    /// orig clades : 123456 | 456 | 123 | 12 | 45
    /// rep 1 clades: 123456 | 456 | 123 | 12 | 45
    /// rep 2 clades: 123456 | 456 | 123 | 13 | 46
    /// rep 3 clades: 123456 | 456 | 123 | 12 | 56
    /// prop.clades : 3 3 3 2 1
    /// ```
    ///
    /// `123456` is the whole tree (the root, which has no support value), so
    /// the four interior branches are `456`, `123`, `12`, `45` with support
    /// `3, 3, 2, 1` out of 3 replicates.
    #[test]
    fn clade_counts_match_ape() {
        // Rows are samples, columns are loci, matching the R fixture.
        let dos = Mat::from_fn(6, 5, |i, j| {
            const DATA: [[f64; 5]; 6] = [
                [0.0, 0.0, 2.0, 2.0, 1.0],
                [0.0, 1.0, 2.0, 1.0, 1.0],
                [1.0, 0.0, 2.0, 2.0, 0.0],
                [2.0, 2.0, 0.0, 0.0, 1.0],
                [2.0, 1.0, 0.0, 0.0, 0.0],
                [2.0, 2.0, 0.0, 1.0, 0.0],
            ];
            DATA[i][j]
        });

        // The original tree, as `hclust(average)` builds it on Euclidean
        // distance. `tip_labels` order is s1..s6, matching R's rows.
        let dist = euclid(&dos);
        let mut tree = upgma(&dist).expect("tree");
        tree.tip_labels = levels(6);

        // The three resamples, as 0-based column indices (R used 1-based).
        let resamples: [&[usize]; 3] = [
            &[0, 0, 1, 2, 3],
            &[1, 1, 2, 3, 4],
            &[0, 2, 2, 4, 4],
        ];

        let target = branch_keys(&tree, true);
        let mut counts: HashMap<Vec<u64>, usize> = HashMap::new();
        for picks in resamples {
            let mut resampled = Mat::zeros(6, picks.len());
            for (dst, &src) in picks.iter().enumerate() {
                for r in 0..6 {
                    resampled.set(r, dst, dos.get(r, src));
                }
            }
            let t = upgma(&euclid(&resampled)).expect("replicate tree");
            for key in branch_keys(&t, true).into_keys() {
                if target.contains_key(&key) {
                    *counts.entry(key).or_insert(0) += 1;
                }
            }
        }

        // Match ape's order, not node-id order. `prop.part` walks the tree
        // depth-first and emits `[whole tree, 456, 123, 12, 45]`; `prop.clades`
        // drops the whole-tree entry, so the reported order is
        // `456, 123, 12, 45`.
        //
        // reed's node ids for this fixture are 7={1,2}, 8={4,5}, 9={4,5,6},
        // 10={1,2,3}, 11=all, so ape's order is ids [9, 10, 7, 8].
        let by_id: HashMap<u32, &Vec<u64>> = target.iter().map(|(k, &v)| (v, k)).collect();
        let got: Vec<usize> = [9u32, 10, 7, 8]
            .iter()
            .map(|id| counts.get(by_id[id]).copied().unwrap_or(0))
            .collect();

        assert_eq!(got, vec![3, 3, 2, 1], "clade counts differ from ape");
    }

    #[test]
    fn too_few_samples_are_rejected() {
        let x = structured(3, 10);
        let tree = tree_of(&x, TreeMethod::Upgma);
        let cfg = BootstrapConfig::default();
        assert!(matches!(
            support(&tree, &x, TreeMethod::Upgma, euclid, &cfg),
            Err(GraphError::TooFewTaxa { need: 4, .. })
        ));
    }

    #[test]
    fn no_sites_are_rejected() {
        let x = Mat::zeros(6, 0);
        let tree = Tree {
            edges: vec![(4, 1), (4, 2), (5, 4), (5, 3)],
            lengths: vec![1.0; 4],
            tip_labels: levels(3),
            node_labels: None,
            rooted: true,
        };
        let cfg = BootstrapConfig::default();
        assert!(support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).is_err());
    }

    #[test]
    fn support_length_matches_internal_node_count() {
        for n in [4usize, 5, 8, 12] {
            let x = structured(n, 24);
            let tree = tree_of(&x, TreeMethod::Upgma);
            let cfg = BootstrapConfig {
                replicates: 5,
                ..Default::default()
            };
            let report = support(&tree, &x, TreeMethod::Upgma, euclid, &cfg).unwrap();
            assert_eq!(
                report.support.len(),
                tree.n_internal(),
                "n={n}: {} vs {}",
                report.support.len(),
                tree.n_internal()
            );
        }
    }
}

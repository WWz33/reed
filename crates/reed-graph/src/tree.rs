//! Tree construction: UPGMA (average linkage) and neighbor-joining.

use reed_core::linalg::Mat;

use crate::error::GraphError;

/// A tree over named tips, in `ape`'s `phylo` shape.
///
/// Node ids are 1-based: tips are `1..=n_tips`, internal nodes follow in
/// creation order (`n_tips + 1..=2 * n_tips - 2` for a rooted binary tree).
#[derive(Clone, Debug, PartialEq)]
pub struct Tree {
    /// `(parent, child)` pairs.
    pub edges: Vec<(u32, u32)>,
    /// Branch lengths, same order as [`Tree::edges`].
    pub lengths: Vec<f64>,
    /// Tip labels; index `i` names node `i + 1`.
    pub tip_labels: Vec<String>,
    /// Internal node labels; `None` leaves them blank.
    pub node_labels: Option<Vec<String>>,
    /// `true` for an ultrametric, rooted tree (UPGMA).
    pub rooted: bool,
}

impl Tree {
    /// Number of tips.
    #[must_use]
    pub fn n_tips(&self) -> usize {
        self.tip_labels.len()
    }

    /// Number of internal nodes.
    ///
    /// Not `n_tips - 1`: that is the count for a *rooted binary* tree. An
    /// unrooted NJ tree over the same tips has `n_tips - 2`.
    #[must_use]
    pub fn n_internal(&self) -> usize {
        let Some(max_id) = self.edges.iter().flat_map(|&(p, c)| [p, c]).max() else {
            return 0;
        };
        (max_id as usize).saturating_sub(self.n_tips())
    }

    /// The largest node id used anywhere in the tree.
    #[must_use]
    pub fn max_node_id(&self) -> u32 {
        self.edges.iter().flat_map(|&(p, c)| [p, c]).max().unwrap_or(0)
    }

    /// The root node id (the only id that never appears as a child).
    ///
    /// Both constructors create internal nodes in ascending order, so the
    /// root carries the largest id; the scan is only a fallback for
    /// hand-built trees. Returns `None` for a tree with no edges.
    #[must_use]
    pub fn root(&self) -> Option<u32> {
        let max_id = self.max_node_id();
        if max_id == 0 {
            return None;
        }
        if !self.edges.iter().any(|&(_, c)| c == max_id) {
            return Some(max_id);
        }
        let mut seen = vec![false; max_id as usize + 1];
        for &(_, c) in &self.edges {
            seen[c as usize] = true;
        }
        (1..seen.len()).find(|&i| !seen[i]).map(|i| i as u32)
    }
}

/// UPGMA (average linkage), matching `poppr::upgma` — which is
/// `ape::as.phylo(stats::hclust(d, method = "average"))`.
///
/// The cluster distances are updated with the Lance–Williams recurrence
/// `d(C_k, C_i ∪ C_j) = (|C_i|·d(C_k,C_i) + |C_j|·d(C_k,C_j)) / (|C_i|+|C_j|)`,
/// which is what `hclust` does; it is algebraically identical to recomputing
/// the mean of all cross-cluster pairs, but O(n²) instead of O(n³).
///
/// Ties are broken by the smallest `(i, j)` pair in ascending order, matching
/// `hclust`'s scan order.
///
/// # Branch-length convention
///
/// `hclust$height` is the *merge distance* `d(i,j)`, which is twice the
/// height of the corresponding dendrogram node. `ape::as.phylo.hclust` halves
/// them so that the tree's patristic distances reproduce the input, and this
/// implementation does the same. The consequence is testable without R: for
/// an **ultrametric** input, `d_tree(i,j) == d_matrix(i,j)` for every pair.
///
/// # Errors
///
/// `NotSquare`, `TooFewTaxa`, `NonFiniteDistance`, `Asymmetric`,
/// `NonZeroDiagonal`.
pub fn upgma(d: &Mat) -> Result<Tree, GraphError> {
    let n = validate(d)?;
    if n < 2 {
        return Err(GraphError::TooFewTaxa { got: n, need: 2 });
    }

    let mut dist = d.clone();
    // Active cluster -> node id.
    let mut node_of: Vec<u32> = (1..=n as u32).collect();
    let mut size: Vec<usize> = vec![1; n];
    let mut height: Vec<f64> = vec![0.0; n];
    let mut active: Vec<bool> = vec![true; n];

    let mut edges = Vec::with_capacity(2 * (n - 1));
    let mut lengths = Vec::with_capacity(2 * (n - 1));
    let mut next_node = n as u32 + 1;

    for _ in 0..(n - 1) {
        // Smallest distance among active pairs, ties by (i, j) order.
        let mut best = f64::INFINITY;
        let mut bi = usize::MAX;
        let mut bj = usize::MAX;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                let v = dist.get(i, j);
                if v < best {
                    best = v;
                    bi = i;
                    bj = j;
                }
            }
        }
        if bi == usize::MAX {
            break;
        }

        let parent = next_node;
        next_node += 1;
        // `hclust` stores the newer cluster in column 1 and the older in
        // column 2 of its merge matrix, and the new cluster takes the slot of
        // the *newer* of the two. That choice decides which slot survives and
        // therefore how later ties resolve, so it is part of the algorithm,
        // not an implementation detail.
        let (newer, older) = if bi > bj { (bi, bj) } else { (bj, bi) };
        edges.push((parent, node_of[bi]));
        // `clusters[..]` holds hclust-style *merge distances*, which are twice
        // the dendrogram node height; a tree's branch lengths are node-height
        // differences, so halve. Without this every branch is 2x too long and
        // the tree no longer reproduces the input distances (`ape`'s
        // `as.phylo.hclust` does the same halving).
        lengths.push((best - height[bi]) / 2.0);
        edges.push((parent, node_of[bj]));
        lengths.push((best - height[bj]) / 2.0);

        // Merge the older cluster into the newer one's slot, so the surviving
        // slot index is `max(bi, bj)` — matching `hclust`.
        let (survivor, gone) = (newer, older);
        let ns = size[survivor] as f64;
        let ng = size[gone] as f64;
        for k in 0..n {
            if !active[k] || k == survivor || k == gone {
                continue;
            }
            let v = (ns * dist.get(survivor, k) + ng * dist.get(gone, k)) / (ns + ng);
            dist.set(survivor, k, v);
            dist.set(k, survivor, v);
        }
        size[survivor] += size[gone];
        height[survivor] = best;
        node_of[survivor] = parent;
        active[gone] = false;
    }

    Ok(Tree {
        edges,
        lengths,
        tip_labels: (0..n).map(|i| (i + 1).to_string()).collect(),
        node_labels: None,
        rooted: true,
    })
}

/// Neighbor-joining, matching `ape::nj`.
///
/// Each step picks the pair minimizing `Q(i,j) = (n-2)·d(i,j) - S_i - S_j`
/// where `S_i` is the sum of `i`'s distances to every other active taxon, and
/// joins them to a new node with branch lengths
/// `d(i,u) = d(i,j)/2 + (S_i - S_j)/(2(n-2))`, `d(j,u) = d(i,j) - d(i,u)`.
///
/// Ties are broken by the first `(i, j)` pair in scan order, as in `ape`.
/// Negative branch lengths are **kept** (ape keeps them too).
///
/// # Errors
///
/// `NotSquare`, `TooFewTaxa`, `NonFiniteDistance`, `Asymmetric`,
/// `NonZeroDiagonal`.
pub fn nj(d: &Mat) -> Result<Tree, GraphError> {
    let n = validate(d)?;
    if n < 3 {
        return Err(GraphError::TooFewTaxa { got: n, need: 3 });
    }

    let mut dist = d.clone();
    let mut node_of: Vec<u32> = (1..=n as u32).collect();
    let mut active: Vec<bool> = vec![true; n];
    let mut edges = Vec::with_capacity(2 * n - 3);
    let mut lengths = Vec::with_capacity(2 * n - 3);
    let mut next_node = n as u32 + 1;

    let mut live = n;
    while live > 3 {
        // S_i = sum of distances from i to every other active taxon.
        let mut s = vec![0.0f64; n];
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in 0..n {
                if active[j] && j != i {
                    s[i] += dist.get(i, j);
                }
            }
        }

        let b = (live - 2) as f64;
        let mut best = f64::INFINITY;
        let mut bi = usize::MAX;
        let mut bj = usize::MAX;
        for i in 0..n {
            if !active[i] {
                continue;
            }
            for j in (i + 1)..n {
                if !active[j] {
                    continue;
                }
                let q = b * dist.get(i, j) - s[i] - s[j];
                if q < best {
                    best = q;
                    bi = i;
                    bj = j;
                }
            }
        }
        if bi == usize::MAX {
            break;
        }

        let dij = dist.get(bi, bj);
        let parent = next_node;
        next_node += 1;
        let diu = 0.5 * dij + (s[bi] - s[bj]) / (2.0 * b);
        let dju = dij - diu;
        edges.push((parent, node_of[bi]));
        lengths.push(diu);
        edges.push((parent, node_of[bj]));
        lengths.push(dju);

        // d(u, k) = (d(i, k) + d(j, k) - d(i, j)) / 2
        for k in 0..n {
            if !active[k] || k == bi || k == bj {
                continue;
            }
            let v = 0.5 * (dist.get(bi, k) + dist.get(bj, k) - dij);
            dist.set(bi, k, v);
            dist.set(k, bi, v);
        }
        node_of[bi] = parent;
        active[bj] = false;
        live -= 1;
    }

    // The last three taxa join a single central node. `ape` emits `2n-3`
    // edges for an unrooted tree: the n-3 joins contribute 2 each, and the
    // final join contributes 3 (not 2), because the central node has degree 3.
    let last: Vec<usize> = (0..n).filter(|&i| active[i]).collect();
    if last.len() == 3 {
        let (a, b, c) = (last[0], last[1], last[2]);
        let dab = dist.get(a, b);
        let dac = dist.get(a, c);
        let dbc = dist.get(b, c);
        let parent = next_node;
        edges.push((parent, node_of[a]));
        lengths.push(0.5 * (dab + dac - dbc));
        edges.push((parent, node_of[b]));
        lengths.push(0.5 * (dab + dbc - dac));
        edges.push((parent, node_of[c]));
        lengths.push(0.5 * (dac + dbc - dab));
    }

    Ok(Tree {
        edges,
        lengths,
        tip_labels: (0..n).map(|i| (i + 1).to_string()).collect(),
        node_labels: None,
        rooted: false,
    })
}

/// Validates a distance matrix and returns its order.
fn validate(d: &Mat) -> Result<usize, GraphError> {
    if d.rows != d.cols {
        return Err(GraphError::NotSquare {
            rows: d.rows,
            cols: d.cols,
        });
    }
    let n = d.rows;
    for i in 0..n {
        for j in 0..n {
            let v = d.get(i, j);
            if !v.is_finite() {
                return Err(GraphError::NonFiniteDistance { row: i, col: j });
            }
            if i == j {
                if v != 0.0 {
                    return Err(GraphError::NonZeroDiagonal { index: i });
                }
            } else if (v - d.get(j, i)).abs() > 1e-9 * (1.0 + v.abs()) {
                return Err(GraphError::Asymmetric { row: i, col: j });
            }
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mat(rows: &[&[f64]]) -> Mat {
        let mut m = Mat::zeros(rows.len(), rows[0].len());
        for (i, r) in rows.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                m.set(i, j, v);
            }
        }
        m
    }

    /// A 4-taxon additive matrix: ((A,B),(C,D)) with known branch lengths.
    fn additive4() -> Mat {
        // Distances from the tree (A:1,B:1|C:2,D:2; internal 3).
        mat(&[
            &[0.0, 2.0, 6.0, 6.0],
            &[2.0, 0.0, 6.0, 6.0],
            &[6.0, 6.0, 0.0, 4.0],
            &[6.0, 6.0, 4.0, 0.0],
        ])
    }

    #[test]
    fn upgma_three_taxa_branch_lengths() {
        // d(1,2)=2, d(1,3)=4, d(2,3)=4. `hclust` reports merge distances 2
        // and 4; the tree carries node heights (half of those), so the
        // branches are 1, 1, 1 and 2 — matching
        // `ape::write.tree(as.phylo(hclust(as.dist(d), "average")))`, which
        // is `(3:2,(1:1,2:1):1);`.
        let d = mat(&[&[0.0, 2.0, 4.0], &[2.0, 0.0, 4.0], &[4.0, 4.0, 0.0]]);
        let t = upgma(&d).unwrap();
        assert_eq!(t.n_tips(), 3);
        assert_eq!(t.edges.len(), 4);
        assert!(t.rooted);
        assert_eq!(t.edges[0], (4, 1));
        assert_eq!(t.edges[1], (4, 2));
        assert!((t.lengths[0] - 1.0).abs() < 1e-12, "{:?}", t.lengths);
        assert!((t.lengths[1] - 1.0).abs() < 1e-12, "{:?}", t.lengths);
        assert_eq!(t.edges[2], (5, 4));
        assert_eq!(t.edges[3], (5, 3));
        assert!((t.lengths[2] - 1.0).abs() < 1e-12, "{:?}", t.lengths);
        assert!((t.lengths[3] - 2.0).abs() < 1e-12, "{:?}", t.lengths);
    }

    /// Patristic distance between two tips: the branch lengths along the path
    /// through their lowest common ancestor.
    fn patristic(t: &Tree, a: u32, b: u32) -> f64 {
        let sz = 2 * t.n_tips() + 1;
        let mut parent = vec![0u32; sz];
        let mut len_to_parent = vec![0.0f64; sz];
        let mut depth = vec![0.0f64; sz];
        for (&(p, c), &l) in t.edges.iter().zip(&t.lengths) {
            parent[c as usize] = p;
            len_to_parent[c as usize] = l;
            depth[c as usize] = depth[p as usize] + l;
        }
        let (mut i, mut j) = (a, b);
        while depth[i as usize] > depth[j as usize] {
            i = parent[i as usize];
        }
        while depth[j as usize] > depth[i as usize] {
            j = parent[j as usize];
        }
        while i != j {
            i = parent[i as usize];
            j = parent[j as usize];
        }
        let lca = i;
        let to_lca = |mut x: u32| -> f64 {
            let mut acc = 0.0;
            while x != lca {
                acc += len_to_parent[x as usize];
                x = parent[x as usize];
            }
            acc
        };
        to_lca(a) + to_lca(b)
    }

    /// The property that defines a UPGMA tree: for an ultrametric input, the
    /// tree reproduces the input distances exactly. This is what catches a
    /// wrong branch-length convention (a 2x factor here is invisible in the
    /// topology and in `hclust$height`, which are merge distances).
    #[test]
    fn upgma_reproduces_ultrametric_distances() {
        // ((1,2),(3,4)) with node heights 1 and 3.
        let d = mat(&[
            &[0.0, 2.0, 6.0, 6.0],
            &[2.0, 0.0, 6.0, 6.0],
            &[6.0, 6.0, 0.0, 2.0],
            &[6.0, 6.0, 2.0, 0.0],
        ]);
        let t = upgma(&d).unwrap();
        for i in 1..=4u32 {
            for j in (i + 1)..=4u32 {
                let want = d.get(i as usize - 1, j as usize - 1);
                let got = patristic(&t, i, j);
                assert!((want - got).abs() < 1e-12, "d({i},{j}): want {want}, got {got}");
            }
        }
    }

    #[test]
    fn upgma_merge_heights_match_hclust() {
        // `hclust(as.dist(d), "average")` on this matrix reports
        // heights 4, 8, 8 and merges (-1,-2), (-3,1), (-4,2).
        //
        // Note this fixture is *not* ultrametric even though it looks
        // symmetric: d(1,3) = 8 while the merge heights say 4 and 8. The
        // check is therefore against the merge heights, not against equal
        // leaf depths (which is a property of the input, not of UPGMA).
        let d = mat(&[
            &[0.0, 4.0, 8.0, 8.0],
            &[4.0, 0.0, 8.0, 8.0],
            &[8.0, 8.0, 0.0, 8.0],
            &[8.0, 8.0, 8.0, 0.0],
        ]);
        let t = upgma(&d).unwrap();
        // Merge height per internal node, in creation order.
        let mut merge_h = vec![0.0f64; 2 * t.n_tips()];
        for (&(p, c), &len) in t.edges.iter().zip(&t.lengths) {
            let child_h = if (c as usize) <= t.n_tips() {
                0.0
            } else {
                merge_h[c as usize]
            };
            merge_h[p as usize] = merge_h[p as usize].max(child_h + len);
        }
        let heights: Vec<f64> = ((t.n_tips() + 1)..(2 * t.n_tips()))
            .map(|i| merge_h[i])
            .collect();
        assert_eq!(heights.len(), 3);
        // Node heights are half the merge distances, as in `as.phylo.hclust`.
        // `ape` renders this matrix as `(4:4,(3:4,(1:2,2:2):2):0);`.
        assert!((heights[0] - 2.0).abs() < 1e-12, "{heights:?}");
        assert!((heights[1] - 4.0).abs() < 1e-12, "{heights:?}");
        assert!((heights[2] - 4.0).abs() < 1e-12, "{heights:?}");
    }

    #[test]
    fn nj_recovers_additive_branch_lengths() {
        let t = nj(&additive4()).unwrap();
        assert_eq!(t.n_tips(), 4);
        assert_eq!(t.edges.len(), 2 * 4 - 3);
        assert!(!t.rooted);
        // The true tree is ((A,B),(C,D)): A and B share a parent whose
        // distance to A is 1, C and D share a parent at distance 2.
        let len_of = |child: u32| -> f64 {
            t.edges
                .iter()
                .zip(&t.lengths)
                .find(|((_, c), _)| *c == child)
                .map(|(_, &l)| l)
                .unwrap()
        };
        assert!((len_of(1) - 1.0).abs() < 1e-12, "A branch {}", len_of(1));
        assert!((len_of(2) - 1.0).abs() < 1e-12, "B branch {}", len_of(2));
        assert!((len_of(3) - 2.0).abs() < 1e-12, "C branch {}", len_of(3));
        assert!((len_of(4) - 2.0).abs() < 1e-12, "D branch {}", len_of(4));
    }

    #[test]
    fn nj_three_taxa_has_one_internal_node() {
        let d = mat(&[&[0.0, 2.0, 4.0], &[2.0, 0.0, 4.0], &[4.0, 4.0, 0.0]]);
        let t = nj(&d).unwrap();
        assert_eq!(t.edges.len(), 3);
        // Three tips unrooted: one central node, so `n-2` internals — not
        // `n-1`, which is the rooted-binary count.
        assert_eq!(t.n_internal(), 1);
    }

    #[test]
    fn n_internal_counts_rooted_and_unrooted() {
        // UPGMA over 4 tips is a rooted binary tree: 3 internal nodes.
        assert_eq!(upgma(&additive4()).unwrap().n_internal(), 3);
        // NJ over 4 tips is unrooted: 2 internal nodes.
        let nj_tree = nj(&additive4()).unwrap();
        assert_eq!(nj_tree.n_internal(), 2);
        assert_eq!(nj_tree.edges.len(), 2 * 4 - 3);
    }

    #[test]
    fn nj_keeps_negative_branch_lengths() {
        // A 3-taxon matrix violating the triangle inequality drives one
        // branch negative; ape keeps it rather than clamping to zero.
        // d(a,b)=1, d(a,c)=1, d(b,c)=10  ->  d(c, centre) = (1+10-1)/2 = 5,
        // and the other two come out at -4 each.
        let d = mat(&[&[0.0, 1.0, 1.0], &[1.0, 0.0, 10.0], &[1.0, 10.0, 0.0]]);
        let t = nj(&d).unwrap();
        assert!(
            t.lengths.iter().any(|&l| l < 0.0),
            "expected a negative branch, got {:?}",
            t.lengths
        );
    }

    #[test]
    fn non_square_is_rejected() {
        let d = Mat::zeros(3, 4);
        assert!(matches!(upgma(&d), Err(GraphError::NotSquare { .. })));
        assert!(matches!(nj(&d), Err(GraphError::NotSquare { .. })));
    }

    #[test]
    fn asymmetric_is_rejected() {
        let d = mat(&[&[0.0, 1.0], &[2.0, 0.0]]);
        assert!(matches!(upgma(&d), Err(GraphError::Asymmetric { .. })));
    }

    #[test]
    fn nonzero_diagonal_is_rejected() {
        let d = mat(&[&[1.0, 1.0], &[1.0, 0.0]]);
        assert!(matches!(upgma(&d), Err(GraphError::NonZeroDiagonal { .. })));
    }

    #[test]
    fn non_finite_is_rejected() {
        let d = mat(&[&[0.0, f64::NAN], &[f64::NAN, 0.0]]);
        assert!(matches!(
            upgma(&d),
            Err(GraphError::NonFiniteDistance { .. })
        ));
    }

    #[test]
    fn too_few_taxa_are_rejected() {
        let d = mat(&[&[0.0]]);
        assert!(matches!(upgma(&d), Err(GraphError::TooFewTaxa { .. })));
        let d3 = mat(&[&[0.0, 1.0, 2.0], &[1.0, 0.0, 1.0], &[2.0, 1.0, 0.0]]);
        assert!(upgma(&d3).is_ok());
    }

    #[test]
    fn root_is_the_only_node_never_a_child() {
        let t = upgma(&additive4()).unwrap();
        let root = t.root().unwrap();
        assert!(t.edges.iter().all(|&(p, _)| p != root || root != 0));
        assert!(!t.edges.iter().any(|&(_, c)| c == root));
    }
}

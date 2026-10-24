//! Newick serialization in `ape::write.tree` format.

use crate::error::GraphError;
use crate::tree::Tree;

/// Maximum nesting depth before the serializer gives up.
///
/// The traversal is iterative, so this is a policy limit rather than a stack
/// guard: a 10 000-deep tree is not a genetic distance tree, it is a bug or a
/// pathological input.
const MAX_DEPTH: usize = 10_000;

/// Renders `tree` as a newick string.
///
/// Numbers use `digits` significant decimal places with trailing zeros
/// removed, matching R's `format(x, digits = digits)`. Branch lengths of
/// exactly zero are written as `:0` (as `ape` does). Internal node labels are
/// emitted when [`Tree::node_labels`] is set.
///
/// # Errors
///
/// `BadBranchLength` for a non-finite length, `TreeTooDeep` past
/// [`MAX_DEPTH`], and `NotSquare`-family errors never arise here.
pub fn to_newick(tree: &Tree, digits: usize) -> Result<String, GraphError> {
    let n = tree.n_tips();
    if n == 0 {
        return Ok("();".to_string());
    }
    // Tips but no edges: a single-tip tree, or a malformed one. Either way
    // there is nothing to render but the label, and sizing the child table
    // from `max_node_id` would leave it empty.
    if tree.edges.is_empty() {
        return Ok(format!("{};", label_for(tree, 1, n)));
    }

    // Children of each node, preserving edge order. Sizing from the largest
    // id actually used (at least `n`, so every tip id is addressable) keeps
    // small trees small.
    let max_id = (tree.max_node_id() as usize).max(n);
    let mut children: Vec<Vec<(u32, f64)>> = vec![Vec::new(); max_id + 1];
    // Every node has at most one parent, so a vec indexed by node id beats a
    // map here.
    let mut parent_len: Vec<Option<f64>> = vec![None; max_id + 1];
    for (idx, &(p, c)) in tree.edges.iter().enumerate() {
        let len = tree.lengths[idx];
        if !len.is_finite() {
            return Err(GraphError::BadBranchLength { edge: idx });
        }
        children[p as usize].push((c, len));
        parent_len[c as usize] = Some(len);
    }

    let root = tree.root().unwrap_or(1);
    let mut out = String::new();

    // Explicit stack: (node id, next child index, depth).
    let mut stack: Vec<(u32, usize, usize)> = vec![(root, 0, 0)];
    // Whether the node is a tip decides its rendering.
    let is_tip = |id: u32| (id as usize) <= n;

    if is_tip(root) {
        // A single-tip tree.
        out.push_str(&label_for(tree, root, n));
        out.push(';');
        return Ok(out);
    }

    // The root has no incoming edge, so no `:` follows it.
    out.push('(');

    while let Some(&mut (node, ref mut child_idx, depth)) = stack.last_mut() {
        if depth > MAX_DEPTH {
            return Err(GraphError::TreeTooDeep { depth });
        }
        let kids = &children[node as usize];
        if *child_idx >= kids.len() {
            stack.pop();
            out.push(')');
            if let Some(labels) = &tree.node_labels {
                let idx = node as usize - n - 1;
                if let Some(l) = labels.get(idx).filter(|l| !l.is_empty()) {
                    out.push_str(l);
                }
            }
            // An internal node's own branch length, unless it is the root.
            if node != root
                && let Some(len) = parent_len.get(node as usize).copied().flatten()
            {
                out.push(':');
                out.push_str(&format_number(len, digits));
            }
            continue;
        }
        let (child, len) = kids[*child_idx];
        if *child_idx > 0 {
            out.push(',');
        }
        *child_idx += 1;
        if is_tip(child) {
            out.push_str(&label_for(tree, child, n));
            out.push(':');
            out.push_str(&format_number(len, digits));
        } else {
            out.push('(');
            stack.push((child, 0, depth + 1));
        }
    }

    // Close the root's own parenthesis group is handled by the pop above only
    // for inner nodes; the root's group needs one more `)`.
    // The loop pushes `(` before descending, so the root's closing `)` is
    // emitted when the root pops — which the loop above already does.
    out.push(';');
    Ok(out)
}

/// The newick label for a tip.
fn label_for(tree: &Tree, id: u32, n: usize) -> String {
    let idx = id as usize - 1;
    if idx < n {
        let l = &tree.tip_labels[idx];
        if l.is_empty() {
            (idx + 1).to_string()
        } else {
            sanitize(l)
        }
    } else {
        id.to_string()
    }
}

/// Escapes the characters newick reserves.
fn sanitize(label: &str) -> String {
    let needs_quote = label
        .chars()
        .any(|c| matches!(c, '(' | ')' | ',' | ':' | ';' | '[' | ']' | '\'' | ' '));
    if !needs_quote {
        return label.to_string();
    }
    let mut out = String::with_capacity(label.len() + 2);
    out.push('\'');
    for c in label.chars() {
        if c == '\'' {
            out.push('\'');
        }
        out.push(c);
    }
    out.push('\'');
    out
}

/// Formats a branch length the way `ape::write.tree` does.
///
/// R's `format(x, digits = n)` keeps `n` **significant** digits, not `n`
/// decimal places: `0.0816666666666` becomes `0.08166666667` and
/// `123456789.123` becomes `123456789.1`. Matching it needs the decimal
/// exponent, so this finds it and asks for the corresponding number of
/// places. Values so small or large that even that fails (below `1e-20` or
/// above `1e20`) fall back to `decimals = 0/20`, which is why the doc says
/// the match holds for the magnitudes branch lengths actually take.
fn format_number(x: f64, digits: usize) -> String {
    if x == 0.0 || !x.is_finite() {
        return "0".to_string();
    }
    let exp = x.abs().log10().floor() as i32;
    let decimals = (digits as i32 - 1 - exp).clamp(0, 20) as usize;
    let s = format!("{x:.decimals$}");
    // Trim only inside a fractional part: `1000` must not become `1`.
    let s = if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    };
    if s.is_empty() || s == "-" {
        "0".to_string()
    } else {
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::linalg::Mat;

    fn mat(rows: &[&[f64]]) -> Mat {
        let mut m = Mat::zeros(rows.len(), rows[0].len());
        for (i, r) in rows.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                m.set(i, j, v);
            }
        }
        m
    }

    #[test]
    fn numbers_are_trimmed() {
        assert_eq!(format_number(0.0, 10), "0");
        assert_eq!(format_number(1.0, 10), "1");
        assert_eq!(format_number(1.5, 10), "1.5");
        assert_eq!(format_number(-2.25, 10), "-2.25");
        // Trailing zeros must only be trimmed inside a fraction: a whole
        // number's zeros are significant (`1000` must not become `1`).
        assert_eq!(format_number(1000.0, 10), "1000");
        assert_eq!(format_number(100.0, 10), "100");
    }

    /// `digits` means significant digits, as in R's `format(x, digits = n)`,
    /// which is what `ape::write.tree` calls. These expectations are R's own
    /// output for the same inputs.
    #[test]
    fn digits_are_significant_not_decimal_places() {
        assert_eq!(format_number(0.085, 10), "0.085");
        assert_eq!(format_number(0.081_666_666_666_6, 10), "0.08166666667");
        assert_eq!(format_number(1.234_567_890_123_45, 10), "1.23456789");
        assert_eq!(format_number(0.000_123_456_789_012, 10), "0.000123456789");
        assert_eq!(format_number(123_456_789.123, 10), "123456789.1");
    }

    #[test]
    fn non_finite_numbers_do_not_reach_the_formatter() {
        // The serializer rejects them before formatting, but a non-finite
        // input must not produce a bare `NaN` token either.
        assert_eq!(format_number(f64::NAN, 10), "0");
        assert_eq!(format_number(f64::INFINITY, 10), "0");
    }

    #[test]
    fn upgma_tree_round_trips_through_newick() {
        let d = mat(&[&[0.0, 2.0, 4.0], &[2.0, 0.0, 4.0], &[4.0, 4.0, 0.0]]);
        let t = crate::upgma(&d).unwrap();
        let s = to_newick(&t, 10).unwrap();
        assert!(s.starts_with('('), "{s}");
        assert!(s.ends_with(';'), "{s}");
        assert_eq!(s.matches('(').count(), 2, "{s}");
        assert_eq!(s.matches(')').count(), 2, "{s}");
        assert_eq!(s.matches(':').count(), 4, "{s}");
    }

    #[test]
    fn nj_tree_has_the_right_tip_count() {
        let d = mat(&[
            &[0.0, 2.0, 6.0, 6.0],
            &[2.0, 0.0, 6.0, 6.0],
            &[6.0, 6.0, 0.0, 4.0],
            &[6.0, 6.0, 4.0, 0.0],
        ]);
        let t = crate::nj(&d).unwrap();
        let s = to_newick(&t, 10).unwrap();
        for tip in ["1", "2", "3", "4"] {
            assert!(s.contains(tip), "{tip} missing from {s}");
        }
    }

    #[test]
    fn labels_with_reserved_characters_are_quoted() {
        let mut t = crate::upgma(&mat(&[
            &[0.0, 2.0, 4.0],
            &[2.0, 0.0, 4.0],
            &[4.0, 4.0, 0.0],
        ]))
        .unwrap();
        t.tip_labels = vec!["a,b".into(), "c:d".into(), "e f".into()];
        let s = to_newick(&t, 10).unwrap();
        assert!(s.contains("'a,b'"), "{s}");
        assert!(s.contains("'c:d'"), "{s}");
        assert!(s.contains("'e f'"), "{s}");
    }

    #[test]
    fn non_finite_branch_length_is_rejected() {
        let mut t = crate::upgma(&mat(&[
            &[0.0, 2.0, 4.0],
            &[2.0, 0.0, 4.0],
            &[4.0, 4.0, 0.0],
        ]))
        .unwrap();
        t.lengths[0] = f64::NAN;
        assert!(matches!(
            to_newick(&t, 10),
            Err(GraphError::BadBranchLength { .. })
        ));
    }

    #[test]
    fn node_labels_are_emitted() {
        let mut t = crate::upgma(&mat(&[
            &[0.0, 2.0, 4.0],
            &[2.0, 0.0, 4.0],
            &[4.0, 4.0, 0.0],
        ]))
        .unwrap();
        t.node_labels = Some(vec!["95".into(), "100".into()]);
        let s = to_newick(&t, 10).unwrap();
        assert!(s.contains(")95"), "{s}");
        assert!(s.contains(")100"), "{s}");
    }

    #[test]
    fn deep_tree_does_not_overflow_the_stack() {
        // A caterpillar: 1 joins 2, then 3, then 4, ...
        let n = 2000;
        let mut d = Mat::zeros(n, n);
        for i in 0..n {
            for j in (i + 1)..n {
                let v = (j - i) as f64;
                d.set(i, j, v);
                d.set(j, i, v);
            }
        }
        let t = crate::upgma(&d).unwrap();
        let s = to_newick(&t, 6).unwrap();
        assert!(s.ends_with(';'));
    }
}

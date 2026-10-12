//! The unified pass0 entry (doc_v6 03 §7): which global quantities the
//! run needs, the mode that collects them, and the analyses that read the
//! matrices back. `--mind`'s keep set, the `--indep-pairwise` site rule
//! and the sketch path's mid stage all read the same pass0 output — the
//! three pre-Y1e normalizations (`MissingJob`, `SketchJob`, `LdPruneJob`)
//! are gone.

pub(super) mod matrix;
pub(super) mod prune;

pub(super) use matrix::{MatrixOutput, SiteMatrix, run_matrix_pass};
pub(super) use prune::{PruneConfig, SiteKeepSet};

use matrix::{ChunkSites, MatrixMode};
use std::path::Path;
use std::sync::Arc;

/// The pass0 mode a spec needs: the matrices are retained when the pruner
/// or the mid stage reads them; the rule bitmaps ride along only for the
/// mid stage; the tail runs inline only when the keep set is final.
pub(super) fn plan_mode(mind: bool, sketch_shape: bool, indep: bool) -> MatrixMode {
    MatrixMode {
        retain: sketch_shape || indep,
        rule_bits: sketch_shape,
        tail_inline: !mind && !sketch_shape,
        // The pruner excludes masked and non-kept samples (doc_v7 D7-4), so
        // the effective mask bitmap rides along whenever it runs.
        mask_plane: indep,
    }
}

/// Packs a per-sample keep vector into the bitmap words the pruner excludes
/// over (doc_v7 D7-4).
pub(super) fn keep_words_of(keep: &[bool]) -> Vec<u64> {
    let mut words = vec![0_u64; keep.len().div_ceil(64).max(1)];
    for (i, &k) in keep.iter().enumerate() {
        if k {
            words[i / 64] |= 1 << (i % 64);
        }
    }
    words
}

/// Resolves the `--mind` keep set from the pass0 counts (doc_v4 02 回退版
/// §3.2): keep iff `(2L − called) · 10⁹ ≤ q · 2L`, with `L` the
/// head-surviving record count and `called` the post-mask called alleles.
/// With no head-surviving record there is no evidence, so every sample
/// stays (`tolerates(0, 0)` is vacuously true). Returns the per-sample
/// keep vector and the number of samples the rule dropped.
pub(super) fn mind_keep(
    rate: &reed_filter::spec::MissingRate,
    called: &[u64],
    n_head: u64,
    keep: &mut [bool],
) -> u64 {
    let total = 2 * n_head;
    let mut n_dropped = 0_u64;
    for (slot, &called) in keep.iter_mut().zip(called) {
        if *slot && !rate.tolerates(total.saturating_sub(called), total) {
            *slot = false;
            n_dropped += 1;
        }
    }
    n_dropped
}

/// The prune configuration from the spec.
pub(super) fn prune_config(prune: &reed_filter::spec::LdPruneSpec) -> PruneConfig {
    PruneConfig {
        threshold: prune.r2,
        window: if prune.kb {
            reed_stats::ld::LdWindow::Bp(prune.window)
        } else {
            reed_stats::ld::LdWindow::Variants(usize::try_from(prune.window).unwrap_or(usize::MAX))
        },
        step: prune.step,
    }
}

/// Builds the global survivor bitmap the post-`--mind` prune filters on
/// from the mid stage's per-chunk decisions.
pub(super) fn survivor_words(
    chunks: &[ChunkSites],
    total_sites: usize,
    decisions: &std::collections::HashMap<u64, Arc<crate::filter::mid_stage::ChunkDecision>>,
) -> Vec<u64> {
    let mut words = vec![0_u64; total_sites.div_ceil(64)];
    for chunk in chunks {
        let Some(decision) = decisions.get(&chunk.id) else {
            continue;
        };
        for s in 0..chunk.n_sites {
            if decision.survivor(s) {
                let global = chunk.start + s;
                words[global / 64] |= 1 << (global % 64);
            }
        }
    }
    words
}

/// The pass0 record-count estimate the `--mind` strategy choice gates on
/// (§11.4: the first chunk's records/byte density, extrapolated to the
/// file size — spec-sanctioned when no index is available).
pub(super) fn estimate_records(input: &Path) -> Result<u64, crate::filter::FilterError> {
    let file_len = std::fs::metadata(input)
        .map_err(crate::filter::FilterError::from)?
        .len();
    let mut opened = crate::input::open(Some(input)).map_err(crate::filter::io_err)?;
    let Some(chunk) = opened.source.next_chunk().map_err(crate::filter::io_err)? else {
        return Ok(0);
    };
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut scratch = Vec::new();
    let view = reed_io::chunk::inflate_chunk(&chunk, &mut inflater, &mut scratch)
        .map_err(crate::filter::io_err)?;
    let body = view.body_bytes(view.text);
    let n_lines = u64::try_from(memchr::memchr_iter(b'\n', body).count()).unwrap_or(0);
    let compressed = u64::try_from(chunk.compressed.len()).unwrap_or(1).max(1);
    if n_lines == 0 {
        return Ok(0);
    }
    Ok(file_len.saturating_mul(n_lines) / compressed)
}

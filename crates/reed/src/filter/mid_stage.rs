//! The sketch path's mid stage (`doc_v2/03-architecture.md` §8, §11.4):
//! per-chunk sketches → per-site counts (keep ∧ planes), tail rules, mask
//! attribution and the survivor's emit data.

use std::sync::Arc;

use reed_core::count::SiteCounts;
use reed_filter::eval::edit_plan::fill_info_set;
use reed_filter::eval::tail::eval_tail;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::vcf::line_edit::InfoValueOut;

use super::sketch::ChunkSketch;

/// One survivor's pass1 emit data (§11.4).
#[derive(Clone, Debug)]
pub(super) struct SiteEmit {
    /// Kept ∧ mask-fired sample bitmap (the EditPlan's `mask_gt`).
    pub mask_words: Vec<u64>,
    /// Recomputed INFO values in fill-tags order.
    pub infos: Vec<(&'static [u8], InfoValueOut)>,
}

/// A chunk's pass1 decisions, produced by the mid stage.
#[derive(Clone, Debug, Default)]
pub(super) struct ChunkDecision {
    /// Records in the chunk (body + trailing seam).
    pub n_records: u64,
    /// Bit per record: passed the head rules.
    pub head_bits: Vec<u64>,
    /// Bit per site (head-passing ordinal): survived the tail rules.
    pub survivor_bits: Vec<u64>,
    /// Per-site emit data (`None` = tail-failed).
    pub survivors: Vec<Option<SiteEmit>>,
    /// The chunk's trailing seam record's emit (the seam site is appended
    /// to this chunk's sketch by the sink).
    pub seam: Option<SiteEmit>,
}

impl ChunkDecision {
    /// Whether record `r` passed the head rules.
    pub(super) fn head_passed(&self, r: usize) -> bool {
        self.head_bits.get(r / 64).is_some_and(|word| word >> (r % 64) & 1 == 1)
    }

    /// Whether site `s` survived the tail rules.
    pub(super) fn survivor(&self, s: usize) -> bool {
        self.survivor_bits.get(s / 64).is_some_and(|word| word >> (s % 64) & 1 == 1)
    }
}

pub(super) fn mid_stage(
    sketches: &[ChunkSketch],
    keep: &[bool],
    plan: &FilterPlan,
    n_samples_out: usize,
    funnel: &mut FunnelPartial,
) -> std::collections::HashMap<u64, Arc<ChunkDecision>> {
    // Clean-padded keep words (the count math anchors on this).
    let w = keep.len().div_ceil(64);
    let mut keep_words = vec![0_u64; w];
    for (i, &k) in keep.iter().enumerate() {
        if k {
            keep_words[i / 64] |= 1 << (i % 64);
        }
    }
    let mut map = std::collections::HashMap::new();
    for sketch in sketches {
        let frag = &sketch.fragment;
        let mut survivor_bits = vec![0_u64; frag.n_sites.div_ceil(64)];
        let mut survivors: Vec<Option<SiteEmit>> = Vec::with_capacity(frag.n_sites);
        for s in 0..frag.n_sites {
            let c = reed_filter::sketch::count_site_sketch(frag, s, &keep_words);
            funnel.masked_gt += c.masked_gt;
            for (r, hits) in
                reed_filter::sketch::rule_hits_of(frag, s, &keep_words).into_iter().enumerate()
            {
                if let Some(slot) = funnel.mask_hits.get_mut(r) {
                    *slot += hits;
                }
            }
            let counts = SiteCounts {
                ac: c.ac,
                an: c.an,
                ns: c.ns,
                n_missing: c.n_missing,
                hom_ref: c.hom_ref,
                het: c.het,
                hom_alt: c.hom_alt,
            };
            let tail_values = reed_filter::expr::SiteValues::from_counts(reed_filter::expr::GtCounts {
                ac: counts.ac,
                an: counts.an,
                n_missing: counts.n_missing,
                n_samples: n_samples_out as u64,
            });
            if eval_tail(&tail_values, &counts, n_samples_out, &plan.tail, &mut funnel.tail_hits) {
                survivor_bits[s / 64] |= 1 << (s % 64);
                let mut infos = Vec::new();
                fill_info_set(&counts, &plan.recompute, &mut infos);
                let mask_words: Vec<u64> = frag
                    .mask_of(s)
                    .iter()
                    .zip(&keep_words)
                    .map(|(m, k)| m & k)
                    .collect();
                survivors.push(Some(SiteEmit { mask_words, infos }));
                funnel.n_output += 1;
            } else {
                survivors.push(None);
            }
        }
        let seam = sketch
            .seam_site
            .and_then(|s| survivors.get(s).cloned())
            .flatten();
        map.insert(
            sketch.id,
            Arc::new(ChunkDecision {
                n_records: frag.n_records,
                head_bits: frag.head_bits.clone(),
                survivor_bits,
                survivors,
                seam,
            }),
        );
    }
    map
}

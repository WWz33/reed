//! The mid stage (`docs/rounds/v2/03-architecture.md` §8, §11.4): the
//! unified pass0's per-contig matrices → per-site counts over the final
//! keep set, tail rules, mask attribution and the survivor's emit data.
//! Before Y1e this read chunk sketches; the site bytes are the same
//! `(4 + R) · w` words, now read back from the matrices via the global
//! site ordinals the pass0 recorded per chunk.

use std::sync::Arc;

use reed_core::count::SiteCounts;
use reed_filter::eval::edit_plan::fill_info_set;
use reed_filter::eval::tail::eval_tail;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::vcf::line_edit::InfoValueOut;

use super::pass0::matrix::{ChunkSites, SiteMatrix};
use reed_engine::budget::ThreadBudget;

/// One survivor's pass1 emit data (§11.4).
#[derive(Clone, Debug)]
pub(super) struct SiteEmit {
    /// Kept ∧ mask-fired sample bitmap (the EditPlan's `mask_gt`).
    pub mask_words: Vec<u64>,
    /// Recomputed INFO values in fill-tags order.
    pub infos: Vec<(&'static [u8], InfoValueOut)>,
    /// Whether the record's FORMAT carried a GT key (doc_v12 D12-2: a
    /// GT-less survivor keeps the input's AC/AN instead of the zeros the
    /// recompute would produce).
    pub has_gt: bool,
}

/// A chunk's pass1 decisions, produced by the mid stage.
#[derive(Clone, Debug, Default)]
pub(super) struct ChunkDecision {
    /// Bit per record: passed the head rules.
    pub head_bits: Vec<u64>,
    /// Bit per site (head-passing ordinal): survived the tail rules.
    pub survivor_bits: Vec<u64>,
    /// Per-site emit data (`None` = tail-failed).
    pub survivors: Vec<Option<SiteEmit>>,
    /// The chunk's trailing seam record's emit (the seam site is appended
    /// to this chunk's stream by the sink).
    pub seam: Option<SiteEmit>,
}

impl ChunkDecision {
    /// Whether record `r` passed the head rules.
    pub(super) fn head_passed(&self, r: usize) -> bool {
        self.head_bits
            .get(r / 64)
            .is_some_and(|word| word >> (r % 64) & 1 == 1)
    }

    /// Whether site `s` survived the tail rules.
    pub(super) fn survivor(&self, s: usize) -> bool {
        self.survivor_bits
            .get(s / 64)
            .is_some_and(|word| word >> (s % 64) & 1 == 1)
    }
}

/// Post-mask per-site counts derived from the raw planes plus the rule
/// bitmaps (the pre-Y1e sketch recount, byte for byte). Identical by
/// construction to a SinglePass recount over `keep ∧ ¬mask` samples.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct PlaneCounts {
    /// ALT allele count.
    pub ac: u64,
    /// Called allele number.
    pub an: u64,
    /// Samples with data (half-calls included).
    pub ns: u64,
    /// Kept samples without data (masked or raw-missing).
    pub n_missing: u64,
    /// Full calls, both alleles REF.
    pub hom_ref: u64,
    /// Full calls, exactly one ALT.
    pub het: u64,
    /// Full calls, both ALT.
    pub hom_alt: u64,
    /// Half-calls among the kept, unmasked samples (exactly one called
    /// slot).
    pub n_half_called: u64,
    /// Called-to-missing transitions this site's mask performed.
    pub masked_gt: u64,
}

/// Counts one site over the kept samples: `planes` is the site's four raw
/// genotype planes, `rules` its per-rule firing bitmaps (their OR is the
/// effective mask), `keep` one word per sample bitmap word.
pub(super) fn count_site_planes(
    planes: [&[u64]; 4],
    rules: &[&[u64]],
    keep: &[u64],
) -> PlaneCounts {
    let mut c = PlaneCounts::default();
    for (i, &keep_word) in keep.iter().enumerate() {
        let mut mask_word = 0_u64;
        for rule in rules {
            mask_word |= rule.get(i).copied().unwrap_or(0);
        }
        let p0a = planes[0].get(i).copied().unwrap_or(0);
        let p0m = planes[1].get(i).copied().unwrap_or(0);
        let p1a = planes[2].get(i).copied().unwrap_or(0);
        let p1m = planes[3].get(i).copied().unwrap_or(0);

        let kept_raw = keep_word & !mask_word;

        // ALT alleles over called slots of kept, unmasked samples.
        c.ac += (kept_raw & p0a).count_ones() as u64;
        c.ac += (kept_raw & p1a).count_ones() as u64;
        // Called slots: everything but the missing bit.
        c.an += 2 * (kept_raw.count_ones() as u64)
            - (kept_raw & p0m).count_ones() as u64
            - (kept_raw & p1m).count_ones() as u64;
        // Samples with data: at least one called slot.
        c.ns += (kept_raw & !(p0m & p1m)).count_ones() as u64;

        // Full calls: both slots called.
        let full = kept_raw & !p0m & !p1m;
        c.het += (full & (p0a ^ p1a)).count_ones() as u64;
        c.hom_alt += (full & p0a & p1a).count_ones() as u64;
        c.hom_ref += (full & !p0a & !p1a).count_ones() as u64;
        // Half-calls: kept and unmasked, exactly one slot missing.
        c.n_half_called += (kept_raw & (p0m ^ p1m)).count_ones() as u64;

        // Masked samples that were called before the mask: the called→
        // missing transitions.
        c.masked_gt += (keep_word & mask_word & !(p0m & p1m)).count_ones() as u64;
    }
    let kept_total: u32 = keep.iter().map(|word| word.count_ones()).sum();
    c.n_missing = u64::from(kept_total) - c.ns;
    c
}

/// Per-rule mask attribution for one site: the raw kept firing count of
/// rule `r` (the caller folds them in funnel order; a sample may fire
/// several rules but SinglePass attributes the FIRST one — mirrored by
/// folding `keep ∧ rule_r ∧ ¬(rule_0 ∨ … ∨ rule_{r-1})`).
pub(super) fn rule_hits_of(planes_rules: &[&[u64]], keep: &[u64]) -> Vec<u64> {
    let w = keep.len();
    let mut or_so_far = vec![0_u64; w];
    let mut hits = vec![0_u64; planes_rules.len()];
    for (r, hit) in hits.iter_mut().enumerate() {
        let rule = planes_rules[r];
        for (i, or_word) in or_so_far.iter_mut().enumerate().take(w) {
            let keep_word = keep[i];
            let first_hit = keep_word & rule[i] & !*or_word;
            *hit += first_hit.count_ones() as u64;
            *or_word |= rule[i];
        }
    }
    hits
}

/// The effective mask of one site: the OR of every rule's firing bitmap
/// (intersected with `keep` for the emit).
fn mask_of(matrix: &SiteMatrix, s: usize, keep: &[u64]) -> Vec<u64> {
    let mut out = vec![0_u64; matrix.w];
    for r in 0..matrix.n_rules {
        for (dst, src) in out.iter_mut().zip(matrix.rule_mask(s, r)) {
            *dst |= src;
        }
    }
    for (dst, k) in out.iter_mut().zip(keep) {
        *dst &= k;
    }
    out
}

/// The mid stage's per-chunk counters (folded in chunk order).
#[derive(Default)]
struct MidPartial {
    masked_gt: u64,
    mask_hits: Vec<u64>,
    tail_hits: Vec<u64>,
    n_output: u64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn mid_stage(
    chunks: &[ChunkSites],
    matrices: &[SiteMatrix],
    contig_base: &[usize],
    keep: &[bool],
    plan: &FilterPlan,
    n_samples_out: usize,
    funnel: &mut FunnelPartial,
    budget: &ThreadBudget,
) -> std::collections::HashMap<u64, Arc<ChunkDecision>> {
    // Clean-padded keep words (the count math anchors on this).
    let w = keep.len().div_ceil(64);
    let mut keep_words = vec![0_u64; w];
    for (i, &k) in keep.iter().enumerate() {
        if k {
            keep_words[i / 64] |= 1 << (i % 64);
        }
    }
    // Independent per-chunk decisions, computed in parallel (middle stage,
    // doc_v3 03 O9); the funnel fold below stays in chunk order, so the
    // result is thread-count independent. The pool comes from the budget,
    // never the rayon global pool (doc_v8 D8-1).
    use rayon::prelude::*;
    let decisions: Vec<(u64, ChunkDecision, MidPartial)> =
        reed_engine::budget::run_auxiliary(budget, |parallel| {
            if parallel {
                chunks
                    .par_iter()
                    .map(|chunk| {
                        decide_chunk(
                            chunk,
                            matrices,
                            contig_base,
                            &keep_words,
                            plan,
                            n_samples_out,
                        )
                    })
                    .collect()
            } else {
                chunks
                    .iter()
                    .map(|chunk| {
                        decide_chunk(
                            chunk,
                            matrices,
                            contig_base,
                            &keep_words,
                            plan,
                            n_samples_out,
                        )
                    })
                    .collect()
            }
        })
        .expect("invariant: the auxiliary pool builds");
    let mut map = std::collections::HashMap::with_capacity(decisions.len());
    for (id, decision, partial) in decisions {
        funnel.masked_gt += partial.masked_gt;
        if funnel.mask_hits.len() < partial.mask_hits.len() {
            funnel.mask_hits.resize(partial.mask_hits.len(), 0);
        }
        for (dst, src) in funnel.mask_hits.iter_mut().zip(&partial.mask_hits) {
            *dst += src;
        }
        if funnel.tail_hits.len() < partial.tail_hits.len() {
            funnel.tail_hits.resize(partial.tail_hits.len(), 0);
        }
        for (dst, src) in funnel.tail_hits.iter_mut().zip(&partial.tail_hits) {
            *dst += src;
        }
        funnel.n_output += partial.n_output;
        map.insert(id, Arc::new(decision));
    }
    map
}

/// One chunk's decision plus its counters (no shared state).
#[allow(clippy::too_many_arguments)]
fn decide_chunk(
    chunk: &ChunkSites,
    matrices: &[SiteMatrix],
    contig_base: &[usize],
    keep_words: &[u64],
    plan: &FilterPlan,
    n_samples_out: usize,
) -> (u64, ChunkDecision, MidPartial) {
    let mut funnel = MidPartial {
        mask_hits: vec![0; plan.mask.len()],
        tail_hits: vec![0; plan.tail.len()],
        ..MidPartial::default()
    };
    let id = chunk.id;
    let decision = decide_chunk_inner(
        chunk,
        matrices,
        contig_base,
        keep_words,
        plan,
        n_samples_out,
        &mut funnel,
    );
    (id, decision, funnel)
}

/// Resolves a global site ordinal to its matrix and local index (the
/// stream is contig-sorted, so the bases partition it).
fn resolve<'a>(
    matrices: &'a [SiteMatrix],
    contig_base: &[usize],
    global: usize,
) -> (&'a SiteMatrix, usize) {
    let c = contig_base.partition_point(|&b| b <= global) - 1;
    let matrix = matrices
        .get(c)
        .expect("invariant: site ordinals stay inside the recorded contigs");
    (matrix, global - contig_base[c])
}

#[allow(clippy::too_many_arguments)]
fn decide_chunk_inner(
    chunk: &ChunkSites,
    matrices: &[SiteMatrix],
    contig_base: &[usize],
    keep_words: &[u64],
    plan: &FilterPlan,
    n_samples_out: usize,
    funnel: &mut MidPartial,
) -> ChunkDecision {
    let mut survivor_bits = vec![0_u64; chunk.n_sites.div_ceil(64)];
    let mut survivors: Vec<Option<SiteEmit>> = Vec::with_capacity(chunk.n_sites);
    // Per-chunk HWE memo (A6-16): the decisions run in parallel, so the
    // cache lives and dies with this chunk.
    let mut hwe = reed_core::math::hwe::HweCache::new();
    for s in 0..chunk.n_sites {
        let (matrix, local) = resolve(matrices, contig_base, chunk.start + s);
        let planes = matrix.planes(local);
        let rules: Vec<&[u64]> = (0..matrix.n_rules)
            .map(|r| matrix.rule_mask(local, r))
            .collect();
        let c = count_site_planes(planes, &rules, keep_words);
        funnel.masked_gt += c.masked_gt;
        for (r, hits) in rule_hits_of(&rules, keep_words).into_iter().enumerate() {
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
            n_half_called: c.n_half_called,
            has_gt: matrix.has_gt(local),
        };
        let tail_values = reed_filter::expr::SiteValues::from_counts(reed_filter::expr::GtCounts {
            ac: counts.ac,
            an: counts.an,
            n_missing: counts.n_missing,
            n_samples: n_samples_out as u64,
        });
        if eval_tail(
            &tail_values,
            &counts,
            n_samples_out,
            &plan.tail,
            &mut funnel.tail_hits,
            &mut hwe,
        ) {
            survivor_bits[s / 64] |= 1 << (s % 64);
            let mut infos = Vec::new();
            fill_info_set(&counts, &plan.recompute, &mut hwe, &mut infos);
            let mask_words = mask_of(matrix, local, keep_words);
            survivors.push(Some(SiteEmit {
                mask_words,
                infos,
                has_gt: counts.has_gt,
            }));
            funnel.n_output += 1;
        } else {
            survivors.push(None);
        }
    }
    let seam = chunk
        .seam_site
        .and_then(|s| survivors.get(s).cloned())
        .flatten();
    ChunkDecision {
        head_bits: chunk.head_bits.clone(),
        survivor_bits,
        survivors,
        seam,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 4 samples, one mask rule. s0 = 0/1, s1 = 1/1, s2 = ./. (masked by
    /// rule 0), s3 = 0/. — the pre-Y1e sketch fragment fixture.
    fn fixture() -> (Vec<Vec<u64>>, Vec<Vec<u64>>) {
        let planes = vec![
            vec![0b0010_u64], // slot0 alt: s1 (s0's slot0 is REF)
            vec![0b0100_u64], // slot0 missing: s2
            vec![0b0011_u64], // slot1 alt: s0, s1
            vec![0b1100_u64], // slot1 missing: s2, s3
        ];
        let rules = vec![vec![0b0100_u64]]; // rule 0 fires on s2
        (planes, rules)
    }

    #[test]
    fn counts_match_hand_computation() {
        let (planes, rules) = fixture();
        let refs: Vec<&[u64]> = planes.iter().map(Vec::as_slice).collect();
        let rule_refs: Vec<&[u64]> = rules.iter().map(Vec::as_slice).collect();
        let keep = [0b1111_u64];
        let c = count_site_planes([refs[0], refs[1], refs[2], refs[3]], &rule_refs, &keep);
        // s0 0/1: ac 1, an 2, het; s1 1/1: ac 2, an 2, hom_alt; s2 masked →
        // missing; s3 0/. half call: ac 0, an 1, ns.
        assert_eq!(c.ac, 3);
        assert_eq!(c.an, 5);
        assert_eq!(c.ns, 3);
        assert_eq!(c.n_missing, 1);
        assert_eq!(c.n_half_called, 1); // s3
        assert_eq!(c.het, 1);
        assert_eq!(c.hom_alt, 1);
        assert_eq!(c.hom_ref, 0);
        // s2 was raw-missing before the mask fired → no transition.
        assert_eq!(c.masked_gt, 0);
    }

    #[test]
    fn keep_restricts_the_denominator() {
        let (planes, rules) = fixture();
        let refs: Vec<&[u64]> = planes.iter().map(Vec::as_slice).collect();
        let rule_refs: Vec<&[u64]> = rules.iter().map(Vec::as_slice).collect();
        let keep = [0b1010_u64];
        let c = count_site_planes([refs[0], refs[1], refs[2], refs[3]], &rule_refs, &keep);
        // s1 1/1 (hom_alt), s3 0/. (half call): ac 2, an 3, ns 2, missing 0.
        assert_eq!(c.ac, 2);
        assert_eq!(c.an, 3);
        assert_eq!(c.ns, 2);
        assert_eq!(c.n_missing, 0);
        assert_eq!(c.n_half_called, 1); // s3
        assert_eq!(c.hom_alt, 1);
    }

    #[test]
    fn masked_called_sample_counts_as_transition() {
        // s0 = 0/1 masked; s1..s3 raw ./. (rule 0 fires on s0 only).
        let planes = [
            vec![0b0001_u64],
            vec![0b1110_u64],
            vec![0b0001_u64],
            vec![0b1110_u64],
        ];
        let rules = [vec![0b0001_u64]];
        let refs: Vec<&[u64]> = planes.iter().map(Vec::as_slice).collect();
        let rule_refs: Vec<&[u64]> = rules.iter().map(Vec::as_slice).collect();
        let keep = [0b1111_u64];
        let c = count_site_planes([refs[0], refs[1], refs[2], refs[3]], &rule_refs, &keep);
        assert_eq!(c.masked_gt, 1);
        assert_eq!(c.n_missing, 4);
        assert_eq!(c.ns, 0);
    }

    #[test]
    fn rule_hits_attribute_the_first_firing_rule() {
        // Two rules both firing on s0; the first owns the attribution.
        let rules = [vec![0b0001_u64], vec![0b0011_u64]];
        let rule_refs: Vec<&[u64]> = rules.iter().map(Vec::as_slice).collect();
        let keep = [0b1111_u64];
        assert_eq!(rule_hits_of(&rule_refs, &keep), vec![1, 1]);
    }
}

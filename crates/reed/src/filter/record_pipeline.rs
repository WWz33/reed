//! The per-record pipeline (`docs/rounds/v2/03-architecture.md` §8): **the single
//! implementation of "head → decode → mask → count → tail"**.
//!
//! Both pass strategies call [`RecordPipeline::process`]; they differ only in
//! what their `body` does with a record (the single pass edits and writes it,
//! the sketch pass fills per-sample planes and rule bitmaps) and in the
//! [`PipelineOptions`] they declare. Nothing in this file knows about
//! strategies, blocks or output.

use reed_core::Bitset;
use reed_core::count::SiteCounts;
use reed_core::genotype::{GenotypePolicies, GtCode, PolicyCounters, STATE_ALT, STATE_MISSING};
use reed_core::math::hwe::HweCache;
use reed_engine::error::EngineError;
use reed_filter::eval::head::{HeadScratch, eval_head};
use reed_filter::eval::tail::eval_tail;
use reed_filter::expr;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::vcf::genotype_scan::{GenotypeScanner, INT_MISSING, ScanOutput};
use reed_io::vcf::info::InfoIndex;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::{ContigTable, SiteView};

use super::job_err;

/// What each strategy asks of the shared stages. Every field is a strategy
/// difference; the stages themselves have no branches on the strategy.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PipelineOptions {
    /// Evaluate the mask for samples outside the keep list (the matrix
    /// pass0 records all samples, the single pass masks kept samples only).
    pub(crate) mask_all_samples: bool,
    /// Resolve CHROM to a header contig id; the index-writing path requires
    /// it (a record on an undeclared contig has no addressable tid).
    pub(crate) require_contig: bool,
    /// Evaluate the tail rules inline with this kept-sample count as the
    /// denominator; `None` defers them to the mid stage (sketch path).
    pub(crate) tail_n_samples: Option<u64>,
    /// Compute [`SiteCounts`] (and the mask attribution/masked-GT counters
    /// that come with it): the single pass needs them for the edit and the
    /// tail; the sketch pass counts from the sketches instead.
    pub(crate) counts: bool,
}

/// One processed record handed to the strategy body.
pub(crate) struct RecordView<'r, 'l> {
    /// The record's line.
    pub(crate) line: &'r LineRef<'l>,
    /// The record's site-level columns.
    pub(crate) site: &'r SiteView<'l>,
    /// The decoded genotypes and FORMAT columns.
    pub(crate) scan: &'r ScanOutput<'r>,
    /// Per-sample counts, present when `PipelineOptions::counts` is set.
    pub(crate) counts: Option<&'r SiteCounts>,
    /// Samples whose mask fired.
    pub(crate) masked: &'r Bitset,
    /// Firing mask rule per sample (`u16::MAX` = no rule fired).
    pub(crate) mask_rule: &'r [u16],
    /// 1-based record position.
    pub(crate) pos1: u64,
    /// Funnel counters: the body adds its own attribution here.
    pub(crate) funnel: &'r mut FunnelPartial,
}

/// Where one record ended up. The default is the head-failed outcome.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Processed {
    /// Whether the head rules passed.
    pub(crate) head_passed: bool,
    /// Whether the record reached the body (head passed and, when the tail
    /// runs inline, the tail passed).
    pub(crate) emitted: bool,
}

/// Per-worker record pipeline: decoder, contig table, head scratch and the
/// reusable mask/count state.
pub(crate) struct RecordPipeline {
    scanner: GenotypeScanner,
    contigs: ContigTable,
    counters: PolicyCounters,
    head: HeadScratch,
    masked: Bitset,
    /// The field rules (expressions evaluated separately), with their
    /// position in the plan's rule list for attribution.
    field_rules: Vec<reed_filter::plan::MaskRule>,
    field_rule_map: Vec<usize>,
    mask_rule: Vec<u16>,
    /// Reused per-sample `FORMAT` value table for `--set-gts .` expressions.
    fmt_storage: Vec<(Vec<u8>, Option<f64>)>,
    /// Memo for the exact HWE test behind the tail rules (A6-16).
    hwe: HweCache,
    counts: SiteCounts,
}

impl RecordPipeline {
    /// Creates the pipeline for `n_samples` sample columns and the plan's
    /// FORMAT columns.
    #[must_use]
    pub(crate) fn new(
        n_samples: usize,
        int_keys: &[&[u8]],
        float_keys: &[&[u8]],
        policies: GenotypePolicies,
        contig_names: &[Vec<u8>],
    ) -> Self {
        Self {
            scanner: GenotypeScanner::new(n_samples, int_keys, float_keys, policies),
            contigs: ContigTable::new(contig_names.iter().cloned()),
            counters: PolicyCounters::default(),
            head: HeadScratch::new(),
            masked: Bitset::new(n_samples),
            field_rules: Vec::new(),
            field_rule_map: Vec::new(),
            mask_rule: vec![u16::MAX; n_samples],
            fmt_storage: Vec::new(),
            hwe: HweCache::new(),
            counts: SiteCounts::default(),
        }
    }

    /// Adapts the pipeline to a plan (mask rule count, `--set-gts .`
    /// expression columns).
    pub(crate) fn resize_for(&mut self, plan: &FilterPlan) {
        self.field_rules = plan
            .mask
            .iter()
            .filter(|r| r.expr.is_none())
            .cloned()
            .collect();
        self.field_rule_map = plan
            .mask
            .iter()
            .enumerate()
            .filter(|(_, r)| r.expr.is_none())
            .map(|(i, _)| i)
            .collect();
        self.fmt_storage = plan
            .fmt_expr_fields
            .iter()
            .map(|f| (f.tag.clone(), None))
            .collect();
    }

    /// Runs head → decode → mask → count → tail for one record and calls
    /// `body` when the strategy gets the record.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records, policy rejections and
    /// output failures reported by `body`.
    pub(crate) fn process<F>(
        &mut self,
        line: &[u8],
        plan: &FilterPlan,
        keep: &[bool],
        opts: PipelineOptions,
        funnel: &mut FunnelPartial,
        body: F,
    ) -> Result<Processed, EngineError>
    where
        F: FnMut(&mut RecordView<'_, '_>) -> Result<(), EngineError>,
    {
        funnel.n_input += 1;
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::job("record line has fewer than 9 columns"));
        };
        let site = if opts.require_contig {
            SiteView::parse_declared(&lr, &mut self.contigs)
        } else {
            SiteView::parse(&lr)
        }
        .map_err(job_err)?;
        let pos1 = site.pos1;

        // ---- head stage (§11.2: order is semantics) ----
        let head_passed = eval_head(&lr, pos1, &plan.head, &mut funnel.head_hits, &mut self.head);
        if !head_passed {
            return Ok(Processed::default());
        }

        // ---- scan stage (the one decoder; policies counted once) ----
        // The scanner *adds* its triggers, so the per-record slot resets first.
        self.counters = PolicyCounters::default();
        let scan = self
            .scanner
            .scan(&lr, &site, &mut self.counters)
            .map_err(job_err)?;
        funnel.policy.add_decode(&self.counters);

        run_stages(
            &lr,
            &site,
            pos1,
            &scan,
            plan,
            keep,
            opts,
            funnel,
            head_passed,
            &self.field_rules,
            &self.field_rule_map,
            &mut self.masked,
            &mut self.mask_rule,
            &mut self.counts,
            &mut self.fmt_storage,
            &mut self.hwe,
            self.head.index_current.then_some(&self.head.index),
            plan.info_trusted,
            body,
        )
    }
}

/// The decoded inputs of a native (BCF) record (doc_v4 O-4): the pipeline
/// skips its own text scanner and consumes these instead.
pub(crate) struct NativeScan<'a> {
    /// One code per sample, in sample order.
    pub(crate) codes: &'a [GtCode],
    /// Integer FORMAT columns, aligned with the plan's integer fields.
    pub(crate) int_columns: &'a [Vec<u32>],
    /// Float FORMAT columns, aligned with the plan's float fields.
    pub(crate) float_columns: &'a [Vec<f32>],
    /// Policy triggers the native decode produced (report parity with the
    /// text path).
    pub(crate) counters: &'a PolicyCounters,
    /// The record's FORMAT carried a GT key (see [`ScanOutput::has_gt`]).
    pub(crate) has_gt: bool,
}

impl RecordPipeline {
    /// Like [`Self::process`], but the record's genotype and FORMAT values
    /// arrive decoded ([`NativeScan`]) instead of being scanned from text.
    /// The head rules and expressions still run on `line` — the BCF decoder
    /// renders only the site columns for them.
    ///
    /// # Errors
    ///
    /// As [`Self::process`].
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn process_scanned<F>(
        &mut self,
        line: &[u8],
        plan: &FilterPlan,
        keep: &[bool],
        opts: PipelineOptions,
        funnel: &mut FunnelPartial,
        native: NativeScan<'_>,
        body: F,
    ) -> Result<Processed, EngineError>
    where
        F: FnMut(&mut RecordView<'_, '_>) -> Result<(), EngineError>,
    {
        funnel.n_input += 1;
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::job("record line has fewer than 9 columns"));
        };
        let site = if opts.require_contig {
            SiteView::parse_declared(&lr, &mut self.contigs)
        } else {
            SiteView::parse(&lr)
        }
        .map_err(job_err)?;
        let pos1 = site.pos1;
        let head_passed = eval_head(&lr, pos1, &plan.head, &mut funnel.head_hits, &mut self.head);
        if !head_passed {
            return Ok(Processed::default());
        }
        funnel.policy.add_decode(native.counters);
        let scan = ScanOutput {
            codes: native.codes,
            int_columns: native.int_columns,
            float_columns: native.float_columns,
            alleles: None,
            has_gt: native.has_gt,
        };
        run_stages(
            &lr,
            &site,
            pos1,
            &scan,
            plan,
            keep,
            opts,
            funnel,
            head_passed,
            &self.field_rules,
            &self.field_rule_map,
            &mut self.masked,
            &mut self.mask_rule,
            &mut self.counts,
            &mut self.fmt_storage,
            &mut self.hwe,
            self.head.index_current.then_some(&self.head.index),
            plan.info_trusted,
            body,
        )
    }
}

/// The native BCF path's reusable decode buffers (doc_v4 O-4), shared by the
/// single-pass job and the pass0 jobs.
pub(crate) struct NativeScratch {
    /// The GT lookup table, built from the plan's policies so the decoded
    /// codes equal the text scanner's.
    pub(crate) table: std::sync::Arc<reed_io::bcf::GtTable>,
    /// The site-only text line the head rules and expressions read.
    pub(crate) site_line: Vec<u8>,
    /// The full rendered record (survivors only).
    pub(crate) full: Vec<u8>,
    /// GT codes, one per sample.
    pub(crate) codes: Vec<GtCode>,
    /// The plan's Integer FORMAT columns.
    pub(crate) int_columns: Vec<Vec<u32>>,
    /// The plan's Float FORMAT columns.
    pub(crate) float_columns: Vec<Vec<f32>>,
    /// Policy triggers from the GT decode.
    pub(crate) counters: PolicyCounters,
    /// Whether the last decoded record's FORMAT carried a GT key.
    pub(crate) has_gt: bool,
}

impl NativeScratch {
    /// Allocates the buffers for a plan's FORMAT keys. The GT table is
    /// shared across workers (A6-10: one 65,536-entry table per run, not
    /// one per worker).
    #[must_use]
    pub(crate) fn new(
        table: std::sync::Arc<reed_io::bcf::GtTable>,
        n_samples: usize,
        int_keys: &[Vec<u8>],
        float_keys: &[Vec<u8>],
    ) -> Self {
        Self {
            table,
            site_line: Vec::with_capacity(1 << 10),
            full: Vec::with_capacity(1 << 14),
            codes: Vec::with_capacity(n_samples),
            int_columns: vec![Vec::with_capacity(n_samples); int_keys.len()],
            float_columns: vec![Vec::with_capacity(n_samples); float_keys.len()],
            counters: PolicyCounters::default(),
            has_gt: false,
        }
    }

    /// Decodes one framed record into the buffers: the site line, the GT
    /// codes and the plan's FORMAT columns. Returns the record's position
    /// (needed by the callers' index keys).
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records and dictionary misses.
    pub(crate) fn decode(
        &mut self,
        dict: &reed_io::bcf::BcfDict,
        bytes: &[u8],
        int_keys: &[Vec<u8>],
        float_keys: &[Vec<u8>],
        n_samples: usize,
    ) -> Result<(), EngineError> {
        let rec = reed_io::bcf::record_ref(bytes).map_err(EngineError::from)?;
        let header = reed_io::bcf::site_header(&rec).map_err(EngineError::from)?;
        self.site_line.clear();
        reed_io::bcf::render_site_line(&rec, dict, &mut self.site_line)
            .map_err(EngineError::from)?;
        let contig = dict
            .contig(header.chrom)
            .ok_or_else(|| EngineError::job("BCF contig index is not in the dictionary"))?;
        let pos1 = u64::try_from(header.pos0).unwrap_or(0) + 1;
        self.counters = PolicyCounters::default();
        self.has_gt = reed_io::bcf::decode_record_gt(
            &rec,
            &header,
            &reed_io::bcf::GtQuery {
                gt_key: reed_io::bcf::format_key(dict, b"GT"),
                n_samples,
                contig_name: contig,
                pos1,
            },
            &self.table,
            &mut self.codes,
            &mut self.counters,
        )
        .map_err(EngineError::from)?;
        let int_keys: Vec<&[u8]> = int_keys.iter().map(Vec::as_slice).collect();
        let float_keys: Vec<&[u8]> = float_keys.iter().map(Vec::as_slice).collect();
        reed_io::bcf::decode_format_columns(
            &rec,
            &header,
            dict,
            &int_keys,
            &float_keys,
            &mut self.int_columns,
            &mut self.float_columns,
        )
        .map_err(EngineError::from)?;
        Ok(())
    }
}

/// The stages after the decoder — mask, count, tail, body — shared by the
/// text path ([`RecordPipeline::process`]) and the native path
/// ([`RecordPipeline::process_scanned`], doc_v4 O-4). The parameters are the
/// pipeline's fields, borrowed disjointly so the scan can hold the decoder.
#[allow(clippy::too_many_arguments)]
fn run_stages<'a, F>(
    lr: &LineRef<'a>,
    site: &SiteView<'a>,
    pos1: u64,
    scan: &ScanOutput<'_>,
    plan: &FilterPlan,
    keep: &[bool],
    opts: PipelineOptions,
    funnel: &mut FunnelPartial,
    head_passed: bool,
    field_rules: &[reed_filter::plan::MaskRule],
    field_rule_map: &[usize],
    masked: &mut Bitset,
    mask_rule: &mut [u16],
    counts: &mut SiteCounts,
    fmt_storage: &mut [(Vec<u8>, Option<f64>)],
    hwe: &mut HweCache,
    info_index: Option<&InfoIndex>,
    info_trusted: bool,
    mut body: F,
) -> Result<Processed, EngineError>
where
    F: FnMut(&mut RecordView<'_, '_>) -> Result<(), EngineError>,
{
    // ---- mask + count stage (A6-17: field rules run columnar — one
    // contiguous pass per rule over the scanned columns; the `--set-gts .`
    // expression path reuses the pipeline's FORMAT table in place, with no
    // per-sample allocation) ----
    masked.clear();
    mask_rule.fill(u16::MAX);
    *counts = SiteCounts::default();
    counts.has_gt = scan.has_gt;
    let evaluated = |i: usize| opts.mask_all_samples || keep.get(i).copied().unwrap_or(true);
    // Site context for the `--set-gts .` expressions (no genotype-derived
    // fields: the plan rejects those, so the counts are not needed).
    let mut expr_site = expr::SiteValues::from_line(lr, pos1, None);
    expr_site.info_index = info_index;
    expr_site.info_trusted = info_trusted;

    // Expression rules first (declaration order among themselves), one
    // sample at a time: the FORMAT slots update in place and the
    // evaluation borrows the table.
    if !fmt_storage.is_empty() {
        for i in 0..scan.codes.len() {
            if !evaluated(i) {
                continue;
            }
            for (slot, field) in fmt_storage.iter_mut().zip(&plan.fmt_expr_fields) {
                slot.1 = fmt_value_at(scan, field, i);
            }
            let mut sample_site = expr_site;
            sample_site.fmt = Some(fmt_storage);
            for (rule_index, rule) in plan.mask.iter().enumerate() {
                let Some(expr) = &rule.expr else {
                    continue;
                };
                let truth = expr::eval(expr, &sample_site);
                let fires = if rule.include { !truth } else { truth };
                if fires {
                    attribute_hit(
                        i,
                        rule_index,
                        keep.get(i).copied().unwrap_or(true),
                        opts.counts,
                        masked,
                        mask_rule,
                        funnel,
                    );
                    break;
                }
            }
        }
    }

    // Field rules columnar: each rule compares its whole scanned column and
    // the masked check preserves the first-hit attribution (the value
    // observation per rule is unchanged — several rules may share one
    // field, DP<5 and DP>100 both read DP).
    for (local, rule) in field_rules.iter().enumerate() {
        let rule_index = field_rule_map[local];
        let mask_on_missing = rule.on_missing == reed_filter::spec::OnMissing::Mask;
        if let Some(j) = rule.float_column {
            let Some(column) = scan.float_columns.get(j) else {
                continue;
            };
            for (i, &x) in column.iter().enumerate() {
                if masked.get(i)
                    || !evaluated(i)
                    || !(x.is_nan() && mask_on_missing)
                    || (!x.is_nan() && !rule.cmp.apply_f64(f64::from(x)))
                {
                    continue;
                }
                attribute_hit(
                    i,
                    rule_index,
                    keep.get(i).copied().unwrap_or(true),
                    opts.counts,
                    masked,
                    mask_rule,
                    funnel,
                );
            }
        } else {
            let Some(column) = scan.int_columns.get(rule.int_column) else {
                continue;
            };
            // Integer thresholds compare in the integer domain (A6-17: no
            // per-element u32→f64 round-trip; the semantics are the same —
            // u32 widens exactly and the op is the same enum).
            match rule.cmp.value {
                reed_filter::spec::Typed::Int(threshold) => {
                    let op = rule.cmp.op;
                    for (i, &x) in column.iter().enumerate() {
                        if masked.get(i)
                            || !evaluated(i)
                            || (x == INT_MISSING && !mask_on_missing)
                            || (x != INT_MISSING && !op.apply_i64(i64::from(x), threshold))
                        {
                            continue;
                        }
                        attribute_hit(
                            i,
                            rule_index,
                            keep.get(i).copied().unwrap_or(true),
                            opts.counts,
                            masked,
                            mask_rule,
                            funnel,
                        );
                    }
                }
                reed_filter::spec::Typed::Float(_) => {
                    for (i, &x) in column.iter().enumerate() {
                        if masked.get(i)
                            || !evaluated(i)
                            || (x == INT_MISSING && !mask_on_missing)
                            || (x != INT_MISSING && !rule.cmp.apply_f64(f64::from(x)))
                        {
                            continue;
                        }
                        attribute_hit(
                            i,
                            rule_index,
                            keep.get(i).copied().unwrap_or(true),
                            opts.counts,
                            masked,
                            mask_rule,
                            funnel,
                        );
                    }
                }
            }
        }
    }

    // `count_site` semantics from the decoded slots: masked calls count as
    // missing (the attribution the single pass reports).
    if opts.counts {
        for (i, &code) in scan.codes.iter().enumerate() {
            if !keep.get(i).copied().unwrap_or(true) {
                continue;
            }
            if !code.is_called_any() {
                counts.n_missing += 1;
            } else if masked.get(i) {
                funnel.masked_gt += 1;
                counts.n_missing += 1;
            } else {
                let (s0, s1) = (code.slot_state(0), code.slot_state(1));
                let n_alt = u64::from(s0 == STATE_ALT) + u64::from(s1 == STATE_ALT);
                let n_miss = u64::from(s0 == STATE_MISSING) + u64::from(s1 == STATE_MISSING);
                counts.ns += 1;
                counts.an += 2 - n_miss;
                counts.ac += n_alt;
                match (n_miss, n_alt) {
                    (0, 0) => counts.hom_ref += 1,
                    (0, 1) => counts.het += 1,
                    (0, _) => counts.hom_alt += 1,
                    // Half-calls feed AN/AC/NS but no genotype class; the
                    // `-g/--genotype` rule reads them as missing (bcftools
                    // `bcf_gt_type` UNKN).
                    _ => counts.n_half_called += 1,
                }
            }
        }
    }

    // ---- tail stage (single pass; the sketch path defers it) ----
    if let Some(n_out) = opts.tail_n_samples.filter(|_| head_passed) {
        let mut values = reed_filter::expr::SiteValues::from_line(
            lr,
            pos1,
            Some(reed_filter::expr::GtCounts {
                ac: counts.ac,
                an: counts.an,
                n_missing: counts.n_missing,
                n_samples: n_out,
            }),
        );
        values.info_index = info_index;
        values.info_trusted = info_trusted;
        if !eval_tail(
            &values,
            counts,
            usize::try_from(n_out).unwrap_or(usize::MAX),
            &plan.tail,
            &mut funnel.tail_hits,
            hwe,
        ) {
            return Ok(Processed {
                head_passed: true,
                emitted: false,
            });
        }
    }

    let counts = opts.counts.then_some(&*counts);
    let mut view = RecordView {
        line: lr,
        site,
        scan,
        counts,
        masked,
        mask_rule,
        pos1,
        funnel,
    };
    body(&mut view)?;
    Ok(Processed {
        head_passed,
        emitted: true,
    })
}

/// One sample's scanned `FORMAT` value for a `--set-gts .` expression field
/// (`None` = missing).
fn fmt_value_at(
    scan: &ScanOutput<'_>,
    field: &reed_filter::plan::FmtExprField,
    i: usize,
) -> Option<f64> {
    match field.float_column {
        Some(j) => scan
            .float_columns
            .get(j)
            .and_then(|column| column.get(i))
            .copied()
            .filter(|x| !x.is_nan())
            .map(f64::from),
        None => scan
            .int_columns
            .get(field.int_column)
            .and_then(|column| column.get(i))
            .copied()
            .filter(|&x| x != INT_MISSING)
            .map(f64::from),
    }
}

/// Attributes a mask hit to sample `i`: the bit, the owning rule and the
/// funnel counter (the first firing rule owns the sample; already-missing
/// samples are attributed too, matching the reported called→missing
/// transitions).
fn attribute_hit(
    i: usize,
    rule_index: usize,
    kept: bool,
    counts_enabled: bool,
    masked: &mut Bitset,
    mask_rule: &mut [u16],
    funnel: &mut FunnelPartial,
) {
    masked.set(i, true);
    if let Ok(slot) = u16::try_from(rule_index) {
        mask_rule[i] = slot;
    }
    if counts_enabled
        && kept
        && let Some(slot) = funnel.mask_hits.get_mut(rule_index)
    {
        *slot += 1;
    }
}

/// Index of `GT` among the colon-separated FORMAT keys.
#[must_use]
pub(super) fn gt_subfield_index(format_text: &[u8]) -> Option<u8> {
    format_text
        .split(|&b| b == b':')
        .position(|k| k == b"GT")
        .and_then(|i| u8::try_from(i).ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_index_found_or_none() {
        assert_eq!(gt_subfield_index(b"GT:DP:GQ"), Some(0));
        assert_eq!(gt_subfield_index(b"DP:GT"), Some(1));
        assert_eq!(gt_subfield_index(b"DP:GQ"), None);
    }
}

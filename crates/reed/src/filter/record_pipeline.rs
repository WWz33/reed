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
use reed_engine::error::EngineError;
use reed_filter::eval::head::{HeadScratch, eval_head};
use reed_filter::eval::mask::eval_mask;
use reed_filter::eval::tail::eval_tail;
use reed_filter::expr::{self, FmtValue};
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::vcf::genotype_scan::{GenotypeScanner, INT_MISSING, ScanOutput};
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
    mask_values: Vec<Option<f64>>,
    /// The field rules (expressions evaluated separately), with their
    /// position in the plan's rule list for attribution.
    field_rules: Vec<reed_filter::plan::MaskRule>,
    field_rule_map: Vec<usize>,
    mask_rule: Vec<u16>,
    /// Reused per-sample `FORMAT` value table for `--set-gts .` expressions.
    fmt_values: Vec<FmtValue<'static>>,
    fmt_storage: Vec<(Vec<u8>, Option<f64>)>,
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
            mask_values: vec![None; 0],
            field_rules: Vec::new(),
            field_rule_map: Vec::new(),
            mask_rule: vec![u16::MAX; n_samples],
            fmt_values: Vec::new(),
            fmt_storage: Vec::new(),
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
        self.mask_values = vec![None; self.field_rules.len()];
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
            &mut self.mask_values,
            &mut self.masked,
            &mut self.mask_rule,
            &mut self.counts,
            &mut self.fmt_storage,
            &mut self.fmt_values,
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
            &mut self.mask_values,
            &mut self.masked,
            &mut self.mask_rule,
            &mut self.counts,
            &mut self.fmt_storage,
            &mut self.fmt_values,
            body,
        )
    }
}

/// The native BCF path's reusable decode buffers (doc_v4 O-4), shared by the
/// single-pass job and the pass0 jobs.
pub(crate) struct NativeScratch {
    /// The GT lookup table, built from the plan's policies so the decoded
    /// codes equal the text scanner's.
    pub(crate) table: reed_io::bcf::GtTable,
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
}

impl NativeScratch {
    /// Allocates the buffers for a plan's FORMAT keys.
    #[must_use]
    pub(crate) fn new(
        policies: GenotypePolicies,
        n_samples: usize,
        int_keys: &[Vec<u8>],
        float_keys: &[Vec<u8>],
    ) -> Self {
        Self {
            table: reed_io::bcf::GtTable::new(policies),
            site_line: Vec::with_capacity(1 << 10),
            full: Vec::with_capacity(1 << 14),
            codes: Vec::with_capacity(n_samples),
            int_columns: vec![Vec::with_capacity(n_samples); int_keys.len()],
            float_columns: vec![Vec::with_capacity(n_samples); float_keys.len()],
            counters: PolicyCounters::default(),
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
        reed_io::bcf::decode_record_gt(
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
    mask_values: &mut [Option<f64>],
    masked: &mut Bitset,
    mask_rule: &mut [u16],
    counts: &mut SiteCounts,
    fmt_storage: &mut [(Vec<u8>, Option<f64>)],
    fmt_values: &mut Vec<FmtValue<'static>>,
    mut body: F,
) -> Result<Processed, EngineError>
where
    F: FnMut(&mut RecordView<'_, '_>) -> Result<(), EngineError>,
{
    // ---- mask + count stage ----
    masked.clear_padding();
    for bit in 0..masked.len() {
        masked.set(bit, false);
    }
    mask_rule.fill(u16::MAX);
    *counts = SiteCounts::default();
    let mask_rules = &plan.mask;
    // Site context for the `--set-gts .` expressions (no genotype-derived
    // fields: the plan rejects those, so the counts are not needed).
    let expr_site = expr::SiteValues::from_line(lr, pos1, None);
    for (i, &code) in scan.codes.iter().enumerate() {
        let kept = keep.get(i).copied().unwrap_or(true);
        if !kept && !opts.mask_all_samples {
            continue;
        }
        // One observation per FIELD rule: rule i reads its own scanned
        // column, so several rules may share one field (DP<5 and
        // DP>100 both read DP). Expression rules (`--set-gts .`) read
        // the per-sample FORMAT table below.
        for (rule, v) in field_rules.iter().zip(mask_values.iter_mut()) {
            *v = match rule.float_column {
                Some(j) => scan
                    .float_columns
                    .get(j)
                    .and_then(|column| column.get(i))
                    .copied()
                    .filter(|x| !x.is_nan())
                    .map(f64::from),
                None => scan
                    .int_columns
                    .get(rule.int_column)
                    .and_then(|column| column.get(i))
                    .copied()
                    .filter(|&x| x != INT_MISSING)
                    .map(f64::from),
            };
        }
        // Expression rules read FORMAT values; field rules read their
        // own column (one observation per rule).
        let mut hit = None;
        if !fmt_storage.is_empty() {
            for (slot, field) in fmt_storage.iter_mut().zip(&plan.fmt_expr_fields) {
                let value = match field.float_column {
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
                };
                slot.1 = value;
            }
            fmt_values.clear();
        }
        for (rule_index, rule) in mask_rules.iter().enumerate() {
            let Some(expr) = &rule.expr else {
                continue;
            };
            let values: Vec<FmtValue<'_>> = fmt_storage
                .iter()
                .map(|(tag, num)| FmtValue { tag, num: *num })
                .collect();
            let mut sample_site = expr_site;
            sample_site.fmt = Some(&values);
            let truth = expr::eval(expr, &sample_site);
            let fires = if rule.include { !truth } else { truth };
            if fires {
                hit = Some(rule_index);
                break;
            }
        }
        if hit.is_none()
            && let Some(local) = eval_mask(mask_values, field_rules)
        {
            hit = field_rule_map.get(local).copied();
        }
        if let Some(rule_index) = hit {
            masked.set(i, true);
            if let Ok(slot) = u16::try_from(rule_index) {
                mask_rule[i] = slot;
            }
            if opts.counts
                && kept
                && let Some(slot) = funnel.mask_hits.get_mut(rule_index)
            {
                // Attribution: the first firing rule owns this sample
                // (already-missing samples are attributed too, matching
                // the reported called→missing transitions).
                *slot += 1;
            }
        }
        if !opts.counts || !kept {
            continue;
        }
        // `count_site` semantics from the decoded slots: masked calls
        // count as missing (the attribution the single pass reports).
        if !code.is_called_any() {
            counts.n_missing += 1;
        } else if hit.is_some() {
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
                // Half-calls feed AN/AC/NS but no genotype class.
                _ => {}
            }
        }
    }

    // ---- tail stage (single pass; the sketch path defers it) ----
    if let Some(n_out) = opts.tail_n_samples.filter(|_| head_passed) {
        let values = reed_filter::expr::SiteValues::from_line(
            lr,
            pos1,
            Some(reed_filter::expr::GtCounts {
                ac: counts.ac,
                an: counts.an,
                n_missing: counts.n_missing,
                n_samples: n_out,
            }),
        );
        if !eval_tail(
            &values,
            counts,
            usize::try_from(n_out).unwrap_or(usize::MAX),
            &plan.tail,
            &mut funnel.tail_hits,
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

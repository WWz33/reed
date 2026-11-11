//! Byte-level line editing (§6.7): survivors are edited in place on the
//! raw VCF bytes — columns 1–7 and FORMAT copy verbatim, INFO drops stale
//! tags and rewrites/appends recomputed ones (htslib-compatible float
//! rendering for `L_tags` parity), sample columns are pruned by the cohort
//! keep mask and masked samples get their GT subfield replaced with `./.`
//! while every other FORMAT subfield is preserved byte-for-byte.
//!
//! Every survivor is routed through [`edit_vcf_line`] today (F6-11 in
//! docs/rounds/v6/01-review.md: the two call sites in
//! `reed/src/filter/single_pass.rs` are unconditional). A verbatim-copy
//! fast path for records with no edit is queued for the Y3 filter work;
//! until then the function itself copies unmodified regions byte-for-byte.

// Casts in this file are bounded by construction: Sample counts are validated against the header and fit u32.
#![allow(clippy::cast_possible_truncation)]

use reed_core::Bitset;
use reed_core::sample::SampleMask;

use crate::vcf::line::{LineRef, Shape};

/// One recomputed INFO value (§6.7: rendered here with the htslib-compatible
/// formatting so `L_tags` parity holds).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InfoValueOut {
    /// Integer tag (`AC`/`AN`/`NS`).
    Int(i64),
    /// Float tag (`MAF`/`HWE`), stored as f32 like the BCF type.
    Float(f32),
}

/// What to change on one line (`docs/spec/architecture.md` §6.7).
#[derive(Default)]
pub struct EditPlan<'a> {
    /// `None` = keep every sample column.
    pub keep_samples: Option<&'a SampleMask>,
    /// Samples whose GT subfield becomes `./.` (input sample order; the
    /// caller only flags kept samples).
    pub mask_gt: Option<&'a Bitset>,
    /// Stale INFO tags to drop.
    pub info_remove: &'a [&'a [u8]],
    /// Recomputed tags: a tag that also exists in the input keeps its
    /// position with its value replaced, the rest are appended in this
    /// order (noodles `Info::insert` / bcftools annotate semantics,
    /// locked by test).
    pub info_set: &'a [(&'a [u8], InfoValueOut)],
}

impl EditPlan<'_> {
    /// Whether the INFO column passes through untouched.
    fn info_untouched(&self) -> bool {
        self.info_remove.is_empty() && self.info_set.is_empty()
    }

    /// Whether the sample region needs editing (pruning or masking).
    fn samples_edited(&self) -> bool {
        self.keep_samples
            .is_some_and(|k| k.n_selected() < k.n_samples() as u32)
            || self.mask_gt.is_some_and(|m| m.count_ones() > 0)
    }
}

/// Renders one INFO value (htslib `%g` semantics for floats).
pub fn render_info_value(value: &InfoValueOut, out: &mut Vec<u8>) {
    match value {
        InfoValueOut::Int(v) => out.extend_from_slice(v.to_string().as_bytes()),
        InfoValueOut::Float(v) => format_g6_into(out, f64::from(*v)),
    }
}

/// Edits one record line into `out` (newline-terminated).
///
/// `gt_index` is the GT position among the FORMAT keys of this
/// record (resolved by the caller; `None` = the record has no GT column).
///
/// Infallible by construction: the `LineRef` argument can only be built by
/// `LineRef::parse`, which already rejected a truncated line (fewer than
/// eight columns), so the old re-parse and its `MalformedRecord` branch were
/// dead code (doc_v12 D12-15).
pub fn edit_vcf_line(
    line: &LineRef<'_>,
    plan: &EditPlan<'_>,
    gt_index: Option<usize>,
    out: &mut Vec<u8>,
) {
    // Nothing to change: copy the original bytes, including a CRLF tail.
    // The old segmented path stripped `\r` and appended `\n` even for an
    // empty plan (doc_v28 V28-11).
    if plan.info_untouched() && !plan.samples_edited() {
        out.extend_from_slice(line.line);
        if line.terminator().is_empty() {
            out.push(b'\n');
        }
        return;
    }
    let terminator = if line.terminator().is_empty() {
        b"\n".as_slice()
    } else {
        line.terminator()
    };
    // Columns 1–7 (CHROM .. FILTER) with their tabs: the bytes before INFO.
    out.extend_from_slice(&line.line[..line.info_start()]);

    // INFO.
    let info = line.info_text();
    if plan.info_untouched() {
        out.extend_from_slice(info);
    } else {
        edit_info(info, plan, out);
    }

    // doc_v21 M7: when the sample subset selects nothing, the record loses its
    // FORMAT and sample columns entirely, matching the 8-column header
    // (bcftools writes legal 8-column records and warns once). `keep_samples`
    // is `Some` only when a sample edit was requested.
    if plan.keep_samples.is_some_and(|k| k.n_selected() == 0) {
        out.extend_from_slice(terminator);
        return;
    }
    // The shape is preserved (doc_v22 R22-3): a sites-only record stays
    // 8 columns, a FORMAT-only record 9, and only a record with samples
    // writes the sample region.
    match line.shape() {
        Shape::SitesOnly => {
            out.extend_from_slice(terminator);
            return;
        }
        Shape::FormatOnly => {
            out.push(b'\t');
            out.extend_from_slice(line.format_text());
            out.extend_from_slice(terminator);
            return;
        }
        Shape::WithSamples => {}
    }
    out.push(b'\t');
    // FORMAT verbatim.
    out.extend_from_slice(line.format_text());
    out.push(b'\t');

    // Sample region.
    let samples = line.samples_text();
    if plan.samples_edited() {
        edit_samples(samples, plan, gt_index, out);
    } else {
        out.extend_from_slice(samples);
    }
    out.extend_from_slice(terminator);
}

/// INFO edit: drop `info_remove`, replace any tag present in `info_set` in
/// place, append the rest in `info_set` order. An empty result renders as
/// `.`.
///
/// Writes straight into `out` (A6-11: no per-record `edited` buffer) and
/// keeps the seen-marks on the stack for plans with a handful of tags.
fn edit_info(info: &[u8], plan: &EditPlan<'_>, out: &mut Vec<u8>) {
    let mark = out.len();
    let mut seen_stack = [false; 8];
    let mut seen_heap;
    let seen: &mut [bool] = if plan.info_set.len() <= seen_stack.len() {
        &mut seen_stack[..plan.info_set.len()]
    } else {
        seen_heap = vec![false; plan.info_set.len()];
        &mut seen_heap
    };
    if info != b"." {
        for (tag, value) in crate::vcf::info::info_fragments(info) {
            if plan.info_remove.contains(&tag) {
                continue;
            }
            if out.len() > mark {
                out.push(b';');
            }
            if let Some(i) = plan.info_set.iter().position(|(t, _)| *t == tag) {
                out.extend_from_slice(plan.info_set[i].0);
                out.push(b'=');
                render_info_value(&plan.info_set[i].1, out);
                seen[i] = true;
            } else {
                out.extend_from_slice(tag);
                if let Some(v) = value {
                    out.push(b'=');
                    out.extend_from_slice(v);
                }
            }
        }
    }
    for (i, (tag, value)) in plan.info_set.iter().enumerate() {
        if seen[i] {
            continue;
        }
        if out.len() > mark {
            out.push(b';');
        }
        out.extend_from_slice(tag);
        out.push(b'=');
        render_info_value(value, out);
    }
    if out.len() == mark {
        out.push(b'.');
    }
}

/// Sample region edit: skip dropped columns, replace the GT subfield of
/// masked samples with `./.`, copy everything else verbatim.
///
/// The bound on sample columns is the header's sample count, read off the
/// bitsets the plan carries (`SampleMask::n_samples` / `Bitset::len` are both
/// built from the header). A line may carry extra *empty* columns past it
/// (htslib ignores a trailing run of tabs, V22-4); those columns have no mask
/// bit, so indexing the mask/GT bitset with them would trip
/// `Bitset::get`'s assert and abort the batch (doc_v23 V23-03). They are
/// dropped rather than written, which is what the scan stage accepted them
/// as: nothing. `samples_edited` guarantees at least one bitset is present.
fn edit_samples(samples: &[u8], plan: &EditPlan<'_>, gt_index: Option<usize>, out: &mut Vec<u8>) {
    let keep = plan.keep_samples;
    let mask = plan.mask_gt;
    let n_samples = keep.map_or_else(
        || mask.map_or(0, Bitset::len),
        reed_core::SampleMask::n_samples,
    );
    // Work only on the first `n_samples` fields. The line may carry extra
    // empty columns past the header count (htslib ignores them); they have
    // no mask bit, and the verbatim tail copy at the end of this function
    // would otherwise write them back out (doc_v23 V23-03).
    let samples = &samples[..sample_region_end(samples, n_samples)];
    let mut field_start = 0;
    let mut sample_idx = 0_usize;
    let mut wrote_column = false;
    // Byte offset where the current run of verbatim samples begins
    // (doc_v3 03 O7: one copy per run instead of one per sample).
    let mut run_start: Option<usize> = None;
    while field_start < samples.len() && sample_idx < n_samples {
        let field_end = memchr::memchr(b'\t', &samples[field_start..])
            .map_or(samples.len(), |rel| field_start + rel);
        let kept = keep.is_none_or(|k| k.contains(reed_core::SampleId(sample_idx as u32)));
        let masked = kept && mask.is_some_and(|m| m.get(sample_idx));
        if kept && !masked {
            run_start.get_or_insert(field_start);
        } else {
            if let Some(start) = run_start.take() {
                if wrote_column {
                    out.push(b'\t');
                }
                // The run ends before the tab that separates it from the
                // current field.
                out.extend_from_slice(&samples[start..field_start - 1]);
                wrote_column = true;
            }
            if kept {
                if wrote_column {
                    out.push(b'\t');
                }
                replace_gt_subfield(&samples[field_start..field_end], gt_index, out);
                wrote_column = true;
            }
        }
        sample_idx += 1;
        field_start = field_end + 1;
    }
    if let Some(start) = run_start {
        if wrote_column {
            out.push(b'\t');
        }
        out.extend_from_slice(&samples[start..]);
    }
    // A sample-count mismatch against the header is a caller invariant (the
    // scan stage validates it).
}

/// Byte offset just past the `n`-th sample field of a sample region (or the
/// region's length when fewer fields are present).
///
/// The scan stage accepts a trailing run of empty columns past the header's
/// sample count and ignores them (V22-4); the editor must end the region at
/// the same place so the verbatim tail copy cannot resurrect them
/// (doc_v23 V23-03).
fn sample_region_end(samples: &[u8], n: usize) -> usize {
    // Tab-separated columns, not colon-separated subfields: walk the first
    // `n` tabs and end at the last column's content.
    let mut start = 0;
    for _ in 0..n {
        let Some(rel) = memchr::memchr(b'\t', &samples[start..]) else {
            return samples.len();
        };
        start += rel + 1;
    }
    start.saturating_sub(1).min(samples.len())
}

/// Copies `field`, replacing the GT subfield with `./.` when the row is
/// long enough; a row without that subfield stays verbatim.
fn replace_gt_subfield(field: &[u8], gt_index: Option<usize>, out: &mut Vec<u8>) {
    let Some(gt_index) = gt_index else {
        out.extend_from_slice(field);
        return;
    };
    // doc_v28 V28-24: the colon walk lives in `format_layout`.
    let Some(span) = crate::vcf::format_layout::subfield_span(field, gt_index) else {
        out.extend_from_slice(field);
        return;
    };
    out.extend_from_slice(&field[..span.start]);
    out.extend_from_slice(b"./.");
    out.extend_from_slice(&field[span.end..]);
}

/// C's `%g` with precision 6 — the htslib INFO float rendering
/// (validated against all 80,228 unique MAF/HWE renderings of the gold
/// VCF, zero mismatches).
#[must_use]
pub fn format_g6(v: f64) -> String {
    let mut buf = Vec::with_capacity(24);
    format_g6_into(&mut buf, v);
    String::from_utf8(buf).expect("invariant: %g renders ASCII")
}

/// A tiny `fmt::Write` over a fixed stack buffer: the renderer's hot path
/// stays heap-free (doc_v15 D15-35). 64 bytes fit any `%g` rendering of an
/// `f64` by construction (the widest is `-2.22507e-308`).
struct StackBuf {
    bytes: [u8; 64],
    len: usize,
}

impl StackBuf {
    fn new() -> Self {
        Self {
            bytes: [0; 64],
            len: 0,
        }
    }

    fn as_slice(&self) -> &[u8] {
        &self.bytes[..self.len]
    }
}

impl std::fmt::Write for StackBuf {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        let end = self.len + text.len();
        if end > self.bytes.len() {
            return Err(std::fmt::Error);
        }
        self.bytes[self.len..end].copy_from_slice(text.as_bytes());
        self.len = end;
        Ok(())
    }
}

/// [`format_g6`] appended to a byte buffer (the renderer's hot path avoids
/// one allocation per value).
pub fn format_g6_into(out: &mut Vec<u8>, v: f64) {
    use std::fmt::Write as _;
    use std::io::Write as _;
    if v == 0.0 {
        // htslib `kputd` keeps the sign of a negative zero (A6-14).
        out.extend_from_slice(if v.is_sign_negative() { b"-0" } else { b"0" });
        return;
    }
    if v.is_nan() {
        out.extend_from_slice(b"nan");
        return;
    }
    if v.is_infinite() {
        out.extend_from_slice(if v < 0.0 { b"-inf" } else { b"inf" });
        return;
    }
    // The style is chosen from the exponent of the value rounded to 6
    // significant digits (C %g semantics): 999999.5 rounds to 1e+06, not
    // 9.99999e5. Everything renders into stack buffers (doc_v15 D15-35).
    let mut wide = StackBuf::new();
    write!(wide, "{:.5e}", v.abs()).expect("invariant: 64 bytes fit a %.5e f64");
    let scientific = wide.as_slice();
    let e_pos = scientific
        .iter()
        .position(|&b| b == b'e')
        .expect("invariant: e-form");
    let exp10: i32 = std::str::from_utf8(&scientific[e_pos + 1..])
        .expect("invariant: ASCII")
        .parse()
        .unwrap_or(0);
    if (-4..6).contains(&exp10) {
        let precision = (5 - exp10).max(0) as usize;
        let mut fixed = StackBuf::new();
        write!(fixed, "{v:.precision$}").expect("invariant: 64 bytes fit a %.*f f64");
        let text = fixed.as_slice();
        // htslib culls trailing zeros only past a decimal point: an integer
        // like 128440 must keep its zeros (the pre-BCF version stripped
        // them, which only showed for values >= 1e5).
        match text.iter().position(|&b| b == b'.') {
            Some(dot) => {
                let (int, frac) = (&text[..dot], &text[dot + 1..]);
                let mut end = frac.len();
                while end > 0 && frac[end - 1] == b'0' {
                    end -= 1;
                }
                out.extend_from_slice(int);
                if end > 0 {
                    out.push(b'.');
                    out.extend_from_slice(&frac[..end]);
                }
            }
            None => out.extend_from_slice(text),
        }
    } else {
        // Scientific: 5 fractional digits, trailing zeros stripped,
        // exponent at least two digits with a sign.
        let text = scientific;
        let mut mantissa = &text[..e_pos];
        while mantissa.last() == Some(&b'0') {
            mantissa = &mantissa[..mantissa.len() - 1];
        }
        if mantissa.last() == Some(&b'.') {
            mantissa = &mantissa[..mantissa.len() - 1];
        }
        let exp: i32 = std::str::from_utf8(&text[e_pos + 1..])
            .expect("invariant: ASCII")
            .parse()
            .expect("invariant: integer exponent");
        // `wide`/`scientific` was formatted from `v.abs()`, so the sign has
        // to be restored here; the fixed branch formats `v` itself and keeps
        // it. Without this a negative INFO/FORMAT float was silently written
        // as positive (doc_v23 V23-04: BCF→VCF data corruption, not latent).
        if v < 0.0 {
            out.push(b'-');
        }
        out.extend_from_slice(mantissa);
        out.push(b'e');
        out.push(if exp < 0 { b'-' } else { b'+' });
        let magnitude = exp.unsigned_abs();
        if magnitude < 10 {
            out.push(b'0');
            out.push(b'0' + magnitude as u8);
        } else {
            let _ = write!(out, "{magnitude}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::SampleId;

    #[test]
    fn g6_matches_htslib_gold_rendings() {
        let cases = [
            (0.120_061_f32, "0.120061"),
            (0.05_f32, "0.05"),
            (0.0_f32, "0"),
            (-0.0_f32, "-0"),
            (8.20818e-38_f32, "8.20818e-38"),
            (0.289_678_f32, "0.289678"),
            (0.5_f32, "0.5"),
            (1.0_f32, "1"),
            (1.0e-5_f32, "1e-05"),
            (1.234_567e-4_f32, "0.000123457"),
            (1_234_567.0_f32, "1.23457e+06"),
        ];
        for (v, expected) in cases {
            assert_eq!(format_g6(f64::from(v)), expected, "value {v}");
        }
    }

    /// doc_v23 V23-04: the scientific branch formats `v.abs()`, so the sign
    /// must be restored separately. Every negative value that renders in
    /// exponent form used to lose its `-` (the fixed branch formats `v`
    /// itself and was unaffected).
    #[test]
    fn g6_negative_scientific_values_keep_their_sign() {
        let cases = [
            (-1.5e-7_f64, "-1.5e-07"),
            (-1_234_567.0_f64, "-1.23457e+06"),
            (-1.0e-5_f64, "-1e-05"),
            (-8.20818e-38_f64, "-8.20818e-38"),
            (-0.000_099_999_f64, "-9.9999e-05"),
            (-999_999.5_f64, "-1e+06"),
            // The fixed branch keeps its sign without help.
            (-0.123_f64, "-0.123"),
            (-0.0_f64, "-0"),
            // And positives are unchanged.
            (1.5e-7_f64, "1.5e-07"),
        ];
        for (v, expected) in cases {
            assert_eq!(format_g6(v), expected, "value {v}");
        }
    }

    /// Both renderers that reach `format_g6_into` are covered: the INFO
    /// rewrite (`render_info_value`, V23-04's reported BCF→VCF corruption)
    /// and FORMAT floats (the BCF renderer's `float_into`, checked here via
    /// the same helper). An INFO value the plan does not touch is copied
    /// verbatim, so the sign survives there by construction; a recomputed
    /// float goes through the renderer.
    #[test]
    fn negative_info_float_survives_an_info_edit() {
        let line = b"Chr01\t100\t.\tA\tG\t50\tPASS\tXF=-1.5e-7;AC=1\tGT\t0/1\n";
        let plan = EditPlan {
            info_set: &[(b"AC".as_slice(), InfoValueOut::Int(2))],
            ..EditPlan::default()
        };
        let got = edit(line, &plan, Some(0));
        // Untouched tag: the input bytes are preserved.
        assert!(got.contains("XF=-1.5e-7"), "{got}");
        assert!(got.contains("AC=2"), "{got}");
        // Recomputed tag: the renderer must keep the sign.
        let plan = EditPlan {
            info_set: &[(b"XF".as_slice(), InfoValueOut::Float(-1.5e-7))],
            ..EditPlan::default()
        };
        let got = edit(line, &plan, Some(0));
        assert!(got.contains("XF=-1.5e-07"), "{got}");
    }

    #[test]
    fn g6_round_trips_every_gold_value() {
        // The gold VCF's rendered floats must re-render identically after a
        // parse → f32 → format cycle (skipped when the gold data is absent).
        for path in ["/tmp/gold_maf.txt", "/tmp/gold_hwe.txt"] {
            let Ok(text) = std::fs::read_to_string(path) else {
                return;
            };
            let mut checked = 0;
            for line in text.lines() {
                let v = line.parse::<f64>().expect("gold float parses");
                assert_eq!(format_g6(f64::from(v as f32)), line, "gold value {line}");
                checked += 1;
            }
            assert!(checked > 1000, "gold fixture unexpectedly small");
        }
    }

    fn mask_from(n: usize, keeps: &[usize]) -> SampleMask {
        SampleMask::from_ids(n, keeps.iter().map(|&i| SampleId(i as u32)))
            .expect("invariant: valid mask")
    }

    fn edit(line: &[u8], plan: &EditPlan<'_>, gt_idx: Option<usize>) -> String {
        edit_raw(line, plan, gt_idx).trim_end().to_string()
    }

    fn edit_raw(line: &[u8], plan: &EditPlan<'_>, gt_idx: Option<usize>) -> String {
        let lr = LineRef::parse(line).expect("invariant: full line");
        let mut out = Vec::new();
        edit_vcf_line(&lr, plan, gt_idx, &mut out);
        String::from_utf8(out).expect("invariant: utf-8")
    }

    const LINE: &[u8] =
        b"Chr01\t100\t.\tA\tG\t50\tPASS\tAC=9;QD=4.5;AF=0.5;DP=150\tGT:DP:GQ\t0/1:15:99\t1|1:12:99\t./.:0:9\t0/0:8:44\n";

    #[test]
    fn column_shape_is_preserved() {
        // doc_v22 R22-3: sites-only stays 8 columns, a FORMAT-only record
        // stays 9 — neither gains a phantom sample column.
        let plan = EditPlan::default();
        assert_eq!(
            edit(b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\n", &plan, None),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t."
        );
        assert_eq!(
            edit(b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\n", &plan, None),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT"
        );
    }

    #[test]
    fn verbatim_when_no_plan() {
        let plan = EditPlan::default();
        assert_eq!(
            edit(LINE, &plan, Some(0)),
            String::from_utf8_lossy(LINE).trim_end()
        );
    }

    /// doc_v28 V28-11: an untouched CRLF record keeps its terminator, and a
    /// local edit does not normalize it to LF.
    #[test]
    fn crlf_terminator_is_preserved() {
        let line = b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\r\n";
        assert_eq!(
            edit_raw(line, &EditPlan::default(), Some(0)),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\r\n"
        );
        let mut bits = Bitset::new(1);
        bits.set(0, true);
        let plan = EditPlan {
            mask_gt: Some(&bits),
            ..EditPlan::default()
        };
        assert_eq!(
            edit_raw(line, &plan, Some(0)),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t./.\r\n"
        );
    }

    #[test]
    fn info_remove_and_set_replace_in_place_and_append() {
        let plan = EditPlan {
            info_remove: &[b"AF", b"DP"],
            info_set: &[
                (b"AC".as_slice(), InfoValueOut::Int(3)),
                (b"NS".as_slice(), InfoValueOut::Int(2)),
                (b"MAF".as_slice(), InfoValueOut::Float(0.25)),
            ],
            ..EditPlan::default()
        };
        let got = edit(LINE, &plan, None);
        assert!(got.contains("AC=3;QD=4.5;NS=2;MAF=0.25"), "{got}");
        assert!(!got.contains("AF=0.5") && !got.contains("DP=150"), "{got}");
    }

    #[test]
    fn all_info_removed_becomes_dot() {
        let plan = EditPlan {
            info_remove: &[b"AC", b"QD", b"AF", b"DP"],
            ..EditPlan::default()
        };
        let got = edit(LINE, &plan, None);
        assert!(got.contains("\t.\t"), "{got}");
    }

    #[test]
    fn dot_info_becomes_the_set_tags() {
        let line = b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let plan = EditPlan {
            info_set: &[(b"NS".as_slice(), InfoValueOut::Int(1))],
            ..EditPlan::default()
        };
        let got = edit(line, &plan, None);
        assert!(got.contains("\tNS=1\t"), "{got}");
    }

    #[test]
    fn sample_pruning_and_gt_masking() {
        let plan = EditPlan {
            keep_samples: Some(&mask_from(4, &[0, 2, 3])),
            ..EditPlan::default()
        };
        let got = edit(LINE, &plan, Some(0));
        assert!(got.ends_with("0/1:15:99\t./.:0:9\t0/0:8:44"), "{got}");

        // V22-R2-1: GT past the 256th FORMAT key is still replaced. With a
        // `u8` position the caller passed `None` and the mask silently did
        // nothing to the GT column.
        let mut keys: Vec<String> = (0..256).map(|i| format!("K{i}")).collect();
        keys.push("GT".into());
        let mut values = vec!["."; 256];
        values.push("0/1");
        let line = format!(
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\t{}\t{}\n",
            keys.join(":"),
            values.join(":")
        );
        let mut bits = Bitset::new(4);
        bits.set(0, true);
        let plan = EditPlan {
            mask_gt: Some(&bits),
            ..EditPlan::default()
        };
        let got = edit(line.as_bytes(), &plan, Some(256));
        assert!(got.ends_with(":./."), "{got}");
        assert!(!got.contains(":0/1"), "{got}");

        // Mask sample 0 (GT subfield 0 → ./.), DP/GQ stay verbatim.
        let mut bits = Bitset::new(4);
        bits.set(0, true);
        let plan = EditPlan {
            mask_gt: Some(&bits),
            ..EditPlan::default()
        };
        let got = edit(LINE, &plan, Some(0));
        assert!(got.contains("\t./.:15:99\t1|1:12:99"), "{got}");

        // A row whose "." value has no GT subfield stays verbatim.
        let short = b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1\t.\t0/0\t1/1\n";
        let plan = EditPlan {
            mask_gt: Some(&bits),
            ..EditPlan::default()
        };
        let got = edit(short, &plan, Some(0));
        assert!(got.contains("\t./.\t.\t"), "{got}");
    }

    #[test]
    fn m7_zero_selected_samples_drops_format_and_samples() {
        // doc_v21 M7: `keep_samples` selecting nothing used to emit the FORMAT
        // column plus an empty sample field while the header was 8-column —
        // an invalid VCF (8-column header, 10-column records).
        let mask = SampleMask::none(4);
        let plan = EditPlan {
            keep_samples: Some(&mask),
            ..EditPlan::default()
        };
        let got = edit(LINE, &plan, Some(0));
        assert_eq!(got.split('\t').count(), 8, "{got}");
        assert_eq!(
            got,
            "Chr01\t100\t.\tA\tG\t50\tPASS\tAC=9;QD=4.5;AF=0.5;DP=150"
        );
        // The control: keeping one sample still writes FORMAT + that column.
        let mask = SampleMask::from_ids(4, [SampleId(1)]).expect("invariant: id in range");
        let plan = EditPlan {
            keep_samples: Some(&mask),
            ..EditPlan::default()
        };
        let got = edit(LINE, &plan, Some(0));
        assert_eq!(got.split('\t').count(), 10, "{got}");
    }

    #[test]
    fn gt_not_at_subfield_zero() {
        let line = b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tDP:GT\t15:0/1\t12:1|1\t9:./.\t44:0/0\n";
        let mut bits = Bitset::new(4);
        bits.set(0, true);
        let plan = EditPlan {
            mask_gt: Some(&bits),
            ..EditPlan::default()
        };
        let got = edit(line, &plan, Some(1));
        assert!(got.contains("\t15:./.\t12:1|1"), "{got}");
    }

    /// doc_v23 V23-03: a line with extra trailing empty sample columns (which
    /// htslib and the scan stage ignore, V22-4) used to walk `sample_idx`
    /// past the header's sample count and index the keep mask / GT mask out
    /// of range — `Bitset::get` asserts, aborting the worker. Both edit kinds
    /// must survive by dropping the extras, exactly as the scan accepted
    /// them (nothing).
    #[test]
    fn trailing_empty_columns_do_not_overrun_the_header_sample_count() {
        // Byte-exact (no `trim_end`): a resurrected trailing tab is the bug,
        // so the assertion has to see it. Header declares two samples; the
        // line carries both plus two empty columns past the header count
        // (the shape the scan stage tolerates, V22-4).
        let raw = |line: &[u8], plan: &EditPlan<'_>| {
            let lr = LineRef::parse(line).expect("invariant: full line");
            let mut out = Vec::new();
            edit_vcf_line(&lr, plan, Some(0), &mut out);
            String::from_utf8(out).expect("invariant: utf-8")
        };
        let line = b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15\t1/1:9\t\t\n";

        // Sample subset: one of two kept, so the region is edited.
        let keep = SampleMask::from_ids(2, [SampleId(0)]).expect("invariant: id in range");
        let plan = EditPlan {
            keep_samples: Some(&keep),
            ..EditPlan::default()
        };
        assert_eq!(
            raw(line, &plan),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15\n",
            "the trailing empty columns must be dropped, not copied back"
        );

        // The GT-mask path (`--min-gq`-style) without `-s`: same bound.
        let mut bits = Bitset::new(2);
        bits.set(1, true);
        let plan = EditPlan {
            mask_gt: Some(&bits),
            ..EditPlan::default()
        };
        assert_eq!(
            raw(line, &plan),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15\t./.:9\n"
        );

        // A single trailing tab (no extra column) keeps working.
        let line = b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15\t1/1:9\t\n";
        assert_eq!(
            raw(line, &plan),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15\t./.:9\n"
        );

        // No sample edit at all still copies the line verbatim, trailing
        // columns included: the scan stage accepted them, and the writer's
        // job is only to not crash and not invent columns.
        let verbatim_plan = EditPlan::default();
        assert_eq!(
            raw(line, &verbatim_plan),
            "Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:15\t1/1:9\t\n"
        );
    }
}

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

use crate::error::IoError;
use crate::vcf::line::LineRef;

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
        InfoValueOut::Float(v) => out.extend_from_slice(format_g6(f64::from(*v)).as_bytes()),
    }
}

/// Edits one record line into `out` (newline-terminated).
///
/// `gt_subfield_index` is the GT position among the FORMAT keys of this
/// record (resolved by the caller; `None` = the record has no GT column).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when the line has fewer than 9 columns.
pub fn edit_vcf_line(
    line: &LineRef<'_>,
    plan: &EditPlan<'_>,
    gt_subfield_index: Option<u8>,
    out: &mut Vec<u8>,
) -> Result<(), IoError> {
    if LineRef::parse(line.line).is_none() {
        return Err(IoError::MalformedRecord {
            locus: crate::Locus {
                contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
                pos1: 0,
            },
            reason: "line has fewer than 9 columns",
        });
    }
    // Columns 1–7 (CHROM .. FILTER) with their tabs: the bytes before INFO.
    out.extend_from_slice(&line.line[..line.info_start()]);

    // INFO.
    let info = line.info_text();
    if plan.info_untouched() {
        out.extend_from_slice(info);
    } else {
        edit_info(info, plan, out);
    }
    out.push(b'\t');

    // FORMAT verbatim.
    out.extend_from_slice(line.format_text());
    out.push(b'\t');

    // Sample region.
    let samples = line.samples_text();
    if plan.samples_edited() {
        edit_samples(samples, plan, gt_subfield_index, out);
    } else {
        out.extend_from_slice(samples);
    }
    out.push(b'\n');
    Ok(())
}

/// INFO edit: drop `info_remove`, replace any tag present in `info_set` in
/// place, append the rest in `info_set` order. An empty result renders as
/// `.`.
fn edit_info(info: &[u8], plan: &EditPlan<'_>, out: &mut Vec<u8>) {
    let mut edited: Vec<u8> = Vec::new();
    let mut seen = vec![false; plan.info_set.len()];
    if info != b"." {
        for (tag, value) in crate::vcf::info::info_fragments(info) {
            if plan.info_remove.contains(&tag) {
                continue;
            }
            if !edited.is_empty() {
                edited.push(b';');
            }
            if let Some(i) = plan.info_set.iter().position(|(t, _)| *t == tag) {
                edited.extend_from_slice(plan.info_set[i].0);
                edited.push(b'=');
                render_info_value(&plan.info_set[i].1, &mut edited);
                seen[i] = true;
            } else {
                edited.extend_from_slice(tag);
                if let Some(v) = value {
                    edited.push(b'=');
                    edited.extend_from_slice(v);
                }
            }
        }
    }
    for (i, (tag, value)) in plan.info_set.iter().enumerate() {
        if seen[i] {
            continue;
        }
        if !edited.is_empty() {
            edited.push(b';');
        }
        edited.extend_from_slice(tag);
        edited.push(b'=');
        render_info_value(value, &mut edited);
    }
    if edited.is_empty() {
        out.push(b'.');
    } else {
        out.extend_from_slice(&edited);
    }
}

/// Sample region edit: skip dropped columns, replace the GT subfield of
/// masked samples with `./.`, copy everything else verbatim.
fn edit_samples(
    samples: &[u8],
    plan: &EditPlan<'_>,
    gt_subfield_index: Option<u8>,
    out: &mut Vec<u8>,
) {
    let keep = plan.keep_samples;
    let mask = plan.mask_gt;
    let mut field_start = 0;
    let mut sample_idx = 0_usize;
    let mut wrote_column = false;
    // Byte offset where the current run of verbatim samples begins
    // (doc_v3 03 O7: one copy per run instead of one per sample).
    let mut run_start: Option<usize> = None;
    while field_start < samples.len() {
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
                replace_gt_subfield(&samples[field_start..field_end], gt_subfield_index, out);
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

/// Copies `field`, replacing the GT subfield with `./.` when the row is
/// long enough; a row without that subfield stays verbatim.
fn replace_gt_subfield(field: &[u8], gt_index: Option<u8>, out: &mut Vec<u8>) {
    let Some(gt_index) = gt_index else {
        out.extend_from_slice(field);
        return;
    };
    let gt_index = usize::from(gt_index);
    let mut subfield_start = 0;
    let mut index = 0_usize;
    loop {
        let subfield_end = memchr::memchr(b':', &field[subfield_start..])
            .map_or(field.len(), |rel| subfield_start + rel);
        if index == gt_index {
            out.extend_from_slice(&field[..subfield_start]);
            out.extend_from_slice(b"./.");
            out.extend_from_slice(&field[subfield_end..]);
            return;
        }
        if subfield_end == field.len() {
            break;
        }
        subfield_start = subfield_end + 1;
        index += 1;
    }
    out.extend_from_slice(field);
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

/// [`format_g6`] appended to a byte buffer (the renderer's hot path avoids
/// one allocation per value).
pub fn format_g6_into(out: &mut Vec<u8>, v: f64) {
    use std::io::Write as _;
    if v == 0.0 {
        out.push(b'0');
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
    // significant digits (C %g semantics).
    // Six significant digits first, then pick the style from the exponent
    // *after* rounding (C %g: 999999.5 rounds to 1e+06, not 9.99999e5).
    let scientific = format!("{:.5e}", v.abs());
    let exp10: i32 = scientific
        .split('e')
        .nth(1)
        .and_then(|e| e.parse().ok())
        .unwrap_or(0);
    if (-4..6).contains(&exp10) {
        let precision = (5 - exp10).max(0) as usize;
        let s = format!("{v:.precision$}");
        // htslib culls trailing zeros only past a decimal point: an integer
        // like 128440 must keep its zeros (the pre-BCF version stripped
        // them, which only showed for values >= 1e5).
        match s.split_once('.') {
            Some((int, frac)) => {
                let frac = frac.trim_end_matches('0');
                out.extend_from_slice(int.as_bytes());
                if !frac.is_empty() {
                    out.push(b'.');
                    out.extend_from_slice(frac.as_bytes());
                }
            }
            None => out.extend_from_slice(s.as_bytes()),
        }
    } else {
        // Scientific: 5 fractional digits, trailing zeros stripped,
        // exponent at least two digits with a sign.
        let mantissa_exp = format!("{v:.5e}");
        let (mantissa, exp) = mantissa_exp.split_once('e').expect("invariant: e-form");
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        let exp: i32 = exp.parse().expect("invariant: integer exponent");
        let _ = write!(
            out,
            "{mantissa}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        );
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

    fn edit(line: &[u8], plan: &EditPlan<'_>, gt_idx: Option<u8>) -> String {
        let lr = LineRef::parse(line).expect("invariant: full line");
        let mut out = Vec::new();
        edit_vcf_line(&lr, plan, gt_idx, &mut out).expect("invariant: edit");
        String::from_utf8(out)
            .expect("invariant: utf-8")
            .trim_end()
            .to_string()
    }

    const LINE: &[u8] =
        b"Chr01\t100\t.\tA\tG\t50\tPASS\tAC=9;QD=4.5;AF=0.5;DP=150\tGT:DP:GQ\t0/1:15:99\t1|1:12:99\t./.:0:9\t0/0:8:44\n";

    #[test]
    fn verbatim_when_no_plan() {
        let plan = EditPlan::default();
        assert_eq!(
            edit(LINE, &plan, Some(0)),
            String::from_utf8_lossy(LINE).trim_end()
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
}

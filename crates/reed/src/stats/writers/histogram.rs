//! `--hist` (doc_v14 02 §3): histogram requests, their mergeable
//! accumulators and the `PREFIX.<m>.hist` TSV renderer.
//!
//! Built-in metrics reuse the QC summary's accumulators (`--hist qual,dp,…`
//! shares the pass with `--summary`, nothing is parsed twice); `gq` and the
//! generic `INFO/TAG`/`FMT/TAG` requests decode their field through the
//! shared scanner. Integer metrics carry exact mean/median/p95
//! (`IntHistogram::quantile`); float metrics use fixed boundaries (B14-2,
//! C14-5).

use std::fmt::Write as _;

use reed_stats::summary::histogram::{IntHistogram, QualDist};
use reed_stats::summary::report::QcReportAcc;


/// The declared type of a generic `--hist` tag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HistKind {
    /// `Type=Integer`.
    Int,
    /// `Type=Float`.
    Float,
}

/// One `--hist` element.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum HistMetric {
    /// QUAL (0.1 bins, SNV and indel separately).
    Qual,
    /// FMT/DP, one observation per genotype.
    Dp,
    /// INFO/DP, one observation per site.
    DpSite,
    /// FMT/GQ, one observation per genotype.
    Gq,
    /// Indel length (insertions positive, deletions negative).
    Idl,
    /// INFO/AC ÷ INFO/AN per ALT (the bcftools AF bins).
    Af,
    /// A generic integer or float tag: `INFO/<TAG>` or `FMT/<TAG>`.
    Tag {
        /// `true` for `FMT/`, `false` for `INFO/`.
        format: bool,
        /// The tag name as spelled.
        tag: String,
        /// The header's declared value kind (the CLI resolves it before the
        /// first record, cli.md §0.7).
        kind: HistKind,
    },
}

impl HistMetric {
    /// Parses one `--hist` element.
    ///
    /// # Errors
    ///
    /// A usage message for an unknown built-in name or a malformed
    /// `INFO/`/`FMT/` prefix.
    pub fn parse(text: &str) -> Result<Self, String> {
        let t = text.trim();
        match t {
            "qual" => Ok(Self::Qual),
            "dp" => Ok(Self::Dp),
            "dp-site" => Ok(Self::DpSite),
            "gq" => Ok(Self::Gq),
            "idl" => Ok(Self::Idl),
            "af" => Ok(Self::Af),
            _ => {
                if let Some(tag) = t.strip_prefix("INFO/") {
                    if tag.is_empty() {
                        return Err(format!("--hist {t:?}: empty tag name"));
                    }
                    Ok(Self::Tag {
                        format: false,
                        tag: tag.to_string(),
                        // The CLI overwrites this after the header check;
                        // a library caller that wants a float tag passes
                        // `HistKind::Float` explicitly.
                        kind: HistKind::Int,
                    })
                } else if let Some(tag) = t.strip_prefix("FMT/") {
                    if tag.is_empty() {
                        return Err(format!("--hist {t:?}: empty tag name"));
                    }
                    Ok(Self::Tag {
                        format: true,
                        tag: tag.to_string(),
                        kind: HistKind::Int,
                    })
                } else {
                    Err(format!(
                        "--hist {t:?}: unknown metric (built-ins: qual, dp, dp-site, gq, idl, \
                         af; tags: INFO/<TAG>, FMT/<TAG>)"
                    ))
                }
            }
        }
    }

    /// The output-file infix: `PREFIX.<infix>.hist`.
    #[must_use]
    pub fn infix(&self) -> String {
        match self {
            Self::Qual => "qual".into(),
            Self::Dp => "dp".into(),
            Self::DpSite => "dp-site".into(),
            Self::Gq => "gq".into(),
            Self::Idl => "idl".into(),
            Self::Af => "af".into(),
            Self::Tag { format, tag, .. } => {
                let scope = if *format { "FMT" } else { "INFO" };
                format!("{scope}_{tag}")
            }
        }
    }

    /// Whether the metric reads FMT/GQ (the scanner must decode it).
    #[must_use]
    pub const fn needs_gq(&self) -> bool {
        matches!(self, Self::Gq)
    }

    /// Whether the metric reads FMT/DP (the scanner must decode the column).
    #[must_use]
    pub const fn needs_dp(&self) -> bool {
        matches!(self, Self::Dp)
    }

    /// Whether the metric reads the QC summary's accumulators (the reports
    /// job runs the QC accumulator in the same pass; doc_v4 O-1).
    #[must_use]
    pub const fn needs_qc(&self) -> bool {
        matches!(self, Self::Qual | Self::Dp | Self::DpSite | Self::Idl | Self::Af)
    }

    /// Whether any metric in `metrics` needs the QC accumulators.
    #[must_use]
    pub fn any_needs_qc(metrics: &[Self]) -> bool {
        metrics.iter().any(Self::needs_qc)
    }

    /// The generic tag this metric reads, if any: `(format, tag, kind)`.
    #[must_use]
    pub fn tag(&self) -> Option<(bool, &str, HistKind)> {
        match self {
            Self::Tag { format, tag, kind } => Some((*format, tag, *kind)),
            _ => None,
        }
    }
}

/// The worker-side accumulators for one `--hist` request. Built-in metrics
/// read the QC summary accumulator (the same pass, doc_v4 O-1); `gq` and the
/// generic tags accumulate here, index-aligned with the request.
pub(crate) struct HistAcc {
    /// `--hist gq`: FMT/GQ, one observation per genotype.
    pub(crate) gq: IntHistogram,
    /// Per request entry: integer observations (`Some` for tag metrics).
    pub(crate) int: Vec<Option<IntHistogram>>,
    /// Per request entry: float observations (`Some` for float tag metrics).
    pub(crate) float: Vec<Option<QualDist>>,
    /// Per request entry: observations that could not be read.
    pub(crate) missing: Vec<u64>,
}

impl HistAcc {
    /// Allocates one slot per requested metric.
    #[must_use]
    pub(crate) fn new(metrics: &[HistMetric]) -> Self {
        let mut acc = Self {
            gq: IntHistogram::new(1_000),
            int: Vec::with_capacity(metrics.len()),
            float: Vec::with_capacity(metrics.len()),
            missing: vec![0; metrics.len()],
        };
        for metric in metrics {
            match metric.tag() {
                Some((_, _, HistKind::Int)) => {
                    acc.int.push(Some(IntHistogram::new(10_000)));
                    acc.float.push(None);
                }
                Some((_, _, HistKind::Float)) => {
                    acc.int.push(None);
                    acc.float.push(Some(QualDist::new()));
                }
                None => {
                    acc.int.push(None);
                    acc.float.push(None);
                }
            }
        }
        acc
    }

    /// Merges another worker's accumulators (chunk partials are independent).
    pub(crate) fn merge(&mut self, other: &Self) {
        self.gq.merge(&other.gq);
        for (dst, src) in self.int.iter_mut().zip(&other.int) {
            if let (Some(dst), Some(src)) = (dst, src) {
                dst.merge(src);
            }
        }
        for (dst, src) in self.float.iter_mut().zip(&other.float) {
            if let (Some(dst), Some(src)) = (dst, src) {
                dst.merge(src);
            }
        }
        for (dst, src) in self.missing.iter_mut().zip(&other.missing) {
            *dst += src;
        }
    }
}

/// One metric's merged output, ready to render.
#[derive(Clone, Debug)]
pub struct HistOut {
    /// The metric (request order).
    pub metric: HistMetric,
    /// Integer bins `(value, count)` when the metric is integer-valued.
    /// Signed: `idl` follows the bcftools IDD convention (deletions
    /// negative, insertions positive).
    pub int_bins: Vec<(i64, u64)>,
    /// Exact mean/median/p95 for integer metrics.
    pub stats: Option<HistStats>,
    /// Float bins `(lower, upper, count)` (0.1 resolution to 100).
    pub float_bins: Vec<(f64, f64, u64)>,
    /// Observations that could not be read.
    pub n_missing: u64,
}

/// Exact summary statistics of an integer histogram.
#[derive(Clone, Copy, Debug)]
pub struct HistStats {
    /// Arithmetic mean.
    pub mean: f64,
    /// Median (the exact quantile).
    pub median: u32,
    /// The 95th percentile (the exact quantile).
    pub p95: u32,
}

/// Merges the per-chunk histogram partials into the request's outputs.
pub(crate) fn merge_hist(
    metrics: &[HistMetric],
    qc: Option<&QcReportAcc>,
    acc: &HistAcc,
) -> Vec<HistOut> {
    let mut out = Vec::with_capacity(metrics.len());
    for (tag_idx, metric) in metrics.iter().enumerate() {
        let mut entry = HistOut {
            metric: metric.clone(),
            int_bins: Vec::new(),
            stats: None,
            float_bins: Vec::new(),
            n_missing: acc.missing.get(tag_idx).copied().unwrap_or(0),
        };
        match metric {
            HistMetric::Qual => {
                if let Some(qc) = qc {
                    let mut bins: Vec<(i64, u64)> = Vec::new();
                    let nbins = qc.qual.snv.nbins().max(qc.qual.indel.nbins());
                    for i in 0..nbins {
                        let (lo, _hi) = QualDist::bin_range(i);
                        let n = qc.qual.snv.count_at(i) + qc.qual.indel.count_at(i);
                        bins.push((i64::from(lo), n));
                    }
                    entry.int_bins = bins;
                }
            }
            HistMetric::Dp => {
                if let Some(qc) = qc {
                    entry.int_bins = qc.dp_genotypes.iter().map(|(v, n)| (i64::from(v), n)).collect();
                    entry.stats = Some(stats_of(&qc.dp_genotypes));
                }
            }
            HistMetric::DpSite => {
                if let Some(qc) = qc {
                    entry.int_bins = qc.dp_sites.iter().map(|(v, n)| (i64::from(v), n)).collect();
                    entry.stats = Some(stats_of(&qc.dp_sites));
                }
            }
            HistMetric::Idl => {
                if let Some(qc) = qc {
                    // bcftools IDD convention: a deletion is `-len`, an
                    // insertion `+len`, one signed length per row
                    // (vcfstats.c's IDD section).
                    let mut bins: Vec<(i64, u64)> = Vec::new();
                    for (len, n) in qc.idd.iter() {
                        bump_signed(&mut bins, -i64::from(len), n);
                    }
                    for (len, n) in qc.idd_ins.iter() {
                        bump_signed(&mut bins, i64::from(len), n);
                    }
                    bins.sort_unstable_by_key(|(v, _)| *v);
                    entry.int_bins = bins;
                }
            }
            HistMetric::Af => {
                if let Some(qc) = qc {
                    let m = qc.af.m_af;
                    for (bin, &snps) in qc.af.snps.iter().enumerate() {
                        let lower = bin as f64 / (m.saturating_sub(2)).max(1) as f64;
                        let upper = (bin + 1) as f64 / (m.saturating_sub(2)).max(1) as f64;
                        entry.float_bins.push((lower, upper, snps));
                    }
                }
            }
            HistMetric::Gq => {
                let h = &acc.gq;
                entry.int_bins = h.iter().map(|(v, n)| (i64::from(v), n)).collect();
                entry.stats = Some(stats_of(h));
            }
            HistMetric::Tag { .. } => {
                if let Some(h) = acc.int.get(tag_idx).and_then(Option::as_ref) {
                    entry.int_bins = h.iter().map(|(v, n)| (i64::from(v), n)).collect();
                    entry.stats = Some(stats_of(h));
                } else if let Some(f) = acc.float.get(tag_idx).and_then(Option::as_ref) {
                    entry.float_bins = float_bins_of(f);
                }
            }
        }
        out.push(entry);
    }
    out
}

/// Adds `n` to the signed bin `value` (the IDD histogram's merge).
fn bump_signed(bins: &mut Vec<(i64, u64)>, value: i64, n: u64) {
    if let Some(slot) = bins.iter_mut().find(|(v, _)| *v == value) {
        slot.1 += n;
    } else {
        bins.push((value, n));
    }
}

fn float_bins_of(dist: &QualDist) -> Vec<(f64, f64, u64)> {
    let mut bins = Vec::new();
    for i in 0..dist.nbins() {
        let (lo, hi) = QualDist::bin_range(i);
        bins.push((f64::from(lo) / 10.0, f64::from(hi) / 10.0, dist.count_at(i)));
    }
    bins
}

fn stats_of(h: &IntHistogram) -> HistStats {
    let total = h.total();
    let mean = if total > 0 {
        h.weighted_sum() as f64 / total as f64
    } else {
        0.0
    };
    HistStats {
        mean,
        median: h.quantile(0.5).unwrap_or(0),
        p95: h.quantile(0.95).unwrap_or(0),
    }
}

/// Renders one metric's TSV (`bin_lower`, `bin_upper`, `count`, `percent`),
/// integer metrics first, then the exact `#mean`/`#median`/`#p95` lines and
/// `#missing` (doc_v14 02 §3.2).
#[must_use]
pub fn render_hist(out: &HistOut) -> String {
    let mut text = String::new();
    text.push_str("bin_lower\tbin_upper\tcount\tpercent\n");
    if out.int_bins.is_empty() && out.float_bins.is_empty() {
        return text;
    }
    if out.int_bins.is_empty() {
        let total: u64 = out.float_bins.iter().map(|(_, _, n)| n).sum();
        for &(lo, hi, n) in &out.float_bins {
            let pct = if total > 0 {
                n as f64 * 100.0 / total as f64
            } else {
                0.0
            };
            let _ = writeln!(text, "{lo:.1}\t{hi:.1}\t{n}\t{pct:.6}");
        }
    } else {
        let total: u64 = out.int_bins.iter().map(|(_, n)| n).sum();
        for &(value, n) in &out.int_bins {
            let pct = if total > 0 {
                n as f64 * 100.0 / total as f64
            } else {
                0.0
            };
            let _ = writeln!(text, "{value}\t{value}\t{n}\t{pct:.6}");
        }
        if let Some(stats) = &out.stats {
            let _ = writeln!(text, "#mean\t{:.6}", stats.mean);
            let _ = writeln!(text, "#median\t{}", stats.median);
            let _ = writeln!(text, "#p95\t{}", stats.p95);
        }
    }
    if out.n_missing > 0 {
        let _ = writeln!(text, "#missing\t{}", out.n_missing);
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metrics_parse_and_name() {
        assert_eq!(HistMetric::parse("qual").expect("builtin"), HistMetric::Qual);
        assert_eq!(HistMetric::parse("dp-site").expect("builtin"), HistMetric::DpSite);
        assert_eq!(
            HistMetric::parse("INFO/DP").expect("info tag"),
            HistMetric::Tag {
                format: false,
                tag: "DP".into(),
                kind: HistKind::Int,
            }
        );
        assert_eq!(
            HistMetric::parse("FMT/GQ").expect("fmt tag").infix(),
            "FMT_GQ"
        );
        assert!(HistMetric::parse("nope").is_err());
        assert!(HistMetric::parse("INFO/").is_err());
    }

    #[test]
    fn integer_histogram_renders_exact_quantiles() {
        let mut h = IntHistogram::new(10);
        for v in [1_u32, 2, 2, 3, 3, 3, 4] {
            h.push(v);
        }
        let out = HistOut {
            metric: HistMetric::Dp,
            int_bins: h.iter().map(|(v, n)| (i64::from(v), n)).collect(),
            stats: Some(stats_of(&h)),
            float_bins: Vec::new(),
            n_missing: 2,
        };
        let text = render_hist(&out);
        assert!(text.contains("2\t2\t2\t28.571429"), "{text}");
        assert!(text.contains("#median\t3"), "{text}");
        assert!(text.contains("#missing\t2"), "{text}");
    }
}

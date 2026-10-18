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
    /// Integer bins `(value, count)` when the metric is integer-valued,
    /// occupied values only. Signed: `idl` follows the bcftools IDD
    /// convention (deletions negative, insertions positive).
    pub int_bins: Vec<(i64, u64)>,
    /// The accumulator's overflow bucket `(max + 1, count)`: observations
    /// above the histogram's range. It is printed as `>max` and kept out of
    /// the exact statistics (its values are unknown).
    pub overflow: Option<(i64, u64)>,
    /// Exact mean/median/p95 over the exact (`int_bins`) observations.
    pub stats: Option<HistStats>,
    /// Float bins `(lower, upper, count)`, occupied only.
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
            overflow: None,
            stats: None,
            float_bins: Vec::new(),
            n_missing: acc.missing.get(tag_idx).copied().unwrap_or(0),
        };
        match metric {
            HistMetric::Qual => {
                if let Some(qc) = qc {
                    // bcftools keys QUAL on `iqual = 1 + round(qual·10)`
                    // (dist.c): a bin index maps to the quality range
                    // `[(i−1)/10, i/10)`. Bin 0 is missing/negative/NaN and
                    // is counted as such, not shown as a `0.0` observation.
                    let mut bins: Vec<(usize, u64)> = Vec::new();
                    for (index, n) in qc.qual.snv.occupied() {
                        bump_index(&mut bins, index, n);
                    }
                    for (index, n) in qc.qual.indel.occupied() {
                        bump_index(&mut bins, index, n);
                    }
                    bins.sort_unstable_by_key(|(index, _)| *index);
                    for (index, n) in bins {
                        if index == 0 {
                            entry.n_missing += n;
                            continue;
                        }
                        let (lo, hi) = QualDist::bin_range(index);
                        entry
                            .float_bins
                            .push((qual_of(lo), qual_of(hi), n));
                    }
                }
            }
            HistMetric::Dp => {
                if let Some(qc) = qc {
                    (entry.int_bins, entry.overflow) = split_int(&qc.dp_genotypes);
                    entry.stats = stats_of_bins(&entry.int_bins);
                }
            }
            HistMetric::DpSite => {
                if let Some(qc) = qc {
                    (entry.int_bins, entry.overflow) = split_int(&qc.dp_sites);
                    entry.stats = stats_of_bins(&entry.int_bins);
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
                    // bcftools's AF section (vcfstats.c): the singletons bin
                    // (index 0) folds into index 1 and the printed frequency
                    // for index i is `(i−1)/(m_af−1)`; the indel column
                    // counts too (the `.hist` has one count column).
                    let m = qc.af.m_af.max(2);
                    let denom = (m - 1) as f64;
                    for i in 1..m {
                        let mut n = qc.af.snps[i] + qc.af.indels[i];
                        if i == 1 {
                            n += qc.af.snps[0] + qc.af.indels[0];
                        }
                        if n == 0 {
                            continue;
                        }
                        entry.float_bins.push((
                            (i as f64 - 1.0) / denom,
                            i as f64 / denom,
                            n,
                        ));
                    }
                }
            }
            HistMetric::Gq => {
                (entry.int_bins, entry.overflow) = split_int(&acc.gq);
                entry.stats = stats_of_bins(&entry.int_bins);
            }
            HistMetric::Tag { .. } => {
                if let Some(h) = acc.int.get(tag_idx).and_then(Option::as_ref) {
                    (entry.int_bins, entry.overflow) = split_int(h);
                    entry.stats = stats_of_bins(&entry.int_bins);
                } else if let Some(f) = acc.float.get(tag_idx).and_then(Option::as_ref) {
                    entry.float_bins = float_bins_of(f);
                }
            }
        }
        out.push(entry);
    }
    out
}

/// The quality of an `iqual` value: `iqual = 1 + qual·10`, so
/// `qual = (iqual − 1)/10`.
fn qual_of(iqual: u32) -> f64 {
    (f64::from(iqual) - 1.0) / 10.0
}

/// Adds `n` to the count at `index` (occupied-bin merges for float dists).
fn bump_index(bins: &mut Vec<(usize, u64)>, index: usize, n: u64) {
    if let Some(slot) = bins.iter_mut().find(|(i, _)| *i == index) {
        slot.1 += n;
    } else {
        bins.push((index, n));
    }
}

/// An integer histogram's occupied exact rows plus its overflow bucket.
type IntSplit = (Vec<(i64, u64)>, Option<(i64, u64)>);

/// Splits an integer accumulator into its occupied exact values and, when
/// non-empty, its overflow bucket (`max + 1`) — the bucket holds values
/// above the histogram's range, so it must not pass as an exact value.
fn split_int(h: &IntHistogram) -> IntSplit {
    let all: Vec<(u32, u64)> = h.iter().collect();
    let last = all.len().saturating_sub(1);
    let mut bins = Vec::new();
    let mut overflow = None;
    for (i, (value, n)) in all.into_iter().enumerate() {
        if n == 0 {
            continue;
        }
        if i == last {
            overflow = Some((i64::from(value), n));
        } else {
            bins.push((i64::from(value), n));
        }
    }
    (bins, overflow)
}

/// Adds `n` to the signed bin `value` (the IDD histogram's merge).
fn bump_signed(bins: &mut Vec<(i64, u64)>, value: i64, n: u64) {
    if let Some(slot) = bins.iter_mut().find(|(v, _)| *v == value) {
        slot.1 += n;
    } else {
        bins.push((value, n));
    }
}

/// The occupied float bins of a distribution, in quality units.
fn float_bins_of(dist: &QualDist) -> Vec<(f64, f64, u64)> {
    let mut bins: Vec<(f64, f64, u64)> = dist
        .occupied()
        .map(|(index, n)| {
            let (lo, hi) = QualDist::bin_range(index);
            (qual_of(lo), qual_of(hi), n)
        })
        .collect();
    bins.sort_by(|a, b| a.0.total_cmp(&b.0));
    bins
}

/// Exact mean/median/p95 over `(value, count)` pairs (the overflow bucket is
/// not one of them: its values are unknown).
fn stats_of_bins(bins: &[(i64, u64)]) -> Option<HistStats> {
    let total: u64 = bins.iter().map(|(_, n)| n).sum();
    if total == 0 {
        return None;
    }
    let weighted: f64 = bins
        .iter()
        .map(|(value, n)| *value as f64 * *n as f64)
        .sum();
    // The quantile rank is `ceil(num·total/den)` (integer arithmetic, so
    // the rank is exact), counted from the smallest value up.
    let quantile = |num: u64, den: u64| -> u32 {
        let target = (num * total).div_ceil(den).clamp(1, total);
        let mut acc = 0_u64;
        for (value, n) in bins {
            acc += n;
            if acc >= target {
                return u32::try_from((*value).max(0)).unwrap_or(0);
            }
        }
        0
    };
    Some(HistStats {
        mean: weighted / total as f64,
        median: quantile(1, 2),
        p95: quantile(95, 100),
    })
}

/// Renders one metric's TSV (`bin_lower`, `bin_upper`, `count`, `percent`),
/// occupied bins only; integer metrics then carry the exact
/// `#mean`/`#median`/`#p95` lines and `#missing` (doc_v14 02 §3.2).
#[must_use]
pub fn render_hist(out: &HistOut) -> String {
    let mut text = String::new();
    text.push_str("bin_lower\tbin_upper\tcount\tpercent\n");
    if out.int_bins.is_empty() && out.float_bins.is_empty() && out.overflow.is_none() {
        if out.n_missing > 0 {
            let _ = writeln!(text, "#missing\t{}", out.n_missing);
        }
        return text;
    }
    if !out.int_bins.is_empty() || out.overflow.is_some() {
        let total: u64 = out
            .int_bins
            .iter()
            .map(|(_, n)| n)
            .chain(out.overflow.iter().map(|(_, n)| n))
            .sum();
        let pct = |n: u64| {
            if total > 0 {
                n as f64 * 100.0 / total as f64
            } else {
                0.0
            }
        };
        for &(value, n) in &out.int_bins {
            let _ = writeln!(text, "{value}\t{value}\t{n}\t{:.6}", pct(n));
        }
        if let Some((value, n)) = out.overflow {
            // Above the accumulator's range: `>max`, not an exact value.
            let max = value - 1;
            let _ = writeln!(text, ">{max}\t>{max}\t{n}\t{:.6}", pct(n));
        }
        if let Some(stats) = &out.stats {
            let _ = writeln!(text, "#mean\t{:.6}", stats.mean);
            let _ = writeln!(text, "#median\t{}", stats.median);
            let _ = writeln!(text, "#p95\t{}", stats.p95);
        }
    } else {
        let total: u64 = out.float_bins.iter().map(|(_, _, n)| n).sum();
        for &(lo, hi, n) in &out.float_bins {
            let pct = if total > 0 {
                n as f64 * 100.0 / total as f64
            } else {
                0.0
            };
            let _ = writeln!(text, "{lo:.6}\t{hi:.6}\t{n}\t{pct:.6}");
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
        let (bins, overflow) = split_int(&h);
        let stats = stats_of_bins(&bins);
        let out = HistOut {
            metric: HistMetric::Dp,
            int_bins: bins,
            overflow,
            stats,
            float_bins: Vec::new(),
            n_missing: 2,
        };
        let text = render_hist(&out);
        assert!(text.contains("2\t2\t2\t28.571429"), "{text}");
        assert!(text.contains("#median\t3"), "{text}");
        assert!(text.contains("#missing\t2"), "{text}");
    }

    /// The overflow bucket is labelled `>max` and kept out of the exact
    /// statistics: two DP values of 600 and 1000 must not read as an exact
    /// `501` with mean 501.
    #[test]
    fn overflow_is_labelled_and_excluded_from_the_stats() {
        let mut h = IntHistogram::new(500);
        h.push(600);
        h.push(1000);
        let (bins, overflow) = split_int(&h);
        let out = HistOut {
            metric: HistMetric::Dp,
            int_bins: bins,
            overflow,
            stats: stats_of_bins(&[]),
            float_bins: Vec::new(),
            n_missing: 0,
        };
        let text = render_hist(&out);
        assert!(text.contains(">500\t>500\t2\t100.000000"), "{text}");
        assert!(!text.contains("#mean"), "no exact stats without exact bins: {text}");
    }

    /// A capped accumulator with in-range data still reports exact stats.
    #[test]
    fn in_range_values_keep_exact_stats() {
        let mut h = IntHistogram::new(500);
        for v in [10_u32, 20, 30, 40] {
            h.push(v);
        }
        let (bins, overflow) = split_int(&h);
        assert!(overflow.is_none());
        let stats = stats_of_bins(&bins).expect("in-range stats");
        assert!((stats.mean - 25.0).abs() < 1e-9);
        assert_eq!(stats.median, 20);
    }

    /// `floats` are labelled in the tag's own units: `iqual = 1 + v·10`, so
    /// `v = (iqual − 1)/10` (the old code divided the iqual index by 10 and
    /// put every value one 0.1 bin high).
    #[test]
    fn float_bins_use_the_value_units() {
        let mut dist = QualDist::new();
        dist.push(Some(1.0));
        dist.push(Some(0.25));
        let bins = float_bins_of(&dist);
        assert_eq!(bins.len(), 2);
        assert_eq!(bins[0], (0.2, 0.3, 1));
        assert_eq!(bins[1], (1.0, 1.1, 1));
    }

    /// The QUAL metric counts iqual 0 (missing/negative/NaN) as missing
    /// rather than as a real `0.0` observation, and labels bins in QUAL
    /// units.
    #[test]
    fn qual_counts_missing_separately() {
        let mut qc = QcReportAcc::new(1);
        qc.qual.snv.push(Some(50.0));
        qc.qual.snv.push(None);
        let acc = HistAcc::new(&[HistMetric::Qual]);
        let out = merge_hist(&[HistMetric::Qual], Some(&qc), &acc);
        let entry = &out[0];
        assert_eq!(entry.float_bins, vec![(50.0, 50.1, 1)]);
        assert_eq!(entry.n_missing, 1);
    }
}

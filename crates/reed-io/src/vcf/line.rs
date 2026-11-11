//! [`LineRef`]: a borrowed VCF record line with the byte offsets of its
//! fixed columns (`docs/spec/architecture.md` §6.7). No allocation, no
//! parsing — offsets only; the field parsers live in [`super::fixed`].

/// Column shape of a record line (doc_v22 R22-3).
///
/// The three shapes are the only legal column layouts; every consumer asks
/// this instead of counting tabs itself, so sites-only input, a FORMAT column
/// without samples, and ordinary sample columns cannot drift apart.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// Eight columns: `CHROM`..`INFO`, no FORMAT or sample region
    /// (bcftools `view -G` output).
    SitesOnly,
    /// Nine columns: a FORMAT column but no sample columns.
    FormatOnly,
    /// Ten or more columns: FORMAT plus at least one sample column.
    WithSamples,
}

/// Why a record's sample columns disagree with the header (doc_v22 R22-3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SampleColumnsError {
    /// Fewer sample columns than the header declares.
    TooFew {
        /// Sample columns present.
        found: usize,
        /// Sample columns the header declares.
        expected: usize,
    },
    /// A non-empty sample column beyond the header's count. Empty columns
    /// after the declared count are trailing tabs and are ignored, matching
    /// htslib.
    TooMany {
        /// One-based position of the first non-empty extra column.
        found: usize,
        /// Sample columns the header declares.
        expected: usize,
    },
}

/// Byte-offset view of one complete VCF record line (terminating `\n`
/// included in `line` but excluded from every field slice).
#[derive(Clone, Copy, Debug)]
pub struct LineRef<'a> {
    /// The full line, `\n` included.
    pub line: &'a [u8],
    /// Positions of the first nine tabs (after CHROM, POS, ID, REF, ALT,
    /// QUAL, FILTER, INFO, FORMAT).
    tabs: [usize; 9],
    /// Tabs actually present: 7 = sites-only (8 columns), 8 = FORMAT but no
    /// sample columns, 9 = FORMAT plus at least one sample column.
    n_tabs: u8,
}

impl<'a> LineRef<'a> {
    /// Parses the tab offsets; `None` when the line has fewer than eight
    /// columns (a truncated line). A sites-only record (`##`-less header
    /// without sample columns, bcftools `view -G`) parses with an empty
    /// [`Self::format_text`] and no sample columns.
    #[must_use]
    pub fn parse(line: &'a [u8]) -> Option<Self> {
        let mut tabs = [0_usize; 9];
        let mut used = 0;
        let mut from = 0;
        while used < 9 {
            let Some(rel) = memchr::memchr(b'\t', &line[from..]) else {
                break;
            };
            from += rel;
            tabs[used] = from;
            from += 1;
            used += 1;
        }
        if used < 7 {
            return None;
        }
        Some(Self {
            line,
            tabs,
            n_tabs: u8::try_from(used).unwrap_or(9),
        })
    }

    fn field(&self, index: usize) -> &'a [u8] {
        let n_tabs = usize::from(self.n_tabs);
        if index > n_tabs {
            return &[];
        }
        let start = if index == 0 {
            0
        } else {
            self.tabs[index - 1] + 1
        };
        // The last present column has no terminating tab; the array's
        // zero-filled tail must not be read as an offset.
        let end = if index < n_tabs {
            self.tabs[index]
        } else {
            self.line.len()
        };
        self.trim_line_end(&self.line[start..end])
    }

    fn trim_line_end(&self, field: &'a [u8]) -> &'a [u8] {
        let no_newline = field.strip_suffix(b"\n").unwrap_or(field);
        no_newline.strip_suffix(b"\r").unwrap_or(no_newline)
    }

    /// The original line terminator (`\n`, `\r\n`, or empty).
    ///
    /// Editors copy it back so an untouched CRLF record is not silently
    /// rewritten to LF (doc_v28 V28-11).
    #[must_use]
    pub fn terminator(&self) -> &'a [u8] {
        let without_lf = self.line.strip_suffix(b"\n").unwrap_or(self.line);
        if without_lf.len() != self.line.len() {
            if without_lf.ends_with(b"\r") {
                &self.line[self.line.len() - 2..]
            } else {
                &self.line[self.line.len() - 1..]
            }
        } else if self.line.ends_with(b"\r") {
            &self.line[self.line.len() - 1..]
        } else {
            &[]
        }
    }

    /// The line without its terminator.
    #[must_use]
    pub fn body(&self) -> &'a [u8] {
        let end = self.line.len() - self.terminator().len();
        &self.line[..end]
    }

    /// CHROM column.
    #[must_use]
    pub fn chrom(&self) -> &'a [u8] {
        self.field(0)
    }

    /// POS column (1-based text).
    #[must_use]
    pub fn pos_text(&self) -> &'a [u8] {
        self.field(1)
    }

    /// ID column.
    #[must_use]
    pub fn id(&self) -> &'a [u8] {
        self.field(2)
    }

    /// REF column.
    #[must_use]
    pub fn reference(&self) -> &'a [u8] {
        self.field(3)
    }

    /// ALT column (comma-separated alleles, unparsed).
    #[must_use]
    pub fn alts(&self) -> &'a [u8] {
        self.field(4)
    }

    /// QUAL column (text).
    #[must_use]
    pub fn qual_text(&self) -> &'a [u8] {
        self.field(5)
    }

    /// FILTER column (text).
    #[must_use]
    pub fn filter_text(&self) -> &'a [u8] {
        self.field(6)
    }

    /// INFO column (text).
    #[must_use]
    pub fn info_text(&self) -> &'a [u8] {
        self.field(7)
    }

    /// FORMAT column (text); empty for a sites-only record.
    #[must_use]
    pub fn format_text(&self) -> &'a [u8] {
        if self.n_tabs < 8 { &[] } else { self.field(8) }
    }

    /// The sample columns region (all text after the FORMAT tab); empty when
    /// the record has no FORMAT/sample columns.
    #[must_use]
    pub fn samples_text(&self) -> &'a [u8] {
        if self.n_tabs < 9 {
            return &[];
        }
        self.trim_line_end(&self.line[self.tabs[8] + 1..])
    }

    /// Number of sample columns (tabs after FORMAT + 1); 0 for a sites-only
    /// record or a FORMAT-only record.
    #[must_use]
    pub fn n_sample_cols(&self) -> usize {
        if self.n_tabs < 9 {
            return 0;
        }
        memchr::memchr_iter(b'\t', self.samples_text()).count() + 1
    }

    /// Byte offset of the first sample column within `line`; `line.len()`
    /// when the record has none.
    #[must_use]
    pub fn samples_start(&self) -> usize {
        if self.n_tabs < 9 {
            self.line.len()
        } else {
            self.tabs[8] + 1
        }
    }

    /// Byte offset of the INFO column text within `line` (columns 1–7 end
    /// here, the FILTER tab included).
    #[must_use]
    pub fn info_start(&self) -> usize {
        self.tabs[6] + 1
    }

    /// The record's column shape (doc_v22 R22-3).
    #[must_use]
    pub fn shape(&self) -> Shape {
        match self.n_tabs {
            0..=7 => Shape::SitesOnly,
            8 => Shape::FormatOnly,
            _ => Shape::WithSamples,
        }
    }

    /// The first required column that is empty (doc_v26 V26-06).
    ///
    /// `ID` and `REF` semantics differ: `ID` may be `.` (the missing-value
    /// symbol), while an empty `ID` is still a column-count error only if the
    /// line is short. Everything checked here is a field that carries data.
    /// `FORMAT` is only checked when the column is present (a sites-only
    /// record has none), and a FORMAT column of `.` is legal — htslib treats
    /// it as "no data" (`vcf_parse_format_empty1`, `vcf.c:3213`).
    #[must_use]
    pub fn first_empty_required(&self) -> Option<EmptyField> {
        if self.chrom().is_empty() {
            return Some(EmptyField::Chrom);
        }
        if self.reference().is_empty() {
            return Some(EmptyField::Ref);
        }
        if self.alts().is_empty() {
            return Some(EmptyField::Alt);
        }
        if self.qual_text().is_empty() {
            return Some(EmptyField::Qual);
        }
        if self.filter_text().is_empty() {
            return Some(EmptyField::Filter);
        }
        if self.info_text().is_empty() {
            return Some(EmptyField::Info);
        }
        if self.n_tabs >= 8 && self.format_text().is_empty() {
            return Some(EmptyField::Format);
        }
        None
    }

    /// Validates the sample columns against a header declaring `expected`
    /// samples (doc_v22 R22-3).
    ///
    /// Empty columns *after* the declared count are trailing tabs and are
    /// ignored (htslib does the same); a non-empty extra column is
    /// [`SampleColumnsError::TooMany`]. An empty column inside the declared
    /// range still counts as a sample — its content is the decoder's
    /// business, not the structural check's.
    ///
    /// # Errors
    ///
    /// [`SampleColumnsError::TooFew`] or [`SampleColumnsError::TooMany`].
    pub fn validate_sample_cols(&self, expected: usize) -> Result<(), SampleColumnsError> {
        if self.shape() != Shape::WithSamples {
            return if expected == 0 {
                Ok(())
            } else {
                Err(SampleColumnsError::TooFew { found: 0, expected })
            };
        }
        let mut seen = 0_usize;
        for (position, field) in self.samples_text().split(|&b| b == b'\t').enumerate() {
            if position < expected {
                seen = position + 1;
            } else if !field.is_empty() {
                return Err(SampleColumnsError::TooMany {
                    found: position + 1,
                    expected,
                });
            }
        }
        if seen < expected {
            Err(SampleColumnsError::TooFew {
                found: seen,
                expected,
            })
        } else {
            Ok(())
        }
    }

    /// The byte range of column `index` in the full `line` (0-based column,
    /// terminator trimmed), or `None` when the line has fewer columns.
    ///
    /// Non-validating and allocation-free (doc_v28 V28-24): a caller that
    /// needs the sample columns *checked* keeps using
    /// [`Self::sample_spans`]; this is the borrowed-range view for callers
    /// that only need one column.
    #[must_use]
    pub fn field_span(&self, index: usize) -> Option<std::ops::Range<usize>> {
        let n_tabs = usize::from(self.n_tabs);
        if index > n_tabs {
            return None;
        }
        let start = if index == 0 {
            0
        } else {
            self.tabs[index - 1] + 1
        };
        let end = if index < n_tabs {
            self.tabs[index]
        } else {
            self.line.len()
        };
        let field = self.trim_line_end(&self.line[start..end]);
        Some(start..start + field.len())
    }

    /// Column `index`'s bytes, or `&[]` when absent — the shape every
    /// consumer of a fixed VCF column already used (chrom, pos, ref, alt…).
    #[must_use]
    pub fn field_at(&self, index: usize) -> &'a [u8] {
        self.field(index)
    }

    /// The byte range of sample column `i` in the full `line`, or `None`
    /// when `i` is past the last column.
    ///
    /// A linear scan from the sample region's start: for a whole-row walk use
    /// [`Self::sample_spans`], which validates the column count as it
    /// records every boundary in one pass.
    #[must_use]
    pub fn sample_span(&self, i: usize) -> Option<std::ops::Range<usize>> {
        let base = self.samples_start();
        if self.shape() != Shape::WithSamples {
            return None;
        }
        let samples = self.samples_text();
        let mut offset = 0_usize;
        for (position, field) in samples.split(|&b| b == b'\t').enumerate() {
            let start = base + offset;
            offset += field.len() + 1;
            if position == i {
                return Some(start..start + field.len());
            }
        }
        None
    }

    /// Sample column `i`'s bytes, or `None` when absent.
    #[must_use]
    pub fn sample_at(&self, i: usize) -> Option<&'a [u8]> {
        self.sample_span(i).map(|span| &self.line[span])
    }

    /// Every sample column as a range, trailing empty columns included, with
    /// no column-count check (doc_v28 V28-24).
    ///
    /// The editor's `-S` path needs exactly the spans of the first `n`
    /// columns it knows exist; validating here would reject a record the scan
    /// stage already accepted.
    pub fn sample_ranges(&self) -> impl Iterator<Item = std::ops::Range<usize>> + 'a {
        let base = self.samples_start();
        let samples = if self.shape() == Shape::WithSamples {
            self.samples_text()
        } else {
            &[][..]
        };
        // An empty region has no columns at all, while `split` would
        // fabricate one empty slice; a non-empty region keeps every column
        // (a trailing tab is a real, empty final column).
        let mut offset = 0_usize;
        let mut done = samples.is_empty();
        std::iter::from_fn(move || {
            if done {
                return None;
            }
            let rest = &samples[offset..];
            let len = memchr::memchr(b'\t', rest).unwrap_or(rest.len());
            let span = base + offset..base + offset + len;
            if len == rest.len() {
                done = true;
            } else {
                offset += len + 1;
            }
            Some(span)
        })
    }

    /// The byte offset where the first `n` sample columns end (one past the
    /// last of their bytes), or the sample region's end when fewer exist —
    /// the prefix length an editor keeps when it rewrites only the first `n`
    /// columns.
    #[must_use]
    pub fn sample_prefix_len(&self, n: usize) -> usize {
        match self.sample_span(n.saturating_sub(1)) {
            Some(span) if n > 0 => span.end,
            _ => self.samples_text().len() + self.samples_start(),
        }
    }

    /// Validates the sample columns **and** records their boundaries in one
    /// pass, appending `(start, end)` byte ranges (relative to the full
    /// `line`) to `out` (doc_v28 V28-17).
    ///
    /// `out` is cleared first and always holds exactly the declared samples
    /// on success. The scanner used to call [`Self::validate_sample_cols`]
    /// and then walk the same sample region again to find its columns; one
    /// pass now serves both, so the region is read once per record.
    ///
    /// # Errors
    ///
    /// [`SampleColumnsError`] exactly as [`Self::validate_sample_cols`];
    /// `out` is left cleared.
    pub fn sample_spans(
        &self,
        expected: usize,
        out: &mut Vec<std::ops::Range<usize>>,
    ) -> Result<(), SampleColumnsError> {
        out.clear();
        if self.shape() != Shape::WithSamples {
            return if expected == 0 {
                Ok(())
            } else {
                Err(SampleColumnsError::TooFew { found: 0, expected })
            };
        }
        let base = self.samples_start();
        let samples = self.samples_text();
        // The same split `validate_sample_cols` walks, but recording each
        // declared column's full-line range while the cursor is there.
        let mut offset = 0_usize;
        for (position, field) in samples.split(|&b| b == b'\t').enumerate() {
            let start = base + offset;
            offset += field.len() + 1;
            if position < expected {
                out.push(start..start + field.len());
            } else if !field.is_empty() {
                out.clear();
                return Err(SampleColumnsError::TooMany {
                    found: position + 1,
                    expected,
                });
            }
        }
        if out.len() < expected {
            // Read the count before clearing: the diagnostic carries it.
            let found = out.len();
            out.clear();
            return Err(SampleColumnsError::TooFew { found, expected });
        }
        Ok(())
    }
}

/// A fixed column that must carry text but is empty (doc_v26 V26-06).
///
/// VCF's missing-value symbol `.` is legal in `ID`, `ALT`, `QUAL`, `FILTER`
/// and `INFO`, but an **empty** field is not: a record like
/// `chr1\t1\t.\t\tG\t.\tPASS\t.` used to classify as a valid record and
/// then flow through filter, stats and convert unchanged, producing output
/// that no other VCF tool accepts. The check is structural only — it does not
/// interpret the column — so every consumer keeps reading the same bytes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmptyField {
    /// The CHROM column is empty.
    Chrom,
    /// The REF column is empty.
    Ref,
    /// The ALT column is empty.
    Alt,
    /// The QUAL column is empty.
    Qual,
    /// The FILTER column is empty.
    Filter,
    /// The INFO column is empty.
    Info,
    /// The FORMAT column is empty (only checked when the column exists).
    Format,
}

impl EmptyField {
    /// The column's name, for diagnostics.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Chrom => "CHROM",
            Self::Ref => "REF",
            Self::Alt => "ALT",
            Self::Qual => "QUAL",
            Self::Filter => "FILTER",
            Self::Info => "INFO",
            Self::Format => "FORMAT",
        }
    }
}

/// The first structural check a line failed (`classify_line`; doc_v25 V25-03).
///
/// Kept as data rather than a message so each consumer can render its own
/// diagnostic (convert's once-per-run warning, stats' fatal job error) without
/// re-deriving which check failed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MalformedLine {
    /// Fewer than eight tab-separated columns.
    FewerThanEightColumns,
    /// The sample columns disagree with the header's declared count.
    SampleColumns(SampleColumnsError),
    /// POS is not a 1-based decimal (`0` included: VCF positions start at 1).
    Position,
    /// A required column is empty (doc_v26 V26-06). `.` is the legal
    /// missing-value symbol and passes; zero bytes does not.
    Empty(EmptyField),
    /// A tag appears more than once in INFO (doc_v26 V26-04). htslib's own
    /// readers take the first match (`bcf_get_info_id`, `vcf.c:6207-6216`),
    /// while reed's index and AC/AN readers disagreed (first vs last), so one
    /// record could filter two different ways. The column is rejected instead
    /// of choosing a convention per consumer — a deliberate divergence from
    /// htslib, which accepts the record silently.
    ///
    /// Also carries the INFO tag-count ceiling (doc_v26 review item 2), which
    /// bounds the duplicate check's quadratic scan.
    InfoColumn(super::info::InfoColumnError),
}

impl MalformedLine {
    /// A short label for warnings and job errors.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::FewerThanEightColumns => "fewer than 8 columns",
            Self::SampleColumns(_) => "a sample-column count differing from the header",
            Self::Position => "an unparsable POS",
            Self::Empty(_) => "an empty required column",
            Self::InfoColumn(super::info::InfoColumnError::Duplicate(_)) => "a duplicate INFO tag",
            Self::InfoColumn(super::info::InfoColumnError::TooManyTags { .. }) => {
                "too many INFO tags"
            }
        }
    }
}

/// How one line of a newline-framed stream classifies (doc_v25 V25-03).
///
/// No `PartialEq`: the [`Record`](Self::Record) variant carries a borrowed
/// [`LineRef`] (offset-only, no value semantics); consumers match on the
/// variant they care about.
#[derive(Clone, Copy, Debug)]
pub enum RecordClass<'a> {
    /// A structurally valid record.
    Record(LineRef<'a>),
    /// Only line terminators (`""`, `"\n"`, `"\r"`, `"\r\n"`): not a record,
    /// not a malformed one either.
    Blank,
    /// A non-blank line that fails a structural check.
    Malformed(MalformedLine),
}

/// The complete lines of a chunk body.
///
/// A body is a run of complete lines ending in `\n` (possibly empty). Splitting
/// on `\n` therefore leaves one empty piece after the final terminator that is
/// *not* a blank line; this iterator drops exactly that piece and keeps genuine
/// blank lines, so a consumer can classify every yielded slice (doc_v25
/// V25-03).
pub fn body_lines(body: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut lines = body.split(|&b| b == b'\n').peekable();
    std::iter::from_fn(move || {
        let line = lines.next()?;
        if line.is_empty() && lines.peek().is_none() {
            return None;
        }
        Some(line)
    })
}

/// Classifies one line for every record consumer (doc_v25 V25-03).
///
/// This is the single definition of "is this line a record": `LineRef` parses
/// it, its sample columns match `n_samples`, and POS is a 1-based decimal.
/// Entrances only choose what to do with [`RecordClass::Blank`] /
/// [`RecordClass::Malformed`]; they never re-derive the judgement.
///
/// `line` may or may not carry its terminator: chunk bodies split on `\n`
/// (no terminator), seam records arrive with one.
#[must_use]
pub fn classify_line(line: &[u8], n_samples: usize) -> RecordClass<'_> {
    classify_line_with_index(line, n_samples, None)
}

/// [`classify_line`], optionally handing back the INFO ranges the column was
/// validated with.
///
/// The filter head needs the same [`InfoIndex`](super::info::InfoIndex) to
/// answer its typed rules and expressions; building it here means one INFO
/// walk per record instead of two (doc_v28 V28-13). The index is left
/// untouched when the line is not a record, so a caller must only reuse it
/// for the [`LineRef`] this call returned.
#[must_use]
pub fn classify_line_with_index<'a>(
    line: &'a [u8],
    n_samples: usize,
    index: Option<&mut super::info::InfoIndex>,
) -> RecordClass<'a> {
    if line.iter().all(|&b| b == b'\n' || b == b'\r') {
        return RecordClass::Blank;
    }
    let Some(lr) = LineRef::parse(line) else {
        return RecordClass::Malformed(MalformedLine::FewerThanEightColumns);
    };
    // Required-field emptiness is judged before the sample-count agreement:
    // a record with no REF is broken regardless of how many sample columns
    // it carries, and reporting it as a column-count mismatch would hide the
    // real defect (doc_v26 V26-06).
    if let Some(field) = lr.first_empty_required() {
        return RecordClass::Malformed(MalformedLine::Empty(field));
    }
    // One walk serves both the duplicate/tag-count rule (doc_v26 V26-04) and
    // the ranges the head reuses; `validate_info` and `InfoIndex::build`
    // share `admit_tag`, so the two entrances cannot disagree.
    let info_check = match index {
        Some(index) => index.build(lr.info_text()),
        None => crate::vcf::info::validate_info(lr.info_text()),
    };
    if let Err(err) = info_check {
        return RecordClass::Malformed(MalformedLine::InfoColumn(err));
    }
    if let Err(err) = lr.validate_sample_cols(n_samples) {
        return RecordClass::Malformed(MalformedLine::SampleColumns(err));
    }
    match super::fixed::parse_pos1(lr.pos_text(), lr.chrom()) {
        Ok(pos1) if pos1 > 0 => RecordClass::Record(lr),
        _ => RecordClass::Malformed(MalformedLine::Position),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn field_access_matches_manual_split() {
        let line = b"chr1\t123\trs1\tA\tG\t99.5\tPASS\tDP=5\tGT:DP\t0/1:9\t1|1:.\n";
        let line_ref = LineRef::parse(line).expect("invariant: 10 columns");
        assert_eq!(line_ref.chrom(), b"chr1");
        assert_eq!(line_ref.pos_text(), b"123");
        assert_eq!(line_ref.id(), b"rs1");
        assert_eq!(line_ref.reference(), b"A");
        assert_eq!(line_ref.alts(), b"G");
        assert_eq!(line_ref.qual_text(), b"99.5");
        assert_eq!(line_ref.filter_text(), b"PASS");
        assert_eq!(line_ref.info_text(), b"DP=5");
        assert_eq!(line_ref.format_text(), b"GT:DP");
        assert_eq!(line_ref.samples_text(), b"0/1:9\t1|1:.");
        assert_eq!(line_ref.n_sample_cols(), 2);
    }

    #[test]
    fn crlf_and_missing_trailing_newline_tolerated() {
        assert_eq!(
            LineRef::parse(b"chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0")
                .expect("invariant")
                .samples_text(),
            b"0/0"
        );
        assert_eq!(
            LineRef::parse(b"chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\r\n")
                .expect("invariant")
                .samples_text(),
            b"0/0"
        );
    }

    #[test]
    fn sites_only_lines_parse_without_format_or_samples() {
        let line = b"chr1\t1\t.\tA\tG\t.\tPASS\tDP=5\n";
        let line_ref = LineRef::parse(line).expect("invariant: 8 columns");
        assert_eq!(line_ref.info_text(), b"DP=5");
        assert_eq!(line_ref.format_text(), b"");
        assert_eq!(line_ref.samples_text(), b"");
        assert_eq!(line_ref.n_sample_cols(), 0);
        assert_eq!(line_ref.samples_start(), line.len());
        // A FORMAT column without samples is 9 columns / 8 tabs.
        let line = b"chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\n";
        let line_ref = LineRef::parse(line).expect("invariant: 9 columns");
        assert_eq!(line_ref.format_text(), b"GT");
        assert_eq!(line_ref.samples_text(), b"");
        assert_eq!(line_ref.n_sample_cols(), 0);
        // Fewer than eight columns is still a truncated record.
        assert!(LineRef::parse(b"chr1\t1\t.\tA\tG\t.\tPASS\n").is_none());
    }

    /// The contract `validate_sample_cols` must satisfy, written from the
    /// naive field split (doc_v22 R22-3).
    fn reference_validation(fields: &[String], expected: usize) -> Result<(), SampleColumnsError> {
        if fields.len() < 10 {
            return if expected == 0 {
                Ok(())
            } else {
                Err(SampleColumnsError::TooFew { found: 0, expected })
            };
        }
        let samples = &fields[9..];
        if let Some(offset) = samples
            .iter()
            .skip(expected)
            .position(|field| !field.is_empty())
        {
            return Err(SampleColumnsError::TooMany {
                found: expected + offset + 1,
                expected,
            });
        }
        let found = expected.min(samples.len());
        if found < expected {
            Err(SampleColumnsError::TooFew { found, expected })
        } else {
            Ok(())
        }
    }

    proptest::proptest! {
        #![proptest_config(proptest::prelude::ProptestConfig {
            cases: 2_000,
            ..proptest::prelude::ProptestConfig::default()
        })]

        /// Every column shape parses deterministically: the accessors agree
        /// with a naive tab split, and no input panics (doc_v22 R22-3).
        #[test]
        fn shape_and_slices_match_a_naive_split(
            fields in proptest::collection::vec(
                proptest::string::string_regex("[A-Za-z0-9._,:=/|.-]{0,10}")
                    .expect("invariant: regex"),
                5..=12,
            ),
            suffix in proptest::prop_oneof![
                proptest::prelude::Just(""),
                proptest::prelude::Just("\n"),
                proptest::prelude::Just("\r\n"),
            ],
        ) {
            let core = fields.join("\t");
            let line = format!("{core}{suffix}");
            let bytes = line.as_bytes();
            let parsed = LineRef::parse(bytes);
            proptest::prop_assert_eq!(parsed.is_some(), fields.len() >= 8);
            let Some(line_ref) = parsed else { return Ok(()) };

            let shape = match fields.len() {
                8 => Shape::SitesOnly,
                9 => Shape::FormatOnly,
                _ => Shape::WithSamples,
            };
            proptest::prop_assert_eq!(line_ref.shape(), shape);
            proptest::prop_assert_eq!(line_ref.chrom(), fields[0].as_bytes());
            proptest::prop_assert_eq!(line_ref.pos_text(), fields[1].as_bytes());
            proptest::prop_assert_eq!(line_ref.id(), fields[2].as_bytes());
            proptest::prop_assert_eq!(line_ref.reference(), fields[3].as_bytes());
            proptest::prop_assert_eq!(line_ref.alts(), fields[4].as_bytes());
            proptest::prop_assert_eq!(line_ref.qual_text(), fields[5].as_bytes());
            proptest::prop_assert_eq!(line_ref.filter_text(), fields[6].as_bytes());
            proptest::prop_assert_eq!(line_ref.info_text(), fields[7].as_bytes());

            let format = fields.get(8).map_or("", String::as_str);
            proptest::prop_assert_eq!(line_ref.format_text(), format.as_bytes());
            let samples = if fields.len() > 9 {
                fields[9..].join("\t")
            } else {
                String::new()
            };
            proptest::prop_assert_eq!(line_ref.samples_text(), samples.as_bytes());
            let sample_cols = fields.len().saturating_sub(9);
            proptest::prop_assert_eq!(line_ref.n_sample_cols(), sample_cols);
            let sample_start = if fields.len() > 9 {
                fields[..9].iter().map(|f| f.len() + 1).sum()
            } else {
                bytes.len()
            };
            proptest::prop_assert_eq!(line_ref.samples_start(), sample_start);
        }

        /// Validation is total and agrees with the reference for every
        /// header sample count, including trailing empty columns.
        #[test]
        fn validate_sample_cols_is_total_and_matches_the_contract(
            fields in proptest::collection::vec(
                proptest::string::string_regex("[A-Za-z0-9._,:=/|.-]{0,6}")
                    .expect("invariant: regex"),
                5..=12,
            ),
            expected in 0_usize..=6,
        ) {
            let line = format!("{}\n", fields.join("\t"));
            let Some(line_ref) = LineRef::parse(line.as_bytes()) else {
                return Ok(());
            };
            let got = line_ref.validate_sample_cols(expected);
            let want = reference_validation(&fields, expected);
            proptest::prop_assert_eq!(got, want);
        }

        /// doc_v28 V28-17: the one-pass `sample_spans` must reach the same
        /// verdicts as `validate_sample_cols`, and on success its ranges must
        /// slice exactly the declared sample columns.
        #[test]
        fn sample_spans_matches_the_validator_and_the_columns(
            fields in proptest::collection::vec(
                proptest::string::string_regex("[A-Za-z0-9._,:=/|.-]{0,6}")
                    .expect("invariant: regex"),
                5..=12,
            ),
            expected in 0_usize..=6,
        ) {
            let line = format!("{}\n", fields.join("\t"));
            let Some(line_ref) = LineRef::parse(line.as_bytes()) else {
                return Ok(());
            };
            let mut spans = Vec::new();
            let got = line_ref.sample_spans(expected, &mut spans);
            let want = reference_validation(&fields, expected);
            proptest::prop_assert_eq!(got, want);
            if got.is_ok() {
                proptest::prop_assert_eq!(spans.len(), expected);
                for (i, span) in spans.iter().enumerate() {
                    let column = std::str::from_utf8(&line.as_bytes()[span.clone()])
                        .expect("invariant: the generator writes ASCII");
                    proptest::prop_assert_eq!(column, fields[9 + i].as_str());
                }
            } else {
                proptest::prop_assert!(spans.is_empty(), "a failed validation leaves no table");
            }
        }
    }

    /// doc_v28 V28-24: the unvalidated accessors must agree with
    /// `sample_spans` and with the fixed-column accessors, including the
    /// missing-column and trailing-empty cases.
    #[test]
    fn borrowed_span_accessors_match_the_validated_walk() {
        let line = b"chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:9\t0/0:4\n";
        let lr = LineRef::parse(line).expect("invariant: 8 fixed columns plus samples");
        // Fixed columns: the whole-line range slices the field's bytes. Only
        // the nine stored tab offsets are addressable (the struct is `Copy`
        // and holds exactly those), so column 9 is the last indexable one.
        for index in 0..10 {
            let span = lr.field_span(index).expect("invariant: column present");
            assert_eq!(&line[span], lr.field_at(index), "column {index}");
        }
        assert_eq!(lr.field_at(0), b"chr1");
        assert_eq!(lr.field_at(3), b"A");
        assert_eq!(lr.field_at(4), b"G");
        // Column 9 starts after the FORMAT tab, so it is the whole sample
        // region: the last column is open-ended by construction.
        assert_eq!(lr.field_at(9), b"0/1:9\t0/0:4");
        // Past the stored columns both views report absence.
        assert_eq!(lr.field_span(10), None);
        assert_eq!(lr.field_at(10), b"");

        // Sample columns: the single-column accessor matches the validated
        // table for every declared column, and past the end is `None`.
        let mut spans = Vec::new();
        lr.sample_spans(2, &mut spans)
            .expect("invariant: two columns");
        assert_eq!(spans.len(), 2);
        for (i, want) in spans.iter().enumerate() {
            let got = lr.sample_span(i).expect("invariant: declared column");
            assert_eq!(&got, want, "sample {i}");
            assert_eq!(lr.sample_at(i), Some(&line[want.clone()]));
        }
        assert_eq!(lr.sample_span(2), None);
        assert_eq!(lr.sample_at(5), None);
        let all: Vec<_> = lr.sample_ranges().collect();
        assert_eq!(all.len(), 2, "exactly the two sample columns");
        assert_eq!(&all[0], &spans[0]);
        assert_eq!(&all[1], &spans[1]);
        // The prefix length ends where the nth column ends.
        assert_eq!(lr.sample_prefix_len(2), spans[1].end);
        assert_eq!(lr.sample_prefix_len(1), spans[0].end);

        // A trailing empty column past the declared count is a column of the
        // region but not a sample: the validated walk accepts it, the
        // iterator sees it, and `.take(n)` drops it — the consumers' rule.
        let padded = b"chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:9\t0/0:4\t\n";
        let lr = LineRef::parse(padded).expect("invariant: trailing empty column");
        let mut spans = Vec::new();
        lr.sample_spans(2, &mut spans)
            .expect("invariant: trailing empties are ignored");
        assert_eq!(spans.len(), 2);
        let all: Vec<_> = lr.sample_ranges().collect();
        assert_eq!(all.len(), 3, "two columns plus the empty tail");
        assert_eq!(all[2].start, all[2].end, "the tail column is empty");

        // A sites-only record has no sample columns at all.
        let sites =
            LineRef::parse(b"chr1\t10\t.\tA\tG\t50\tPASS\t.\n").expect("invariant: 8 columns");
        assert_eq!(sites.sample_span(0), None);
        assert_eq!(sites.sample_ranges().count(), 0);
        assert_eq!(sites.sample_prefix_len(3), sites.samples_start());
    }

    /// Every terminator-only form is Blank, not malformed: a blank line may
    /// straddle a chunk seam (`"\n"`) or be a CRLF blank (`"\r\n"`), and no
    /// consumer should count a terminator as a record (doc_v25 V25-03).
    #[test]
    fn terminator_only_lines_are_blank() {
        for line in [b"".as_slice(), b"\n", b"\r", b"\r\n", b"\n\r", b"\r\r\n\n"] {
            assert!(
                matches!(classify_line(line, 1), RecordClass::Blank),
                "{:?}",
                String::from_utf8_lossy(line)
            );
        }
    }

    /// The classifier names the failed check and returns the parsed line
    /// otherwise; a record without its terminator (body) and with one (seam)
    /// classify the same.
    #[test]
    fn classify_line_names_the_failed_check() {
        let record = b"chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        assert!(matches!(
            classify_line(record, 1),
            RecordClass::Record(lr) if lr.samples_text() == b"0/1"
        ));
        assert!(matches!(
            classify_line(&record[..record.len() - 1], 1),
            RecordClass::Record(_)
        ));
        assert!(matches!(
            classify_line(b"chr1\t1\t.\tA\tG\t50\tPASS\n", 1),
            RecordClass::Malformed(MalformedLine::FewerThanEightColumns)
        ));
        assert!(matches!(
            classify_line(record, 2),
            RecordClass::Malformed(MalformedLine::SampleColumns(SampleColumnsError::TooFew {
                found: 1,
                expected: 2,
            }))
        ));
        assert!(matches!(
            classify_line(b"chr1\tx\t.\tA\tG\t50\tPASS\t.\n", 0),
            RecordClass::Malformed(MalformedLine::Position)
        ));
        assert!(matches!(
            classify_line(b"chr1\t0\t.\tA\tG\t50\tPASS\t.\n", 0),
            RecordClass::Malformed(MalformedLine::Position)
        ));
        // Whitespace is not blank: it is a malformed line.
        assert!(matches!(
            classify_line(b" ", 0),
            RecordClass::Malformed(MalformedLine::FewerThanEightColumns)
        ));
    }

    /// doc_v26 V26-06: required columns must carry text. `.` is the legal
    /// missing-value symbol and passes; zero bytes does not. The record that
    /// motivated the check is the empty REF below, which every other VCF tool
    /// rejects and reed used to pass through unchanged.
    #[test]
    fn empty_required_columns_are_malformed() {
        let cases: &[(&[u8], EmptyField)] = &[
            (b"\t1\t.\tA\tG\t50\tPASS\t.\n", EmptyField::Chrom),
            (b"chr1\t1\t.\t\tG\t50\tPASS\t.\n", EmptyField::Ref),
            (b"chr1\t1\t.\tA\t\t50\tPASS\t.\n", EmptyField::Alt),
            (b"chr1\t1\t.\tA\tG\t\tPASS\t.\n", EmptyField::Qual),
            (b"chr1\t1\t.\tA\tG\t50\t\t.\n", EmptyField::Filter),
            (b"chr1\t1\t.\tA\tG\t50\tPASS\t\n", EmptyField::Info),
            (
                b"chr1\t1\t.\tA\tG\t50\tPASS\t.\t\t0/1\n",
                EmptyField::Format,
            ),
        ];
        for (line, field) in cases {
            assert!(
                matches!(
                    classify_line(line, 1),
                    RecordClass::Malformed(MalformedLine::Empty(seen)) if seen == *field
                ),
                "{:?} must report {:?}",
                String::from_utf8_lossy(line),
                field
            );
        }
        // The legal missing-value symbols all pass.
        assert!(matches!(
            classify_line(b"chr1\t1\t.\tA\t.\t.\t.\t.\tGT\t0/1\n", 1),
            RecordClass::Record(_)
        ));
        // A sites-only record has no FORMAT column to be empty.
        assert!(matches!(
            classify_line(b"chr1\t1\t.\tA\tG\t.\tPASS\t.\n", 0),
            RecordClass::Record(_)
        ));
    }

    /// doc_v26 V26-04: a duplicated INFO tag is rejected by the classifier,
    /// so every consumer sees the same judgement instead of the index's
    /// first match and the AC/AN reader's last one disagreeing.
    #[test]
    fn duplicate_info_tags_are_malformed() {
        let line = b"chr1\t1\t.\tA\tG\t50\tPASS\tAC=1;AC=3;AN=4\tGT\t0/1\n";
        match classify_line(line, 1) {
            RecordClass::Malformed(MalformedLine::InfoColumn(
                crate::vcf::info::InfoColumnError::Duplicate(key),
            )) => {
                assert_eq!(key.key_len, 2);
            }
            other => panic!("expected a duplicate-INFO rejection, got {other:?}"),
        }
    }

    /// `body_lines` drops only the piece after the body's final terminator;
    /// genuine blank lines survive so the classifier can see them.
    #[test]
    fn body_lines_keeps_genuine_blanks_and_drops_the_artifact() {
        let cases: &[(&[u8], &[&[u8]])] = &[
            (b"", &[]),
            (b"a\n", &[b"a"]),
            (b"a\n\n", &[b"a", b""]),
            (b"\n", &[b""]),
            (b"a\n\nb\n", &[b"a", b"", b"b"]),
            (b"a", &[b"a"]),
        ];
        for (body, want) in cases {
            let got: Vec<&[u8]> = body_lines(body).collect();
            assert_eq!(got, *want, "body {:?}", String::from_utf8_lossy(body));
        }
    }

    /// The labels are stable: they are the strings consumers put in
    /// warnings and job errors.
    #[test]
    fn malformed_labels_name_the_check() {
        assert_eq!(
            MalformedLine::FewerThanEightColumns.label(),
            "fewer than 8 columns"
        );
        assert_eq!(
            MalformedLine::SampleColumns(SampleColumnsError::TooFew {
                found: 0,
                expected: 1
            })
            .label(),
            "a sample-column count differing from the header"
        );
        assert_eq!(MalformedLine::Position.label(), "an unparsable POS");
    }
}

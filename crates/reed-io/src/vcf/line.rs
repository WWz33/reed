//! [`LineRef`]: a borrowed VCF record line with the byte offsets of its
//! fixed columns (`docs/spec/architecture.md` §6.7). No allocation, no
//! parsing — offsets only; the field parsers live in [`super::fixed`].

/// Byte-offset view of one complete VCF record line (terminating `\n`
/// included in `line` but excluded from every field slice).
#[derive(Clone, Copy, Debug)]
pub struct LineRef<'a> {
    /// The full line, `\n` included.
    pub line: &'a [u8],
    /// Positions of the first nine tabs (after CHROM, POS, ID, REF, ALT,
    /// QUAL, FILTER, INFO, FORMAT).
    tabs: [usize; 9],
}

impl<'a> LineRef<'a> {
    /// Parses the tab offsets; `None` when the line has fewer than nine
    /// tabs (sites-only VCF or a truncated line).
    #[must_use]
    pub fn parse(line: &'a [u8]) -> Option<Self> {
        let mut tabs = [0_usize; 9];
        let mut used = 0;
        let mut from = 0;
        while used < 9 {
            let rel = memchr::memchr(b'\t', &line[from..])?;
            from += rel;
            tabs[used] = from;
            from += 1;
            used += 1;
        }
        Some(Self { line, tabs })
    }

    fn field(&self, index: usize) -> &'a [u8] {
        let start = if index == 0 {
            0
        } else {
            self.tabs[index - 1] + 1
        };
        self.trim_line_end(&self.line[start..self.tabs[index]])
    }

    fn trim_line_end(&self, field: &'a [u8]) -> &'a [u8] {
        let no_newline = field.strip_suffix(b"\n").unwrap_or(field);
        no_newline.strip_suffix(b"\r").unwrap_or(no_newline)
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

    /// FORMAT column (text).
    #[must_use]
    pub fn format_text(&self) -> &'a [u8] {
        self.field(8)
    }

    /// The sample columns region (all text after the FORMAT tab).
    #[must_use]
    pub fn samples_text(&self) -> &'a [u8] {
        self.trim_line_end(&self.line[self.tabs[8] + 1..])
    }

    /// Number of sample columns (tabs after FORMAT + 1).
    #[must_use]
    pub fn n_sample_cols(&self) -> usize {
        memchr::memchr_iter(b'\t', self.samples_text()).count() + 1
    }

    /// Byte offset of the first sample column within `line`.
    #[must_use]
    pub fn samples_start(&self) -> usize {
        self.tabs[8] + 1
    }

    /// Byte offset of the INFO column text within `line` (columns 1–7 end
    /// here, the FILTER tab included).
    #[must_use]
    pub fn info_start(&self) -> usize {
        self.tabs[6] + 1
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
    fn sites_only_lines_do_not_parse() {
        assert!(LineRef::parse(b"chr1\t1\t.\tA\tG\t.\tPASS\t.\n").is_none());
    }
}

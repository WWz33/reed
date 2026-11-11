//! The VCF header as an opaque value (doc_v26 V26-10).
//!
//! [`HeaderFacts`](super::HeaderFacts) used to carry a public
//! `noodles::vcf::Header` field, which put a third-party type in `reed-io`'s
//! public API (and, through the facade's re-exports, in `reed`'s). A noodles
//! major upgrade then became a breaking change for every downstream user of
//! this crate, and *only* this crate's format needs ever used the type.
//!
//! [`VcfHeader`] owns the parsed header and exposes the operations reed's
//! commands actually perform:
//!
//! * read the sample/contig/INFO/FORMAT dictionaries (already flattened in
//!   [`HeaderFacts`](super::HeaderFacts));
//! * prune the sample list for a subset (`-S`, `--keep`-style masks);
//! * add the INFO and FILTER definitions a recompute or a `PASS` record
//!   needs;
//! * serialize to VCF text or to a BCF file header;
//! * build the BCF encoders.
//!
//! Nothing else escapes. A caller cannot name the inner type, cannot mutate
//! an arbitrary header line, and a noodles upgrade no longer changes this
//! crate's surface.

use crate::error::IoError;

/// The header of an opened VCF/BCF stream.
///
/// Constructed by [`super::HeaderFacts::parse`] (per open) or
/// [`VcfHeader::parse`]; both reject a header noodles cannot parse.
#[derive(Clone, Debug)]
pub struct VcfHeader {
    /// The parsed header. Private: the boundary is the point of this type
    /// (doc_v26 V26-10).
    inner: noodles::vcf::Header,
}

/// One INFO definition to add when a command recomputes a tag the input
/// header does not declare.
///
/// The `Number`/`Type` mapping is the caller's (it is semantics, not format),
/// but the spelling of the header line is this module's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InfoDefinition {
    /// The tag name (e.g. `AC`).
    pub id: &'static str,
    /// `Number` spelling: `A`, `R`, `G`, a fixed count, or `.`.
    pub number: InfoNumber,
    /// Value type.
    pub ty: InfoType,
    /// The `Description` text (quoted by this module).
    pub description: &'static str,
}

/// The `Number` column of an INFO definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoNumber {
    /// One value per ALT allele.
    A,
    /// One value per allele (REF included).
    R,
    /// One value per genotype.
    G,
    /// A fixed count.
    Fixed(u32),
    /// Unbounded (`.`).
    Unknown,
}

/// The `Type` column of an INFO definition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InfoType {
    /// Integer.
    Integer,
    /// Float.
    Float,
    /// Flag.
    Flag,
    /// String.
    String,
}

impl VcfHeader {
    /// Parses header text (through the `#CHROM` line).
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when the text is not UTF-8 or noodles rejects it.
    pub fn parse(text: &[u8]) -> Result<Self, IoError> {
        let s = std::str::from_utf8(text)
            .map_err(|_| IoError::Header("header text is not UTF-8".into()))?;
        let inner: noodles::vcf::Header = s
            .parse()
            .map_err(|err| IoError::Header(format!("noodles rejected the header: {err}")))?;
        Ok(Self { inner })
    }

    /// The parsed header, for crate-internal consumers only.
    pub(crate) fn noodles(&self) -> &noodles::vcf::Header {
        &self.inner
    }

    /// Sample names in header order.
    #[must_use]
    pub fn sample_names(&self) -> Vec<String> {
        self.inner
            .sample_names()
            .iter()
            .map(ToString::to_string)
            .collect()
    }

    /// Number of declared samples.
    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.inner.sample_names().len()
    }

    /// Contig names in header order.
    #[must_use]
    pub fn contig_names(&self) -> Vec<String> {
        self.inner
            .contigs()
            .keys()
            .map(ToString::to_string)
            .collect()
    }

    /// Contig lengths, parallel to [`Self::contig_names`].
    #[must_use]
    pub fn contig_lengths(&self) -> Vec<Option<u64>> {
        self.inner
            .contigs()
            .values()
            .map(|c| c.length().map(|l| l as u64))
            .collect()
    }

    /// INFO tag names in header order.
    #[must_use]
    pub fn info_keys(&self) -> Vec<String> {
        self.inner.infos().keys().map(ToString::to_string).collect()
    }

    /// FORMAT tag names in header order.
    #[must_use]
    pub fn format_keys(&self) -> Vec<String> {
        self.inner
            .formats()
            .keys()
            .map(ToString::to_string)
            .collect()
    }

    /// Whether the header declares the INFO tag `id`.
    #[must_use]
    pub fn declares_info(&self, id: &str) -> bool {
        self.inner.infos().contains_key(id)
    }

    /// Whether the header declares the FILTER `id`.
    #[must_use]
    pub fn declares_filter(&self, id: &str) -> bool {
        self.inner.filters().contains_key(id)
    }

    /// Contig names in **BCF dictionary order** (the space a BCF record's
    /// `tid` lives in), with gaps filled by the empty name.
    ///
    /// This is the order noodles' `StringMaps` (and therefore the native
    /// encoder) assigns, honouring an explicit `IDX=` and falling back to the
    /// next free slot exactly as [`crate::bcf::BcfDict`] does. It differs from
    /// header line order whenever the header pins `IDX=` (doc_v28 V28-28).
    ///
    /// A contig whose `IDX=` leaves a hole produces an empty name at that
    /// position: the dictionary is sparse by definition, and a caller that
    /// needs a dense table (an index reference list) must reject the holes
    /// rather than silently shift every later entry.
    #[must_use]
    pub fn contig_names_in_dict_order(&self) -> Vec<String> {
        let mut names: Vec<Option<&str>> = Vec::new();
        for (name, map) in self.inner.contigs() {
            let idx = map.idx().unwrap_or_else(|| {
                // No explicit IDX: htslib takes the next free slot.
                let mut next = 0;
                while names.get(next).is_some_and(Option::is_some) {
                    next += 1;
                }
                next
            });
            if idx >= names.len() {
                names.resize(idx + 1, None);
            }
            names[idx] = Some(name.as_str());
        }
        names
            .into_iter()
            .map(|name| name.unwrap_or_default().to_string())
            .collect()
    }

    /// A copy with the declared samples pruned to `keep` (a full-length
    /// mask, one entry per declared sample, in header order).
    ///
    /// # Panics
    ///
    /// Panics when `keep` is not exactly `n_samples` long — a caller
    /// invariant (doc_v25 V25-07b: the old short-mask default silently kept
    /// the tail).
    #[must_use]
    pub fn with_samples_pruned(&self, keep: &[bool]) -> Self {
        assert_eq!(
            keep.len(),
            self.inner.sample_names().len(),
            "invariant: keep mask length differs from the header's sample count"
        );
        let mut out = self.inner.clone();
        let kept: Vec<String> = self
            .inner
            .sample_names()
            .iter()
            .map(ToString::to_string)
            .zip(keep)
            .filter(|(_, keep)| **keep)
            .map(|(name, _)| name)
            .collect();
        let names = out.sample_names_mut();
        names.clear();
        for name in kept {
            names.insert(name);
        }
        Self { inner: out }
    }

    /// A copy with `definitions` added, skipping tags the header already
    /// declares (a redefinition is never written).
    #[must_use]
    pub fn with_info_definitions(&self, definitions: &[InfoDefinition]) -> Self {
        use noodles::vcf::header::record::value::Map;
        use noodles::vcf::header::record::value::map::Info;
        let mut out = self.inner.clone();
        for def in definitions {
            if out.infos().contains_key(def.id) {
                continue;
            }
            let number = match def.number {
                InfoNumber::A => noodles::vcf::header::record::value::map::info::Number::A,
                InfoNumber::R => noodles::vcf::header::record::value::map::info::Number::R,
                InfoNumber::G => noodles::vcf::header::record::value::map::info::Number::G,
                InfoNumber::Fixed(n) => {
                    noodles::vcf::header::record::value::map::info::Number::Count(n as usize)
                }
                InfoNumber::Unknown => {
                    noodles::vcf::header::record::value::map::info::Number::Unknown
                }
            };
            let ty = match def.ty {
                InfoType::Integer => noodles::vcf::header::record::value::map::info::Type::Integer,
                InfoType::Float => noodles::vcf::header::record::value::map::info::Type::Float,
                InfoType::Flag => noodles::vcf::header::record::value::map::info::Type::Flag,
                InfoType::String => noodles::vcf::header::record::value::map::info::Type::String,
            };
            out.infos_mut().insert(
                def.id.to_string(),
                Map::<Info>::new(number, ty, def.description),
            );
        }
        Self { inner: out }
    }

    /// A copy declaring `PASS` when the header does not.
    ///
    /// htslib declares it on read and write ("the filter PASS must appear
    /// first in the dictionary"); without it a `FILTER=PASS` record fails the
    /// BCF encoder's dictionary lookup and the written BCF header lacks the
    /// line htslib assumes at index 0 (doc_v10 N7/D10-7).
    #[must_use]
    pub fn with_pass_filter(&self) -> Self {
        use noodles::vcf::header::record::value::Map;
        use noodles::vcf::header::record::value::map::Filter;
        let mut out = self.inner.clone();
        if !out.filters().contains_key("PASS") {
            out.filters_mut()
                .insert(String::from("PASS"), Map::<Filter>::pass());
        }
        Self { inner: out }
    }

    /// A copy with the meta-information line `key=value` inserted.
    ///
    /// `key` is a bare agreed key such as `reed_command` or `reference`; the
    /// value must not be empty (a VCF meta value may not be, and noodles
    /// rejects it).
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when `key` is not a valid agreed key or the value
    /// is rejected.
    pub fn with_meta_line(&self, key: &str, value: &str) -> Result<Self, IoError> {
        let mut out = self.inner.clone();
        out.insert(
            key.parse()
                .map_err(|e| IoError::Header(format!("header insert {key}: {e}")))?,
            value.into(),
        )
        .map_err(|e| IoError::Header(format!("header insert {key}: {e}")))?;
        Ok(Self { inner: out })
    }

    /// The header as VCF text (through the `#CHROM` line).
    #[must_use]
    pub fn to_text(&self) -> Vec<u8> {
        let mut out = Vec::new();
        // noodles' writer produces exactly the stored text; it fails only on
        // I/O to a `Vec`, which cannot fail.
        if noodles::vcf::io::Writer::new(&mut out)
            .write_header(&self.inner)
            .is_err()
        {
            debug_assert!(false, "invariant: writing a header to a Vec cannot fail");
        }
        out
    }

    /// The BCF file header (magic + length + NUL-terminated VCF text).
    ///
    /// # Errors
    ///
    /// [`IoError`] when the header cannot be encoded.
    pub fn to_bcf_file_header(&self) -> Result<Vec<u8>, IoError> {
        crate::bcf::bcf_file_header(&self.inner)
    }

    /// Builds the native BCF record encoder for this header.
    ///
    /// # Errors
    ///
    /// [`IoError`] when the header's maps cannot be built.
    pub fn native_encoder(&self) -> Result<crate::bcf::NativeEncoder, IoError> {
        crate::bcf::NativeEncoder::new(&self.inner)
    }

    /// Builds the noodles-backed BCF encoder (the differential path) for
    /// this header, with its header already written.
    ///
    /// # Errors
    ///
    /// [`IoError`] when the header cannot be encoded.
    pub fn noodles_encoder(&self) -> Result<crate::bcf::BcfEncoder, IoError> {
        let mut encoder = crate::bcf::BcfEncoder::new();
        encoder.write_header(&self.inner)?;
        Ok(encoder)
    }

    /// Encodes one VCF record line to BCF bytes through the noodles path —
    /// the differential "oracle A/B" route used when a run opts into it.
    ///
    /// # Errors
    ///
    /// [`IoError`] when the line cannot be encoded.
    pub fn encode_record_noodles(
        &self,
        encoder: &mut crate::bcf::BcfEncoder,
        line: &[u8],
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        encoder.encode_record(&self.inner, line, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header_text(samples: &[&str]) -> Vec<u8> {
        let mut text = String::from(
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
             ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        let mut cols: Vec<&str> = Vec::new();
        let _ = &mut cols;
        for s in samples {
            text.push('\t');
            text.push_str(s);
        }
        text.push('\n');
        text.into_bytes()
    }

    #[test]
    fn facts_round_trip_through_the_opaque_type() {
        let header = VcfHeader::parse(&header_text(&["s1", "s2"])).expect("invariant: header");
        assert_eq!(header.n_samples(), 2);
        assert_eq!(header.sample_names(), ["s1", "s2"]);
        assert_eq!(header.contig_names(), ["chr1"]);
        assert_eq!(header.contig_lengths(), [Some(1000)]);
        assert_eq!(header.info_keys(), ["DP"]);
        assert_eq!(header.format_keys(), ["GT"]);
        assert!(header.declares_info("DP"));
        assert!(!header.declares_info("AC"));
        // Serializing reproduces the same lines (noodles writes ##contig
        // after ##INFO/##FORMAT, which is the ordering every reed output
        // already had: the text path used the same writer before this type).
        let mut got: Vec<String> = String::from_utf8(header.to_text())
            .expect("invariant: utf8")
            .lines()
            .map(str::to_string)
            .collect();
        let mut want: Vec<String> = String::from_utf8(header_text(&["s1", "s2"]))
            .expect("invariant: utf8")
            .lines()
            .map(str::to_string)
            .collect();
        got.sort();
        want.sort();
        assert_eq!(got, want);
    }

    #[test]
    fn pruning_and_definition_additions_are_immutable() {
        let header =
            VcfHeader::parse(&header_text(&["s1", "s2", "s3"])).expect("invariant: header");
        let pruned = header.with_samples_pruned(&[true, false, true]);
        assert_eq!(pruned.sample_names(), ["s1", "s3"]);
        assert_eq!(
            header.sample_names(),
            ["s1", "s2", "s3"],
            "the original is untouched"
        );

        let added = header.with_info_definitions(&[
            InfoDefinition {
                id: "AC",
                number: InfoNumber::A,
                ty: InfoType::Integer,
                description: "Allele count",
            },
            InfoDefinition {
                id: "DP",
                number: InfoNumber::Fixed(1),
                ty: InfoType::Integer,
                description: "should not replace the existing definition",
            },
        ]);
        assert!(added.declares_info("AC"));
        let text = String::from_utf8(added.to_text()).expect("invariant: utf8");
        assert!(
            text.contains("##INFO=<ID=AC,Number=A,Type=Integer"),
            "{text}"
        );
        assert_eq!(
            text.matches("##INFO=<ID=DP").count(),
            1,
            "an existing definition is not duplicated"
        );
    }

    #[test]
    fn pass_is_declared_once() {
        let header = VcfHeader::parse(&header_text(&[])).expect("invariant: header");
        assert!(!header.declares_filter("PASS"));
        let with_pass = header.with_pass_filter();
        assert!(with_pass.declares_filter("PASS"));
        assert!(with_pass.with_pass_filter().declares_filter("PASS"));
    }

    /// doc_v28 V28-28: a BCF record's `tid` is the dictionary index, which an
    /// explicit `IDX=` reorders. For a header that declares `chr2` first but
    /// pins it to `IDX=1`, the dictionary order is `[chr1, chr2]` — the exact
    /// inverse of line order. Writers that build their index reference table
    /// from line order then key every chunk at the wrong reference sequence,
    /// and htslib returns no records at all for either contig.
    #[test]
    fn contig_names_in_dict_order_honours_explicit_idx() {
        let text = b"##fileformat=VCFv4.3\n\
                     ##contig=<ID=chr2,length=200,IDX=1>\n\
                     ##contig=<ID=chr1,length=100,IDX=0>\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header = VcfHeader::parse(text).expect("invariant: header");
        assert_eq!(header.contig_names(), ["chr2", "chr1"], "line order");
        assert_eq!(
            header.contig_names_in_dict_order(),
            ["chr1", "chr2"],
            "dictionary order (record tid space)"
        );
    }

    /// Without `IDX=` the two orders coincide, and a gap left by a sparse
    /// `IDX=` keeps its position instead of shifting later entries down.
    #[test]
    fn contig_names_in_dict_order_fills_sparse_and_identity_cases() {
        let plain = b"##fileformat=VCFv4.3\n\
                      ##contig=<ID=chr1,length=10>\n\
                      ##contig=<ID=chr2,length=20>\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header = VcfHeader::parse(plain).expect("invariant: header");
        assert_eq!(header.contig_names_in_dict_order(), ["chr1", "chr2"]);

        // IDX=2 leaves slots 0 and 1 empty; the name must stay at index 2.
        let sparse = b"##fileformat=VCFv4.3\n\
                       ##contig=<ID=chrX,length=10,IDX=2>\n\
                       #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header = VcfHeader::parse(sparse).expect("invariant: header");
        assert_eq!(header.contig_names_in_dict_order(), ["", "", "chrX"]);
    }
}

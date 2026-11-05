//! VCF header reading: the header text is parsed once per open and exposed
//! through [`HeaderFacts`] plus the opaque [`VcfHeader`] (doc_v26 V26-10).
//!
//! The parsed third-party type never appears in this crate's public API:
//! [`VcfHeader`] owns it and exposes the operations reed's commands perform
//! (prune samples, add definitions, serialize, build encoders). A noodles
//! upgrade is therefore a change confined to `reed-io`.

mod vcf_header;

pub use vcf_header::{InfoDefinition, InfoNumber, InfoType, VcfHeader};

use crate::error::IoError;

/// Facts extracted from the VCF header.
#[derive(Clone, Debug)]
pub struct HeaderFacts {
    /// Full header text (through the `#CHROM` line, newline included).
    pub text: Vec<u8>,
    /// Sample names in header order.
    pub sample_names: Vec<String>,
    /// Contig names in header order (the source of `ContigId`).
    pub contig_names: Vec<String>,
    /// Contig declared lengths, parallel to `contig_names` (`None` when the
    /// header declares no length).
    pub contig_lengths: Vec<Option<u64>>,
    /// INFO definitions in header order (the source of `TagId`).
    pub info_keys: Vec<String>,
    /// FORMAT definitions in header order (the source of `TagId`).
    pub format_keys: Vec<String>,
    /// The parsed header, as an opaque value (doc_v26 V26-10: this used to be
    /// a public `noodles::vcf::Header` field).
    pub header: VcfHeader,
}

/// Builds the plan compiler's tag table from a parsed header (doc_v3 04 §C3):
/// `reed-filter` only ever sees this view.
#[must_use]
pub fn tag_table(header: &VcfHeader) -> reed_core::tag_table::TagTable {
    use reed_core::tag_table::{TagNumber, TagType};
    let noodles = header.noodles();
    let number = |number: noodles::vcf::header::record::value::map::info::Number| match number {
        noodles::vcf::header::record::value::map::info::Number::Count(n) => {
            TagNumber::Count(u32::try_from(n).unwrap_or(u32::MAX))
        }
        noodles::vcf::header::record::value::map::info::Number::A => TagNumber::A,
        noodles::vcf::header::record::value::map::info::Number::R => TagNumber::R,
        noodles::vcf::header::record::value::map::info::Number::G => TagNumber::G,
        noodles::vcf::header::record::value::map::info::Number::Unknown => TagNumber::Unknown,
    };
    let info_type = |ty: noodles::vcf::header::record::value::map::info::Type| match ty {
        noodles::vcf::header::record::value::map::info::Type::Integer => TagType::Integer,
        noodles::vcf::header::record::value::map::info::Type::Float => TagType::Float,
        noodles::vcf::header::record::value::map::info::Type::Flag => TagType::Flag,
        noodles::vcf::header::record::value::map::info::Type::String
        | noodles::vcf::header::record::value::map::info::Type::Character => TagType::Character,
    };
    let format_type = |ty: noodles::vcf::header::record::value::map::format::Type| match ty {
        noodles::vcf::header::record::value::map::format::Type::Integer => TagType::Integer,
        noodles::vcf::header::record::value::map::format::Type::Float => TagType::Float,
        _ => TagType::Character,
    };
    let mut table = reed_core::tag_table::TagTable::new();
    for (key, map) in noodles.infos() {
        table.add_info(key.as_str(), info_type(map.ty()), number(map.number()));
    }
    for (name, _) in noodles.contigs() {
        table.add_contig(name.as_str());
    }
    for (key, map) in noodles.formats() {
        let format_number = match map.number() {
            noodles::vcf::header::record::value::map::format::Number::Count(n) => {
                TagNumber::Count(u32::try_from(n).unwrap_or(u32::MAX))
            }
            noodles::vcf::header::record::value::map::format::Number::AlternateBases => {
                TagNumber::A
            }
            noodles::vcf::header::record::value::map::format::Number::ReferenceAlternateBases => {
                TagNumber::R
            }
            noodles::vcf::header::record::value::map::format::Number::Samples => TagNumber::G,
            _ => TagNumber::Unknown,
        };
        table.add_format(key.as_str(), format_type(map.ty()), format_number);
    }
    table
}

impl HeaderFacts {
    /// Parses header text (accumulated by the chunk reader's header scan).
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when noodles rejects the header.
    pub fn parse(text: &[u8]) -> Result<Self, IoError> {
        let header = VcfHeader::parse(text)?;
        Ok(Self {
            text: text.to_vec(),
            sample_names: header.sample_names(),
            contig_names: header.contig_names(),
            contig_lengths: header.contig_lengths(),
            info_keys: header.info_keys(),
            format_keys: header.format_keys(),
            header,
        })
    }

    /// Number of samples.
    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_minimal_header() {
        let text = b"##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
                      ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
                      ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
                      #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";
        let facts = HeaderFacts::parse(text).expect("invariant: valid header");
        assert_eq!(facts.n_samples(), 2);
        assert_eq!(facts.sample_names, ["s1", "s2"]);
        assert_eq!(facts.contig_names, ["chr1"]);
        assert_eq!(facts.info_keys, ["DP"]);
        assert_eq!(facts.format_keys, ["GT"]);
        // Same lines, in noodles' serialization order (see the round-trip
        // test in `vcf_header`).
        let mut got: Vec<String> = String::from_utf8(facts.header.to_text())
            .expect("invariant: utf8")
            .lines()
            .map(str::to_string)
            .collect();
        let mut want: Vec<String> = String::from_utf8(text.to_vec())
            .expect("invariant: utf8")
            .lines()
            .map(str::to_string)
            .collect();
        got.sort();
        want.sort();
        assert_eq!(got, want);
    }

    /// doc_v23 V23-42: pruning keeps the declared order and drops exactly
    /// the false entries.
    #[test]
    fn pruned_header_preserves_order_and_selection() {
        let text = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\
                      s1\ts2\ts3\ts4\n";
        let facts = HeaderFacts::parse(text).expect("invariant: valid header");
        let out = facts.header.with_samples_pruned(&[true, false, true, true]);
        assert_eq!(out.sample_names(), ["s1", "s3", "s4"]);
    }

    /// doc_v25 V25-07b: the keep mask must be full length; the old
    /// short-mask default kept the tail and hid a caller bug.
    #[test]
    fn pruned_header_requires_a_full_length_keep() {
        let text = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\
                      s1\ts2\ts3\n";
        let facts = HeaderFacts::parse(text).expect("invariant: valid header");
        let out = facts.header.with_samples_pruned(&[false, true, true]);
        assert_eq!(out.sample_names(), ["s2", "s3"]);
        let short = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let _ = facts.header.with_samples_pruned(&[false]);
        }));
        assert!(short.is_err(), "a short keep mask must panic");
    }
}

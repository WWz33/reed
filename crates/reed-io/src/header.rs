//! VCF header reading: the header text is parsed by noodles (per
//! `docs/spec/architecture.md` §6.4), and [`HeaderFacts`] exposes what the
//! field plans and the output layer need.

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
    /// The parsed noodles header (one construction per open; drivers never
    /// re-open the file just to parse the header).
    pub header: noodles::vcf::Header,
}

/// Builds the plan compiler's tag table from a parsed header (doc_v3 04 §C3):
/// `reed-filter` only ever sees this view.
#[must_use]
pub fn tag_table(header: &noodles::vcf::Header) -> reed_core::tag_table::TagTable {
    use noodles::vcf::header::record::value::map::{format, info};
    use reed_core::tag_table::{TagNumber, TagType};
    let number = |number: info::Number| match number {
        info::Number::Count(n) => TagNumber::Count(u32::try_from(n).unwrap_or(u32::MAX)),
        info::Number::A => TagNumber::A,
        info::Number::R => TagNumber::R,
        info::Number::G => TagNumber::G,
        info::Number::Unknown => TagNumber::Unknown,
    };
    let info_type = |ty: info::Type| match ty {
        info::Type::Integer => TagType::Integer,
        info::Type::Float => TagType::Float,
        info::Type::Flag => TagType::Flag,
        info::Type::String | info::Type::Character => TagType::Character,
    };
    let format_type = |ty: format::Type| match ty {
        format::Type::Integer => TagType::Integer,
        format::Type::Float => TagType::Float,
        _ => TagType::Character,
    };
    let mut table = reed_core::tag_table::TagTable::new();
    for (key, map) in header.infos() {
        table.add_info(key.as_str(), info_type(map.ty()), number(map.number()));
    }
    for (key, map) in header.formats() {
        let format_number = match map.number() {
            format::Number::Count(n) => TagNumber::Count(u32::try_from(n).unwrap_or(u32::MAX)),
            format::Number::AlternateBases => TagNumber::A,
            format::Number::ReferenceAlternateBases => TagNumber::R,
            format::Number::Samples => TagNumber::G,
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
        let s = std::str::from_utf8(text)
            .map_err(|_| IoError::Header("header text is not UTF-8".into()))?;
        let header: noodles::vcf::Header = s
            .parse()
            .map_err(|err| IoError::Header(format!("noodles rejected the header: {err}")))?;
        Ok(Self {
            text: text.to_vec(),
            sample_names: header
                .sample_names()
                .iter()
                .map(ToString::to_string)
                .collect(),
            contig_names: header
                .contigs()
                .iter()
                .map(|(name, _)| name.to_string())
                .collect(),
            contig_lengths: header
                .contigs()
                .iter()
                .map(|(_, contig)| contig.length().map(|l| l as u64))
                .collect(),
            info_keys: header.infos().keys().map(|key| key.to_string()).collect(),
            format_keys: header.formats().keys().map(|key| key.to_string()).collect(),
            header,
        })
    }

    /// Number of samples.
    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }

    /// The output header for a sample subset: the declared samples pruned to
    /// the kept ones, everything else untouched.
    #[must_use]
    pub fn pruned_header(&self, keep: &[bool]) -> noodles::vcf::Header {
        let mut out = self.header.clone();
        let kept: Vec<String> = self
            .sample_names
            .iter()
            .enumerate()
            .filter(|(i, _)| keep.get(*i).copied().unwrap_or(true))
            .map(|(_, name)| name.clone())
            .collect();
        let names = out.sample_names_mut();
        names.clear();
        for name in kept {
            names.insert(name);
        }
        out
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
    }
}

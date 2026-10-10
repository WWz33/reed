//! Index writing (`docs/rounds/v2/03-architecture.md` §5): CSI (default, min_shift
//! 14 / depth 6 so every VCF-legal position is covered) or TBI (positions
//! ≥ 2^29 are refused).

use std::path::{Path, PathBuf};

use crate::error::IoError;

/// TBI coordinates are limited to 29 bits (htslib); anything wider needs CSI.
pub const TBI_MAX_POS1: u64 = (1 << 29) - 1;

/// Which index flavor to write. CSI is the default (bcftools parity).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WriteIndex {
    /// No index.
    None,
    /// Coordinate-sorted index (`.csi`), the default.
    #[default]
    Csi,
    /// Tabix index (`.tbi`); positions ≥ 2^29 are refused.
    Tbi,
}

impl WriteIndex {
    /// Parses the optional `--write-index[=csi|tbi]` value: absent = no
    /// index (bcftools semantics — the index is opt-in via `-W`), `-W` alone
    /// or `-W=csi` = CSI, `-W=tbi` = TBI.
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] on an unknown flavor name.
    pub fn parse(value: Option<&str>) -> Result<Self, IoError> {
        match value {
            None => Ok(Self::None),
            Some("") | Some("csi") => Ok(Self::Csi),
            Some("tbi") => Ok(Self::Tbi),
            Some(other) => Err(IoError::Io(format!(
                "unknown index flavor {other:?} (expected csi or tbi)"
            ))),
        }
    }

    /// Whether an index is requested at all.
    #[must_use]
    pub const fn is_some(self) -> bool {
        !matches!(self, Self::None)
    }

    /// The index file suffix.
    #[must_use]
    pub const fn suffix(self) -> &'static str {
        match self {
            Self::Tbi => ".tbi",
            Self::Csi | Self::None => ".csi",
        }
    }

    /// The index path next to `output`.
    #[must_use]
    pub fn path_for(self, output: &Path) -> PathBuf {
        let mut s = output.as_os_str().to_os_string();
        s.push(self.suffix());
        PathBuf::from(s)
    }
}

/// The index writer, per flavor. CSI is keyed by header contig id; TBI maps
/// names itself and is limited to 29-bit coordinates.
pub enum Indexer {
    /// Coordinate-sorted index (contig-id keyed, with the contig count).
    Csi(
        Box<
            noodles::csi::binning_index::indexer::Indexer<
                noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex,
            >,
        >,
        usize,
    ),
    /// Tabix index (name keyed).
    Tbi(Box<noodles::tabix::index::Indexer>),
}

impl Indexer {
    /// Builds an indexer for `flavor` with `contig_names` in header order.
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] when the contig count overflows CSI's ids.
    pub fn new(flavor: WriteIndex, contig_names: &[Vec<u8>]) -> Result<Self, IoError> {
        use noodles::csi::binning_index::index::header::ReferenceSequenceNames;
        let names: ReferenceSequenceNames = contig_names
            .iter()
            .map(|v| bstr::BString::from(v.as_slice()))
            .collect();
        let header = noodles::csi::binning_index::index::Header::builder()
            .set_format(noodles::csi::binning_index::index::header::Format::Vcf)
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'#')
            .set_reference_sequence_names(names)
            .build();
        match flavor {
            WriteIndex::Tbi => {
                let mut indexer = noodles::tabix::index::Indexer::default();
                indexer.set_header(header);
                Ok(Self::Tbi(Box::new(indexer)))
            }
            WriteIndex::Csi | WriteIndex::None => {
                use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
                let indexer =
                    noodles::csi::binning_index::indexer::Indexer::<BinnedIndex>::new(14, 6)
                        .expect("invariant: min_shift/depth are valid")
                        .set_header(header);
                Ok(Self::Csi(Box::new(indexer), contig_names.len()))
            }
        }
    }

    /// Adds one record's chunk `[voff_start, voff_end)`.
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] for out-of-range positions or an over-long TBI
    /// coordinate.
    pub fn add(
        &mut self,
        name: &str,
        contig: reed_core::ContigId,
        pos1: u64,
        ref_len: u32,
        voff_start: u64,
        voff_end: u64,
    ) -> Result<(), IoError> {
        use noodles::bgzf::VirtualPosition;
        use noodles::core::Position as CorePosition;
        use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
        let chunk = Chunk::new(
            VirtualPosition::from(voff_start),
            VirtualPosition::from(voff_end),
        );
        let start = usize::try_from(pos1)
            .ok()
            .and_then(|p| CorePosition::try_from(p).ok())
            .ok_or_else(|| IoError::Io(format!("index position {pos1} out of range")))?;
        let end_pos1 = pos1.saturating_add(u64::from(ref_len.saturating_sub(1)));
        let end = usize::try_from(end_pos1)
            .ok()
            .and_then(|p| CorePosition::try_from(p).ok())
            .ok_or_else(|| IoError::Io(format!("index end {end_pos1} out of range")))?;
        match self {
            Self::Csi(indexer, _) => indexer
                .add_record(
                    Some((
                        usize::try_from(contig.get()).expect("invariant: contig id fits usize"),
                        start,
                        end,
                        true,
                    )),
                    chunk,
                )
                .map_err(|e| IoError::Io(format!("index add_record: {e}"))),
            Self::Tbi(indexer) => {
                if end_pos1 > TBI_MAX_POS1 {
                    return Err(IoError::Io(format!(
                        "{end_pos1}: TBI coordinates are limited to {TBI_MAX_POS1}; use --write-index=csi"
                    )));
                }
                indexer
                    .add_record(name, start, end, chunk)
                    .map_err(|e| IoError::Io(format!("index add_record: {e}")))
            }
        }
    }

    /// Writes the index to `path`.
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] on write failure.
    pub fn write(self, path: &Path) -> Result<(), IoError> {
        match self {
            Self::Csi(indexer, n_ref) => noodles::csi::fs::write(path, &indexer.build(n_ref))
                .map_err(|e| IoError::Io(e.to_string())),
            Self::Tbi(indexer) => noodles::tabix::fs::write(path, &indexer.build())
                .map_err(|e| IoError::Io(e.to_string())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flavor_parsing_and_suffixes() {
        // No `-W` at all = no index (bcftools); `-W`/`-W=csi` selects CSI.
        assert_eq!(
            WriteIndex::parse(None).expect("invariant: default"),
            WriteIndex::None
        );
        assert_eq!(
            WriteIndex::parse(Some("csi")).expect("invariant: csi"),
            WriteIndex::Csi
        );
        assert_eq!(
            WriteIndex::parse(Some("")).expect("invariant: bare -W"),
            WriteIndex::Csi
        );
        assert_eq!(
            WriteIndex::parse(Some("tbi")).expect("invariant: tbi"),
            WriteIndex::Tbi
        );
        assert!(WriteIndex::parse(Some("nope")).is_err());
        assert_eq!(WriteIndex::Tbi.suffix(), ".tbi");
        assert_eq!(WriteIndex::Csi.suffix(), ".csi");
        assert_eq!(
            WriteIndex::Tbi.path_for(Path::new("/tmp/out.vcf.gz")),
            Path::new("/tmp/out.vcf.gz.tbi")
        );
        assert!(!WriteIndex::None.is_some());
    }
}

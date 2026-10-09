//! Contig identifiers.

/// Numeric contig id: the index of the contig in the header's `##contig`
/// order — the tid space of the CSI/TBI indexes reed writes and queries
/// (`docs/rounds/v2/03-architecture.md` §5).
///
/// A record on a contig the header does not declare has no id; readers that
/// address the index must reject it rather than invent one.
pub type ContigId = u32;

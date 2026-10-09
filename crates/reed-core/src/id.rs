//! Strong identifier newtypes (`docs/spec/architecture.md` §5.1).
//!
//! Every identifier is a dense index into a table that lives elsewhere
//! (`ContigTable`, `SampleTable`, `GroupTable`, the header tag dictionary, the
//! chunk stream). Public APIs must not accept bare integers where one of these
//! is meant; hot loops may unbox with `.get()`.

use std::fmt;

macro_rules! id_newtype {
    ($(#[$meta:meta])* $name:ident($inner:ty)) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(pub $inner);

        impl $name {
            /// Unboxes the identifier for hot loops that index arrays directly.
            #[must_use]
            pub const fn get(self) -> $inner {
                self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}", self.0)
            }
        }
    };
}

id_newtype! {
    /// Index of a contig in header `##contig` order (target semantics,
    /// superseding the first-appearance ids of the legacy window path).
    ContigId(u32)
}
id_newtype! {
    /// Index of a sample in the header sample-column order.
    SampleId(u32)
}
id_newtype! {
    /// Index of a group in [`crate::group::GroupTable`] order.
    GroupId(u16)
}
id_newtype! {
    /// Index of an INFO/FORMAT tag in header-definition order.
    TagId(u32)
}
id_newtype! {
    /// Input-order index of a chunk; the only basis of reduction order.
    ChunkId(u64)
}

/// Position of a record inside the chunk stream: chunk index plus the record's
/// local index within the chunk. Sketches and survivor bitmaps are aligned by
/// this pair.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RecordOrdinal {
    /// Input-order chunk index.
    pub chunk: ChunkId,
    /// Record index within the chunk.
    pub local: u32,
}

impl fmt::Display for RecordOrdinal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "chunk {} record {}", self.chunk, self.local)
    }
}

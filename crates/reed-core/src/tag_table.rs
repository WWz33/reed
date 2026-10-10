//! The declared-tag table (doc_v3 04 §C3): the plan compiler's only view of
//! the header, so `reed-filter` stays free of I/O dependencies and the plan
//! can be tested without a VCF header. `reed-io` builds it from the parsed
//! header.

use std::collections::HashMap;

/// A tag's declared value type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TagType {
    /// Integer.
    Integer,
    /// Float.
    Float,
    /// Character (FORMAT GT strings and `Type=String` in general).
    Character,
    /// INFO flag (`Number=0`).
    Flag,
}

impl TagType {
    /// Whether numeric rules accept this type.
    #[must_use]
    pub const fn is_numeric(self) -> bool {
        matches!(self, Self::Integer | Self::Float)
    }
}

/// A tag's declared cardinality.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TagNumber {
    /// A fixed count per record/sample.
    Count(u32),
    /// One value per ALT allele.
    A,
    /// One value per allele (REF included).
    R,
    /// One value per genotype.
    G,
    /// `.` — unknown.
    Unknown,
}

/// One declared tag.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TagDef {
    /// Tag name.
    pub name: String,
    /// Declared type.
    pub ty: TagType,
    /// Declared number.
    pub number: TagNumber,
}

/// The header's declared INFO and FORMAT tags.
#[derive(Clone, Debug, Default)]
pub struct TagTable {
    info: HashMap<String, TagDef>,
    format: HashMap<String, TagDef>,
}

impl TagTable {
    /// An empty table (tests and hand-built plans).
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Declares an INFO tag.
    pub fn add_info(&mut self, name: &str, ty: TagType, number: TagNumber) {
        self.info.insert(
            name.to_string(),
            TagDef {
                name: name.to_string(),
                ty,
                number,
            },
        );
    }

    /// Declares a FORMAT tag.
    pub fn add_format(&mut self, name: &str, ty: TagType, number: TagNumber) {
        self.format.insert(
            name.to_string(),
            TagDef {
                name: name.to_string(),
                ty,
                number,
            },
        );
    }

    /// The INFO tag's definition.
    #[must_use]
    pub fn info(&self, name: &str) -> Option<&TagDef> {
        self.info.get(name)
    }

    /// The FORMAT tag's definition.
    #[must_use]
    pub fn format(&self, name: &str) -> Option<&TagDef> {
        self.format.get(name)
    }

    /// INFO tag names, sorted by name.
    #[must_use]
    pub fn info_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.info.keys().cloned().collect();
        names.sort();
        names
    }

    /// FORMAT tag names, sorted by name.
    #[must_use]
    pub fn format_names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.format.keys().cloned().collect();
        names.sort();
        names
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookups_and_names() {
        let mut table = TagTable::new();
        table.add_info("AC", TagType::Integer, TagNumber::A);
        table.add_format("GT", TagType::Character, TagNumber::Count(1));
        assert_eq!(table.info("AC").expect("AC").ty, TagType::Integer);
        assert!(table.format("AC").is_none());
        assert_eq!(table.format_names(), vec!["GT".to_string()]);
        assert!(TagType::Float.is_numeric());
        assert!(!TagType::Character.is_numeric());
    }
}

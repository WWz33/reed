//! Provenance (§11.8): the canonical spec JSON (keys sorted, floats
//! shortest round-trip), its sha256, and the `##reed_*` header lines that
//! replace the v0.2 Debug-output hash (A19).

use sha2::{Digest, Sha256};

use crate::spec::FilterSpecV2;
use crate::FilterError;

/// Provenance of one run.
#[derive(Clone, Debug)]
pub struct Provenance {
    /// `##reed_version=`.
    pub version: String,
    /// `##reed_command=`.
    pub command: String,
    /// `##reed_filterCommand=` (the canonical flag spelling of the spec,
    /// `FilterSpecV2::to_args`; 02-cli.md §1).
    pub filter_command: String,
    /// `##reed_filter_spec=` (canonical JSON).
    pub canonical_json: String,
    /// `##reed_filter_sha256=`.
    pub sha256: String,
}

impl Provenance {
    /// Builds the provenance: canonical JSON of the spec, hashed together
    /// with the reed version.
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] when the spec fails to serialize (never in
    /// practice: every field is JSON-native).
    pub fn new(spec: &FilterSpecV2, command_line: &str) -> Result<Self, FilterError> {
        let value = serde_json::to_value(spec)
            .map_err(|e| FilterError::Spec(format!("spec serialization failed: {e}")))?;
        // serde_json's Map is a BTreeMap unless `preserve_order` is on, so
        // key order is sorted; floats use the shortest round-trip form.
        let canonical = value.to_string();
        let mut hasher = Sha256::new();
        hasher.update(canonical.as_bytes());
        hasher.update(env!("CARGO_PKG_VERSION").as_bytes());
        let sha256 = format!("{:x}", hasher.finalize());
        Ok(Self {
            version: env!("CARGO_PKG_VERSION").to_string(),
            command: command_line.to_string(),
            filter_command: spec.to_args().join(" "),
            canonical_json: canonical,
            sha256,
        })
    }

    /// The `##reed_*` header lines in output order.
    #[must_use]
    pub fn header_lines(&self, ascertainment: &str) -> Vec<(String, String)> {
        vec![
            (String::from("reed_version"), self.version.clone()),
            (String::from("reed_command"), self.command.clone()),
            (String::from("reed_filterCommand"), self.filter_command.clone()),
            (String::from("reed_filter_spec"), self.canonical_json.clone()),
            (String::from("reed_filter_sha256"), self.sha256.clone()),
            (String::from("reed_ascertainment"), ascertainment.to_string()),
        ]
    }
}

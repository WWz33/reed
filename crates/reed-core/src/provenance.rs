//! Run provenance (doc_v25 V25-08c): the recorded command line, the reed
//! version, and — when the caller has a typed spec — its canonical JSON,
//! SHA-256 and flag rendering.
//!
//! The type lives in core so every command can attach the same header
//! instead of re-spelling `##reed_command` per writer (filter's
//! `Provenance` used to be the only implementation).

use sha2::{Digest, Sha256};

/// Provenance of one run.
#[derive(Clone, Debug)]
pub struct Provenance {
    /// `##reed_version=`.
    pub version: String,
    /// `##reed_command=` (the shell-quoted argv).
    pub command: String,
    /// `##reed_filterCommand=` (the canonical flag spelling of a typed
    /// spec); `None` for commands without one.
    pub spec_command: Option<String>,
    /// `##reed_filter_spec=` (canonical JSON); `None` without a spec.
    pub canonical_json: Option<String>,
    /// `##reed_filter_sha256=`.
    pub sha256: Option<String>,
}

impl Provenance {
    /// Provenance for a run without a typed spec: version + command.
    #[must_use]
    pub fn new(command_line: &str) -> Self {
        Self {
            version: env!("CARGO_PKG_VERSION").to_string(),
            command: command_line.to_string(),
            spec_command: None,
            canonical_json: None,
            sha256: None,
        }
    }

    /// Adds a typed spec: its canonical JSON (keys sorted, floats shortest
    /// round-trip), the hash of JSON + version, and the spec's flag
    /// rendering.
    ///
    /// # Errors
    ///
    /// [`serde_json::Error`] when the spec is not JSON-serializable (never in
    /// practice for the typed specs).
    pub fn with_spec<S: serde::Serialize>(
        mut self,
        spec: &S,
        spec_command: String,
    ) -> Result<Self, serde_json::Error> {
        let value = serde_json::to_value(spec)?;
        // serde_json's Map is a BTreeMap unless `preserve_order` is on, so
        // key order is sorted; floats use the shortest round-trip form.
        let canonical = value.to_string();
        let mut hasher = Sha256::new();
        hasher.update(canonical.as_bytes());
        hasher.update(env!("CARGO_PKG_VERSION").as_bytes());
        self.sha256 = Some(format!("{:x}", hasher.finalize()));
        self.spec_command = Some(spec_command);
        self.canonical_json = Some(canonical);
        Ok(self)
    }

    /// The `##reed_*` header lines in output order.
    #[must_use]
    pub fn header_lines(&self) -> Vec<(String, String)> {
        let mut lines = vec![
            (String::from("reed_version"), self.version.clone()),
            (String::from("reed_command"), self.command.clone()),
        ];
        if let (Some(command), Some(json), Some(hash)) =
            (&self.spec_command, &self.canonical_json, &self.sha256)
        {
            lines.push((String::from("reed_filterCommand"), command.clone()));
            lines.push((String::from("reed_filter_spec"), json.clone()));
            lines.push((String::from("reed_filter_sha256"), hash.clone()));
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(serde::Serialize)]
    struct Spec {
        b: u8,
        a: &'static str,
    }

    #[test]
    fn a_spec_contributes_its_canonical_json_and_hash() {
        let p = Provenance::new("reed filter -i 'AC>1'")
            .with_spec(&Spec { b: 2, a: "x" }, "-i AC>1".into())
            .expect("serializable");
        assert_eq!(p.canonical_json.as_deref(), Some(r#"{"a":"x","b":2}"#));
        assert_eq!(p.sha256.as_ref().map(String::len), Some(64));
        let lines = p.header_lines();
        let keys: Vec<&str> = lines.iter().map(|(k, _)| k.as_str()).collect();
        assert_eq!(
            keys,
            [
                "reed_version",
                "reed_command",
                "reed_filterCommand",
                "reed_filter_spec",
                "reed_filter_sha256",
            ]
        );
    }

    #[test]
    fn without_a_spec_only_the_common_lines_appear() {
        let p = Provenance::new("reed pca in.vcf");
        let lines = p.header_lines();
        let keys: Vec<&str> = lines.iter().map(|(k, _)| k.as_str()).collect();
        assert_eq!(keys, ["reed_version", "reed_command"]);
    }
}

//! The CLI parameter owners (doc_v25 V25-06).
//!
//! Every parameter that has caused repeated fixes lives here once: the
//! `--threads`/`--memory` domain is `reed_engine::budget`'s value parsers,
//! and `-O` is parsed here for every subcommand. The acceptance gate
//! (`ci/check-structure.sh`) allows the `-O` letter/level parsing only in
//! this module.

use reed::Error;

/// One parsed `-O` spelling: the type letter and a 1–9 compression level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct OutputType {
    /// The type letter (`v`, `z`, `b` or `u`).
    pub kind: char,
    /// Compression level (1–9; 6 when the spelling carried none).
    pub level: u8,
}

impl OutputType {
    /// The BCF level (`b` keeps the parsed level, `u` is level 0).
    #[must_use]
    pub(crate) fn bcf_level(self) -> Option<u8> {
        match self.kind {
            'b' => Some(self.level),
            'u' => Some(0),
            _ => None,
        }
    }
}

/// Parses one `-O` spelling: `v|z|b|u` with an optional 1–9 level.
///
/// `spellings` names the accepted kinds in the error text; a command that
/// accepts a subset checks [`OutputType::kind`] afterwards (its own message
/// explains which subset it writes).
///
/// # Errors
///
/// An empty spelling, an unknown type letter, or a level outside 1–9.
pub(crate) fn parse_output_type(spec: &str, spellings: &str) -> Result<OutputType, Error> {
    if spec.is_empty() {
        return Err(Error::usage(format!(
            "invalid -O '': expected {spellings}"
        )));
    }
    let kind = spec.chars().next().expect("invariant: non-empty");
    let rest = &spec[kind.len_utf8()..];
    let level = if rest.is_empty() {
        6
    } else {
        let parsed = rest
            .parse::<u8>()
            .map_err(|_| Error::usage(format!("invalid -O level in {spec:?} (want 1-9)")))?;
        if !(1..=9).contains(&parsed) {
            return Err(Error::usage(format!(
                "invalid -O level in {spec:?} (want 1-9)"
            )));
        }
        parsed
    };
    if !matches!(kind, 'v' | 'z' | 'b' | 'u') {
        return Err(Error::usage(format!(
            "unknown -O type {kind:?} (expected {spellings})"
        )));
    }
    Ok(OutputType { kind, level })
}

/// Text-only `-O` (`v` or `z`): the rename and hapmap→VCF spelling.
///
/// # Errors
///
/// [`parse_output_type`]'s errors, or the command's `unsupported` message
/// when the kind is `b`/`u`.
pub(crate) fn parse_text_output_type(
    spec: &str,
    spellings: &str,
    unsupported: &str,
) -> Result<(bool, u8), Error> {
    let parsed = parse_output_type(spec, spellings)?;
    match parsed.kind {
        'v' => Ok((false, parsed.level)),
        'z' => Ok((true, parsed.level)),
        _ => Err(Error::usage(unsupported.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed::ErrorKind;

    #[test]
    fn output_type_parses_and_validates_once() {
        let z = parse_output_type("z", "v, z, b or u").expect("z");
        assert_eq!((z.kind, z.level), ('z', 6));
        let z9 = parse_output_type("z9", "v, z, b or u").expect("z9");
        assert_eq!(z9.level, 9);
        let u = parse_output_type("u", "v, z, b or u").expect("u");
        assert_eq!(u.kind, 'u');
        assert_eq!(u.bcf_level(), Some(0));

        for bad in ["", "z0", "z10", "x"] {
            let err = parse_output_type(bad, "v, z, b or u").expect_err(bad);
            assert_eq!(err.kind(), ErrorKind::Usage, "{bad}: {err}");
        }
    }

    #[test]
    fn text_output_rejects_bcf() {
        let (compressed, level) =
            parse_text_output_type("z3", "v or z", "text only").expect("z3");
        assert!(compressed);
        assert_eq!(level, 3);
        let err = parse_text_output_type("b", "v or z", "text only").expect_err("b");
        assert_eq!(err.kind(), ErrorKind::Usage);
        assert!(err.to_string().contains("text only"), "{err}");
    }
}

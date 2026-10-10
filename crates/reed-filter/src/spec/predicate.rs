//! Comparison primitives (§11.1): a `keep_if`/`when` string is one fixed
//! binary comparison — not an expression language (DESIGN §16 verdict).

use serde::{Deserialize, Serialize};

use crate::FilterError;

/// Comparison operator.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Op {
    /// `>=`
    #[serde(rename = ">=")]
    Ge,
    /// `>`
    #[serde(rename = ">")]
    Gt,
    /// `<=`
    #[serde(rename = "<=")]
    Le,
    /// `<`
    #[serde(rename = "<")]
    Lt,
    /// `==`
    #[serde(rename = "==")]
    Eq,
    /// `!=`
    #[serde(rename = "!=")]
    Ne,
}

impl Op {
    /// Parses an operator prefix off `text`, returning the rest.
    pub fn parse_from(text: &str) -> Result<(Self, &str), FilterError> {
        for (prefix, op) in [
            (">=", Op::Ge),
            ("<=", Op::Le),
            ("!=", Op::Ne),
            ("==", Op::Eq),
            (">", Op::Gt),
            ("<", Op::Lt),
        ] {
            if let Some(rest) = text.strip_prefix(prefix) {
                return Ok((op, rest));
            }
        }
        Err(FilterError::Spec(format!(
            "comparison must start with one of >= <= == != > <: {text:?}"
        )))
    }

    /// Applies the comparison to two f64 values.
    ///
    /// Exact by design: the threshold and the tag value are the numbers the
    /// user and the file provide; an epsilon would change the boundary.
    #[must_use]
    #[allow(clippy::float_cmp)]
    pub fn apply_f64(self, a: f64, b: f64) -> bool {
        match self {
            Self::Ge => a >= b,
            Self::Gt => a > b,
            Self::Le => a <= b,
            Self::Lt => a < b,
            Self::Eq => a == b,
            Self::Ne => a != b,
        }
    }

    /// Applies the comparison to two i64 values (exact for Integer tags).
    #[must_use]
    pub fn apply_i64(self, a: i64, b: i64) -> bool {
        match self {
            Self::Ge => a >= b,
            Self::Gt => a > b,
            Self::Le => a <= b,
            Self::Lt => a < b,
            Self::Eq => a == b,
            Self::Ne => a != b,
        }
    }
}

/// The typed right-hand side of a comparison.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum Typed {
    /// Integer value (header Type=Integer).
    Int(i64),
    /// Floating value (header Type=Float).
    Float(f64),
}

impl Typed {
    /// Parses a scalar from text (integer when it has no `.`/`e`).
    pub fn parse(text: &str) -> Result<Self, FilterError> {
        let text = text.trim();
        if text.contains('.') || text.contains('e') || text.contains('E') {
            text.parse::<f64>()
                .map(Typed::Float)
                .map_err(|_| FilterError::Spec(format!("not a number: {text:?}")))
        } else {
            match text.parse::<i64>() {
                Ok(v) => Ok(Typed::Int(v)),
                Err(_) => text
                    .parse::<f64>()
                    .map(Typed::Float)
                    .map_err(|_| FilterError::Spec(format!("not a number: {text:?}"))),
            }
        }
    }
}

/// A parsed fixed-format comparison.
#[derive(Clone, Debug, PartialEq)]
pub struct Cmp {
    /// The operator.
    pub op: Op,
    /// The right-hand side.
    pub value: Typed,
}

impl Cmp {
    /// Parses `">= 4.0"` / `"<= 50"` / `"< 20"`.
    pub fn parse(text: &str) -> Result<Self, FilterError> {
        let (op, rest) = Op::parse_from(text)?;
        let value = Typed::parse(rest)?;
        Ok(Self { op, value })
    }

    /// Applies to an integer observation.
    #[must_use]
    pub fn apply_i64(&self, observed: i64) -> bool {
        match self.value {
            Typed::Int(v) => self.op.apply_i64(observed, v),
            Typed::Float(v) => self.op.apply_f64(observed as f64, v),
        }
    }

    /// Applies to a float observation.
    #[must_use]
    pub fn apply_f64(&self, observed: f64) -> bool {
        match self.value {
            Typed::Int(v) => self.op.apply_f64(observed, v as f64),
            Typed::Float(v) => self.op.apply_f64(observed, v),
        }
    }
}

/// Behaviour when the compared tag/value is missing.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum OnMissing {
    /// The predicate passes (bcftools missing-comparison semantics); the
    /// default for INFO rules.
    #[default]
    #[serde(rename = "pass")]
    Pass,
    /// The predicate fails.
    #[serde(rename = "fail")]
    Fail,
    /// The sample is not masked (the default for genotype mask rules).
    #[serde(rename = "keep")]
    Keep,
    /// The sample is masked.
    #[serde(rename = "mask")]
    Mask,
}

/// Value quantifier for `Number=A` tags over multiallelic records.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Quantifier {
    /// At least one ALT value must satisfy the comparison.
    #[default]
    #[serde(rename = "any")]
    Any,
    /// Every ALT value must satisfy the comparison.
    #[serde(rename = "all")]
    All,
}

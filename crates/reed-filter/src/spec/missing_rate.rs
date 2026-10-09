//! The missing-rate threshold type (`doc_v4/02-cli-semantics.md` 回退版 §1,
//! §4).
//!
//! `--geno X` / `--mind X` are upper bounds on the **missing** fraction
//! (`doc_v4` 02 §1: the quantity to suppress is written `≤`, like `--max-maf`
//! and `-Q`). The threshold is stored as the rational `q / 10^9`, so the
//! boundary comparison is an exact integer one — the CLI path performs no
//! `1 − X` arithmetic, and the site/sample rules count **alleles** (a half
//! call contributes one called allele, a masked call none).

use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::FilterError;

/// The denominator of the rational representation (10^9 covers every
/// decimal a user writes).
pub const MISSING_RATE_DENOMINATOR: u64 = 1_000_000_000;

/// A "keep if the missing fraction is at most X" threshold, stored as
/// `q / MISSING_RATE_DENOMINATOR`.
///
/// [`Self::tolerates`] decides `missing · 10^9 ≤ q · total` in `u128`, which
/// is exact; `parse` quantises the decimal to `q = round(X · 10^9)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MissingRate {
    q: u64,
}

impl MissingRate {
    /// Parses a decimal in `[0, 1]` into the exact rational form.
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] when the text is not a number or the value is
    /// outside `[0, 1]`.
    pub fn parse(text: &str) -> Result<Self, FilterError> {
        let value: f64 = text
            .trim()
            .parse()
            .map_err(|_| FilterError::Spec(format!("missing rate is not a number: {text:?}")))?;
        Self::from_f64(value)
    }

    /// Quantises an already-parsed value (validating the range).
    ///
    /// # Errors
    ///
    /// [`FilterError::Spec`] when the value is not finite or outside `[0, 1]`.
    pub fn from_f64(value: f64) -> Result<Self, FilterError> {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(FilterError::Spec(format!(
                "missing rate out of range [0,1]: {value}"
            )));
        }
        // `value · 10^9` is within ±1e-7 of an integer for any decimal the
        // user writes; the rounding is what makes the form exact.
        let q = (value * MISSING_RATE_DENOMINATOR as f64).round();
        Ok(Self { q: q as u64 })
    }

    /// The numerator (`q` in `q / 10^9`).
    #[must_use]
    pub const fn q(self) -> u64 {
        self.q
    }

    /// Whether `missing` out of `total` satisfies the threshold:
    /// `missing · 10^9 ≤ q · total` (exact `u128` arithmetic).
    ///
    /// `total == 0` is vacuously true (`0 ≤ 0`).
    #[must_use]
    pub fn tolerates(self, missing: u64, total: u64) -> bool {
        u128::from(missing) * u128::from(MISSING_RATE_DENOMINATOR)
            <= u128::from(self.q) * u128::from(total)
    }

    /// `1 − X` in the same exact form: the value the withdrawn keep-direction
    /// flags translate to (`--max-missing 0.9` → `--geno 0.1`).
    #[must_use]
    pub const fn complement(self) -> Self {
        Self { q: MISSING_RATE_DENOMINATOR - self.q }
    }

    /// The shortest decimal that parses back to the same threshold
    /// (round-trip for `to_args`, CLI and canonical JSON).
    #[must_use]
    pub fn render(self) -> String {
        let int = self.q / MISSING_RATE_DENOMINATOR;
        let mut frac = format!("{:09}", self.q % MISSING_RATE_DENOMINATOR);
        while frac.ends_with('0') {
            frac.pop();
        }
        if frac.is_empty() {
            int.to_string()
        } else {
            format!("{int}.{frac}")
        }
    }
}

impl std::fmt::Display for MissingRate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.render())
    }
}

impl Serialize for MissingRate {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        // A JSON number keeps the canonical spec readable; deserialisation
        // re-quantises to the same `q`.
        serializer.serialize_f64(self.render().parse().expect("invariant: rendered as a decimal"))
    }
}

impl<'de> Deserialize<'de> for MissingRate {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct RateVisitor;
        impl Visitor<'_> for RateVisitor {
            type Value = MissingRate;

            fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str("a missing rate (number or decimal string) in [0,1]")
            }

            fn visit_f64<E: de::Error>(self, v: f64) -> Result<MissingRate, E> {
                MissingRate::from_f64(v).map_err(E::custom)
            }

            fn visit_i64<E: de::Error>(self, v: i64) -> Result<MissingRate, E> {
                self.visit_f64(v as f64)
            }

            fn visit_u64<E: de::Error>(self, v: u64) -> Result<MissingRate, E> {
                self.visit_f64(v as f64)
            }

            fn visit_str<E: de::Error>(self, v: &str) -> Result<MissingRate, E> {
                MissingRate::parse(v).map_err(E::custom)
            }
        }
        deserializer.deserialize_any(RateVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v4 02 回退版 §7.6: `X = 0.2` keeps exactly 2 of 10 missing
    /// (equality) and drops 3; `X = 0.3` keeps 3 (0.3 has no exact binary
    /// form — the rational comparison does not care). The site rule passes
    /// `total = 2N` (allele slots), the sample rule `2L`.
    #[test]
    fn exact_boundary_decisions() {
        let rate = MissingRate::parse("0.2").expect("valid");
        assert!(rate.tolerates(2, 10), "2/10 = 0.2 keeps (equality)");
        assert!(!rate.tolerates(3, 10), "3/10 > 0.2 drops");
        let rate = MissingRate::parse("0.3").expect("valid");
        assert!(rate.tolerates(3, 10), "3/10 ≤ 0.3 keeps");
        assert!(!rate.tolerates(4, 10), "4/10 > 0.3 drops");
        // The extremes.
        assert!(MissingRate::parse("1").expect("valid").tolerates(20, 20));
        assert!(MissingRate::parse("0").expect("valid").tolerates(0, 20));
        assert!(!MissingRate::parse("0").expect("valid").tolerates(1, 20));
    }

    #[test]
    fn render_round_trips_and_is_shortest() {
        for text in ["0", "1", "0.9", "0.8", "0.7", "0.25", "0.123456789", "1.0"] {
            let rate = MissingRate::parse(text).expect("valid");
            let rendered = rate.render();
            assert_eq!(
                MissingRate::parse(&rendered).expect("renders back"),
                rate,
                "{text} -> {rendered}"
            );
        }
        assert_eq!(MissingRate::parse("0.20").expect("valid").render(), "0.2");
        assert_eq!(MissingRate::parse("1.0").expect("valid").render(), "1");
    }

    /// The withdrawn flags' error text needs `1 − X` exactly:
    /// `--max-missing 0.9` → `--geno 0.1`, not `0.09999999999999998`.
    #[test]
    fn complement_is_exact_decimal() {
        assert_eq!(MissingRate::parse("0.9").expect("valid").complement().render(), "0.1");
        assert_eq!(MissingRate::parse("0.8").expect("valid").complement().render(), "0.2");
        assert_eq!(MissingRate::parse("0").expect("valid").complement().render(), "1");
    }

    #[test]
    fn out_of_range_is_rejected() {
        assert!(MissingRate::parse("1.5").is_err());
        assert!(MissingRate::parse("-0.1").is_err());
        assert!(MissingRate::parse("nan").is_err());
        assert!(MissingRate::parse("x").is_err());
    }

    #[test]
    fn serde_round_trips() {
        let rate = MissingRate::parse("0.1").expect("valid");
        let json = serde_json::to_string(&rate).expect("serialize");
        assert_eq!(json, "0.1");
        assert_eq!(serde_json::from_str::<MissingRate>(&json).expect("parse"), rate);
        assert_eq!(
            serde_json::from_str::<MissingRate>("\"0.85\"").expect("parse"),
            MissingRate::parse("0.85").expect("valid")
        );
    }
}

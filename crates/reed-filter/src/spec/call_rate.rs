//! The call-rate threshold type (`doc_v4/02-cli-semantics.md` §1, §4).
//!
//! Every call-rate flag is a **keep condition** stated in the direction of
//! vcftools' `--max-missing`: the value is the fraction that must be called
//! (`--max-missing 0.8` keeps sites with call rate ≥ 0.8). The threshold is
//! stored as the rational `q / 10^9` so the boundary comparison is an exact
//! integer one: `1 − 0.8` is not representable in binary floating point, and
//! a decision written as `missing/N ≤ 1 − X` would drop the boundary site
//! that both vcftools and plink2 keep.

use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::FilterError;

/// The denominator of the rational representation (10^9 covers every
/// decimal a user writes).
pub const CALL_RATE_DENOMINATOR: u64 = 1_000_000_000;

/// A "keep if the called fraction is at least X" threshold, stored as
/// `q / CALL_RATE_DENOMINATOR`.
///
/// `parse` quantises the decimal to `q = round(X · 10^9)`; [`Self::keeps`]
/// then decides `called · 10^9 ≥ q · total` in `u128`, which is exact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CallRate {
    q: u64,
}

impl CallRate {
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
            .map_err(|_| FilterError::Spec(format!("call rate is not a number: {text:?}")))?;
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
                "call rate out of range [0,1]: {value}"
            )));
        }
        // `value · 10^9` is within ±1e-7 of an integer for any decimal the
        // user writes; the rounding is what makes the form exact.
        let q = (value * CALL_RATE_DENOMINATOR as f64).round();
        Ok(Self { q: q as u64 })
    }

    /// The numerator (`q` in `q / 10^9`).
    #[must_use]
    pub const fn q(self) -> u64 {
        self.q
    }

    /// Whether `called` out of `total` satisfies the threshold:
    /// `called · 10^9 ≥ q · total` (exact `u128` arithmetic).
    ///
    /// `total == 0` is vacuously true (`0 ≥ 0`).
    #[must_use]
    pub fn keeps(self, called: u64, total: u64) -> bool {
        u128::from(called) * u128::from(CALL_RATE_DENOMINATOR)
            >= u128::from(self.q) * u128::from(total)
    }

    /// `1 − X` in the same exact form: the value a removed drop-direction
    /// flag (`--geno 0.1`) translates to (`--max-missing 0.9`).
    #[must_use]
    pub const fn complement(self) -> Self {
        Self { q: CALL_RATE_DENOMINATOR - self.q }
    }

    /// The shortest decimal that parses back to the same threshold
    /// (round-trip for `to_args`, CLI and canonical JSON).
    #[must_use]
    pub fn render(self) -> String {
        let int = self.q / CALL_RATE_DENOMINATOR;
        let mut frac = format!("{:09}", self.q % CALL_RATE_DENOMINATOR);
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

impl std::fmt::Display for CallRate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.render())
    }
}

impl Serialize for CallRate {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        // A JSON number keeps the canonical spec readable; deserialisation
        // re-quantises to the same `q`.
        serializer.serialize_f64(self.render().parse().expect("invariant: rendered as a decimal"))
    }
}

impl<'de> Deserialize<'de> for CallRate {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct RateVisitor;
        impl Visitor<'_> for RateVisitor {
            type Value = CallRate;

            fn expecting(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.write_str("a call rate (number or decimal string) in [0,1]")
            }

            fn visit_f64<E: de::Error>(self, v: f64) -> Result<CallRate, E> {
                CallRate::from_f64(v).map_err(E::custom)
            }

            fn visit_i64<E: de::Error>(self, v: i64) -> Result<CallRate, E> {
                self.visit_f64(v as f64)
            }

            fn visit_u64<E: de::Error>(self, v: u64) -> Result<CallRate, E> {
                self.visit_f64(v as f64)
            }

            fn visit_str<E: de::Error>(self, v: &str) -> Result<CallRate, E> {
                CallRate::parse(v).map_err(E::custom)
            }
        }
        deserializer.deserialize_any(RateVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_boundary_decisions() {
        // 02 §7.6: X·2N is an integer with N = 10, X = 0.8, missing = 2.
        // The float form `missing/N ≤ 1 − X` would drop this site
        // (1 − 0.8 = 0.19999999999999996 < 0.2).
        let rate = CallRate::parse("0.8").expect("valid");
        assert!(rate.keeps(18, 20), "18/20 ≥ 0.8 keeps");
        // X = 0.7 with 1 − X = 0.30000000000000004: missing = 3 keeps.
        let rate = CallRate::parse("0.7").expect("valid");
        assert!(rate.keeps(17, 20), "17/20 ≥ 0.7 keeps");
        assert!(!rate.keeps(13, 20), "13/20 < 0.7 drops");
        // The extremes.
        assert!(CallRate::parse("0").expect("valid").keeps(0, 20));
        assert!(CallRate::parse("1").expect("valid").keeps(20, 20));
        assert!(!CallRate::parse("1").expect("valid").keeps(19, 20));
    }

    #[test]
    fn render_round_trips_and_is_shortest() {
        for text in ["0", "1", "0.9", "0.8", "0.7", "0.25", "0.123456789", "1.0"] {
            let rate = CallRate::parse(text).expect("valid");
            let rendered = rate.render();
            assert_eq!(
                CallRate::parse(&rendered).expect("renders back"),
                rate,
                "{text} -> {rendered}"
            );
        }
        assert_eq!(CallRate::parse("0.90").expect("valid").render(), "0.9");
        assert_eq!(CallRate::parse("1.0").expect("valid").render(), "1");
        assert_eq!(CallRate::parse("0.1").expect("valid").render(), "0.1");
    }

    #[test]
    fn complement_is_exact_decimal() {
        assert_eq!(CallRate::parse("0.1").expect("valid").complement().render(), "0.9");
        assert_eq!(CallRate::parse("0.8").expect("valid").complement().render(), "0.2");
        assert_eq!(CallRate::parse("0").expect("valid").complement().render(), "1");
    }

    #[test]
    fn out_of_range_is_rejected() {
        assert!(CallRate::parse("1.5").is_err());
        assert!(CallRate::parse("-0.1").is_err());
        assert!(CallRate::parse("nan").is_err());
        assert!(CallRate::parse("x").is_err());
    }

    #[test]
    fn serde_round_trips() {
        let rate = CallRate::parse("0.9").expect("valid");
        let json = serde_json::to_string(&rate).expect("serialize");
        assert_eq!(json, "0.9");
        assert_eq!(serde_json::from_str::<CallRate>(&json).expect("parse"), rate);
        assert_eq!(serde_json::from_str::<CallRate>("\"0.85\"").expect("parse"), CallRate::parse("0.85").expect("valid"));
    }
}

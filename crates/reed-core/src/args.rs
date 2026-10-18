//! The `--flag value` renderer shared by the two `to_args` text contracts
//! (`FilterSpec::to_args`, `PopgenRequest::to_args`). Kept here (L0) so both
//! renderers use one implementation instead of near-copies.

/// Appends a flag plus its value to `args`.
pub fn push_flag(args: &mut Vec<String>, name: &str, value: impl std::fmt::Display) {
    args.push(name.to_string());
    args.push(value.to_string());
}

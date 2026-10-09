//! INFO fragment iteration (`docs/spec/architecture.md` §6.4).
//!
//! Typed access by `TagId` lands with the header tag dictionary (M5 filter
//! head predicates); for now this is the zero-copy fragment view.

/// Iterates INFO fragments as `(key, value)`; a Flag yields `value = None`.
/// Empty INFO (`.`) yields nothing.
pub fn info_fragments(info: &[u8]) -> impl Iterator<Item = (&[u8], Option<&[u8]>)> {
    info.split(|&b| b == b';')
        .filter(|f| !f.is_empty() && f != b".")
        .map(|fragment| match memchr::memchr(b'=', fragment) {
            Some(eq) => (&fragment[..eq], Some(&fragment[eq + 1..])),
            None => (fragment, None),
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fragments_split_and_flags() {
        let fragments: Vec<_> = info_fragments(b"DP=5;SOMETHING;MQ=58.2").collect();
        assert_eq!(fragments[0], (b"DP".as_slice(), Some(b"5".as_slice())));
        assert_eq!(fragments[1], (b"SOMETHING".as_slice(), None));
        assert_eq!(fragments[2], (b"MQ".as_slice(), Some(b"58.2".as_slice())));
    }

    #[test]
    fn missing_info_is_empty() {
        assert_eq!(info_fragments(b".").count(), 0);
        assert_eq!(info_fragments(b"").count(), 0);
    }
}

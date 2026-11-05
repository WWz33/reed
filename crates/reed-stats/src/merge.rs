//! The merge-shape invariant (doc_v25 V25-05).
//!
//! Every accumulator's `merge` is component-wise addition over vectors whose
//! shapes are fixed by the request (projection size, sample count, bin count,
//! region count). A plain `zip` over `iter_mut()` silently truncates a
//! longer right-hand side, so a caller bug (two partials built for different
//! shapes) used to produce a short result with wrong numbers in release
//! builds — V23-28 proved the failure mode is real. `merge_zip_eq` makes
//! shape equality a checked precondition on every build.

/// Merges `src` into `dst` component-wise, panicking when the shapes differ.
///
/// `ctx` names the accumulator in the panic message (release builds check
/// too: a wrong merge must never silently drop tail components).
pub fn merge_zip_eq<T>(ctx: &'static str, dst: &mut [T], src: &[T], add: impl Fn(&mut T, &T)) {
    assert_eq!(
        dst.len(),
        src.len(),
        "{ctx}: merge operands have different shapes ({} vs {})",
        dst.len(),
        src.len()
    );
    for (d, s) in dst.iter_mut().zip(src) {
        add(d, s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shapes_must_match_in_every_build() {
        let mut dst = vec![1_u64, 2];
        merge_zip_eq("test dst", &mut dst, &[10, 20], |d, s| *d += s);
        assert_eq!(dst, [11, 22]);

        let err = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut short = vec![1_u64];
            merge_zip_eq("projection size", &mut short, &[1, 2], |d, s| *d += s);
        }))
        .expect_err("a short shape must panic");
        let msg = err
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or_default();
        assert!(msg.contains("projection size"), "{msg}");
        assert!(msg.contains("1 vs 2"), "{msg}");
    }
}

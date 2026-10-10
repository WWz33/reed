//! Fixed-length bit set over `u64` words (`docs/spec/architecture.md` §5.3).
//!
//! Supersedes the legacy `Bits` of the frozen 0.2 oracle with word-level
//! operations so
//! counting kernels can combine masks and planes without per-bit loops.
//!
//! Invariants: trailing padding bits (indices `len..words.len()*64`) are
//! always zero, in every public view of the words.

/// A dense, fixed-length bit vector with word-level operations.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Bitset {
    words: Box<[u64]>,
    len: usize,
}

impl Bitset {
    /// Creates a zeroed bit set of `len` bits.
    #[must_use]
    pub fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)].into_boxed_slice(),
            len,
        }
    }

    /// Number of bits.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Whether there are no bits.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Number of `u64` words backing this bit set.
    #[must_use]
    pub const fn n_words(&self) -> usize {
        self.len.div_ceil(64)
    }

    /// Whether the bit set is all zeros.
    #[must_use]
    pub fn is_clear(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }

    /// Returns bit `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= len` — an out-of-range index is a caller invariant
    /// violation, not an input-data condition.
    #[must_use]
    pub fn get(&self, i: usize) -> bool {
        assert!(i < self.len, "invariant: bit index out of range");
        self.words[i / 64] >> (i % 64) & 1 != 0
    }

    /// Sets bit `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= len`.
    pub fn set(&mut self, i: usize, value: bool) {
        assert!(i < self.len, "invariant: bit index out of range");
        let mask = 1u64 << (i % 64);
        if value {
            self.words[i / 64] |= mask;
        } else {
            self.words[i / 64] &= !mask;
        }
    }

    /// Read-only view of the backing words (padding bits are zero).
    #[must_use]
    pub fn words(&self) -> &[u64] {
        &self.words
    }

    /// Mutable view of the backing words.
    ///
    /// The caller must not set padding bits (indices `>= len`); call
    /// [`Bitset::clear_padding`] if there is any doubt.
    pub fn words_mut(&mut self) -> &mut [u64] {
        &mut self.words
    }

    /// Clears every bit (word-at-a-time; padding stays zero per the
    /// invariant).
    pub fn clear(&mut self) {
        self.words.fill(0);
    }

    /// Forces trailing padding bits to zero, restoring the type invariant.
    pub fn clear_padding(&mut self) {
        let used = self.len % 64;
        if used != 0 {
            let last = self.words.len() - 1;
            self.words[last] &= (1 << used) - 1;
        }
    }

    /// Number of set bits.
    #[must_use]
    pub fn count_ones(&self) -> u64 {
        self.words.iter().map(|w| u64::from(w.count_ones())).sum()
    }

    /// Word-wise AND; errors on length mismatch.
    pub fn and(&self, other: &Bitset) -> Result<Bitset, crate::error::CoreError> {
        self.zip(other, |a, b| a & b)
    }

    fn zip(
        &self,
        other: &Bitset,
        f: impl Fn(u64, u64) -> u64,
    ) -> Result<Bitset, crate::error::CoreError> {
        if self.len != other.len {
            return Err(crate::error::CoreError::LengthMismatch {
                expected: self.len,
                got: other.len,
            });
        }
        let words = self
            .words
            .iter()
            .zip(other.words.iter())
            .map(|(&a, &b)| f(a, b))
            .collect::<Vec<_>>()
            .into_boxed_slice();
        Ok(Bitset {
            words,
            len: self.len,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_get_roundtrip_at_word_boundaries() {
        let mut b = Bitset::new(130);
        assert_eq!(b.n_words(), 3);
        b.set(0, true);
        b.set(64, true);
        b.set(129, true);
        assert!(b.get(0) && b.get(64) && b.get(129));
        assert!(!b.get(1) && !b.get(128));
        assert_eq!(b.count_ones(), 3);
        b.set(64, false);
        assert_eq!(b.count_ones(), 2);
    }

    #[test]
    fn padding_bits_stay_zero() {
        let mut b = Bitset::new(65);
        b.set(0, true);
        b.set(64, true);
        // The last word holds one used bit and 63 padding bits.
        assert_eq!(b.words()[1], 1);
        b.words_mut()[1] = u64::MAX;
        b.clear_padding();
        assert_eq!(b.words()[1], 1);
    }

    #[test]
    fn word_ops_match_scalar_truth() {
        let len = 100;
        let mut a = Bitset::new(len);
        let mut b = Bitset::new(len);
        for i in (0..len).step_by(3) {
            a.set(i, true);
        }
        for i in (0..len).step_by(5) {
            b.set(i, true);
        }
        let both = a.and(&b).unwrap();
        for i in 0..len {
            assert_eq!(both.get(i), a.get(i) && b.get(i));
        }
        // Multiples of 15 below 100, including 0.
        assert_eq!(both.count_ones(), 7);
    }

    #[test]
    fn length_mismatch_is_an_error() {
        let a = Bitset::new(64);
        let b = Bitset::new(65);
        assert!(a.and(&b).is_err());
    }
}

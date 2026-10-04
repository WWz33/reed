//! Minimal bit vector over `u64` words.

/// A dense bit vector.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Bits {
    words: Vec<u64>,
    len: usize,
}

impl Bits {
    /// Creates a zeroed bit vector of `len` bits.
    pub fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
            len,
        }
    }

    /// Creates a bit vector from booleans.
    pub fn from_bools(values: &[bool]) -> Self {
        let mut bits = Self::new(values.len());
        for (i, &v) in values.iter().enumerate() {
            bits.set(i, v);
        }
        bits
    }

    /// Number of bits.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether there are no bits.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns bit `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= len`.
    pub fn get(&self, i: usize) -> bool {
        assert!(i < self.len, "bit index out of range");
        self.words[i / 64] >> (i % 64) & 1 != 0
    }

    /// Sets bit `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= len`.
    pub fn set(&mut self, i: usize, value: bool) {
        assert!(i < self.len, "bit index out of range");
        let mask = 1u64 << (i % 64);
        if value {
            self.words[i / 64] |= mask;
        } else {
            self.words[i / 64] &= !mask;
        }
    }

    /// Number of set bits.
    pub fn count_ones(&self) -> u64 {
        self.words.iter().map(|w| u64::from(w.count_ones())).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_get_count() {
        let mut b = Bits::new(130);
        assert_eq!(b.len(), 130);
        assert!(!b.get(0));
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
    fn from_bools_roundtrip() {
        let vals = [true, false, true, true, false, false, true];
        let b = Bits::from_bools(&vals);
        for (i, &v) in vals.iter().enumerate() {
            assert_eq!(b.get(i), v);
        }
    }
}

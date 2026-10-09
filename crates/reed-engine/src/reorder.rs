//! Bounded reorder buffer (`docs/spec/architecture.md` §7.4).
//!
//! In-flight chunk ids form a contiguous window `[expected, expected +
//! capacity)`, so a fixed ring indexed by `id % capacity` suffices — the
//! unbounded `BTreeMap` of the legacy engine is gone (issue A14).

// Casts in this file are bounded by construction: Ring slots are indexed by chunk ids, which fit `usize` on every supported target.
#![allow(clippy::cast_possible_truncation)]

use reed_core::ChunkId;

/// Fixed-capacity ring delivering results in chunk order.
#[derive(Debug)]
pub struct ReorderBuffer<T> {
    ring: Vec<Option<T>>,
    expected: u64,
}

impl<T> ReorderBuffer<T> {
    /// A buffer holding up to `capacity` in-flight results.
    ///
    /// # Panics
    ///
    /// Panics when `capacity == 0`.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        assert!(capacity > 0, "invariant: reorder capacity is positive");
        Self {
            ring: (0..capacity).map(|_| None).collect(),
            expected: 0,
        }
    }

    /// Capacity (the ticket count).
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.ring.len()
    }

    /// The next chunk id to be delivered.
    #[must_use]
    pub const fn expected(&self) -> u64 {
        self.expected
    }

    /// Whether the next expected result has already arrived.
    #[must_use]
    pub fn ready(&self) -> bool {
        self.ring[(self.expected as usize) % self.ring.len()].is_some()
    }

    /// Inserts a result for `id`.
    ///
    /// # Panics
    ///
    /// Panics when `id` is outside the window `[expected, expected +
    /// capacity)` or the slot is occupied — both are engine invariants.
    pub fn insert(&mut self, id: ChunkId, value: T) {
        let capacity = self.ring.len() as u64;
        let offset = id.get().wrapping_sub(self.expected);
        assert!(
            offset < capacity,
            "invariant: chunk {} outside the reorder window [{}..{})",
            id.get(),
            self.expected,
            self.expected + capacity
        );
        let capacity = self.ring.len();
        let slot = &mut self.ring[(id.get() as usize) % capacity];
        assert!(slot.is_none(), "invariant: duplicate delivery of chunk {}", id.get());
        *slot = Some(value);
    }

    /// Takes the next in-order result, if it has arrived.
    pub fn take_ready(&mut self) -> Option<(ChunkId, T)> {
        let capacity = self.ring.len();
        let slot = self.ring[(self.expected as usize) % capacity].take()?;
        let id = ChunkId(self.expected);
        self.expected += 1;
        Some((id, slot))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delivers_in_order_despite_arrival_order() {
        let mut ring = ReorderBuffer::new(4);
        ring.insert(ChunkId(2), "c");
        ring.insert(ChunkId(1), "b");
        // Chunks 1 and 2 arrived; chunk 0 is still outstanding.
        assert!(!ring.ready());
        ring.insert(ChunkId(0), "a");
        assert!(ring.ready());
        assert_eq!(ring.take_ready().unwrap().0.get(), 0);
        assert_eq!(ring.take_ready().unwrap().0.get(), 1);
        assert_eq!(ring.take_ready().unwrap().0.get(), 2);
        assert!(!ring.ready());
    }

    #[test]
    fn wraps_the_ring_window() {
        let mut ring = ReorderBuffer::new(2);
        ring.insert(ChunkId(0), "a");
        ring.insert(ChunkId(1), "b");
        assert_eq!(ring.take_ready().unwrap(), (ChunkId(0), "a"));
        // Chunk 2 reuses slot 0.
        ring.insert(ChunkId(2), "c");
        assert_eq!(ring.take_ready().unwrap(), (ChunkId(1), "b"));
        assert_eq!(ring.take_ready().unwrap(), (ChunkId(2), "c"));
    }

    #[test]
    #[should_panic(expected = "outside the reorder window")]
    fn rejects_ids_beyond_the_window() {
        let mut ring = ReorderBuffer::new(2);
        ring.insert(ChunkId(5), "far");
    }

    #[test]
    #[should_panic(expected = "duplicate delivery")]
    fn rejects_duplicates() {
        let mut ring = ReorderBuffer::new(2);
        ring.insert(ChunkId(0), "a");
        ring.insert(ChunkId(0), "a again");
    }
}

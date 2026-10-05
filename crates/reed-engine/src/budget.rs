//! Resource budgets (`docs/architecture-spec.md` §7.3, §7.7).

use std::path::Path;

use crate::error::EngineError;

/// In-flight chunks per worker (§7.4): the ticket window capacity is
/// `TICKETS_PER_WORKER × workers`.
pub const TICKETS_PER_WORKER: usize = 2;

/// The single compute-resource knob: total threads reed may use.
///
/// Layout: 1 reader + 1 reducer/writer + `workers()` compute workers.
/// `threads <= 2` selects the sequential mode (§7.3), whose output must be
/// byte-identical to the parallel mode (locked by tests).
#[derive(Clone, Copy, Debug)]
pub struct ThreadBudget {
    /// Total threads (default: all available cores).
    pub threads: usize,
    /// In-memory upper bound for in-flight chunk data (§7.7).
    pub mem: MemBudget,
}

impl ThreadBudget {
    /// Budget using all available cores and the default memory budget.
    #[must_use]
    pub fn all() -> Self {
        Self {
            threads: std::thread::available_parallelism()
                .map_or(1, std::num::NonZeroUsize::get),
            mem: MemBudget::default_half_memory(),
        }
    }

    /// Compute-worker count (the rest of the budget goes to reader/reducer).
    #[must_use]
    pub fn workers(&self) -> usize {
        self.threads.saturating_sub(2).max(1)
    }

    /// Ticket-window capacity: the hard upper bound on in-flight chunks.
    #[must_use]
    pub fn tickets(&self) -> usize {
        TICKETS_PER_WORKER * self.workers()
    }

    /// Whether the run takes the single-threaded sequential path.
    #[must_use]
    pub const fn is_sequential(&self) -> bool {
        self.threads <= 2
    }
}

impl Default for ThreadBudget {
    fn default() -> Self {
        Self::all()
    }
}

/// In-memory upper bound for in-flight chunk data (§7.7, §2 principle 11).
#[derive(Clone, Copy, Debug)]
pub struct MemBudget {
    /// Maximum in-flight compressed bytes.
    pub bytes: u64,
}

impl MemBudget {
    /// An explicit byte budget.
    #[must_use]
    pub const fn from_bytes(bytes: u64) -> Self {
        Self { bytes }
    }

    /// Default: 50% of physical memory (Linux; other platforms: unlimited).
    #[must_use]
    pub fn default_half_memory() -> Self {
        Self {
            bytes: physical_memory_bytes().map_or(u64::MAX, |total| total / 2),
        }
    }

    /// Checks an observed in-flight peak against the budget.
    ///
    /// # Errors
    ///
    /// [`EngineError::MemBudgetExceeded`] when `peak_bytes` exceeds the
    /// budget.
    pub fn check(&self, peak_bytes: u64) -> Result<(), EngineError> {
        if peak_bytes > self.bytes {
            return Err(EngineError::MemBudgetExceeded {
                peak_bytes,
                budget_bytes: self.bytes,
            });
        }
        Ok(())
    }
}

/// Physical memory in bytes (Linux `/proc/meminfo`), `None` if unknown.
#[must_use]
fn physical_memory_bytes() -> Option<u64> {
    let info = std::fs::read_to_string(Path::new("/proc/meminfo")).ok()?;
    for line in info.lines() {
        if let Some(rest) = line.strip_prefix("MemTotal:") {
            let kib: u64 = rest.trim().trim_end_matches(" kB").trim().parse().ok()?;
            return Some(kib * 1024);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn threads(n: usize) -> ThreadBudget {
        ThreadBudget {
            threads: n,
            mem: MemBudget::from_bytes(u64::MAX),
        }
    }

    #[test]
    fn thread_layout_follows_the_spec() {
        assert!(threads(2).is_sequential());
        assert!(threads(1).is_sequential());
        assert_eq!(threads(2).workers(), 1);
        assert_eq!(threads(56).workers(), 54);
        assert_eq!(threads(56).tickets(), 108);
    }

    #[test]
    fn mem_budget_rejects_overages() {
        let budget = MemBudget::from_bytes(100);
        assert!(budget.check(99).is_ok());
        let err = budget.check(101).unwrap_err();
        assert!(matches!(err, EngineError::MemBudgetExceeded { .. }));
    }
}

#![allow(unsafe_code)] // the counting GlobalAlloc below is test-only

//! Allocation-count regression for the per-record INFO validation
//! (doc_v26 review item 5).
//!
//! `validate_info` runs once per record from `classify_line`. The first
//! version built a `Vec<&[u8]>` for the seen keys, so every record with INFO
//! paid a heap allocation; its doc comment claimed "allocation-free", which
//! was false. This test counts allocations across many calls and requires
//! zero after the first, so a reintroduced `Vec` fails here.
//!
//! The counting allocator is global to this test binary; this file contains
//! exactly one test so no sibling test can perturb the counter (the same
//! pattern as `reed/tests/alloc_per_record.rs`).

use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

use reed_io::vcf::info::validate_info;

struct CountingAlloc;

static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOC: CountingAlloc = CountingAlloc;

#[test]
fn validate_info_is_allocation_free() {
    // Warm up anything the first call might lazily initialise (formatting,
    // panic infrastructure) so the measured window is steady state.
    let info = b"DP=5;AC=1,2;AN=4;SOMETHING;MQ=58.2;HWE=0.1";
    let duplicate = b"DP=5;DP=6";
    validate_info(info).expect("invariant: unique tags");
    let _ = validate_info(duplicate);

    let before = ALLOCATIONS.load(Ordering::Relaxed);
    for _ in 0..2000 {
        validate_info(info).expect("invariant: unique tags");
        assert!(validate_info(duplicate).is_err());
    }
    let after = ALLOCATIONS.load(Ordering::Relaxed);
    assert_eq!(
        after - before,
        0,
        "validate_info must not allocate; saw {} allocations over 4000 calls",
        after - before
    );
}

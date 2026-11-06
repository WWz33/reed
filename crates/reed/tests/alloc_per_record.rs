#![allow(unsafe_code)] // the counting GlobalAlloc below is test-only

//! Allocation-count regression for the filter's per-record path (doc_v24
//! V24-07).
//!
//! The test measures the *slope* — allocations of a 2N-record run minus an
//! N-record run, divided by N — so one-time setup (plan, buffers, threads)
//! cancels out. A per-record `Vec` reintroduced anywhere in the hot path
//! shows up as a slope jump.
//!
//! The counting allocator is global to this test binary; this file contains
//! exactly one test so no sibling test can perturb the counter.

use std::alloc::{GlobalAlloc, Layout, System};
use std::fmt::Write as _;
use std::sync::atomic::{AtomicUsize, Ordering};

use reed::{FilterOptions, run_filter};
use reed::spec::FilterSpec;

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
static ALLOCATOR: CountingAlloc = CountingAlloc;

/// Maximum allocations per record for the single-pass filter. V24-07
/// measured 0.001/record (one allocation across the 2,000-record delta)
/// after hoisting the `auto` list and borrowing the stale-tag list; 0.1
/// leaves room for allocator/test noise while still catching a
/// per-record `Vec`.
const MAX_ALLOCS_PER_RECORD: f64 = 0.1;

fn write_fixture(path: &std::path::Path, n_records: usize) {
    let mut text = String::from(
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=chr1,length=10000000>\n\
         ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
         ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
         ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3\n",
    );
    for i in 0..n_records {
        let pos = i + 1;
        let _ = writeln!(
            text,
            "chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT:DP:AD\t0/1:20:10,10\t0/0:30:30,0\t\
             1/1:25:0,25\t0/1:18:9,9"
        );
    }
    std::fs::write(path, text).expect("invariant: write fixture");
}

/// Runs one plain-VCF filter and returns the allocation count.
fn allocations_for(dir: &std::path::Path, tag: &str, n_records: usize) -> usize {
    let input = dir.join(format!("in-{tag}.vcf"));
    let output = dir.join(format!("out-{tag}.vcf"));
    write_fixture(&input, n_records);
    let spec = FilterSpec::default();
    let opts = FilterOptions {
        threads: 1,
        output: Some(output),
        no_version: true,
        no_command: true,
        ..FilterOptions::default()
    };
    // Warm the allocator's own state once per path so the first run's lazy
    // initialisation does not count as per-record work.
    let before = ALLOCATIONS.load(Ordering::Relaxed);
    run_filter(Some(&input), &spec, &opts).expect("invariant: filter runs");
    ALLOCATIONS.load(Ordering::Relaxed) - before
}

#[test]
fn single_pass_allocations_per_record_are_bounded() {
    let dir = std::env::temp_dir().join(format!("reed-v24-alloc-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let n = 2_000_usize;
    let small = allocations_for(&dir, "small", n);
    let large = allocations_for(&dir, "large", 2 * n);
    let per_record = (large.saturating_sub(small)) as f64 / n as f64;
    eprintln!(
        "allocations: small={small} large={large} slope={per_record:.3}/record (limit {MAX_ALLOCS_PER_RECORD})"
    );
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        per_record <= MAX_ALLOCS_PER_RECORD,
        "{per_record:.3} allocations/record exceeds the {MAX_ALLOCS_PER_RECORD} budget"
    );
}

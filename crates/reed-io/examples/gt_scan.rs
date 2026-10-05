//! Parallel GT-only scan over a BGZF VCF: the BAPP end-to-end driver used
//! for the M2 acceptance gate (`docs/architecture-spec.md` §19: full-file
//! GT-only scan ≤ 20 s on 56 threads).
//!
//! This is also the prototype of the M3 engine skeleton: reader → bounded
//! chunk channel → workers (inflate + scan fused) → ordered reducer that
//! chain-stitches seam records.
//!
//! Usage: `gt_scan <input.vcf.gz> [--threads N]` (default: all cores).

use std::io::BufReader;
use std::time::Instant;

use crossbeam_channel::bounded;
use reed_core::genotype::{GenotypePolicies, PlaneBlock, PolicyCounters};
use reed_io::chunk::bgzf::BgzfChunkReader;
use reed_io::chunk::inflate_chunk;
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::sample_scan::SampleScanner;
use reed_io::Locus;

#[allow(clippy::too_many_lines)]
fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let Some(input) = args.next() else {
        eprintln!("usage: gt_scan <input.vcf.gz> [--threads N]");
        std::process::exit(2);
    };
    let threads = match (args.next().as_deref(), args.next().as_deref()) {
        (Some("--threads"), Some(n)) => n.parse::<usize>()?,
        _ => std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
    };
    let n_workers = threads.saturating_sub(2).max(1);

    let t_open = Instant::now();
    let file = std::fs::File::open(&input)?;
    let mut reader = BgzfChunkReader::new(BufReader::with_capacity(1 << 20, file));

    // Header: serial inflate of the leading blocks (BAPP step 1).
    let mut header_inflater = reed_io::bgzf::inflate::Inflater::new();
    let header = reader.scan_header(&mut header_inflater)?;
    let n_samples = header.n_samples;
    let t_header = t_open.elapsed();
    eprintln!(
        "header: {} bytes, n_samples={} ({:.2?})",
        header.text.len(),
        n_samples,
        t_header
    );

    let policies = GenotypePolicies::default();
    let (chunk_tx, chunk_rx) = bounded::<reed_io::chunk::ChunkBytes>(2 * n_workers);
    let (result_tx, result_rx) =
        bounded::<ChunkResult>(2 * n_workers);

    let t_scan = Instant::now();
    std::thread::scope(|scope| -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
        // Reader: cut chunks at block boundaries (serial, header-parse only).
        scope.spawn({
            let chunk_tx = chunk_tx.clone();
            move || {
                while let Ok(Some(chunk)) = reader.next_chunk() {
                    if chunk_tx.send(chunk).is_err() {
                        break;
                    }
                }
            }
        });
        drop(chunk_tx);

        // Workers: inflate + fixed columns + GT planes, fused per chunk.
        for _ in 0..n_workers {
            scope.spawn({
                let chunk_rx = chunk_rx.clone();
                let result_tx = result_tx.clone();
                move || {
                    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
                    let mut scratch = Vec::new();
                    let mut scanner = SampleScanner::new(&[]);
                    let mut planes = PlaneBlock::new(0, n_samples, false);
                    let mut counters = PolicyCounters::default();
                    while let Ok(chunk) = chunk_rx.recv() {
                        let result = scan_chunk(
                            &chunk,
                            &mut inflater,
                            &mut scratch,
                            &mut scanner,
                            &mut planes,
                            &mut counters,
                            &policies,
                        );
                        let out = match result {
                            Ok(out) => out,
                            Err(err) => {
                                let _ = result_tx.send(ChunkResult::Failed(err));
                                return;
                            }
                        };
                        if result_tx.send(ChunkResult::Done(out)).is_err() {
                            return;
                        }
                    }
                }
            });
        }
        drop(result_tx);

        // Reducer (this thread): collect, then walk chunks in id order and
        // chain-stitch seam fragments (BAPP steps 4–5).
        let mut partials = Vec::new();
        for result in result_rx {
            match result {
                ChunkResult::Done(partial) => partials.push(partial),
                ChunkResult::Failed(err) => return Err(Box::new(err)),
            }
        }
        partials.sort_by_key(|p| p.chunk_id);
        let mut n_records = 0_u64;
        let mut fragment: Option<Vec<u8>> = None;
        let mut seam_scanner = SampleScanner::new(&[]);
        let mut seam_planes = PlaneBlock::new(1, n_samples, false);
        let mut seam_counters = PolicyCounters::default();
        for partial in &partials {
            if let Some(head) = &partial.head {
                let mut seam = fragment.take().unwrap_or_default();
                seam.extend_from_slice(head);
                if seam.ends_with(b"\n") {
                    process_seam_line(
                        &seam,
                        &mut seam_scanner,
                        &mut seam_planes,
                        &mut seam_counters,
                        &policies,
                    )?;
                    n_records += 1;
                } else {
                    fragment = Some(seam);
                }
            }
            n_records += partial.n_records;
            fragment = match (fragment.take(), &partial.tail) {
                (Some(mut acc), Some(tail)) => {
                    acc.extend_from_slice(tail);
                    Some(acc)
                }
                (None, Some(tail)) => Some(tail.clone()),
                (acc, None) => acc,
            };
        }
        if let Some(leftover) = fragment {
            let has_record_bytes = leftover.iter().any(|&b| b != b'\n');
            if has_record_bytes {
                return Err(format!("truncated final record: {} bytes", leftover.len()).into());
            }
        }
        let elapsed = t_scan.elapsed();
        eprintln!(
            "scan: n_records={} n_samples={} workers={} wall={:.2?} ({:.1} Mrec/min)",
            n_records,
            n_samples,
            n_workers,
            elapsed,
            n_records as f64 / 1e6 / elapsed.as_secs_f64() * 60.0
        );
        Ok(())
    })?;

    Ok(())
}

struct ChunkPartial {
    chunk_id: reed_core::ChunkId,
    n_records: u64,
    head: Option<Vec<u8>>,
    tail: Option<Vec<u8>>,
}

enum ChunkResult {
    Done(ChunkPartial),
    Failed(reed_io::IoError),
}

#[allow(clippy::type_complexity)]
fn scan_chunk(
    chunk: &reed_io::chunk::ChunkBytes,
    inflater: &mut reed_io::bgzf::inflate::Inflater,
    scratch: &mut Vec<u8>,
    scanner: &mut SampleScanner,
    planes: &mut PlaneBlock,
    counters: &mut PolicyCounters,
    policies: &GenotypePolicies,
) -> Result<ChunkPartial, reed_io::IoError> {
    let view = inflate_chunk(chunk, inflater, scratch)?;
    planes.resize_sites(2_000);
    let mut site_idx = 0_usize;
    let mut n_records = 0_u64;
    for line_bytes in view.body_bytes(view.text).split(|&b| b == b'\n') {
        if line_bytes.is_empty() {
            continue;
        }
        let line = LineRef::parse(line_bytes).ok_or_else(|| reed_io::IoError::MalformedRecord {
            locus: Locus { contig_name: "?".into(), pos1: 0 },
            reason: "fewer than 9 tabs in a record line",
        })?;
        let locus = reed_io::Locus {
            contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
            pos1: 0,
        };
        let pos1 = parse_pos1(line.pos_text(), &locus)?;
        if site_idx >= planes.n_sites() {
            planes.resize_sites(planes.n_sites() * 2 + 1);
        }
        scanner.scan_record(&line, site_idx, pos1, planes, &mut [], policies, counters)?;
        site_idx += 1;
        n_records += 1;
    }
    Ok(ChunkPartial {
        chunk_id: chunk.id,
        n_records,
        head: view.head_bytes(view.text).map(Vec::from),
        tail: view.tail_bytes(view.text).map(Vec::from),
    })
}

fn process_seam_line(
    seam: &[u8],
    scanner: &mut SampleScanner,
    planes: &mut PlaneBlock,
    counters: &mut PolicyCounters,
    policies: &GenotypePolicies,
) -> Result<(), reed_io::IoError> {
    planes.resize_sites(1);
    let line = LineRef::parse(seam).ok_or_else(|| reed_io::IoError::MalformedRecord {
        locus: reed_io::Locus { contig_name: "seam".into(), pos1: 0 },
        reason: "fewer than 9 tabs in a seam record",
    })?;
    let locus = reed_io::Locus {
        contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
        pos1: 0,
    };
    let pos1 = parse_pos1(line.pos_text(), &locus)?;
    scanner.scan_record(&line, 0, pos1, planes, &mut [], policies, counters)
}

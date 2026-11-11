//! Shared integration-test fixtures.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

/// A BCF whose first record is complete and whose second record is cut by
/// three bytes, then wrapped in a valid BGZF stream with an EOF marker.
///
/// The compressed input is structurally valid; only the BCF tail is
/// truncated. This is the shape that used to publish filter/index/report
/// output before `job.finish`/`FrameChain` reported the truncation
/// (doc_v28 V28-06).
pub(crate) fn truncated_bcf(dir: &Path, name: &str) -> PathBuf {
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("../reed-io/tests/fixtures/bcf/fix.bcf");
    let bytes = std::fs::read(src).expect("invariant: BCF fixture");
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut inflated = Vec::new();
    let mut rest = bytes.as_slice();
    while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
        inflater
            .inflate_block(&rest[..total], &mut inflated)
            .expect("invariant: fixture block");
        rest = &rest[total..];
    }
    let header_len = u32::from_le_bytes(
        inflated[5..9]
            .try_into()
            .expect("invariant: BCF length prefix"),
    ) as usize;
    let mut records = &inflated[9 + header_len..];
    let first = reed_io::bcf::record_len(records).expect("invariant: first record");
    records = &records[first..];
    let second = reed_io::bcf::record_len(records).expect("invariant: second record");
    inflated.truncate(9 + header_len + first + second - 3);

    let mut deflater =
        reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level 6 is valid");
    let mut out = Vec::new();
    for part in inflated.chunks(13) {
        deflater
            .deflate_block(part, &mut out)
            .expect("invariant: block fits");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);

    std::fs::create_dir_all(dir).expect("invariant: temp dir");
    let path = dir.join(name);
    std::fs::write(&path, out).expect("invariant: write truncated BCF");
    path
}

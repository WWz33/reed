//! Streams a BCF file's records as VCF text (development harness for the
//! bcftools byte-parity check).
//!
//! Usage: bcf_dump FILE [N]

use std::io::{Read, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let path = args.next().ok_or("usage: bcf_dump FILE [N]")?;
    let n: Option<usize> = args.next().map(|v| v.parse().expect("N"));

    let file = std::fs::File::open(&path)?;
    let mut reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();

    // pending: inflated bytes not yet consumed; off: cursor into pending.
    let mut pending: Vec<u8> = Vec::new();
    let mut off = 0_usize;
    let mut dict: Option<reed_io::bcf::BcfDict> = None;
    let mut rendered = 0_usize;
    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::with_capacity(1 << 20, stdout.lock());
    let mut line = Vec::with_capacity(1 << 16);

    let mut header = [0_u8; 18];
    'stream: loop {
        if off > 0 {
            pending.drain(..off);
            off = 0;
        }
        if reader.read_exact(&mut header).is_err() {
            break;
        }
        let total = reed_io::bgzf::block::parse_block_size(&header).ok_or("not BGZF")?;
        let mut block = header.to_vec();
        block.resize(total, 0);
        reader.read_exact(&mut block[18..])?;
        inflater.inflate_block(&block, &mut pending)?;

        if dict.is_none() {
            if pending.len() < 9 {
                continue;
            }
            let l_text =
                u32::from_le_bytes([pending[5], pending[6], pending[7], pending[8]]) as usize;
            if pending.len() < 9 + l_text {
                continue;
            }
            let text = pending[9..9 + l_text].to_vec();
            dict = Some(reed_io::bcf::BcfDict::from_header(&text)?);
            off = 9 + l_text;
        }
        let dict = dict.as_ref().expect("invariant: header parsed");

        while let Some(len) = reed_io::bcf::record_len(&pending[off..]) {
            let rec = reed_io::bcf::record_ref(&pending[off..off + len])?;
            line.clear();
            reed_io::bcf::render_record(&rec, dict, &mut line)?;
            out.write_all(&line)?;
            off += len;
            rendered += 1;
            if n.is_some_and(|n| rendered >= n) {
                break 'stream;
            }
        }
    }
    out.flush()?;
    eprintln!("{rendered} records rendered");
    Ok(())
}

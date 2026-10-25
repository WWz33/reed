//! Drain-time renderers for the convert targets (doc_v13 §3.6, §7.6):
//! PHYLIP/FASTA/NEXUS/binary-NEXUS, the `--012` triple, the TASSEL hapmap
//! and the reference-driven MSA.
//!
//! Every row is written by streaming the merged per-sample column blocks
//! (no full-matrix copy beyond the columns themselves), and every byte is a
//! pure function of the merged partials, so output is thread-count
//! independent.

// Casts in this file are bounded by construction: sample indices and
// reference offsets are far below u32::MAX (the addressable sequence space).
#![allow(clippy::cast_possible_truncation)]

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;

use std::collections::HashMap;
use std::sync::Arc;

use super::PhyloFormats;
use super::job::{MsaFrameReader, MsaOutput, RowsBlob, SpilledCount, SpilledMsa, SpilledPhylo};

/// `PREFIX` + suffix appended to the full path (no extension juggling).
fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_os_string();
    s.push(suffix);
    PathBuf::from(s)
}

/// Opens a buffered writer for `path`.
fn writer(path: &Path) -> Result<BufWriter<File>, EngineError> {
    let file = File::create(path).map_err(EngineError::Io)?;
    Ok(BufWriter::new(file))
}

/// The row order: the outgroup first (when named), then the header order.
fn row_order(names: &[String], outgroup: Option<&str>) -> Vec<usize> {
    let mut order: Vec<usize> = (0..names.len()).collect();
    if let Some(outgroup) = outgroup
        && let Some(pos) = names.iter().position(|n| n == outgroup)
    {
        order.retain(|&i| i != pos);
        order.insert(0, pos);
    }
    order
}

/// The name column width (`max_name_len + 3`).
fn name_width(names: &[String]) -> usize {
    names.iter().map(String::len).max().unwrap_or(0) + 3
}

/// Writes one padded name + sequence row (PHYLIP/NEXUS layout).
fn write_padded_row(
    out: &mut impl Write,
    name: &str,
    width: usize,
    seq: &RowsBlob,
) -> Result<(), EngineError> {
    out.write_all(name.as_bytes()).map_err(EngineError::Io)?;
    let pad = width.saturating_sub(name.len());
    const SPACES: &[u8; 64] = &[b' '; 64];
    let mut left = pad;
    while left > 0 {
        let chunk = left.min(SPACES.len());
        out.write_all(&SPACES[..chunk]).map_err(EngineError::Io)?;
        left -= chunk;
    }
    seq.write_to(out)?;
    out.write_all(b"\n").map_err(EngineError::Io)?;
    Ok(())
}

/// Writes the phylogenetic alignment files for `plan.formats` plus the
/// binary NEXUS and `used_sites` outputs (v13 §3.6).
///
/// # Errors
///
/// I/O failures.
pub(crate) fn write_phylo(
    prefix: &Path,
    plan: &super::job::PhyloPlan,
    names: &[String],
    partial: Option<&SpilledPhylo>,
) -> Result<Vec<PathBuf>, EngineError> {
    let n = names.len();
    let order = row_order(names, plan.options.outgroup.as_deref());
    let width = name_width(names);
    let empty;
    let partial = match partial {
        Some(p) => p,
        None => {
            empty = SpilledPhylo {
                n_sites: 0,
                n_bin_sites: 0,
                cols: (0..n).map(|_| RowsBlob::default()).collect(),
                bin_cols: (0..n).map(|_| RowsBlob::default()).collect(),
                used: RowsBlob::default(),
            };
            &empty
        }
    };
    let mut files = Vec::new();
    let formats: PhyloFormats = plan.formats;

    if formats.phylip {
        let path = with_suffix(prefix, ".phy");
        let mut out = writer(&path)?;
        writeln!(out, "{} {}", n, partial.n_sites).map_err(EngineError::Io)?;
        for &i in &order {
            write_padded_row(&mut out, &names[i], width, &partial.cols[i])?;
        }
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    if formats.fasta {
        let path = with_suffix(prefix, ".fasta");
        let mut out = writer(&path)?;
        for &i in &order {
            out.write_all(b">").map_err(EngineError::Io)?;
            out.write_all(names[i].as_bytes())
                .map_err(EngineError::Io)?;
            out.write_all(b"\n").map_err(EngineError::Io)?;
            partial.cols[i].write_to(&mut out)?;
            out.write_all(b"\n").map_err(EngineError::Io)?;
        }
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    if formats.nexus {
        let path = with_suffix(prefix, ".nexus");
        let mut out = writer(&path)?;
        write_nexus_header(&mut out, n, partial.n_sites, b"DNA", b"N").map_err(EngineError::Io)?;
        for &i in &order {
            write_padded_row(&mut out, &names[i], width, &partial.cols[i])?;
        }
        out.write_all(b";\nEND;\n").map_err(EngineError::Io)?;
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    if plan.bin {
        let path = with_suffix(prefix, ".bin.nexus");
        let mut out = writer(&path)?;
        write_nexus_header(&mut out, n, partial.n_bin_sites, b"SNP", b"?")
            .map_err(EngineError::Io)?;
        for &i in &order {
            write_padded_row(&mut out, &names[i], width, &partial.bin_cols[i])?;
        }
        out.write_all(b";\nEND;\n").map_err(EngineError::Io)?;
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    if plan.options.write_used_sites {
        let path = with_suffix(prefix, ".used_sites.tsv");
        let mut out = writer(&path)?;
        out.write_all(b"#CHROM\tPOS\tNUM_SAMPLES\n")
            .map_err(EngineError::Io)?;
        partial.used.write_to(&mut out)?;
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    Ok(files)
}

/// The NEXUS header (vcf2phylip's exact spelling).
fn write_nexus_header(
    out: &mut impl Write,
    n: usize,
    nchar: u64,
    datatype: &[u8],
    missing: &[u8],
) -> std::io::Result<()> {
    out.write_all(b"#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=")?;
    write!(out, "{n}")?;
    out.write_all(b" NCHAR=")?;
    write!(out, "{nchar}")?;
    out.write_all(b";\n\tFORMAT DATATYPE=")?;
    out.write_all(datatype)?;
    out.write_all(b" MISSING=")?;
    out.write_all(missing)?;
    out.write_all(b" GAP=- ;\nMATRIX\n")
}

/// Writes the `--012` triple (vcftools layout, v13 §3.8): `PREFIX.012`,
/// `PREFIX.012.pos`, `PREFIX.012.indv`. The matrix row's leading index is
/// the sample's *original* header index (vcftools writes `ui`).
///
/// # Errors
///
/// I/O failures.
pub(crate) fn write_count_matrix(
    prefix: &Path,
    names: &[String],
    keep: &[bool],
    partial: Option<&SpilledCount>,
) -> Result<Vec<PathBuf>, EngineError> {
    let original: Vec<u32> = keep
        .iter()
        .enumerate()
        .filter_map(|(i, &k)| k.then_some(i as u32))
        .collect();
    let idx_path = with_suffix(prefix, ".012.indv");
    let mut indv = writer(&idx_path)?;
    for name in names {
        indv.write_all(name.as_bytes()).map_err(EngineError::Io)?;
        indv.write_all(b"\n").map_err(EngineError::Io)?;
    }
    indv.flush().map_err(EngineError::Io)?;

    let pos_path = with_suffix(prefix, ".012.pos");
    let mut pos = writer(&pos_path)?;
    if let Some(partial) = partial {
        partial.pos.write_to(&mut pos)?;
    }
    pos.flush().map_err(EngineError::Io)?;

    let matrix_path = with_suffix(prefix, ".012");
    let mut matrix = writer(&matrix_path)?;
    for (row, name_index) in original.iter().enumerate() {
        write!(matrix, "{name_index}").map_err(EngineError::Io)?;
        if let Some(partial) = partial {
            partial.cells[row].visit(|chunk| {
                for &byte in chunk {
                    matrix.write_all(b"\t")?;
                    matrix.write_all(match byte as i8 {
                        -1 => &b"-1"[..],
                        0 => &b"0"[..],
                        1 => &b"1"[..],
                        _ => &b"2"[..],
                    })?;
                }
                Ok(())
            })?;
        }
        matrix.write_all(b"\n").map_err(EngineError::Io)?;
    }
    matrix.flush().map_err(EngineError::Io)?;
    Ok(vec![matrix_path, pos_path, idx_path])
}

/// Writes the TASSEL hapmap `PREFIX.hmp.txt` (v13 §3.9): the 11 fixed
/// columns plus one column per kept sample, in header order.
///
/// # Errors
///
/// I/O failures.
pub(crate) fn write_hapmap(
    prefix: &Path,
    names: &[String],
    rows: Option<&RowsBlob>,
) -> Result<Vec<PathBuf>, EngineError> {
    let path = with_suffix(prefix, ".hmp.txt");
    let mut out = writer(&path)?;
    out.write_all(
        b"rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode",
    )
    .map_err(EngineError::Io)?;
    for name in names {
        out.write_all(b"\t").map_err(EngineError::Io)?;
        out.write_all(name.as_bytes()).map_err(EngineError::Io)?;
    }
    out.write_all(b"\n").map_err(EngineError::Io)?;
    if let Some(rows) = rows {
        rows.write_to(&mut out)?;
    }
    out.flush().map_err(EngineError::Io)?;
    Ok(vec![path])
}

/// Writes the MSA outputs (v13 §3.10): one FASTA per region (or per
/// reference contig), each record `>contig[:start-end]_sample`.
///
/// # Errors
///
/// I/O failures and a reference/output mismatch (an internal invariant).
pub(crate) fn write_msa(
    prefix: &Path,
    reference: &HashMap<Vec<u8>, Arc<[u8]>>,
    outputs: &[MsaOutput],
    names: &[String],
    partial: Option<&SpilledMsa>,
) -> Result<Vec<PathBuf>, EngineError> {
    let mut files = Vec::new();
    let empty: Vec<(Vec<u8>, RowsBlob)> = Vec::new();
    let runs: &[(Vec<u8>, RowsBlob)] = match partial {
        Some(p) => &p.runs,
        None => &empty,
    };
    for output in outputs {
        let Some(sequence) = reference.get(&output.contig) else {
            return Err(EngineError::job(format!(
                "MSA output for contig {} is not in the reference",
                String::from_utf8_lossy(&output.contig)
            )));
        };
        let start0 = output.start1 - 1;
        let end0 = output.end1;
        let blobs: Vec<&RowsBlob> = runs
            .iter()
            .filter(|(name, _)| *name == output.contig)
            .map(|(_, blob)| blob)
            .collect();
        let name_suffix = if output.explicit {
            format!("_{}-{}", output.start1, output.end1)
        } else {
            String::new()
        };
        let contig_text = String::from_utf8_lossy(&output.contig).into_owned();
        let path = with_suffix(prefix, &format!("_{contig_text}{name_suffix}.fasta"));
        let mut out = writer(&path)?;
        // One pass per sample: the frames are re-walked for every sample, so
        // no region-sized buffer is ever materialized (doc_v15 D15-7; each
        // pass re-reads the spill, which the frame reader supports by
        // rewinding).
        let mut frame = Vec::new();
        let mut chosen = Vec::new();
        for (sample, name) in names.iter().enumerate() {
            out.write_all(b">").map_err(EngineError::Io)?;
            out.write_all(contig_text.as_bytes())
                .map_err(EngineError::Io)?;
            out.write_all(format!(":{}-{}_", output.start1, output.end1).as_bytes())
                .map_err(EngineError::Io)?;
            out.write_all(name.as_bytes()).map_err(EngineError::Io)?;
            out.write_all(b"\n").map_err(EngineError::Io)?;

            let mut cursor = start0;
            let mut pending: Option<u64> = None;
            let mut last_pos: Option<u64> = None;
            for blob in &blobs {
                let mut reader = MsaFrameReader::open(blob, names.len())?;
                while let Some(pos0) = reader.next_frame_for_sample(sample, &mut frame)? {
                    // doc_v20 review of R12: the sortedness check must see
                    // every frame of the contig, including those outside the
                    // emitted window — otherwise a regression before
                    // `start0` passes and the contract depends on which
                    // region was requested.
                    if let Some(prev) = last_pos
                        && pos0 < prev
                    {
                        return Err(EngineError::job(format!(
                            "contig {} site {} follows {}: input is not contig-sorted",
                            String::from_utf8_lossy(&output.contig),
                            pos0 + 1,
                            prev + 1
                        )));
                    }
                    last_pos = Some(pos0);
                    if pos0 < start0 {
                        continue;
                    }
                    if pos0 >= end0 {
                        break;
                    }
                    // doc_v20 R12: frames must strictly ascend per contig
                    // across blobs (a revisited contig — chr1…, chr2…, chr1…
                    // — arrives as a later blob with lower positions). The
                    // old walk clamped the cursor and re-emitted the region,
                    // silently misaligning the FASTA. Only an exact duplicate
                    // of the pending site keeps the last-record-wins rule;
                    // a position below the newest seen site — pending or
                    // already flushed — is unsorted input and refuses the
                    // run, the same contract as the stats sink.
                    if let Some(prev) = pending
                        && pos0 < prev
                    {
                        return Err(EngineError::job(format!(
                            "contig {} site {} follows {}: input is not contig-sorted",
                            String::from_utf8_lossy(&output.contig),
                            pos0 + 1,
                            prev + 1
                        )));
                    }
                    if pending == Some(pos0) {
                        // Duplicate positions: the last record wins
                        // (vcf2msa's per-record overwrite).
                        chosen.clear();
                        chosen.extend_from_slice(&frame);
                        continue;
                    }
                    flush_msa_site(&mut out, sequence, &mut cursor, pending, &mut chosen)?;
                    chosen.extend_from_slice(&frame);
                    pending = Some(pos0);
                }
            }
            flush_msa_site(&mut out, sequence, &mut cursor, pending, &mut chosen)?;
            out.write_all(&sequence[cursor as usize..end0 as usize])
                .map_err(EngineError::Io)?;
            out.write_all(b"\n").map_err(EngineError::Io)?;
        }
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    Ok(files)
}

/// The MSA walk's flush step: writes the reference run between `cursor` and
/// the pending site, then the site's bytes.
fn flush_msa_site(
    out: &mut impl Write,
    sequence: &[u8],
    cursor: &mut u64,
    pending: Option<u64>,
    chosen: &mut Vec<u8>,
) -> Result<(), EngineError> {
    if let Some(pos0) = pending {
        // doc_v20 review of R12: the regressing-frame guard makes
        // `cursor <= pos0` an invariant, so the old `min` clamp was dead
        // code that would hide a future violation.
        debug_assert!(*cursor <= pos0, "MSA cursor {cursor} ran past site {pos0}");
        out.write_all(&sequence[*cursor as usize..pos0 as usize])
            .map_err(EngineError::Io)?;
        out.write_all(chosen).map_err(EngineError::Io)?;
        *cursor = pos0 + 1;
        chosen.clear();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::job::{MsaOutput, RowsBlob, SpilledMsa};

    /// One frame for a single sample, one byte wide (the on-disk layout the
    /// sink writes: `pos0: u64 LE || width: u32 LE || vals`).
    fn frame(pos0: u64, val: u8) -> RowsBlob {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&pos0.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        bytes.push(val);
        RowsBlob::Memory(bytes)
    }

    /// The `write_msa` arguments for one c1 output over a 20 bp reference.
    struct MsaCase {
        reference: HashMap<Vec<u8>, Arc<[u8]>>,
        outputs: Vec<MsaOutput>,
        names: Vec<String>,
        partial: SpilledMsa,
    }

    fn msa_case(runs: Vec<(Vec<u8>, RowsBlob)>) -> MsaCase {
        let mut reference = HashMap::new();
        reference.insert(b"c1".to_vec(), Arc::from(&b"AAAAAAAAAAAAAAAAAAAA"[..]));
        MsaCase {
            reference,
            outputs: vec![MsaOutput {
                contig: b"c1".to_vec(),
                start1: 1,
                end1: 20,
                explicit: false,
            }],
            names: vec!["s0".to_string()],
            partial: SpilledMsa { runs },
        }
    }

    /// doc_v20 R12: a revisited contig (chr1…, chr2…, chr1…) arrives as a
    /// second run blob with regressing positions. The walk used to clamp the
    /// cursor (`from = cursor.min(pos0)`) and re-emit the region, silently
    /// writing a misaligned FASTA; it must refuse the frames instead — the
    /// same "input is not contig-sorted" contract the stats sink enforces.
    #[test]
    fn r12_write_msa_refuses_regressing_frames() {
        let d = std::env::temp_dir().join(format!("reed-r12-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("invariant: temp dir");
        let case = msa_case(vec![
            (b"c1".to_vec(), frame(5, b'G')),
            (b"c1".to_vec(), frame(10, b'G')),
            (b"c1".to_vec(), frame(3, b'T')), // revisited contig
        ]);
        let err = write_msa(
            &d.join("msa"),
            &case.reference,
            &case.outputs,
            &case.names,
            Some(&case.partial),
        )
        .expect_err("regressing frames must be refused");
        assert!(err.to_string().contains("not contig-sorted"), "{err}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// doc_v20 review of R12: a regression that happens before the emitted
    /// window must still refuse — the contract is "unsorted input is an
    /// error", not "unsorted inside the requested region".
    #[test]
    fn r12_regression_outside_the_window_is_still_refused() {
        let d = std::env::temp_dir().join(format!("reed-r12-win-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("invariant: temp dir");
        let mut case = msa_case(vec![
            (b"c1".to_vec(), frame(5, b'G')),
            (b"c1".to_vec(), frame(15, b'G')),
            (b"c1".to_vec(), frame(7, b'T')), // regression below start0
        ]);
        // Emit only 11..20, so both 5 and 7 sit outside the window and the
        // regression is invisible to the old in-window check.
        case.outputs[0].start1 = 11;
        let err = write_msa(
            &d.join("msa"),
            &case.reference,
            &case.outputs,
            &case.names,
            Some(&case.partial),
        )
        .expect_err("a regression outside the window must still be refused");
        assert!(err.to_string().contains("not contig-sorted"), "{err}");
        let _ = std::fs::remove_dir_all(&d);
    }

    /// The legal seam: a duplicate position across blob boundaries keeps the
    /// last-record-wins rule, and strictly ascending frames across blobs are
    /// fine (sorted multi-chunk input produces one blob per chunk).
    #[test]
    fn r12_write_msa_keeps_seam_duplicates_and_ascending_blobs() {
        let d = std::env::temp_dir().join(format!("reed-r12-legal-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).expect("invariant: temp dir");
        let case = msa_case(vec![
            (b"c1".to_vec(), frame(5, b'G')),
            (b"c1".to_vec(), frame(10, b'G')),
            (b"c1".to_vec(), frame(10, b'C')), // seam duplicate: last wins
            (b"c1".to_vec(), frame(15, b'T')),
        ]);
        let files = write_msa(
            &d.join("msa"),
            &case.reference,
            &case.outputs,
            &case.names,
            Some(&case.partial),
        )
        .expect("sorted frames across blobs are legal");
        let text = std::fs::read_to_string(&files[0]).expect("invariant: fasta");
        let seq: Vec<&str> = text.lines().filter(|l| !l.starts_with('>')).collect();
        assert_eq!(seq, ["AAAAAGAAAACAAAATAAAA"], "one base per site: {text}");
        let _ = std::fs::remove_dir_all(&d);
    }
}

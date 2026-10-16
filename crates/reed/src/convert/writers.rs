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
use super::job::{CountPartial, HmpPartial, MsaOutput, MsaSite, PhyloPartial};

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
    seq: &[u8],
) -> std::io::Result<()> {
    out.write_all(name.as_bytes())?;
    let pad = width.saturating_sub(name.len());
    const SPACES: &[u8; 64] = &[b' '; 64];
    let mut left = pad;
    while left > 0 {
        let chunk = left.min(SPACES.len());
        out.write_all(&SPACES[..chunk])?;
        left -= chunk;
    }
    out.write_all(seq)?;
    out.write_all(b"\n")
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
    partial: Option<&PhyloPartial>,
) -> Result<Vec<PathBuf>, EngineError> {
    let n = names.len();
    let order = row_order(names, plan.options.outgroup.as_deref());
    let width = name_width(names);
    let empty;
    let partial = match partial {
        Some(p) => p,
        None => {
            empty = PhyloPartial {
                n_sites: 0,
                n_bin_sites: 0,
                cols: vec![Vec::new(); n],
                bin_cols: vec![Vec::new(); n],
                used: Vec::new(),
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
            write_padded_row(&mut out, &names[i], width, &partial.cols[i])
                .map_err(EngineError::Io)?;
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
            out.write_all(&partial.cols[i]).map_err(EngineError::Io)?;
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
            write_padded_row(&mut out, &names[i], width, &partial.cols[i])
                .map_err(EngineError::Io)?;
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
            write_padded_row(&mut out, &names[i], width, &partial.bin_cols[i])
                .map_err(EngineError::Io)?;
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
        out.write_all(&partial.used).map_err(EngineError::Io)?;
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
    partial: Option<&CountPartial>,
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
        pos.write_all(&partial.pos).map_err(EngineError::Io)?;
    }
    pos.flush().map_err(EngineError::Io)?;

    let matrix_path = with_suffix(prefix, ".012");
    let mut matrix = writer(&matrix_path)?;
    for (row, name_index) in original.iter().enumerate() {
        write!(matrix, "{name_index}").map_err(EngineError::Io)?;
        if let Some(partial) = partial {
            for cell in &partial.cells[row] {
                matrix.write_all(b"\t").map_err(EngineError::Io)?;
                match cell {
                    -1 => matrix.write_all(b"-1").map_err(EngineError::Io)?,
                    0 => matrix.write_all(b"0").map_err(EngineError::Io)?,
                    1 => matrix.write_all(b"1").map_err(EngineError::Io)?,
                    _ => matrix.write_all(b"2").map_err(EngineError::Io)?,
                }
            }
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
    partial: Option<&HmpPartial>,
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
    if let Some(partial) = partial {
        out.write_all(&partial.rows).map_err(EngineError::Io)?;
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
    partial: Option<&super::job::MsaPartial>,
) -> Result<Vec<PathBuf>, EngineError> {
    let mut files = Vec::new();
    let empty_runs: Vec<super::job::MsaRun> = Vec::new();
    let runs: &[super::job::MsaRun] = match partial {
        Some(p) => &p.runs,
        None => &empty_runs,
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
        // The sites of this contig, in position order.
        let sites: Vec<&MsaSite> = runs
            .iter()
            .filter(|r| r.contig == output.contig)
            .flat_map(|r| r.sites.iter())
            .filter(|s| s.pos0 >= start0 && s.pos0 < end0)
            .collect();
        let mut buffers: Vec<Vec<u8>> = vec![Vec::new(); names.len()];
        let mut site_idx = 0_usize;
        for pos0 in start0..end0 {
            if site_idx < sites.len() && sites[site_idx].pos0 == pos0 {
                // Duplicate positions: the last record wins (vcf2msa's
                // per-record overwrite).
                let mut chosen = sites[site_idx];
                while site_idx < sites.len() && sites[site_idx].pos0 == pos0 {
                    chosen = sites[site_idx];
                    site_idx += 1;
                }
                let width = usize::from(chosen.width);
                for (s, buffer) in buffers.iter_mut().enumerate() {
                    buffer.extend_from_slice(&chosen.vals[s * width..(s + 1) * width]);
                }
            } else {
                let base = sequence[pos0 as usize];
                for buffer in &mut buffers {
                    buffer.push(base);
                }
            }
        }
        let name_suffix = if output.explicit {
            format!("_{}-{}", output.start1, output.end1)
        } else {
            String::new()
        };
        let contig_text = String::from_utf8_lossy(&output.contig).into_owned();
        let path = with_suffix(prefix, &format!("_{contig_text}{name_suffix}.fasta"));
        let mut out = writer(&path)?;
        for (sample, buffer) in buffers.iter().enumerate() {
            out.write_all(b">").map_err(EngineError::Io)?;
            out.write_all(contig_text.as_bytes())
                .map_err(EngineError::Io)?;
            out.write_all(format!(":{}-{}_", output.start1, output.end1).as_bytes())
                .map_err(EngineError::Io)?;
            out.write_all(names[sample].as_bytes())
                .map_err(EngineError::Io)?;
            out.write_all(b"\n").map_err(EngineError::Io)?;
            out.write_all(buffer).map_err(EngineError::Io)?;
            out.write_all(b"\n").map_err(EngineError::Io)?;
        }
        out.flush().map_err(EngineError::Io)?;
        files.push(path);
    }
    Ok(files)
}

//! Allele-table and MSA helpers for the convert job
//! (split out of `job/mod.rs`, doc_v30 P4).

use super::convert_job::Malformed;
use super::{Arc, DrawKey, GtSlots, HashMap, LineRef, MsaPlan, PathBuf, ReferenceMap, iupac};

/// The raw ALT allele count (`.` counts zero; else commas + 1).
pub(super) fn raw_alt_count(lr: &LineRef<'_>) -> u32 {
    let alts = lr.alts();
    if alts == b"." {
        return 0;
    }
    memchr::memchr_iter(b',', alts).count() as u32 + 1
}

/// vcf2phylip's SNP gate on the processed allele table: REF (raw) is one
/// character and every ALT allele is one character.
pub(super) fn is_snp(lr: &LineRef<'_>, ranges: &[(u32, u32)]) -> bool {
    if lr.reference().len() != 1 {
        return false;
    }
    ranges[1..].iter().all(|&(_, len)| len == 1)
}

/// The `--012` cell for one sample (v13 §3.8, B13-8): `0`/`1`/`2` for a full
/// diploid call, `-1` for missing, half calls and haploids.
pub(super) fn cell_012(gt: &GtSlots) -> Option<i8> {
    if gt.bad || gt.overflow || gt.n == 0 {
        return None;
    }
    if gt.n != 2 {
        // Haploid (and any other ploidy): reed emits -1; vcftools emits
        // 1/0/1 (B13-8).
        return Some(-1);
    }
    match (gt.slots[0], gt.slots[1]) {
        (Some(a), Some(b)) => {
            if a > 1 || b > 1 {
                return None;
            }
            Some(match (a, b) {
                (0, 0) => 0,
                (1, 1) => 2,
                _ => 1,
            })
        }
        // Missing or half call.
        _ => Some(-1),
    }
}

/// Resolves one sample's MSA value (vcf2msa `genotype_resolve`): single
/// characters join into an IUPAC letter; indel alleles use the literal text
/// (the shortest on conflict, ties lexicographic); a fully missing genotype
/// yields the REF allele (vcf2msa's step 4).
pub(super) fn resolve_msa_value(
    gt: &GtSlots,
    ranges: &[(u32, u32)],
    flat: &[u8],
    out: &mut Vec<u8>,
) -> Result<(), Malformed> {
    let allele = |index: u32| -> Option<&[u8]> {
        let &(start, len) = ranges.get(index as usize)?;
        Some(&flat[start as usize..(start + len) as usize])
    };
    let mut called = 0;
    let mut var: Vec<u8> = Vec::new();
    let mut indel: Vec<Vec<u8>> = Vec::new();
    for slot in gt.slots.iter().take(usize::from(gt.n)) {
        let Some(index) = *slot else { continue };
        let text = allele(index).ok_or(Malformed)?;
        called += 1;
        if text.len() == 1 {
            if !var.contains(&text[0]) {
                var.push(text[0]);
            }
        } else {
            let value = text.to_vec();
            if !indel.contains(&value) {
                indel.push(value);
            }
        }
    }
    if called == 0 {
        // No called allele: the REF string (vcf2msa step 4).
        out.extend_from_slice(allele(0).ok_or(Malformed)?);
        return Ok(());
    }
    if !var.is_empty() {
        if var.len() == 1 {
            out.push(var[0]);
            return Ok(());
        }
        var.sort_unstable();
        var.dedup();
        let letter = iupac::ambig_lookup(&var).ok_or(Malformed)?;
        out.push(letter);
        return Ok(());
    }
    // Indel-only: the shortest literal, ties lexicographic.
    indel.sort();
    indel.dedup();
    let shortest = indel.iter().min_by_key(|v| v.len()).ok_or(Malformed)?;
    out.extend_from_slice(shortest);
    Ok(())
}

/// The deterministic `--resolve-iupac` draw (splitmix64 over
/// `(chunk, line, sample, seed)`, v13 §7.4.3).
pub(super) fn resolve_pick(key: DrawKey, sample: u32, n: usize, seed: u64) -> usize {
    let mut x = seed
        ^ key.chunk.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ u64::from(key.line).wrapping_mul(0xC2B2_AE3D_27D4_EB4F)
        ^ u64::from(sample).wrapping_mul(0x1656_67B1_9E37_79F9);
    x ^= x >> 30;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    (x % n as u64) as usize
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// One MSA output unit (a region or a whole reference contig).
pub(crate) struct MsaOutput {
    pub(crate) contig: Vec<u8>,
    pub(crate) start1: u64,
    pub(crate) end1: u64,
    /// The file name carries the coordinates only for explicit regions.
    pub(crate) explicit: bool,
}

/// Reads a FASTA into per-contig sequences (the MSA reference).
///
/// Streams the file line by line: `fs::read` held the whole reference (a
/// human genome is ~3 GB) and then copied every contig into an `Arc<[u8]>`,
/// so the peak was twice the file plus all the sequence data (doc_v23
/// V23-11). Only the current contig's sequence is buffered twice (the
/// `Vec` and the `Arc` it is moved into), and the `Arc` is the only copy
/// that outlives a record.
pub(super) fn read_reference(path: &PathBuf) -> Result<Arc<ReferenceMap>, crate::Error> {
    use std::io::BufRead as _;
    let file = std::fs::File::open(path).map_err(|e| {
        crate::Error::usage(format!("cannot read --reference {}: {e}", path.display()))
    })?;
    let reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut map: HashMap<Vec<u8>, Arc<[u8]>> = HashMap::new();
    let mut name: Option<Vec<u8>> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.split(b'\n') {
        let line = line.map_err(|e| {
            crate::Error::usage(format!("cannot read --reference {}: {e}", path.display()))
        })?;
        let line = line.strip_suffix(b"\r").unwrap_or(&line);
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(n) = name.take() {
                let duplicate = map
                    .insert(n.clone(), Arc::from(std::mem::take(&mut seq)))
                    .is_some();
                if duplicate {
                    return Err(crate::Error::usage(format!(
                        "--reference {}: duplicate FASTA record name {}",
                        path.display(),
                        String::from_utf8_lossy(&n)
                    )));
                }
            }
            let token = header
                .split(|&b| b == b' ' || b == b'\t')
                .next()
                .unwrap_or(header);
            if token.is_empty() {
                return Err(crate::Error::usage("--reference: empty FASTA record name"));
            }
            name = Some(token.to_vec());
        } else if !line.is_empty() {
            if name.is_none() {
                return Err(crate::Error::usage(
                    "--reference: FASTA sequence before the first `>` header",
                ));
            }
            seq.extend_from_slice(line);
        }
    }
    if let Some(n) = name.take() {
        let duplicate = map.insert(n.clone(), Arc::from(seq)).is_some();
        if duplicate {
            return Err(crate::Error::usage(format!(
                "--reference {}: duplicate FASTA record name {}",
                path.display(),
                String::from_utf8_lossy(&n)
            )));
        }
    }
    if map.is_empty() {
        return Err(crate::Error::usage("--reference: no FASTA records"));
    }
    Ok(Arc::new(map))
}

/// Builds the MSA plan and output list from the reference and the request
/// regions (whole contigs when none were given).
///
/// The spellings are parsed against the reference FASTA's own contig names
/// (`Region::parse_in`, doc_v30 R04): a contig literally named `chr1:100`
/// used to be split into contig `chr1` and position 100 by the header-less
/// heuristic.
pub(super) fn build_msa(
    reference: &Arc<ReferenceMap>,
    region_specs: &[String],
    flank: u64,
) -> Result<(MsaPlan, Vec<MsaOutput>), crate::Error> {
    let mut predicates: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    let mut outputs: Vec<MsaOutput> = Vec::new();
    if region_specs.is_empty() {
        // Deterministic file order: the reference sorted by name.
        let mut names: Vec<&Vec<u8>> = reference.keys().collect();
        names.sort();
        for name in names {
            let len = reference[name].len() as u64;
            predicates.insert(name.clone(), vec![(0, len)]);
            outputs.push(MsaOutput {
                contig: name.clone(),
                start1: 1,
                end1: len,
                explicit: false,
            });
        }
    } else {
        let contigs: Vec<String> = reference
            .keys()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect();
        for spec in region_specs {
            let region = crate::spec::Region::parse_in(spec, &contigs)
                .map_err(|e| crate::Error::usage(format!("invalid --region {spec:?}: {e}")))?;
            if region.start == Some(0) || region.end == Some(0) {
                return Err(crate::Error::usage(format!(
                    "--region {spec:?} uses a 0 bound; positions are 1-based"
                )));
            }
            if let (Some(start), Some(end)) = (region.start, region.end)
                && start > end
            {
                return Err(crate::Error::usage(format!(
                    "--region {spec:?} has a start above its end"
                )));
            }
            let contig = region.contig.as_bytes().to_vec();
            let Some(sequence) = reference.get(&contig) else {
                return Err(crate::Error::usage(format!(
                    "--region names contig {} which the reference FASTA does not contain",
                    region.contig
                )));
            };
            let len = sequence.len() as u64;
            let start1 = region.start.unwrap_or(1).max(1);
            let end1 = region.end.unwrap_or(len).min(len);
            if start1 > len || end1 < start1 {
                return Err(crate::Error::usage(format!(
                    "--region {}:{start1} exceeds the reference ({len} bp)",
                    region.contig
                )));
            }
            let start0 = (start1 - 1).saturating_sub(flank);
            let end0 = (end1 + flank).min(len);
            predicates
                .entry(contig.clone())
                .or_default()
                .push((start0, end0));
            outputs.push(MsaOutput {
                contig,
                start1: start0 + 1,
                end1: end0,
                explicit: true,
            });
        }
    }
    Ok((
        MsaPlan {
            regions: Arc::new(predicates),
            reference: Arc::clone(reference),
        },
        outputs,
    ))
}

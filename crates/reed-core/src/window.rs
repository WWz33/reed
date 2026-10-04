//! The window buffer: L2 materialisation of decoded genotypes.
//!
//! Rules encoded here (see `reed/DESIGN.md`):
//!
//! - `contig_len` is a first-class requirement; zero is a construction error.
//! - v0.1 is diploid and biallelic: any allele index above 1 is rejected, never
//!   silently dropped; any genotype that is not diploid is rejected by the reader.
//! - Missingness is stored both in the genotype sentinel (`-1`) and a parallel
//!   bit vector; the constructor enforces they agree.
//! - Phasing is per genotype (a bit), matching the observed reality that GT and
//!   physical phasing (PGT/PID) are separate.

use std::fmt;

use crate::bits::Bits;
use crate::manifest::SiteCompleteness;

/// Index of a contig within one stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ContigId(pub u32);

/// Violations of the window-buffer invariants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WindowError {
    ZeroContigLength,
    EmptyWindow { start: u64, end: u64 },
    PositionCount { expected: usize, got: usize },
    GenotypeCount { expected: usize, got: usize },
    MissingCount { expected: usize, got: usize },
    PhasedCount { expected: usize, got: usize },
    PositionOutOfWindow { index: usize, pos: u32, start: u64, end: u64 },
    MissingBitMismatch { index: usize },
    /// v0.1 rejects allele indices above 1 (biallelic) and below -1.
    AlleleOutOfRange { index: usize, allele: i8 },
    /// Window span extends past the declared contig length.
    EndBeyondContigLength { end: u64, contig_len: u64 },
    /// Positions must be strictly ascending.
    PositionNotAscending { index: usize, prev: u32, pos: u32 },
}

impl fmt::Display for WindowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WindowError::ZeroContigLength => {
                write!(f, "contig length is required (zero not allowed)")
            }
            WindowError::EmptyWindow { start, end } => {
                write!(f, "empty window [{start}, {end})")
            }
            WindowError::PositionCount { expected, got } => {
                write!(f, "expected {expected} positions, got {got}")
            }
            WindowError::GenotypeCount { expected, got } => {
                write!(f, "expected {expected} genotype alleles, got {got}")
            }
            WindowError::MissingCount { expected, got } => {
                write!(f, "expected {expected} missing bits, got {got}")
            }
            WindowError::PhasedCount { expected, got } => {
                write!(f, "expected {expected} phased bits, got {got}")
            }
            WindowError::PositionOutOfWindow { index, pos, start, end } => write!(
                f,
                "position {pos} (site {index}) outside window [{start}, {end})"
            ),
            WindowError::MissingBitMismatch { index } => write!(
                f,
                "missing bit for genotype {index} disagrees with the -1 allele sentinel"
            ),
            WindowError::AlleleOutOfRange { index, allele } => write!(
                f,
                "allele index {allele} at genotype {index} is out of range (v0.1: biallelic, -1 = missing)"
            ),
            WindowError::EndBeyondContigLength { end, contig_len } => {
                write!(f, "window end {end} exceeds the contig length {contig_len}")
            }
            WindowError::PositionNotAscending { index, prev, pos } => write!(
                f,
                "positions must be strictly ascending: {prev} followed by {pos} (site {index})"
            ),
        }
    }
}

impl std::error::Error for WindowError {}

/// Allele tallies at one site over a sample subset, per-allele deletion.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AlleleCounts {
    pub ref_count: u64,
    pub alt_count: u64,
    /// Number of missing alleles among the requested samples.
    pub missing: u64,
}

/// Per-site diploid genotype category counts over a sample subset.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GenotypeCounts {
    pub hom_ref: u64,
    pub het: u64,
    pub hom_alt: u64,
    /// Individuals with both alleles called.
    pub n_called: u64,
}

/// Decoded genotypes for one window of one contig.
#[derive(Debug, Clone)]
pub struct WindowBuffer {
    pub manifest: SiteCompleteness,
    pub contig: ContigId,
    pub contig_len: u32,
    /// 0-based, half-open `[start, end)`; `end` never exceeds `contig_len`.
    pub start: u64,
    pub end: u64,
    pub n_samples: u32,
    /// 0-based positions, ascending.
    pub pos: Vec<u32>,
    /// Whether the site is segregating within the sample set (any allele > 0
    /// among non-missing calls). ALT present in the record with all-0/0 calls is
    /// *not* variant under this definition.
    pub is_variant: Bits,
    /// `n_sites * n_samples * 2` allele indices in `{-1, 0, 1}` (-1 = missing).
    pub genotypes: Vec<i8>,
    /// `n_sites * n_samples` missing bits (true iff either allele is `-1`).
    pub missing: Bits,
    /// `n_sites * n_samples` phasing bits (per genotype).
    pub phased: Bits,
}

impl WindowBuffer {
    /// Builds a buffer, validating every invariant.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        manifest: SiteCompleteness,
        contig: ContigId,
        contig_len: u32,
        start: u64,
        end: u64,
        n_samples: u32,
        pos: Vec<u32>,
        is_variant: Bits,
        genotypes: Vec<i8>,
        missing: Bits,
        phased: Bits,
    ) -> Result<Self, WindowError> {
        if contig_len == 0 {
            return Err(WindowError::ZeroContigLength);
        }
        if end <= start {
            return Err(WindowError::EmptyWindow { start, end });
        }
        if end > contig_len as u64 {
            return Err(WindowError::EndBeyondContigLength {
                end,
                contig_len: contig_len as u64,
            });
        }
        let n_sites = pos.len();
        if is_variant.len() != n_sites {
            return Err(WindowError::PositionCount {
                expected: n_sites,
                got: is_variant.len(),
            });
        }
        let n_gt = n_sites * n_samples as usize * 2;
        if genotypes.len() != n_gt {
            return Err(WindowError::GenotypeCount {
                expected: n_gt,
                got: genotypes.len(),
            });
        }
        let n_calls = n_sites * n_samples as usize;
        if missing.len() != n_calls {
            return Err(WindowError::MissingCount {
                expected: n_calls,
                got: missing.len(),
            });
        }
        if phased.len() != n_calls {
            return Err(WindowError::PhasedCount {
                expected: n_calls,
                got: phased.len(),
            });
        }
        for (index, &p) in pos.iter().enumerate() {
            if (p as u64) < start || (p as u64) >= end {
                return Err(WindowError::PositionOutOfWindow {
                    index,
                    pos: p,
                    start,
                    end,
                });
            }
            if index > 0 && p <= pos[index - 1] {
                return Err(WindowError::PositionNotAscending {
                    index,
                    prev: pos[index - 1],
                    pos: p,
                });
            }
        }
        for index in 0..n_calls {
            let a = genotypes[2 * index];
            let b = genotypes[2 * index + 1];
            if !(-1..=1).contains(&a) {
                return Err(WindowError::AlleleOutOfRange { index, allele: a });
            }
            if !(-1..=1).contains(&b) {
                return Err(WindowError::AlleleOutOfRange { index, allele: b });
            }
            let is_missing = a == -1 || b == -1;
            if missing.get(index) != is_missing {
                return Err(WindowError::MissingBitMismatch { index });
            }
        }
        Ok(Self {
            manifest,
            contig,
            contig_len,
            start,
            end,
            n_samples,
            pos,
            is_variant,
            genotypes,
            missing,
            phased,
        })
    }

    /// Number of sites in the window.
    pub fn n_sites(&self) -> usize {
        self.pos.len()
    }

    /// Window span in base pairs.
    pub fn window_bp(&self) -> u64 {
        self.end - self.start
    }

    /// Allele indices of one sample's genotype.
    pub fn genotype(&self, site: usize, sample: usize) -> (i8, i8) {
        let i = site * self.n_samples as usize + sample;
        (self.genotypes[2 * i], self.genotypes[2 * i + 1])
    }

    /// Allele tallies at a site over the given samples, per-allele deletion.
    pub fn allele_counts(&self, site: usize, samples: &[u32]) -> AlleleCounts {
        let n_samples = self.n_samples as usize;
        let mut counts = AlleleCounts::default();
        for &s in samples {
            let i = site * n_samples + s as usize;
            for k in 0..2 {
                match self.genotypes[2 * i + k] {
                    -1 => counts.missing += 1,
                    0 => counts.ref_count += 1,
                    _ => counts.alt_count += 1,
                }
            }
        }
        counts
    }

    /// Genotype category counts at a site over the given samples; individuals
    /// with any missing allele are excluded.
    pub fn genotype_counts(&self, site: usize, samples: &[u32]) -> GenotypeCounts {
        let mut counts = GenotypeCounts::default();
        for &s in samples {
            match self.genotype(site, s as usize) {
                (0, 0) => counts.hom_ref += 1,
                (1, 1) => counts.hom_alt += 1,
                (0, 1) | (1, 0) => counts.het += 1,
                _ => continue,
            }
            counts.n_called += 1;
        }
        counts
    }
}

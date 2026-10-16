//! `reed convert` (docs/rounds/v13/02-architecture.md): the multi-target
//! conversion driver — phylogenetic matrices (vcf2phylip semantics), the
//! genotype-count matrix (`--012`, vcftools), TASSEL hapmap, reference-driven
//! MSA (vcf2msa semantics) and the projected-SFS export (mirroring
//! `reed stats --sfs`).
//!
//! Every target folds over the same streamed chunk pipeline (§7): one line
//! parse per record, per-target column/row folding, order-preserving additive
//! partials, and a drain-time render. Output bytes are independent of the
//! thread count (each partial merges in chunk order).

mod hapmap;
mod job;
pub(crate) mod writers;

pub use job::{ConvertCounters, ConvertReport, run_convert, run_convert_opened};

use std::path::PathBuf;

use crate::stats::request::{Populations, SfsRequest};

/// One conversion target (v13 §2.1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConvertTarget {
    /// The projected SFS + its dadi/fastsimcoal2 export; semantics identical
    /// to `reed stats --sfs` (v13 §4).
    Sfs {
        /// The spectra request (`--sfs [folded] --proj …`).
        request: SfsRequest,
        /// `--sfs-int`: round the spectra to integers.
        int_dtype: bool,
    },
    /// The nucleotide alignment formats sharing one column model (v13 §2.3).
    Phylo {
        /// Which of PHYLIP/FASTA/NEXUS to write (at least one).
        formats: PhyloFormats,
        /// The alignment-level options.
        options: PhyloOptions,
    },
    /// The binary NEXUS stream (SNAPP; v13 §3.4).
    NexusBin,
    /// The 0/1/2 genotype-count matrix (vcftools `--012`, v13 §3.8).
    CountMatrix,
    /// The TASSEL hapmap `.hmp.txt` (v13 §3.9).
    Hapmap,
    /// The reference-driven whole-site MSA (vcf2msa semantics, v13 §3.10).
    Msa {
        /// The reference genome FASTA (`--reference`).
        reference: PathBuf,
        /// The output regions (empty = every reference contig).
        regions: Vec<MsaRegion>,
        /// Bases added to both sides of each region (`--flank`).
        flank: u64,
        /// Collapse sample names at the first occurrence of this string
        /// (`--collapse-suffix`); `None` keeps names verbatim.
        collapse_suffix: Option<String>,
    },
    /// The hapmap→VCF reverse conversion (v13 §3.12). Exclusive with every
    /// VCF-input target: the input must be a hapmap (`.hmp.txt`).
    Vcf {
        /// `-Oz` (BGZF) instead of `-Ov` (plain text).
        compressed: bool,
        /// The deflate level (`-Oz<level>`, 1–9).
        level: u8,
        /// The output file (`-o`); `None` writes the VCF to stdout.
        output: Option<PathBuf>,
        /// The `##reed_command` value (captured argv); empty = omitted.
        command_line: String,
        /// `--no-version`: drop the `##reed_*` provenance lines.
        no_version: bool,
    },
}

/// The nucleotide alignment formats sharing one column model (§2.3).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PhyloFormats {
    /// `--phylip`.
    pub phylip: bool,
    /// `--fasta`.
    pub fasta: bool,
    /// `--nexus`.
    pub nexus: bool,
}

impl PhyloFormats {
    /// At least one nucleotide format selected.
    #[must_use]
    pub const fn any(self) -> bool {
        self.phylip || self.fasta || self.nexus
    }
}

/// The alignment-level options (`Phylo` and `NexusBin` share them, v13 §3.5).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PhyloOptions {
    /// `--outgroup NAME`: this sample is written first in every alignment.
    pub outgroup: Option<String>,
    /// `--resolve-iupac`: resolve ambiguous genotypes to one allele
    /// (deterministic, v13 §3.3).
    pub resolve_iupac: bool,
    /// `--write-used-sites`: write `PREFIX.used_sites.tsv`.
    pub write_used_sites: bool,
}

/// One `--msa` output region (1-based inclusive bounds; `None` = open end).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MsaRegion {
    /// Contig name (must be in the reference FASTA).
    pub contig: String,
    /// Inclusive 1-based lower bound (`None` = contig start).
    pub start1: Option<u64>,
    /// Inclusive 1-based upper bound (`None` = contig end).
    pub end1: Option<u64>,
}

/// One conversion run: the targets, the sample subset and the output prefix.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConvertRequest {
    /// The targets to produce (validated non-empty by the driver).
    pub targets: Vec<ConvertTarget>,
    /// Output prefix: `PREFIX.phy`, `PREFIX_sfs/…`, … (the CLI resolves the
    /// default from the input file name; v13 §3.6).
    pub prefix: PathBuf,
    /// The population assignment for the SFS target (`--pop`; `All` = one
    /// group, v13 §1.2).
    pub populations: Populations,
    /// `--force-samples`: group-file rows naming absent samples only warn.
    pub force_samples: bool,
    /// Static sample subset (`-s/-S/--keep/--remove`), header order.
    pub keep_samples: Option<Vec<bool>>,
    /// The explicit `--min-samples-locus N` (the phylo default is 4 when
    /// unset; `--012`/`--hapmap` apply the gate only when it is given,
    /// v13 §3.1/§3.8).
    pub min_samples_locus: Option<u32>,
    /// The `--resolve-iupac` draw seed (default 0; fixed per the determinism
    /// contract, v13 §3.3).
    pub seed: u64,
}

impl Default for ConvertRequest {
    fn default() -> Self {
        Self {
            targets: Vec::new(),
            prefix: PathBuf::new(),
            populations: Populations::All,
            force_samples: false,
            keep_samples: None,
            min_samples_locus: None,
            seed: 0,
        }
    }
}

impl ConvertRequest {
    /// Renders the request back to CLI spellings (`to_args` is the text
    /// contract the CLI contract tests round-trip, v11 N3).
    ///
    /// # Errors
    ///
    /// A request that cannot be spelled (a target missing its required
    /// companion, e.g. `Msa` without a reference path).
    pub fn to_args(&self) -> Result<Vec<String>, crate::Error> {
        let mut args: Vec<String> = Vec::new();
        for target in &self.targets {
            match target {
                ConvertTarget::Sfs { request, int_dtype } => {
                    args.push("--sfs".into());
                    args.push(if request.folded { "folded" } else { "unfolded" }.into());
                    if !request.projections.is_empty() {
                        args.push("--proj".into());
                        args.push(
                            request
                                .projections
                                .iter()
                                .map(u32::to_string)
                                .collect::<Vec<_>>()
                                .join(","),
                        );
                    }
                    if *int_dtype {
                        args.push("--sfs-int".into());
                    }
                }
                ConvertTarget::Phylo { formats, options } => {
                    if formats.phylip {
                        args.push("--phylip".into());
                    }
                    if formats.fasta {
                        args.push("--fasta".into());
                    }
                    if formats.nexus {
                        args.push("--nexus".into());
                    }
                    if formats.any() {
                        if let Some(outgroup) = &options.outgroup {
                            args.push("--outgroup".into());
                            args.push(outgroup.clone());
                        }
                        if options.resolve_iupac {
                            args.push("--resolve-iupac".into());
                        }
                        if options.write_used_sites {
                            args.push("--write-used-sites".into());
                        }
                    }
                }
                ConvertTarget::Vcf {
                    compressed,
                    level,
                    output,
                    no_version,
                    ..
                } => {
                    args.push("--vcf".into());
                    args.push("-O".into());
                    args.push(if *compressed {
                        format!("z{level}")
                    } else {
                        "v".into()
                    });
                    if let Some(path) = output {
                        let text = path
                            .to_str()
                            .ok_or_else(|| crate::Error::usage("output path is not UTF-8"))?;
                        args.push("-o".into());
                        args.push(text.into());
                    }
                    if *no_version {
                        args.push("--no-version".into());
                    }
                }
                ConvertTarget::NexusBin => args.push("--nexus-binary".into()),
                ConvertTarget::CountMatrix => args.push("--012".into()),
                ConvertTarget::Hapmap => args.push("--hapmap".into()),
                ConvertTarget::Msa {
                    reference,
                    regions,
                    flank,
                    collapse_suffix,
                } => {
                    args.push("--msa".into());
                    args.push("--reference".into());
                    let path = reference
                        .to_str()
                        .ok_or_else(|| crate::Error::usage("reference path is not UTF-8"))?;
                    args.push(path.into());
                    for region in regions {
                        use std::fmt::Write as _;
                        let mut spelling = region.contig.clone();
                        match (region.start1, region.end1) {
                            (Some(s), Some(e)) => {
                                let _ = write!(spelling, ":{s}-{e}");
                            }
                            (Some(s), None) => {
                                let _ = write!(spelling, ":{s}-");
                            }
                            (None, Some(e)) => {
                                let _ = write!(spelling, ":-{e}");
                            }
                            (None, None) => {}
                        }
                        args.push("--region".into());
                        args.push(spelling);
                    }
                    if *flank > 0 {
                        args.push("--flank".into());
                        args.push(flank.to_string());
                    }
                    if let Some(suffix) = collapse_suffix {
                        args.push("--collapse-suffix".into());
                        args.push(suffix.clone());
                    }
                }
            }
        }
        if let Some(gate) = self.min_samples_locus {
            args.push("--min-samples-locus".into());
            args.push(gate.to_string());
        }
        Ok(args)
    }
}

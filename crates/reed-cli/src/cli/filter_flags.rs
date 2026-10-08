//! Shared filtering flag groups (04-cleanup §2.1): `reed filter` and
//! `reed stats` accept the same region and sample selection spellings, so
//! the definitions live here once.

use clap::Args;

/// Region selection (`doc_v2/02-cli.md` §2): indexed and streamed regions.
#[derive(Args, Clone, Debug, Default)]
pub struct RegionArgs {
    /// Regions to keep, comma-separated (`Chr01`, `Chr01:100-200`);
    /// needs an index next to the input (bcftools `-r`)
    #[arg(short = 'r', long = "regions", value_name = "LIST")]
    pub regions: Option<String>,
    /// Regions to keep, one per line; needs an index (bcftools `-R`)
    #[arg(short = 'R', long = "regions-file", value_name = "FILE")]
    pub regions_file: Option<String>,
    /// Regions to keep, streamed without an index (bcftools `-t`)
    #[arg(short = 't', long = "targets", value_name = "LIST")]
    pub targets: Option<String>,
    /// Regions to keep, one per line (bcftools `-T`); `^` prefix excludes
    #[arg(short = 'T', long = "targets-file", value_name = "FILE")]
    pub targets_file: Option<String>,
}

/// Sample selection (`doc_v2/02-cli.md` §3): keep and remove lists.
#[derive(Args, Clone, Debug, Default)]
pub struct SampleArgs {
    /// Samples to keep: comma-separated names; `^` prefix excludes
    /// (`^A,B` = keep all but A and B; bcftools `-s`)
    #[arg(short = 's', long = "samples", value_name = "[^]LIST")]
    pub samples: Vec<String>,
    /// Samples to keep, one name per line (bcftools `-S`)
    #[arg(short = 'S', long = "samples-file", value_name = "FILE")]
    pub samples_file: Option<String>,
    /// Samples to keep, one name per line (plink2/vcftools spelling)
    #[arg(long, value_name = "FILE")]
    pub keep: Option<String>,
    /// Samples to drop, one name per line
    #[arg(long, value_name = "FILE")]
    pub remove: Option<String>,
}

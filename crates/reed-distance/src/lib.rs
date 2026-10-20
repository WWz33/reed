//! Individual-by-individual matrices for VCF data (doc_v17).
//!
//! Every estimator here answers "how similar are individuals i and j", but
//! they answer it differently, and the differences are what make cross-tool
//! comparisons fail. Each function therefore documents which reference
//! implementation it reproduces and cites the source line it was derived from
//! (see `docs/rounds/v17/01-research.md`).
//!
//! | matrix | reference | output |
//! |---|---|---|
//! | [`allele_ct_distance`] | plink `--distance` | `.dist` |
//! | [`ibs`] / [`one_minus_ibs`] | plink `--distance ibs` / `1-ibs` | `.mibs` / `.mdist` |
//! | [`grm`] / [`covariance`] | plink2 `--make-rel [cov]`, GCTA | `.rel` |
//! | [`king`] | plink2 `--make-king` | `.king` |
//! | [`pi_hat`] | plink `--genome` | `.genome` |
//! | [`phi`] | vcftools `--relatedness2` | `.relatedness2` |
//! | [`ajk`] | vcftools `--relatedness` | `.relatedness` |
//! | [`discordance`] | bcftools `gtcheck` | `.gtcheck` |

#![allow(clippy::needless_range_loop)]
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::suboptimal_flops)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::needless_borrow)]
#![allow(clippy::if_same_then_else)]

pub mod grm;
pub mod matrices;
pub mod output;
pub mod standard;

pub use matrices::{
    ajk, allele_ct_distance, discordance, grm, ibs, king, one_minus_ibs, phi, pi_hat,
    AjkOptions, DiscordanceMode, DiscordanceOptions, KingCounts, MatrixKind, MatrixResult,
    PiHatOptions, PiHatResult,
};
pub use output::{write_matrix, MatrixFormat};
pub use grm::{Grm, MissingPolicy};
pub use standard::{standardize, standardize_binom2, Stand, Standardized};

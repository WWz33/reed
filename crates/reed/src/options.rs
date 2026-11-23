//! The typed option vocabularies shared by the facade's requests (doc_v30
//! P1-6): a finite set of CLI spellings is a type here, not a `String` that
//! is only checked at run time.
//!
//! Each enum implements [`std::str::FromStr`] with `Err = String` (the error
//! text the CLI prints) and [`std::fmt::Display`] (the canonical spelling the
//! CLI and the provenance header use). The lower crates keep their own
//! enums; the facade converts at the boundary so those types do not become
//! part of the public surface.

use std::fmt;
use std::str::FromStr;

use reed_distance::matrices::MatrixKind as LowerMatrixKind;
use reed_distance::output::MatrixFormat as LowerMatrixFormat;

/// Genotype standardization: `binom2` divides by `sqrt(2p(1-p))`, `binom`
/// by `sqrt(p(1-p))`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Stand {
    /// `binom2` (GCTA/plink2/flashpca).
    #[default]
    Binom2,
    /// `binom` (older EIGENSOFT).
    Binom,
}

impl Stand {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "binom2" => Ok(Self::Binom2),
            "binom" => Ok(Self::Binom),
            other => Err(format!("unknown standardization {other:?} (binom2|binom)")),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Binom2 => "binom2",
            Self::Binom => "binom",
        }
    }

    /// The lower crate's enum for the solve.
    pub(crate) const fn lower(self) -> reed_reduce::pca::Stand {
        match self {
            Self::Binom2 => reed_reduce::pca::Stand::Binom2,
            Self::Binom => reed_reduce::pca::Stand::Binom,
        }
    }

    /// The distance crate's enum for the estimator.
    pub(crate) const fn distance(self) -> reed_distance::standard::Stand {
        match self {
            Self::Binom2 => reed_distance::standard::Stand::Binom2,
            Self::Binom => reed_distance::standard::Stand::Binom,
        }
    }
}

impl fmt::Display for Stand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for Stand {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

/// Missing-genotype policy for the matrix estimators.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MissingPolicy {
    /// Each pair uses its own denominator (plink2/GCTA default).
    #[default]
    Pairwise,
    /// Missing genotypes are mean-imputed (flashpca/plink2 meanimpute).
    MeanImpute,
}

impl MissingPolicy {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "pairwise" => Ok(Self::Pairwise),
            "meanimpute" | "mean-impute" => Ok(Self::MeanImpute),
            other => Err(format!(
                "unknown missing-genotype policy {other:?} (pairwise|meanimpute)"
            )),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Pairwise => "pairwise",
            Self::MeanImpute => "meanimpute",
        }
    }

    pub(crate) const fn distance(self) -> reed_distance::matrices::MissingPolicy {
        match self {
            Self::Pairwise => reed_distance::matrices::MissingPolicy::Pairwise,
            Self::MeanImpute => reed_distance::matrices::MissingPolicy::MeanImpute,
        }
    }
}

impl fmt::Display for MissingPolicy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for MissingPolicy {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

/// The PCA solver request.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PcaSolver {
    /// Pick dense when the dense path's peak fits the memory budget.
    #[default]
    Auto,
    /// Materialize the `N x N` matrix and extract the full spectrum.
    Dense,
    /// Matrix-free block Lanczos.
    Stream,
}

impl PcaSolver {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "auto" => Ok(Self::Auto),
            "dense" => Ok(Self::Dense),
            "stream" => Ok(Self::Stream),
            other => Err(format!("unknown PCA solver {other:?} (auto|dense|stream)")),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Dense => "dense",
            Self::Stream => "stream",
        }
    }
}

impl fmt::Display for PcaSolver {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for PcaSolver {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

/// How `--xval` summarizes per-repetition accuracy.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum XvalResultMode {
    /// Mean over every held-out individual.
    #[default]
    Overall,
    /// Mean of the per-group means.
    GroupMean,
}

impl XvalResultMode {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "overall" => Ok(Self::Overall),
            "group-mean" => Ok(Self::GroupMean),
            other => Err(format!(
                "unknown --xval-result {other:?} (overall|group-mean)"
            )),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Overall => "overall",
            Self::GroupMean => "group-mean",
        }
    }

    pub(crate) const fn lower(self) -> reed_reduce::xval::XvalResult {
        match self {
            Self::Overall => reed_reduce::xval::XvalResult::Overall,
            Self::GroupMean => reed_reduce::xval::XvalResult::GroupMean,
        }
    }
}

impl fmt::Display for XvalResultMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for XvalResultMode {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

/// Which matrix or per-pair table `reed distance` computes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DistanceKind {
    /// plink `--distance` (allele-count units).
    AlleleCountDistance,
    /// plink `--distance ibs` (identity by state).
    Ibs,
    /// plink `--distance 1-ibs` (genomic distance).
    OneMinusIbs,
    /// plink2 `--make-rel` / GCTA (variance-standardized relationship).
    #[default]
    Grm,
    /// plink2 `--make-rel cov`.
    Covariance,
    /// plink2 `--make-king`.
    King,
    /// plink `--genome` (PI_HAT).
    PiHat,
    /// vcftools `--relatedness2` (Manichaikul phi).
    Phi,
    /// vcftools `--relatedness` (Yang Ajk).
    Ajk,
    /// bcftools `gtcheck` (discordance).
    Discordance,
}

impl DistanceKind {
    /// Parses the CLI spelling (including the reference-tool aliases).
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "allele-ct" | "dist" => Ok(Self::AlleleCountDistance),
            "ibs" => Ok(Self::Ibs),
            "1-ibs" => Ok(Self::OneMinusIbs),
            "grm" | "rel" => Ok(Self::Grm),
            "cov" => Ok(Self::Covariance),
            "king" => Ok(Self::King),
            "pihat" | "genome" => Ok(Self::PiHat),
            "phi" | "relatedness2" => Ok(Self::Phi),
            "ajk" | "relatedness" => Ok(Self::Ajk),
            "discordance" | "gtcheck" => Ok(Self::Discordance),
            other => Err(format!(
                "unknown --type {other:?} \
                 (allele-ct|ibs|1-ibs|grm|cov|king|pihat|phi|ajk|discordance)"
            )),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::AlleleCountDistance => "allele-ct",
            Self::Ibs => "ibs",
            Self::OneMinusIbs => "1-ibs",
            Self::Grm => "grm",
            Self::Covariance => "cov",
            Self::King => "king",
            Self::PiHat => "pihat",
            Self::Phi => "phi",
            Self::Ajk => "ajk",
            Self::Discordance => "discordance",
        }
    }

    /// Whether this kind is written as a matrix (plus `.id`) rather than a
    /// per-pair table.
    #[must_use]
    pub const fn is_matrix(self) -> bool {
        matches!(
            self,
            Self::AlleleCountDistance
                | Self::Ibs
                | Self::OneMinusIbs
                | Self::Grm
                | Self::Covariance
                | Self::King
        )
    }

    /// The output suffix the reference tool uses.
    #[must_use]
    pub fn suffix(self) -> &'static str {
        self.lower().suffix()
    }

    pub(crate) const fn lower(self) -> LowerMatrixKind {
        match self {
            Self::AlleleCountDistance => LowerMatrixKind::AlleleCountDistance,
            Self::Ibs => LowerMatrixKind::Ibs,
            Self::OneMinusIbs => LowerMatrixKind::OneMinusIbs,
            Self::Grm => LowerMatrixKind::Grm,
            Self::Covariance => LowerMatrixKind::Covariance,
            Self::King => LowerMatrixKind::King,
            Self::PiHat => LowerMatrixKind::PiHat,
            Self::Phi => LowerMatrixKind::Phi,
            Self::Ajk => LowerMatrixKind::Ajk,
            Self::Discordance => LowerMatrixKind::Discordance,
        }
    }
}

impl fmt::Display for DistanceKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for DistanceKind {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

/// Output layout for the matrix kinds.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DistanceFormat {
    /// Lower triangle, row-major (plink's default).
    #[default]
    Triangle,
    /// Full symmetric square matrix.
    Square,
    /// Full square with the upper-right triangle zeroed (plink `square0`).
    Square0,
    /// Binary double-precision square matrix (R `readBin`).
    Bin,
}

impl DistanceFormat {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "triangle" => Ok(Self::Triangle),
            "square" => Ok(Self::Square),
            "square0" => Ok(Self::Square0),
            "bin" => Ok(Self::Bin),
            other => Err(format!(
                "unknown matrix format {other:?} (triangle|square|square0|bin)"
            )),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Triangle => "triangle",
            Self::Square => "square",
            Self::Square0 => "square0",
            Self::Bin => "bin",
        }
    }

    pub(crate) const fn lower(self) -> LowerMatrixFormat {
        match self {
            Self::Triangle => LowerMatrixFormat::Triangle,
            Self::Square => LowerMatrixFormat::Square,
            Self::Square0 => LowerMatrixFormat::Square0,
            Self::Bin => LowerMatrixFormat::Bin,
        }
    }
}

impl fmt::Display for DistanceFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for DistanceFormat {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

/// The tree construction `reed distance --tree` uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TreeMethod {
    /// Average linkage (rooted; support by clade).
    #[default]
    Upgma,
    /// Neighbor-joining (unrooted; support by bipartition).
    Nj,
}

impl TreeMethod {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "upgma" => Ok(Self::Upgma),
            "nj" => Ok(Self::Nj),
            other => Err(format!("unknown --tree {other:?} (upgma|nj)")),
        }
    }

    /// The canonical CLI spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Upgma => "upgma",
            Self::Nj => "nj",
        }
    }

    pub(crate) const fn lower(self) -> reed_graph::TreeMethod {
        match self {
            Self::Upgma => reed_graph::TreeMethod::Upgma,
            Self::Nj => reed_graph::TreeMethod::Nj,
        }
    }
}

impl fmt::Display for TreeMethod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

impl FromStr for TreeMethod {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        Self::parse(text.trim())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_finite_vocabulary_round_trips_through_its_cli_spelling() {
        for text in ["binom2", "binom"] {
            let parsed: Stand = text.parse().expect("invariant: spelling");
            assert_eq!(parsed.to_string(), text);
        }
        for text in ["pairwise", "meanimpute"] {
            let parsed: MissingPolicy = text.parse().expect("invariant: spelling");
            assert_eq!(parsed.to_string(), text);
        }
        for text in ["auto", "dense", "stream"] {
            let parsed: PcaSolver = text.parse().expect("invariant: spelling");
            assert_eq!(parsed.to_string(), text);
        }
        for text in ["overall", "group-mean"] {
            let parsed: XvalResultMode = text.parse().expect("invariant: spelling");
            assert_eq!(parsed.to_string(), text);
        }
        for text in ["allele-ct", "dist"] {
            assert_eq!(
                DistanceKind::parse(text).expect("invariant: spelling"),
                DistanceKind::AlleleCountDistance
            );
        }
        assert!(DistanceKind::parse("angle-ct").is_err());
        for text in ["triangle", "square", "square0", "bin"] {
            let parsed: DistanceFormat = text.parse().expect("invariant: spelling");
            assert_eq!(parsed.to_string(), text);
        }
        for text in ["upgma", "nj"] {
            let parsed: TreeMethod = text.parse().expect("invariant: spelling");
            assert_eq!(parsed.to_string(), text);
        }
    }
}

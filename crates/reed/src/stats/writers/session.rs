//! Streaming window output sessions (doc_v27 V27-18).
//!
//! The CLI writes each completed row to a sibling temporary file immediately;
//! the final paths are published only after the driver's input-integrity
//! checks pass. This is the same per-file atomic contract as V26-08, but
//! avoids retaining every `WindowRow` until the run ends.

use std::fmt::Write as _;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use reed_engine::error::EngineError;
use reed_io::write::atomic::AtomicFile;
use reed_stats::request::StatKind;

use crate::stats::{InvariantMode, WindowLayoutOptions, WindowRow};

/// Output-session metadata that only the CLI supplies. Library requests
/// leave it unset and keep their full `WindowRow` accumulators instead.
#[derive(Clone, Debug)]
pub struct WindowOutputConfig {
    /// `-o PREFIX` for the window-layout files.
    pub prefix: PathBuf,
    /// Command line recorded in the run log.
    pub command: String,
    /// Package version recorded in the run log.
    pub version: String,
    /// Whether an SFS request shares this run.
    pub sfs: bool,
    /// Whether f-statistics share this run.
    pub fstats: bool,
    /// Whether LD shares this run.
    pub ld: bool,
}

struct Pending {
    path: PathBuf,
    writer: std::io::BufWriter<std::fs::File>,
    atomic: Option<AtomicFile>,
}

impl Pending {
    fn create(path: PathBuf, header: &str) -> Result<Self, EngineError> {
        let (atomic, file) = AtomicFile::create(&path).map_err(EngineError::Io)?;
        let mut writer = std::io::BufWriter::with_capacity(1 << 16, file);
        writer
            .write_all(header.as_bytes())
            .map_err(EngineError::Io)?;
        Ok(Self {
            path,
            writer,
            atomic: Some(atomic),
        })
    }

    fn line(&mut self, text: &str) -> Result<(), EngineError> {
        self.writer
            .write_all(text.as_bytes())
            .map_err(EngineError::Io)
    }

    fn finish(mut self) -> Result<PathBuf, EngineError> {
        self.writer.flush().map_err(EngineError::Io)?;
        drop(self.writer);
        self.atomic
            .take()
            .expect("invariant: uncommitted output")
            .commit()
            .map_err(EngineError::Io)?;
        Ok(self.path)
    }
}

/// The set of files requested by the current `kinds`. Files are created at
/// run start, including the empty-result case, and committed at the end.
pub struct WindowOutputSession {
    kinds: Vec<StatKind>,
    group_names: Vec<String>,
    variant_only: bool,
    pi: Option<Pending>,
    dxy: Option<Pending>,
    theta: Option<Pending>,
    tajima: Option<Pending>,
    fst: Option<Pending>,
    fst_global: Option<Pending>,
    pbs: Vec<(String, Pending)>,
    command: String,
    version: String,
    sfs: bool,
    fstats: bool,
    ld: bool,
}

impl WindowOutputSession {
    /// Creates the requested files and writes their headers.
    pub fn new(
        config: &WindowOutputConfig,
        kinds: &[StatKind],
        group_names: &[String],
        layout: Option<&WindowLayoutOptions>,
    ) -> Result<Self, EngineError> {
        let prefix = &config.prefix;
        let has = |predicate: fn(&StatKind) -> bool| kinds.iter().any(predicate);
        let pi = has(|k| matches!(k, StatKind::Pi { .. }))
            .then(|| {
                Pending::create(
                    with_suffix(prefix, "_pi.txt"),
                    "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_pi\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n",
                )
            })
            .transpose()?;
        let dxy = has(|k| matches!(k, StatKind::Dxy { .. }))
            .then(|| {
                Pending::create(
                    with_suffix(prefix, "_dxy.txt"),
                    "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_dxy\tno_sites\tcount_diffs\tcount_comparisons\tcount_missing\n",
                )
            })
            .transpose()?;
        let theta = has(|k| matches!(k, StatKind::ThetaW { .. }))
            .then(|| {
                Pending::create(
                    with_suffix(prefix, "_watterson_theta.txt"),
                    "pop\tchromosome\twindow_pos_1\twindow_pos_2\tavg_watterson_theta\tno_sites\traw_watterson_theta\tno_var_sites\tweighted_no_sites\n",
                )
            })
            .transpose()?;
        let tajima = has(|k| matches!(k, StatKind::TajimaD { .. }))
            .then(|| {
                Pending::create(
                    with_suffix(prefix, "_tajima_d.txt"),
                    "pop\tchromosome\twindow_pos_1\twindow_pos_2\ttajima_d\tno_sites\traw_pi\traw_watterson_theta\ttajima_d_stdev\n",
                )
            })
            .transpose()?;
        let fst_method = kinds.iter().find_map(|kind| match kind {
            StatKind::Fst { kind, .. } => Some(match kind {
                reed_stats::fst::FstKind::Hudson => "hudson",
                reed_stats::fst::FstKind::WeirCockerham84 => "wc",
            }),
            _ => None,
        });
        let fst = fst_method
            .map(|method| {
                Pending::create(
                    with_suffix(prefix, "_fst.txt"),
                    &format!(
                        "pop1\tpop2\tchromosome\twindow_pos_1\twindow_pos_2\tavg_{method}_fst\tno_snps\n"
                    ),
                )
            })
            .transpose()?;
        let fst_global = has(|k| matches!(k, StatKind::FstGlobal))
            .then(|| {
                Pending::create(
                    with_suffix(prefix, "_fst_global.txt"),
                    "pops\tchromosome\twindow_pos_1\twindow_pos_2\tavg_wc_fst\tno_snps\n",
                )
            })
            .transpose()?;
        let mut pbs = Vec::new();
        let mut stems = Vec::new();
        for kind in kinds {
            let StatKind::Pbs { a, b, c, .. } = kind else {
                continue;
            };
            let name = |i: usize| group_names.get(i).map_or("?", String::as_str);
            let (na, nb, nc) = (name(*a), name(*b), name(*c));
            let base = pbs_stem(na, nb, nc);
            let mut stem = base.clone();
            let mut n = 2;
            while stems.contains(&stem) {
                stem = format!("{base}-{n}");
                n += 1;
            }
            stems.push(stem.clone());
            let suffix = if kinds
                .iter()
                .filter(|k| matches!(k, StatKind::Pbs { .. }))
                .count()
                == 1
            {
                "_pbs.txt".to_string()
            } else {
                format!("_pbs.{stem}.txt")
            };
            let header = format!(
                "{na}\t{nb}\t{nc}\tchromosome\twindow_pos_1\twindow_pos_2\tpbs\tno_snps\tfst_12\tfst_13\tfst_23\n"
            );
            pbs.push((
                stem,
                Pending::create(with_suffix(prefix, &suffix), &header)?,
            ));
        }
        Ok(Self {
            kinds: kinds.to_vec(),
            group_names: group_names.to_vec(),
            variant_only: layout
                .is_some_and(|layout| layout.invariant == InvariantMode::VariantOnly),
            pi,
            dxy,
            theta,
            tajima,
            fst,
            fst_global,
            pbs,
            command: config.command.clone(),
            version: config.version.clone(),
            sfs: config.sfs,
            fstats: config.fstats,
            ld: config.ld,
        })
    }

    /// Writes one completed window/contig/genome row.
    pub fn write_row(&mut self, row: &WindowRow) -> Result<(), EngineError> {
        let mut pbs_index = 0;
        for (index, kind) in self.kinds.iter().enumerate() {
            let Some(output) = row.outputs.get(index) else {
                continue;
            };
            let name = |i: usize| self.group_names.get(i).map_or("?", String::as_str);
            match kind {
                StatKind::Pi { group, .. } => {
                    let Some((diffs, comparable, missing)) = output.pair_counts() else {
                        continue;
                    };
                    let n_sites = output.n_present().unwrap_or(0);
                    let value = value_or_na(output.value(), n_sites == 0 && !self.variant_only);
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        name(*group),
                        row.contig,
                        row.start1,
                        row.end,
                        value,
                        n_sites,
                        diffs,
                        comparable,
                        missing
                    );
                    self.pi.as_mut().expect("invariant: pi file").line(&line)?;
                }
                StatKind::Dxy { a, b, .. } => {
                    let Some((diffs, comparable, missing)) = output.pair_counts() else {
                        continue;
                    };
                    let n_sites = output.n_present().unwrap_or(0);
                    let value = value_or_na(output.value(), n_sites == 0 && !self.variant_only);
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        name(*a),
                        name(*b),
                        row.contig,
                        row.start1,
                        row.end,
                        value,
                        n_sites,
                        diffs,
                        comparable,
                        missing
                    );
                    self.dxy
                        .as_mut()
                        .expect("invariant: dxy file")
                        .line(&line)?;
                }
                StatKind::ThetaW { group, .. } => {
                    let Some((raw, no_sites, no_var, weighted)) = output.theta_columns() else {
                        continue;
                    };
                    let avg = value_or_na(output.value(), no_sites == 0 && !self.variant_only);
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        name(*group),
                        row.contig,
                        row.start1,
                        row.end,
                        avg,
                        no_sites,
                        raw,
                        no_var,
                        weighted
                    );
                    self.theta
                        .as_mut()
                        .expect("invariant: theta file")
                        .line(&line)?;
                }
                StatKind::TajimaD { group, .. } => {
                    let Some((raw_pi, raw_theta, stdev, no_sites)) = output.tajima_columns() else {
                        continue;
                    };
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        name(*group),
                        row.contig,
                        row.start1,
                        row.end,
                        cell(output.value()),
                        no_sites,
                        cell(raw_pi),
                        cell(raw_theta),
                        cell(stdev)
                    );
                    self.tajima
                        .as_mut()
                        .expect("invariant: tajima file")
                        .line(&line)?;
                }
                StatKind::Fst { a, b, .. } => {
                    let Some((value, no_snps)) = output.fst_columns() else {
                        continue;
                    };
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        name(*a),
                        name(*b),
                        row.contig,
                        row.start1,
                        row.end,
                        cell(value),
                        no_snps
                    );
                    self.fst
                        .as_mut()
                        .expect("invariant: fst file")
                        .line(&line)?;
                }
                StatKind::FstGlobal => {
                    let Some((value, no_snps)) = output.fst_columns() else {
                        continue;
                    };
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        self.group_names.join(","),
                        row.contig,
                        row.start1,
                        row.end,
                        cell(value),
                        no_snps
                    );
                    self.fst_global
                        .as_mut()
                        .expect("invariant: global fst file")
                        .line(&line)?;
                }
                StatKind::Pbs { a, b, c, .. } => {
                    let Some((value, no_snps, (f_ab, f_ac, f_bc))) = output.pbs_columns() else {
                        continue;
                    };
                    let line = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        name(*a),
                        name(*b),
                        name(*c),
                        row.contig,
                        row.start1,
                        row.end,
                        cell(value),
                        no_snps,
                        cell(f_ab),
                        cell(f_ac),
                        cell(f_bc)
                    );
                    self.pbs[pbs_index].1.line(&line)?;
                    pbs_index += 1;
                }
                StatKind::Sfs { .. } => {}
            }
        }
        Ok(())
    }

    /// Writes the run log and publishes every file. Called only after all
    /// input-integrity checks have passed; dropping the session instead
    /// removes all temporary files.
    #[allow(clippy::too_many_arguments)]
    pub fn finish(
        mut self,
        n_windows: u64,
        n_records: u64,
        n_invariant_sites: u64,
        n_multiallelic_skipped: u64,
        invariant: InvariantMode,
        prefix: &Path,
        sfs_report: Option<&crate::stats::SfsReport>,
    ) -> Result<Vec<PathBuf>, EngineError> {
        let mut files = Vec::new();
        for pending in [self.pi.take(), self.dxy.take(), self.theta.take()]
            .into_iter()
            .flatten()
        {
            files.push(pending.finish()?);
        }
        for (_, pending) in self.pbs.drain(..) {
            files.push(pending.finish()?);
        }
        if let Some(pending) = self.fst.take() {
            files.push(pending.finish()?);
        }
        if let Some(pending) = self.tajima.take() {
            files.push(pending.finish()?);
        }
        if let Some(pending) = self.fst_global.take() {
            files.push(pending.finish()?);
        }
        let methods = crate::stats::writers::windows::methods_text(
            &self.kinds,
            self.sfs,
            self.fstats,
            self.ld,
        );
        let mut log = format!(
            "reed {}\ncommand: {}\ninvariant-sites: {:?}\nwindows: {}\nrecords: {}\n\
             invariant sites seen: {}\nmultiallelic skipped: {}\n{methods}",
            self.version,
            self.command,
            invariant,
            n_windows,
            n_records,
            n_invariant_sites,
            n_multiallelic_skipped
        );
        if let Some(sfs) = sfs_report {
            let _ = writeln!(
                log,
                "sfs sites absorbed: {}\nsfs skipped (m < n): {}",
                sfs.n_sites.iter().sum::<u64>(),
                sfs.n_skipped_m_lt_n.iter().sum::<u64>()
            );
        }
        let log_path = with_suffix(prefix, ".log");
        reed_io::write::atomic::write_atomic(&log_path, log).map_err(EngineError::Io)?;
        files.push(log_path);
        Ok(files)
    }
}

fn value_or_na(value: Option<f64>, empty: bool) -> String {
    if empty { "NA".to_string() } else { cell(value) }
}

fn cell(value: Option<f64>) -> String {
    value.map_or_else(|| "NA".to_string(), |v| format!("{v}"))
}

fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

fn pbs_stem(a: &str, b: &str, c: &str) -> String {
    let mut stem = String::with_capacity(a.len() + b.len() + c.len() + 2);
    for (index, name) in [a, b, c].iter().enumerate() {
        if index > 0 {
            stem.push('-');
        }
        for ch in name.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | '-') {
                stem.push(ch);
            } else {
                stem.push('_');
            }
        }
    }
    stem
}

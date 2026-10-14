//! `StatsRequest::to_args`: the text contract (doc_v5 08 §6.2) for the
//! statistics flags — render a request back to the CLI spelling clap
//! accepts. The name-based spellings (`--pbs A,B,C`, `--f2 A,B`) need the
//! population assignment order, so unlike `FilterSpec::to_args` this takes
//! the group names.
//!
//! Not every request has a lossless flag spelling (the wide-table kind set,
//! per-bp kinds, a `Regions` grid); those return a usage error instead of
//! panicking (doc_v11 N3), so `to_args` composes with the crate's own
//! `default_kinds`/`window_kinds` builders.

use reed_stats::fst::FstKind;
use reed_stats::pbs::PbsKind;
use reed_stats::request::StatKind;

use reed_core::args::push_flag as flag;
use reed_core::grid::WindowGrid;

use crate::stats::request::{FStatTest, InvariantMode, Levels, Populations, StatsRequest};

impl StatsRequest {
    /// Renders the windowed form of this request (`layout: Some`) to the
    /// flags `reed stats` accepts; `group_names` is the assignment order
    /// (`--pop` resolved), required by the name-based spellings.
    ///
    /// Not rendered (no lossless flag spelling): `keep_samples` (a resolved
    /// mask), `regions`/`targets`/`targets_exclude`, the header-merged
    /// `contig_lengths`, a `Regions` grid (the grid carries no source
    /// path), and the wide-table kind set (`layout: None`).
    ///
    /// # Errors
    ///
    /// A usage error for a request the CLI cannot produce: a kind with no
    /// flag spelling (`per_bp` variants, `StatKind::Sfs`, mixed Fst or Pbs
    /// flavours within one request), an `LdRequest` with `decay_max_bp:
    /// None` or `decay_bin_bp != 1000`, or a `Regions` grid.
    pub fn to_args(&self, group_names: &[String]) -> Result<Vec<String>, crate::error::Error> {
        let name = |i: usize| -> String {
            group_names
                .get(i)
                .map_or_else(|| format!("group{i}"), Clone::clone)
        };
        let mut args: Vec<String> = Vec::new();
        if let Populations::File(path) = &self.populations {
            flag(&mut args, "--pop", path.display());
        }
        // ---- kinds: one scan, aggregate the flag family ----
        let mut has_pi = false;
        let mut has_dxy = false;
        let mut has_theta = false;
        let mut has_tajima = false;
        let mut fst_kind: Option<FstKind> = None;
        let mut fst_global = false;
        let mut pbs_n1 = false;
        let mut triples: Vec<String> = Vec::new();
        for kind in &self.kinds {
            match kind {
                StatKind::Pi { per_bp, .. } => {
                    if *per_bp {
                        return Err(crate::error::Error::usage(
                            "no --flag spelling for per-bp Pi",
                        ));
                    }
                    has_pi = true;
                }
                StatKind::Dxy { per_bp, .. } => {
                    if *per_bp {
                        return Err(crate::error::Error::usage(
                            "no --flag spelling for per-bp Dxy",
                        ));
                    }
                    has_dxy = true;
                }
                StatKind::ThetaW { per_bp, .. } => {
                    if *per_bp {
                        return Err(crate::error::Error::usage(
                            "no --flag spelling for per-bp ThetaW",
                        ));
                    }
                    has_theta = true;
                }
                StatKind::TajimaD { .. } => {
                    has_tajima = true;
                }
                StatKind::Fst { kind, .. } => {
                    if fst_kind.is_some_and(|seen| seen != *kind) {
                        return Err(crate::error::Error::usage(
                            "mixed Fst flavours have no flag spelling",
                        ));
                    }
                    fst_kind = Some(*kind);
                }
                StatKind::FstGlobal => fst_global = true,
                StatKind::Pbs { a, b, c, kind } => {
                    if pbs_n1 && *kind != PbsKind::PbsN1 {
                        return Err(crate::error::Error::usage(
                            "mixed Pbs flavours have no flag spelling",
                        ));
                    }
                    pbs_n1 |= *kind == PbsKind::PbsN1;
                    triples.push(format!("{},{},{}", name(*a), name(*b), name(*c)));
                }
                StatKind::Sfs { .. } => {
                    return Err(crate::error::Error::usage(
                        "no --sfs-kind flag spelling (use SfsRequest)",
                    ));
                }
            }
        }
        if has_pi {
            args.push("--pi".to_string());
        }
        if has_dxy {
            args.push("--dxy".to_string());
        }
        if has_theta {
            args.push("--theta".to_string());
        }
        if has_tajima {
            args.push("--tajima-d".to_string());
        }
        if let Some(kind) = fst_kind {
            flag(
                &mut args,
                "--fst",
                match kind {
                    FstKind::WeirCockerham84 => "wc84",
                    FstKind::Hudson => "hudson",
                },
            );
        } else if fst_global {
            // clap requires --fst alongside --fst-global; the value is
            // consumed by the global form (it yields no pairwise Fst kind),
            // so the round trip stays closed.
            flag(&mut args, "--fst", "wc84");
        }
        if fst_global {
            args.push("--fst-global".to_string());
        }
        for triple in &triples {
            flag(&mut args, "--pbs", triple);
        }
        if pbs_n1 {
            args.push("--pbs-normed".to_string());
        }
        // ---- sfs ----
        if let Some(sfs) = &self.sfs {
            flag(
                &mut args,
                "--sfs",
                if sfs.folded { "folded" } else { "unfolded" },
            );
            let projections: Vec<String> = sfs.projections.iter().map(u32::to_string).collect();
            flag(&mut args, "--proj", projections.join(","));
        }
        // ---- f-statistics ----
        if let Some(fstats) = &self.fstats {
            for test in &fstats.tests {
                let (flag_name, value): (&str, String) = match test {
                    FStatTest::F2 { a, b } => ("--f2", format!("{},{}", name(*a), name(*b))),
                    FStatTest::F3 { c, a, b } => {
                        ("--f3", format!("{};{},{}", name(*c), name(*a), name(*b)))
                    }
                    FStatTest::F4 { a, b, c, d } => (
                        "--f4",
                        format!("{},{};{},{}", name(*a), name(*b), name(*c), name(*d)),
                    ),
                    FStatTest::Dstat { a, b, c, d } => (
                        "--dstat",
                        format!("{},{};{},{}", name(*a), name(*b), name(*c), name(*d)),
                    ),
                    // Mirror of `build_fstats_request`: `A,B;C,D;X` →
                    // a = parts[0], o = parts[1], c = parts[2], b = parts[3],
                    // x = parts[4] (doc_v7 D7-8).
                    FStatTest::F4Ratio { o, a, x, b, c } => (
                        "--f4-ratio",
                        format!(
                            "{},{};{},{};{}",
                            name(*a),
                            name(*o),
                            name(*c),
                            name(*b),
                            name(*x)
                        ),
                    ),
                };
                flag(&mut args, flag_name, value);
            }
            flag(&mut args, "--block-size", fstats.block_bp);
        }
        // ---- ld ----
        if let Some(ld) = &self.ld {
            let Some(max) = ld.decay_max_bp else {
                return Err(crate::error::Error::usage(
                    "--ld-decay-less LD has no flag spelling",
                ));
            };
            if ld.decay_bin_bp != 1000 {
                return Err(crate::error::Error::usage(
                    "decay_bin_bp has no flag spelling",
                ));
            }
            flag(&mut args, "--ld-decay", max);
            if ld.pairs {
                args.push("--ld-pairs".to_string());
            }
        }
        // ---- grid (Bp only; Regions carries no source path) ----
        match &self.grid {
            WindowGrid::Bp { size, step } => {
                flag(&mut args, "--window", size);
                if step != size {
                    flag(&mut args, "--step", step);
                }
            }
            WindowGrid::Regions(_) => {
                return Err(crate::error::Error::usage(
                    "a Regions grid carries no --bed-file path to render",
                ));
            }
            WindowGrid::Sites { .. } => {
                return Err(crate::error::Error::usage(
                    "site-count windows are not implemented",
                ));
            }
        }
        // ---- levels / misc ----
        let default_levels = Levels::default();
        if self.levels != default_levels {
            let mut parts: Vec<&str> = Vec::new();
            if self.levels.window {
                parts.push("window");
            }
            if self.levels.contig {
                parts.push("contig");
            }
            if self.levels.genome {
                parts.push("genome");
            }
            flag(&mut args, "--levels", parts.join(","));
        }
        if self.skip_empty {
            args.push("--skip-empty".to_string());
        }
        if let Some(layout) = &self.layout {
            match layout.invariant {
                InvariantMode::VariantOnly => args.push("--variant-only".to_string()),
                InvariantMode::Bypass => args.push("--bypass-invariant-check".to_string()),
                InvariantMode::Require => {}
            }
        }
        if self.force_samples {
            args.push("--force-samples".to_string());
        }
        if let Some(spec) = &self.site {
            args.extend(spec.to_args());
        }
        Ok(args)
    }
}

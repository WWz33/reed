//! `reed distance`: the individual-by-individual matrix family (doc_v17).

use std::io::Write as _;

use reed::budget::ThreadBudget;
use reed_distance::matrices::{
    ajk, allele_ct_distance, discordance, grm, ibs, king, one_minus_ibs, phi, pi_hat,
    AjkOptions, DistanceOptions, DiscordanceMode, DiscordanceOptions, MatrixKind, MatrixResult,
    MissingPolicy, PiHatOptions,
};
use reed_distance::output::{is_matrix_kind, write_ids, write_matrix};
use reed_distance::{standardize, MatrixFormat, Stand};

use crate::cli::DistanceCli;

/// Which matrix the user asked for.
fn parse_kind(s: &str) -> Result<MatrixKind, String> {
    Ok(match s {
        "allele-ct" | "dist" => MatrixKind::AlleleCountDistance,
        "ibs" => MatrixKind::Ibs,
        "1-ibs" => MatrixKind::OneMinusIbs,
        "grm" | "rel" => MatrixKind::Grm,
        "cov" => MatrixKind::Covariance,
        "king" => MatrixKind::King,
        "pihat" | "genome" => MatrixKind::PiHat,
        "phi" | "relatedness2" => MatrixKind::Phi,
        "ajk" | "relatedness" => MatrixKind::Ajk,
        "discordance" | "gtcheck" => MatrixKind::Discordance,
        other => {
            return Err(format!(
                "unknown --type {other:?} \
                 (allele-ct|ibs|1-ibs|grm|cov|king|pihat|phi|ajk|discordance)"
            ));
        }
    })
}

/// Runs `reed distance`.
///
/// # Errors
///
/// Usage errors for bad flags, and I/O or engine failures.
pub(crate) fn run_distance_cli(args: DistanceCli) -> Result<(), reed::Error> {
    let kind = parse_kind(&args.matrix_type).map_err(reed::Error::usage)?;
    let format = MatrixFormat::parse(&args.format).map_err(reed::Error::usage)?;
    let stand = Stand::parse(&args.stand).map_err(reed::Error::usage)?;
    let Some(prefix) = args.output.as_deref() else {
        return Err(reed::Error::usage(
            "-o PREFIX is required: distance writes PREFIX.<suffix>",
        ));
    };

    let path = crate::units::input_path(&args.input);
    let opened = reed::open_regions(path.as_deref(), &[])?;
    let sample_names = opened.facts().sample_names.clone();
    if sample_names.is_empty() {
        return Err(reed::Error::usage(
            "the header declares no samples; nothing to compute",
        ));
    }
    let keep = super::common::parse_sample_keep(&args.samples, &sample_names)
        .map_err(reed::Error::usage)?;

    let threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .expect("invariant: parallelism")
            .get()
    });
    let budget = ThreadBudget {
        threads,
        mem: args
            .memory
            .map_or_else(reed::budget::MemBudget::default_half_memory, |m| {
                reed::budget::MemBudget::from_bytes(m.saturating_mul(1 << 20))
            }),
    };

    let (x, ids) = reed::read_genotype_matrix(opened, keep.as_deref(), &budget)?;
    if x.rows == 0 || x.cols == 0 {
        return Err(reed::Error::usage("no biallelic sites with a diploid GT"));
    }
    let st = standardize(&x, stand);

    let (mut result, aux) = match kind {
        MatrixKind::AlleleCountDistance => {
            (allele_ct_distance(&st, &DistanceOptions::default()), None)
        }
        MatrixKind::Ibs => (ibs(&st, &DistanceOptions::default()), None),
        MatrixKind::OneMinusIbs => (one_minus_ibs(&st, &DistanceOptions::default()), None),
        MatrixKind::Grm | MatrixKind::Covariance => {
            let policy = if args.mean_impute {
                MissingPolicy::MeanImpute
            } else {
                MissingPolicy::Pairwise
            };
            (grm(&st, policy, kind == MatrixKind::Grm), None)
        }
        MatrixKind::King => {
            let (m, _counts) = king(&st);
            (m, None)
        }
        MatrixKind::PiHat => {
            let opts = PiHatOptions { unbounded: args.unbounded, nudge: args.nudge };
            let (m, _res) = pi_hat(&st, &opts);
            (m, None)
        }
        MatrixKind::Phi => {
            let (m, a) = phi(&st);
            (m, Some(a))
        }
        MatrixKind::Ajk => (ajk(&st, &AjkOptions::default()), None),
        MatrixKind::Discordance => {
            let opts = DiscordanceOptions {
                mode: if args.mismatch_count {
                    DiscordanceMode::Count
                } else {
                    DiscordanceMode::Probability
                },
                error_probability: args.error_probability,
            };
            let (m, a) = discordance(&st, &opts);
            (m, Some(a))
        }
    };
    result.ids.clone_from(&ids);
    // `pi_hat` fills its own aux (the Z triplet); only overwrite when the
    // estimator returned one explicitly.
    if aux.is_some() {
        result.aux = aux;
    }

    let suffix = kind.suffix();
    let out_path = std::path::PathBuf::from(format!("{prefix}.{suffix}"));

    if is_matrix_kind(kind) {
        write_matrix(&out_path, &result, format).map_err(reed::Error::from)?;
        write_ids(
            &std::path::PathBuf::from(format!("{prefix}.{suffix}.id")),
            &ids,
        )
        .map_err(reed::Error::from)?;
    } else {
        write_table(&out_path, &result)?;
    }

    eprintln!(
        "reed distance: {} samples, {} sites, type={}, format={}",
        x.rows,
        x.cols,
        args.matrix_type,
        format.name()
    );
    Ok(())
}

/// Writes the per-pair table formats (`.genome`, `.relatedness2`,
/// `.relatedness`, `.gtcheck`) with the reference tools' exact headers.
fn write_table(path: &std::path::Path, m: &MatrixResult) -> Result<(), reed::Error> {
    let mut f = std::fs::File::create(path).map_err(reed::Error::from)?;
    let n = m.n;
    match m.kind {
        MatrixKind::Phi => {
            writeln!(f, "INDV1\tINDV2\tN_AaAa\tN_AAaa\tN1_Aa\tN2_Aa\tRELATEDNESS_PHI")
                .map_err(reed::Error::from)?;
            let aux = m.aux.as_ref();
            // vcftools emits the full ordered matrix, both directions.
            for i in 0..n {
                for j in 0..n {
                    let hethet = aux.map_or(0, |a| a.n_hethet[i * n + j]);
                    let opp = aux.map_or(0, |a| a.n_opposite_hom[i * n + j]);
                    // N_Aa[i] is the marginal het count; recover it from the
                    // diagonal of n_hethet is wrong, so recompute from phi's
                    // denominator is not possible here. Write the pair counts
                    // that vcftools reports and the phi value.
                    let phi_v = m.values[i * n + j];
                    let n1 = aux.map_or(0, |a| a.n_het.get(i).copied().unwrap_or(0));
                    let n2 = aux.map_or(0, |a| a.n_het.get(j).copied().unwrap_or(0));
                    writeln!(
                        f,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{:.9}",
                        m.ids[i], m.ids[j], hethet, opp, n1, n2, phi_v
                    )
                    .map_err(reed::Error::from)?;
                }
            }
        }
        MatrixKind::Ajk => {
            writeln!(f, "INDV1\tINDV2\tRELATEDNESS_AJK").map_err(reed::Error::from)?;
            for i in 0..n {
                for j in i..n {
                    writeln!(f, "{}\t{}\t{:.9}", m.ids[i], m.ids[j], m.values[i * n + j])
                        .map_err(reed::Error::from)?;
                }
            }
        }
        MatrixKind::Discordance => {
            writeln!(
                f,
                "#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t\
                 [5]Average -log P(HWE)\t[6]Number of sites compared\t\
                 [7]Number of matching genotypes"
            )
            .map_err(reed::Error::from)?;
            let aux = m.aux.as_ref();
            for i in 0..n {
                for j in 0..i {
                    let cmp = aux.map_or(0, |a| a.n_compared[i * n + j]);
                    let mt = aux.map_or(0, |a| a.n_matching[i * n + j]);
                    writeln!(
                        f,
                        "DCv2\t{}\t{}\t{:.9e}\t0\t{}\t{}",
                        m.ids[i], m.ids[j], m.values[i * n + j], cmp, mt
                    )
                    .map_err(reed::Error::from)?;
                }
            }
        }
        MatrixKind::PiHat => {
            writeln!(
                f,
                "FID1\tIID1\tFID2\tIID2\tRT\tEZ\tZ0\tZ1\tZ2\tPI_HAT\tPHE\tDST\tPPC\tRATIO"
            )
            .map_err(reed::Error::from)?;
            let aux = m.aux.as_ref();
            // plink writes FID and IID as separate columns; reed's ids are
            // the combined "FID_IID" sample names, so split on the last '_'
            // when it is present.
            let split = |id: &str| -> (String, String) {
                match id.rsplit_once('_') {
                    Some((fid, iid)) => (fid.to_string(), iid.to_string()),
                    None => (id.to_string(), id.to_string()),
                }
            };
            for i in 0..n {
                for j in 0..i {
                    let (f1, i1) = split(&m.ids[i]);
                    let (f2, i2) = split(&m.ids[j]);
                    let z0 = aux.map_or(f64::NAN, |a| a.z0[i * n + j]);
                    let z1 = aux.map_or(f64::NAN, |a| a.z1[i * n + j]);
                    let z2 = aux.map_or(f64::NAN, |a| a.z2[i * n + j]);
                    writeln!(
                        f,
                        "{f1}\t{i1}\t{f2}\t{i2}\tOT\tNA\t{z0:.4}\t{z1:.4}\t{z2:.4}\t{:.4}\t-9\t0\t0\t0",
                        m.values[i * n + j]
                    )
                    .map_err(reed::Error::from)?;
                }
            }
        }
        _ => unreachable!("is_matrix_kind covers the rest"),
    }
    Ok(())
}

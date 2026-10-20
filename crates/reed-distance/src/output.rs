//! Writing a matrix in the formats the reference tools use.

use std::io::Write;

use crate::matrices::{MatrixKind, MatrixResult};

/// Output layout, mirroring plink's modifiers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixFormat {
    /// Lower triangle, row-major: row 1 holds `d(1,2)`, row 2 holds
    /// `d(1,3) d(2,3)`, and so on. This is plink's default.
    Triangle,
    /// Full symmetric square matrix.
    Square,
    /// Full square with the upper-right triangle zeroed (plink `square0`).
    Square0,
    /// Binary double-precision square matrix, no header (R `readBin`).
    Bin,
}

impl MatrixFormat {
    /// Parses a CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "triangle" => Ok(MatrixFormat::Triangle),
            "square" => Ok(MatrixFormat::Square),
            "square0" => Ok(MatrixFormat::Square0),
            "bin" => Ok(MatrixFormat::Bin),
            other => Err(format!(
                "unknown matrix format {other:?} (triangle|square|square0|bin)"
            )),
        }
    }

    /// The CLI spelling.
    pub fn name(self) -> &'static str {
        match self {
            MatrixFormat::Triangle => "triangle",
            MatrixFormat::Square => "square",
            MatrixFormat::Square0 => "square0",
            MatrixFormat::Bin => "bin",
        }
    }
}

/// Writes `m` to `path` in `format`.
///
/// # Errors
///
/// I/O failures.
pub fn write_matrix(
    path: &std::path::Path,
    m: &MatrixResult,
    format: MatrixFormat,
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    match format {
        MatrixFormat::Bin => write_bin(&mut f, m),
        MatrixFormat::Triangle => write_triangle(&mut f, m),
        MatrixFormat::Square => write_square(&mut f, m, false),
        MatrixFormat::Square0 => write_square(&mut f, m, true),
    }
}

fn write_bin(f: &mut std::fs::File, m: &MatrixResult) -> std::io::Result<()> {
    let mut buf = Vec::with_capacity(m.values.len() * 8);
    for v in &m.values {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    f.write_all(&buf)
}

/// plink's triangle layout **includes the diagonal**: row `i` holds
/// `d(i,0) .. d(i,i)`, so row 0 has one value and row `n-1` has `n`.
fn write_triangle(f: &mut std::fs::File, m: &MatrixResult) -> std::io::Result<()> {
    for i in 0..m.n {
        for j in 0..=i {
            if j > 0 {
                write!(f, "\t")?;
            }
            write!(f, "{:.9}", m.values[i * m.n + j])?;
        }
        writeln!(f)?;
    }
    Ok(())
}

fn write_square(f: &mut std::fs::File, m: &MatrixResult, zero_upper: bool) -> std::io::Result<()> {
    for i in 0..m.n {
        for j in 0..m.n {
            if j > 0 {
                write!(f, "\t")?;
            }
            let v = if zero_upper && j > i {
                0.0
            } else {
                m.values[i * m.n + j]
            };
            write!(f, "{v:.9}")?;
        }
        writeln!(f)?;
    }
    Ok(())
}

/// Writes the `.id` sidecar (two columns, `FID IID`; reed uses the sample
/// name for both, matching `reed pca`).
///
/// # Errors
///
/// I/O failures.
pub fn write_ids(path: &std::path::Path, ids: &[String]) -> std::io::Result<()> {
    let mut f = std::fs::File::create(path)?;
    for id in ids {
        writeln!(f, "{id}\t{id}")?;
    }
    Ok(())
}

/// Whether a kind is written as a matrix (triangle/square/bin) rather than as
/// a per-pair table with a header.
pub fn is_matrix_kind(kind: MatrixKind) -> bool {
    matches!(
        kind,
        MatrixKind::AlleleCountDistance
            | MatrixKind::Ibs
            | MatrixKind::OneMinusIbs
            | MatrixKind::Grm
            | MatrixKind::Covariance
            | MatrixKind::King
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrices::{MatrixKind, MatrixResult};

    fn sample() -> MatrixResult {
        let mut m = MatrixResult {
            values: vec![1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0],
            n: 3,
            kind: MatrixKind::Grm,
            ids: vec!["a".into(), "b".into(), "c".into()],
            aux: None,
        };
        m.values[0] = 1.0;
        m
    }

    fn tmp(name: &str) -> std::path::PathBuf {
        let p = std::env::temp_dir().join(format!("reed-distance-{name}-{}", std::process::id()));
        let _ = std::fs::remove_file(&p);
        p
    }

    #[test]
    fn triangle_has_one_fewer_entry_per_row() {
        let p = tmp("tri");
        write_matrix(&p, &sample(), MatrixFormat::Triangle).expect("write");
        let text = std::fs::read_to_string(&p).expect("read");
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 3);
        // plink's triangle includes the diagonal.
        assert_eq!(lines[0].split('\t').count(), 1); // d(0,0)
        assert_eq!(lines[1].split('\t').count(), 2); // d(1,0) d(1,1)
        assert_eq!(lines[2].split('\t').count(), 3); // d(2,0) d(2,1) d(2,2)
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn square0_zeroes_the_upper_triangle() {
        let p = tmp("sq0");
        write_matrix(&p, &sample(), MatrixFormat::Square0).expect("write");
        let text = std::fs::read_to_string(&p).expect("read");
        let first: Vec<&str> = text.lines().next().unwrap().split('\t').collect();
        assert_eq!(first[0], "1.000000000");
        assert_eq!(first[1], "0.000000000");
        assert_eq!(first[2], "0.000000000");
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn bin_is_n_squared_doubles() {
        let p = tmp("bin");
        write_matrix(&p, &sample(), MatrixFormat::Bin).expect("write");
        let bytes = std::fs::read(&p).expect("read");
        assert_eq!(bytes.len(), 9 * 8);
        let first = f64::from_le_bytes(bytes[0..8].try_into().unwrap());
        assert!((first - 1.0).abs() < 1e-12);
        let _ = std::fs::remove_file(&p);
    }
}

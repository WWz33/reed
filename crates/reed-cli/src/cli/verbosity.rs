//! The shared `--verbosity` contract (doc_v4 03-architecture.md §3.3, 01
//! F4-3). One enum for every subcommand; the driver returns a report and the
//! CLI decides what to print.

/// Verbosity levels.
///
/// | level | behaviour |
/// |---|---|
/// | 0 `silent` | nothing on stderr |
/// | 1 `summary` (default) | one summary line (+ the funnel table for `filter`) |
/// | 2 `timing` | + the `--timing` breakdown on stderr |
/// | 3 `debug` | + the `RUST_LOG=debug` tracing stream |
///
/// Levels above 3 are rejected at parse time: nothing is implemented for
/// them, and accepting them would be a silent no-op.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Verbosity {
    /// 0: no stderr output.
    Silent,
    /// 1: the summary line.
    #[default]
    Summary,
    /// 2: summary + phase timing.
    Timing,
    /// 3: timing + tracing.
    Debug,
}

impl Verbosity {
    /// Parses the `--verbosity` value; the byte above 3 is rejected.
    ///
    /// # Errors
    ///
    /// A message naming the supported range (clap turns it into exit code 2).
    pub(crate) fn parse(text: &str) -> Result<Self, String> {
        let level: u8 = text.trim().parse().map_err(|_| {
            format!("invalid --verbosity {text:?}: expected an integer 0-3")
        })?;
        match level {
            0 => Ok(Self::Silent),
            1 => Ok(Self::Summary),
            2 => Ok(Self::Timing),
            3 => Ok(Self::Debug),
            other => Err(format!(
                "unsupported --verbosity {other}: levels are 0 silent, 1 summary, \
                 2 stage timing, 3 debug"
            )),
        }
    }

    /// Whether the summary line is printed.
    pub(crate) fn shows_summary(self) -> bool {
        self >= Self::Summary
    }

    /// Whether the timing breakdown is printed.
    pub(crate) fn shows_timing(self) -> bool {
        self >= Self::Timing
    }

    /// Initializes the tracing subscriber. Level 0 installs nothing; levels
    /// 1–2 print the backend's `warn`+ events (08 §6.1: the backend only
    /// emits `tracing` events, the front end decides whether to show them);
    /// level 3 uses `RUST_LOG` (default `debug`).
    pub(crate) fn init_tracing(self) {
        if self < Self::Summary {
            return;
        }
        use tracing_subscriber::EnvFilter;
        let filter = if self >= Self::Debug {
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("debug"))
        } else {
            EnvFilter::new("warn")
        };
        let _ = tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_writer(std::io::stderr)
            .try_init();
    }
}

impl std::fmt::Display for Verbosity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let level = match self {
            Self::Silent => 0,
            Self::Summary => 1,
            Self::Timing => 2,
            Self::Debug => 3,
        };
        write!(f, "{level}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levels_map_and_above_three_is_rejected() {
        assert_eq!(Verbosity::parse("0").expect("0"), Verbosity::Silent);
        assert_eq!(Verbosity::parse("3").expect("3"), Verbosity::Debug);
        let err = Verbosity::parse("4").expect_err("4 is not implemented");
        assert!(err.contains("0 silent"), "{err}");
        assert!(Verbosity::Summary.shows_summary());
        assert!(!Verbosity::Summary.shows_timing());
        assert!(Verbosity::Timing.shows_timing());
        assert!(Verbosity::Silent < Verbosity::Summary);
    }
}

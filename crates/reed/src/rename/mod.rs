//! `reed rename` (doc_v14 02 §4): the identifier editor. It touches only the
//! header lines, the CHROM column, the ID column and the sample names —
//! genotypes are never decoded (D14-4), which is the structural guard this
//! module is built around: there is no `GenotypeScanner`/`GtCode` here.
//!
//! Evaluation order (02 §4.3): header replacement → contig rename → ID
//! operations → sample rename → header meta edits. VCF input takes the
//! byte-edit path (records pass through except for the edited columns);
//! BCF input re-encodes the header block (and the contig dictionary when
//! contigs are renamed).

use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::chunk::ChunkBytes;
use reed_io::vcf::line::LineRef;

use crate::input::OpenedInput;

/// The identifier edits (doc_v14 02 §4.2), applied in the fixed order of
/// §4.3. Every field is empty unless the operator asked for it.
#[derive(Clone, Debug, Default)]
pub struct RenameSpec {
    /// `-h/--header FILE`: replace the whole header (must keep the samples).
    pub header: Option<PathBuf>,
    /// `-f/--fai FILE`: set `##contig` lengths from a `.fai`.
    pub fai: Option<PathBuf>,
    /// `--add-meta LINE` (repeatable): append `##key=value` lines.
    pub add_meta: Vec<String>,
    /// `--remove-meta KEY` (repeatable): drop `##KEY=` lines.
    pub remove_meta: Vec<String>,
    /// `--rename-chrs FILE`: `old<TAB>new`, applied to `##contig` and CHROM.
    pub rename_chrs: Option<PathBuf>,
    /// `--chr-style`: the naming convention to normalize to.
    pub chr_style: Option<ChrStyle>,
    /// `--set-id PATTERN`: rewrite every ID (`%CHROM %POS %REF %ALT`).
    pub set_id: Option<String>,
    /// `--set-missing-var-ids PATTERN`: rewrite only `.` IDs.
    pub set_missing_var_ids: Option<String>,
    /// `--update-name FILE`: `old<TAB>new` ID mapping.
    pub update_name: Option<PathBuf>,
    /// `--idspace-to SEP`: replace the ID separator.
    pub idspace_to: Option<String>,
    /// `--rm-dup-ids`: suffix repeated IDs with `_1`, `_2`, ….
    pub rm_dup_ids: bool,
    /// `--clear-ids`: set every ID to `.`.
    pub clear_ids: bool,
    /// `-s/--samples PAIRS`: `old new,old2 new2`.
    pub samples: Option<String>,
    /// `-N/--samples-file FILE`: `old<TAB>new` (or `old new`) pairs.
    pub samples_file: Option<PathBuf>,
    /// `--replace-sample-name REGEX`: rewrite names matching the pattern.
    pub replace_sample_name: Option<String>,
}

/// The chromosome naming convention `--chr-style` normalizes to (B14-4: a
/// finite three-state enum, not an arbitrary mapping).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChrStyle {
    /// No `chr` prefix (`1`, `MT`).
    Plain,
    /// `chr` prefix (`chr1`, `chrM`).
    Ucsc,
    /// `chr` prefix with mitochondrial as `chrM` (`mt-is-m`).
    MtIsM,
}

impl ChrStyle {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// A usage message for an unknown style.
    pub fn parse(text: &str) -> Result<Self, String> {
        match text {
            "plain" => Ok(Self::Plain),
            "ucsc" => Ok(Self::Ucsc),
            "mt-is-m" => Ok(Self::MtIsM),
            other => Err(format!(
                "unknown --chr-style {other:?} (expected plain, ucsc or mt-is-m)"
            )),
        }
    }

    /// Applies the convention to one contig name.
    #[must_use]
    pub fn apply(self, name: &str) -> String {
        let bare = name.strip_prefix("chr").unwrap_or(name);
        match self {
            Self::Plain => bare.to_string(),
            Self::Ucsc => format!("chr{bare}"),
            Self::MtIsM => {
                let mito = matches!(bare, "M" | "MT" | "m" | "mt");
                if mito {
                    "chrM".to_string()
                } else {
                    format!("chr{bare}")
                }
            }
        }
    }
}

/// The output target (the `-O`/`-o` flags, filter's spellings).
#[derive(Clone, Debug)]
pub struct RenameOutput {
    /// Output path; `None` = stdout (plain VCF text).
    pub path: Option<PathBuf>,
    /// BGZF-compress the output.
    pub compressed: bool,
    /// BGZF compression level.
    pub level: u8,
    /// The `##reed_command` line (captured argv).
    pub command_line: String,
    /// Suppress the `##reed_*` provenance lines.
    pub no_version: bool,
    /// Suppress the `##reed_command` line only (`--no-command`).
    pub no_command: bool,
    /// `-W[=csi|tbi]`: write an index next to the output (filter's
    /// semantics; omitted = no index).
    pub write_index: Option<crate::filter::WriteIndex>,
}

/// What the run did: the files written and the per-plane edit counts.
#[derive(Clone, Debug, Default)]
pub struct RenameReport {
    /// Files written (empty for stdout).
    pub files: Vec<PathBuf>,
    /// Records rewritten.
    pub n_records: u64,
    /// Contig names renamed (header + records).
    pub n_contigs: u64,
    /// Sample names renamed.
    pub n_samples: u64,
    /// IDs rewritten (any ID plane).
    pub n_ids: u64,
    /// IDs that collided and got a `_N` suffix (`--rm-dup-ids`).
    pub n_dup_ids: u64,
    /// A pre-existing index next to the output that renaming made stale;
    /// the CLI prints the warning (the backend never prints, 08 §6.4.2).
    pub stale_index: Option<PathBuf>,
}

/// The compiled, per-record plan (doc_v14 02 §4.3): the ID plane resolved
/// once, so the worker does table lookups and a template render per record.
struct CompiledPlan {
    /// Contig name → header index (the block writer's `ContigId`), keyed by
    /// the *post-rename* name so records resolve after the CHROM edit.
    contigs: HashMap<Vec<u8>, reed_core::ContigId>,
    chrs: HashMap<Vec<u8>, Vec<u8>>,
    chr_style: Option<ChrStyle>,
    set_id: Option<IdTemplate>,
    set_missing_var_ids: Option<IdTemplate>,
    /// ID mapping with the slot of each key in `update_keys` (the hit flags
    /// are indexed by slot).
    update_name: HashMap<Vec<u8>, (Vec<u8>, usize)>,
    /// The `--update-name` keys in file order and one hit flag each; a key
    /// that never matched an ID is a typo and fails the run (02 §4.3: a
    /// misspelt name must not be skipped silently). The flags are shared by
    /// the workers, so no per-chunk merge is needed.
    update_keys: Vec<Vec<u8>>,
    update_hits: Vec<AtomicBool>,
    idspace_to: Option<Vec<u8>>,
    rm_dup_ids: bool,
    clear_ids: bool,
}

/// A parsed `%CHROM:%POS:%REF:%ALT` template: literal and placeholder runs.
struct IdTemplate {
    parts: Vec<TemplatePart>,
}

enum TemplatePart {
    Literal(Vec<u8>),
    Chrom,
    Pos,
    Ref,
    Alt,
}

impl CompiledPlan {
    /// The header index of a (possibly renamed) contig name.
    fn contig_id(&self, name: &[u8]) -> Option<reed_core::ContigId> {
        self.contigs.get(name).copied()
    }
}

impl IdTemplate {
    /// Parses a template; `%%` is a literal `%`.
    ///
    /// # Errors
    ///
    /// A usage message for an unknown placeholder.
    fn parse(text: &str) -> Result<Self, String> {
        let bytes = text.as_bytes();
        let mut parts = Vec::new();
        let mut literal = Vec::new();
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] != b'%' {
                literal.push(bytes[i]);
                i += 1;
                continue;
            }
            if i + 1 < bytes.len() && bytes[i + 1] == b'%' {
                literal.push(b'%');
                i += 2;
                continue;
            }
            if !literal.is_empty() {
                parts.push(TemplatePart::Literal(std::mem::take(&mut literal)));
            }
            let rest = &text[i + 1..];
            let (part, used) = if rest.starts_with("CHROM") {
                (TemplatePart::Chrom, 5)
            } else if rest.starts_with("POS") {
                (TemplatePart::Pos, 3)
            } else if rest.starts_with("REF") {
                (TemplatePart::Ref, 3)
            } else if rest.starts_with("ALT") {
                (TemplatePart::Alt, 3)
            } else {
                return Err(format!(
                    "--set-id: unknown placeholder in {text:?} (supported: %CHROM %POS %REF \
                     %ALT, %% for a literal %)"
                ));
            };
            parts.push(part);
            i += 1 + used;
        }
        if !literal.is_empty() {
            parts.push(TemplatePart::Literal(literal));
        }
        Ok(Self { parts })
    }

    /// Renders the ID for one record; `chrom` is the *post-rename* CHROM
    /// (the ID plane runs after the contig plane, 02 §4.3).
    fn render(&self, chrom: &[u8], lr: &LineRef<'_>, pos1: u64) -> Vec<u8> {
        let mut out = Vec::new();
        for part in &self.parts {
            match part {
                TemplatePart::Literal(text) => out.extend_from_slice(text),
                TemplatePart::Chrom => out.extend_from_slice(chrom),
                TemplatePart::Pos => out.extend_from_slice(pos1.to_string().as_bytes()),
                TemplatePart::Ref => out.extend_from_slice(lr.reference()),
                TemplatePart::Alt => out.extend_from_slice(lr.alts()),
            }
        }
        out
    }
}

/// Reads a two-column `old<TAB>new` mapping file (whitespace also accepted).
///
/// # Errors
///
/// A usage message for an unreadable or malformed file.
fn read_pairs(path: &Path, flag: &str) -> Result<Vec<(String, String)>, crate::Error> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| crate::Error::usage(format!("{flag}: cannot read {}: {e}", path.display())))?;
    let mut pairs = Vec::new();
    for (line_no, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut cols = line.split_whitespace();
        let (Some(old), Some(new)) = (cols.next(), cols.next()) else {
            return Err(crate::Error::usage(format!(
                "{flag}: {}:{}: expected `old<TAB>new`, got {line:?}",
                path.display(),
                line_no + 1
            )));
        };
        pairs.push((old.to_string(), new.to_string()));
    }
    Ok(pairs)
}

/// Parses `-s "old new,old2 new2"` into rename pairs.
///
/// # Errors
///
/// A usage message for a pair without exactly two names.
fn parse_inline_pairs(text: &str) -> Result<Vec<(String, String)>, crate::Error> {
    let mut pairs = Vec::new();
    for element in text.split(',').map(str::trim).filter(|e| !e.is_empty()) {
        let names: Vec<&str> = element.split_whitespace().collect();
        if names.len() != 2 {
            return Err(crate::Error::usage(format!(
                "-s wants `old new` pairs (comma-separated), got {element:?}"
            )));
        }
        pairs.push((names[0].to_string(), names[1].to_string()));
    }
    Ok(pairs)
}

/// The `--replace-sample-name` rule: `pre*` strips a prefix, `*suf` a
/// suffix, and any other pattern removes the text wherever it occurs. No
/// regex engine: a `$1`-style substitution is not supported, and the
/// workspace takes no dependency for one flag (B14-4's spirit). A pattern
/// carrying regex metacharacters is a usage error rather than a silent
/// no-op: the flag's help and cli.md §12.4 document the literal rule.
enum NameRewrite {
    Prefix(String),
    Suffix(String),
    Remove(String),
}

impl NameRewrite {
    /// Regex metacharacters the literal rule cannot honour; refusing them
    /// keeps `s[0-9]` or `^s` from looking like a successful rewrite.
    const METACHARS: &'static [char] = &['^', '$', '[', ']', '(', ')', '|', '+', '?', '{', '}', '\\'];

    /// Parses the CLI pattern.
    ///
    /// # Errors
    ///
    /// A usage message when the pattern uses regex syntax.
    fn parse(pattern: &str) -> Result<Self, crate::Error> {
        let body = pattern.strip_suffix('*').or_else(|| pattern.strip_prefix('*')).unwrap_or(pattern);
        if let Some(bad) = body.chars().find(|c| Self::METACHARS.contains(c)) {
            return Err(crate::Error::usage(format!(
                "--replace-sample-name: {bad:?} is regex syntax; the pattern is literal \
                 (`pre*` strips a prefix, `*suf` a suffix, anything else is deleted)"
            )));
        }
        Ok(if let Some(prefix) = pattern.strip_suffix('*') {
            Self::Prefix(prefix.to_string())
        } else if let Some(suffix) = pattern.strip_prefix('*') {
            Self::Suffix(suffix.to_string())
        } else {
            Self::Remove(pattern.to_string())
        })
    }

    /// Applies the rule to one name.
    fn apply(&self, name: &str) -> String {
        match self {
            Self::Prefix(prefix) => name.strip_prefix(prefix).unwrap_or(name).to_string(),
            Self::Suffix(suffix) => name.strip_suffix(suffix).unwrap_or(name).to_string(),
            Self::Remove(text) => name.replace(text.as_str(), ""),
        }
    }
}

/// Runs `reed rename`.
///
/// # Errors
///
/// Usage errors from the flag validation (02 §4.3) and I/O failures.
pub fn run_rename(
    input: Option<&Path>,
    spec: &RenameSpec,
    output: &RenameOutput,
    budget: &ThreadBudget,
) -> Result<RenameReport, crate::Error> {
    let opened = crate::input::open(input).map_err(EngineError::from)?;
    run_rename_opened(opened, spec, output, budget)
}

/// Runs `reed rename` over an already-opened input (the CLI opens once so
/// stdin streams).
///
/// # Errors
///
/// As [`run_rename`].
pub fn run_rename_opened(
    opened: OpenedInput,
    spec: &RenameSpec,
    output: &RenameOutput,
    budget: &ThreadBudget,
) -> Result<RenameReport, crate::Error> {
    validate(spec)?;
    let facts = opened.facts().clone();
    if opened.bcf_dict().is_some() {
        return Err(crate::Error::usage(
            "reed rename: BCF input is not supported yet; convert to VCF first \
             (`reed filter -Ov` or `bcftools view`)",
        ));
    }
    let sample_names = facts.sample_names.clone();
    let n_samples = sample_names.len();
    // A sites-only VCF (8 columns) has no parsable record form for this
    // editor: the CHROM/ID planes would silently do nothing while the run
    // "succeeded". Refuse instead (the header planes alone still work).
    let record_planes = spec.rename_chrs.is_some()
        || spec.chr_style.is_some()
        || spec.set_id.is_some()
        || spec.set_missing_var_ids.is_some()
        || spec.update_name.is_some()
        || spec.idspace_to.is_some()
        || spec.rm_dup_ids
        || spec.clear_ids;
    if n_samples == 0 && record_planes {
        return Err(crate::Error::usage(
            "the input has no sample columns (a sites-only VCF); the CHROM/ID planes need \
             9-column records — use the header planes only, or add the samples first",
        ));
    }

    // ---- header: the replacement, the contig edits and the meta edits ----
    let mut header_text = String::from_utf8_lossy(&facts.text).into_owned();
    if let Some(path) = &spec.header {
        let text = std::fs::read_to_string(path).map_err(|e| {
            crate::Error::usage(format!("--header: cannot read {}: {e}", path.display()))
        })?;
        // The replacement must keep every sample (02 §4.3: a header without
        // the samples would silently drop columns).
        let new_samples = header_samples(&text)?;
        if new_samples != sample_names {
            return Err(crate::Error::usage(format!(
                "--header: the replacement names {} samples but the input has {n_samples}; \
                 the sample set must not change (use -s/-N to rename)",
                new_samples.len()
            )));
        }
        header_text = text;
    }
    let mut contig_map: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    if let Some(path) = &spec.rename_chrs {
        for (old, new) in read_pairs(path, "--rename-chrs")? {
            if let Some(contig) = facts.contig_names.iter().find(|c| **c == old) {
                contig_map.insert(contig.as_bytes().to_vec(), new.into_bytes());
            } else {
                return Err(crate::Error::usage(format!(
                    "--rename-chrs: {old:?} is not a contig of the input"
                )));
            }
        }
    }
    if let Some(style) = spec.chr_style {
        // 02 §4.3 step 2: the explicit rename and the convention stack, so
        // the style applies to the *mapped* name (the two used to fight,
        // with the style silently discarding the rename).
        for name in &facts.contig_names {
            let base = contig_map.get(name.as_bytes()).map_or_else(
                || name.clone(),
                |new| String::from_utf8_lossy(new).into_owned(),
            );
            let renamed = style.apply(&base);
            if renamed.as_bytes() != base.as_bytes() {
                contig_map.insert(name.as_bytes().to_vec(), renamed.into_bytes());
            }
        }
    }
    // Contigs whose name actually changed (a style that maps a contig onto
    // itself is not an edit).
    let n_contigs = u64::try_from(
        contig_map
            .iter()
            .filter(|(old, new)| old.as_slice() != new.as_slice())
            .count(),
    )
    .unwrap_or(u64::MAX);
    if !contig_map.is_empty() {
        header_text = rewrite_contig_lines(&header_text, &contig_map);
    }
    if let Some(path) = &spec.fai {
        let lengths = read_pairs(path, "--fai")?;
        header_text = apply_fai(&header_text, &lengths);
    }
    for key in &spec.remove_meta {
        header_text = remove_meta_lines(&header_text, key);
    }
    for line in &spec.add_meta {
        let line = line.strip_prefix("##").unwrap_or(line);
        if !header_text.ends_with('\n') {
            header_text.push('\n');
        }
        header_text.push_str("##");
        header_text.push_str(line);
        header_text.push('\n');
    }

    // ---- samples ----
    let mut sample_map: HashMap<String, String> = HashMap::new();
    if let Some(inline) = &spec.samples {
        for (old, new) in parse_inline_pairs(inline)? {
            add_sample_pair(&sample_names, &mut sample_map, &old, &new)?;
        }
    }
    if let Some(path) = &spec.samples_file {
        for (old, new) in read_pairs(path, "-N")? {
            add_sample_pair(&sample_names, &mut sample_map, &old, &new)?;
        }
    }
    if let Some(pattern) = &spec.replace_sample_name {
        let rewrite = NameRewrite::parse(pattern)?;
        for name in &sample_names {
            let new = rewrite.apply(name);
            if new != *name {
                add_sample_pair(&sample_names, &mut sample_map, name, &new)?;
            }
        }
    }
    if !sample_map.is_empty() {
        header_text = rewrite_sample_line(&header_text, &sample_map);
    }
    let n_samples_renamed = u64::try_from(sample_map.len()).unwrap_or(u64::MAX);

    // ---- provenance ----
    if !output.no_version {
        let mut prov = String::new();
        {
            use std::fmt::Write as _;
            let _ = writeln!(prov, "##reed_version={}", env!("CARGO_PKG_VERSION"));
            if !output.command_line.is_empty() && !output.no_command {
                let _ = writeln!(prov, "##reed_command={}", output.command_line);
            }
            let mut planes = Vec::new();
            if spec.header.is_some() {
                planes.push("header");
            }
            if n_contigs > 0 {
                planes.push("contigs");
            }
            if spec.set_id.is_some() || spec.set_missing_var_ids.is_some() {
                planes.push("ids");
            }
            if n_samples_renamed > 0 {
                planes.push("samples");
            }
            let _ = writeln!(prov, "##reed_rename={}", planes.join(","));
        }
        header_text = insert_before_chrom(&header_text, &prov);
    }

    // ---- per-record plan ----
    // The *written* header is the truth for the index and for the record
    // keys: `-h` may replace the contig set/order, and the renames above
    // change the names. Deriving from `facts` here made `-h` + `-W` build an
    // index that described the input's contigs (silently wrong queries).
    let contig_names: Vec<Vec<u8>> = header_contig_names(&header_text);
    let mut contigs: HashMap<Vec<u8>, reed_core::ContigId> = HashMap::new();
    for (idx, name) in contig_names.iter().enumerate() {
        contigs.insert(
            name.clone(),
            reed_core::ContigId(u32::try_from(idx).unwrap_or(u32::MAX)),
        );
    }
    // `--update-name`: the pairs in file order (the hit flags are indexed by
    // slot, so a typo can be named at the end of the run).
    let update_pairs: Vec<(Vec<u8>, Vec<u8>)> = match &spec.update_name {
        Some(path) => read_pairs(path, "--update-name")?
            .into_iter()
            .map(|(old, new)| (old.into_bytes(), new.into_bytes()))
            .collect(),
        None => Vec::new(),
    };
    let update_keys: Vec<Vec<u8>> = update_pairs.iter().map(|(old, _)| old.clone()).collect();
    let update_name: HashMap<Vec<u8>, (Vec<u8>, usize)> = update_pairs
        .into_iter()
        .enumerate()
        .map(|(slot, (old, new))| (old, (new, slot)))
        .collect();
    let plan = Arc::new(CompiledPlan {
        contigs,
        chrs: contig_map,
        chr_style: spec.chr_style,
        set_id: spec
            .set_id
            .as_deref()
            .map(IdTemplate::parse)
            .transpose()
            .map_err(crate::Error::usage)?,
        set_missing_var_ids: spec
            .set_missing_var_ids
            .as_deref()
            .map(IdTemplate::parse)
            .transpose()
            .map_err(crate::Error::usage)?,
        update_name,
        update_keys: update_keys.clone(),
        update_hits: update_keys.iter().map(|_| AtomicBool::new(false)).collect(),
        idspace_to: spec.idspace_to.as_ref().map(|s| s.as_bytes().to_vec()),
        rm_dup_ids: spec.rm_dup_ids,
        clear_ids: spec.clear_ids,
    });
    if plan.set_id.is_some() && plan.set_missing_var_ids.is_some() {
        return Err(crate::Error::usage(
            "--set-id and --set-missing-var-ids are mutually exclusive",
        ));
    }
    if plan.clear_ids
        && (plan.set_id.is_some()
            || !plan.update_name.is_empty()
            || plan.idspace_to.is_some()
            || plan.rm_dup_ids)
    {
        return Err(crate::Error::usage(
            "--clear-ids cannot be combined with the other ID flags",
        ));
    }
    if plan.set_id.is_some() && plan.idspace_to.is_some() {
        return Err(crate::Error::usage(
            "--idspace-to cannot be combined with --set-id (write the separator in the template)",
        ));
    }

    // ---- write ----
    let writer: Box<dyn std::io::Write + Send> = match &output.path {
        None => Box::new(std::io::BufWriter::with_capacity(
            1 << 20,
            std::io::stdout(),
        )),
        Some(p) => Box::new(std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(p).map_err(|e| {
                crate::Error::usage(format!("cannot create {}: {e}", p.display()))
            })?,
        )),
    };
    // `-W`: an index needs BGZF framing and a real file; the contig names
    // are the *written* header's (derived above), so the index matches the
    // file it ships with.
    let contig_names: Arc<[Vec<u8>]> = Arc::from(contig_names.into_boxed_slice());
    let indexable = output.compressed && output.path.is_some() && output.write_index.is_some();
    if indexable && contig_names.is_empty() {
        return Err(crate::Error::usage(
            "-W needs the header's ##contig lines and this header has none; add them with \
             -f FILE or -h FILE",
        ));
    }
    let need_key = indexable;
    let indexer = if indexable {
        let flavor = output.write_index.expect("invariant: checked above");
        Some(
            reed_io::write::index::Indexer::new(flavor, &contig_names)
                .map_err(EngineError::from)?,
        )
    } else {
        None
    };
    let index_path = output
        .path
        .as_ref()
        .filter(|_| indexable)
        .map(|p| {
            output
                .write_index
                .expect("invariant: checked above")
                .path_for(p)
        });
    // A pre-existing index next to the output would describe the old file;
    // the CLI reports it (the backend never prints, 08 §6.4.2).
    let stale_index = match (&output.path, &index_path) {
        (Some(path), None) => [".csi", ".tbi"].iter().find_map(|suffix| {
            let mut sidecar = path.as_os_str().to_os_string();
            sidecar.push(suffix);
            let sidecar = PathBuf::from(sidecar);
            sidecar.exists().then_some(sidecar)
        }),
        _ => None,
    };
    let mut out = reed_io::write::ordered_output::OrderedOutput::new(
        writer,
        output.compressed,
        output.level,
        indexer,
        index_path,
        contig_names,
    );
    out.write_header(header_text.as_bytes())
        .map_err(EngineError::from)?;

    let job = RenameJob {
        plan: Arc::clone(&plan),
        level: output.level,
        compressed: output.compressed,
        need_key,
    };
    let report_slot: Arc<std::sync::Mutex<Option<RenameReport>>> =
        Arc::new(std::sync::Mutex::new(None));
    let sink = RenameSink {
        out,
        plan: Arc::clone(&plan),
        seen_ids: std::collections::HashSet::new(),
        next_suffix: std::collections::HashMap::new(),
        dedup_writer: plan.rm_dup_ids.then(|| {
            reed_io::write::block_writer::BlockWriter::new(output.level, output.compressed)
                .expect("invariant: the level comes from the CLI clamp")
        }),
        report: RenameReport {
            files: output.path.iter().cloned().collect(),
            n_contigs,
            n_samples: n_samples_renamed,
            stale_index,
            ..RenameReport::default()
        },
        slot: Arc::clone(&report_slot),
    };
    let abort = reed_engine::abort::abort_for(budget);
    let mut opened = opened;
    let run_report = {
        let source = opened.source_mut();
        run(source, &job, sink, budget, &abort)?
    };
    let mut report = report_slot
        .lock()
        .expect("invariant: rename slot not poisoned")
        .take()
        .unwrap_or_default();
    // The sink's `finish` already counts the records; `sink_report` carries
    // the same number (the engine returns it).
    report.n_records = run_report.sink_report;
    // Every `--update-name` key must have matched an ID: a typo would
    // otherwise rename nothing and still exit 0 (02 §4.3). Only meaningful
    // when the ID plane actually consulted the map.
    if plan.set_id.is_none()
        && !plan.clear_ids
        && let Some((slot, key)) = plan
            .update_keys
            .iter()
            .enumerate()
            .find(|(slot, _)| !plan.update_hits[*slot].load(Ordering::Relaxed))
    {
        return Err(crate::Error::usage(format!(
            "--update-name: {:?} is not an ID of the input (key #{}); names are not \
             skipped silently",
            String::from_utf8_lossy(key),
            slot + 1
        )));
    }
    Ok(report)
}

/// The flag validation that runs before any I/O (02 §4.3).
fn validate(spec: &RenameSpec) -> Result<(), crate::Error> {
    let any = spec.header.is_some()
        || spec.fai.is_some()
        || !spec.add_meta.is_empty()
        || !spec.remove_meta.is_empty()
        || spec.rename_chrs.is_some()
        || spec.chr_style.is_some()
        || spec.set_id.is_some()
        || spec.set_missing_var_ids.is_some()
        || spec.update_name.is_some()
        || spec.idspace_to.is_some()
        || spec.rm_dup_ids
        || spec.clear_ids
        || spec.samples.is_some()
        || spec.samples_file.is_some()
        || spec.replace_sample_name.is_some();
    if !any {
        return Err(crate::Error::usage(
            "no rename requested: pass one of -h/-f/--add-meta/--remove-meta/--rename-chrs/\
             --chr-style/--set-id/--set-missing-var-ids/--update-name/--idspace-to/\
             --rm-dup-ids/--clear-ids/-s/-N/--replace-sample-name",
        ));
    }
    Ok(())
}

/// Adds one sample rename pair, refusing an unknown old name (02 §4.3: a
/// typo must not be skipped silently) or a collision with a kept name.
fn add_sample_pair(
    sample_names: &[String],
    map: &mut HashMap<String, String>,
    old: &str,
    new: &str,
) -> Result<(), crate::Error> {
    if !sample_names.iter().any(|n| n == old) {
        return Err(crate::Error::usage(format!(
            "sample {old:?} is not in the header"
        )));
    }
    if new.is_empty() {
        return Err(crate::Error::usage(format!(
            "sample {old:?}: the new name is empty"
        )));
    }
    // The final names must stay unique: a rename onto another sample's name
    // would create duplicate columns.
    let renamed: Vec<&str> = sample_names
        .iter()
        .map(|name| {
            if name == old {
                new
            } else {
                map.get(name).map_or(name.as_str(), String::as_str)
            }
        })
        .collect();
    if renamed.iter().filter(|n| **n == new).count() > 1 {
        return Err(crate::Error::usage(format!(
            "renaming {old:?} to {new:?} would duplicate a sample name"
        )));
    }
    map.insert(old.to_string(), new.to_string());
    Ok(())
}

/// The sample names of a `#CHROM` line in `header` text.
fn header_samples(header: &str) -> Result<Vec<String>, crate::Error> {
    for line in header.lines() {
        if let Some(rest) = line.strip_prefix("#CHROM") {
            let cols = rest.split('\t').skip(9);
            return Ok(cols.map(str::to_string).collect());
        }
    }
    Err(crate::Error::usage(
        "--header: the replacement has no #CHROM line",
    ))
}

/// Rewrites the `##contig` lines' IDs per `map`.
fn rewrite_contig_lines(header: &str, map: &HashMap<Vec<u8>, Vec<u8>>) -> String {
    let mut out = String::with_capacity(header.len());
    for line in header.split_inclusive('\n') {
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let id_end = rest.find(',').or_else(|| rest.find('>')).unwrap_or(rest.len());
            let id = &rest[..id_end];
            if let Some(new) = map.get(id.as_bytes()) {
                out.push_str("##contig=<ID=");
                out.push_str(&String::from_utf8_lossy(new));
                out.push_str(&rest[id_end..]);
                continue;
            }
        }
        out.push_str(line);
    }
    out
}

/// Applies `--fai` lengths to the `##contig` lines (adding a line when the
/// contig has none).
fn apply_fai(header: &str, lengths: &[(String, String)]) -> String {
    let mut out = String::with_capacity(header.len());
    let mut seen: Vec<&str> = Vec::new();
    for line in header.split_inclusive('\n') {
        let mut replaced = false;
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let id_end = rest.find(',').or_else(|| rest.find('>')).unwrap_or(rest.len());
            let id = &rest[..id_end];
            if let Some((_, len)) = lengths.iter().find(|(name, _)| name == id) {
                seen.push(id);
                // Replace an existing length= or append one.
                if let Some(start) = rest.find("length=") {
                    let value_start = start + "length=".len();
                    let value_end = rest[value_start..]
                        .find([',', '>'])
                        .map_or(rest.len(), |i| value_start + i);
                    out.push_str("##contig=<ID=");
                    out.push_str(&rest[..value_start]);
                    out.push_str(len);
                    out.push_str(&rest[value_end..]);
                } else {
                    let close = rest.find('>').unwrap_or(rest.len());
                    out.push_str("##contig=<ID=");
                    out.push_str(&rest[..close]);
                    let _ = write!(out, ",length={len}");
                    out.push_str(&rest[close..]);
                }
                replaced = true;
            }
        }
        if !replaced {
            out.push_str(line);
        }
    }
    // Contigs the header lacks get a fresh line, *before* `#CHROM`: a `##`
    // line after it is not a valid VCF header (this used to append at the
    // end of the whole header text).
    let mut missing = String::new();
    for (name, len) in lengths {
        if !seen.iter().any(|s| *s == name) {
            let _ = writeln!(missing, "##contig=<ID={name},length={len}>");
        }
    }
    if missing.is_empty() {
        out
    } else {
        insert_before_chrom(&out, &missing)
    }
}

/// Drops every `##KEY=` line.
fn remove_meta_lines(header: &str, key: &str) -> String {
    let prefix = format!("##{key}=");
    header
        .split_inclusive('\n')
        .filter(|line| !line.starts_with(&prefix))
        .collect()
}

/// Inserts provenance text just before the `#CHROM` line.
fn insert_before_chrom(header: &str, text: &str) -> String {
    match header.find("#CHROM") {
        Some(pos) => {
            let mut out = String::with_capacity(header.len() + text.len());
            out.push_str(&header[..pos]);
            out.push_str(text);
            out.push_str(&header[pos..]);
            out
        }
        None => {
            let mut out = header.to_string();
            out.push_str(text);
            out
        }
    }
}

/// Splits a line into its body and its terminator (`\n`, `\r\n`, `\r` or
/// none), so the last column of a CRLF header is not read with a trailing
/// `\r` (which made sample renames silently miss).
fn split_line_end(line: &str) -> (&str, &str) {
    let body = line.trim_end_matches(['\n', '\r']);
    (&line[..body.len()], &line[body.len()..])
}

/// Rewrites the `#CHROM` line's sample columns per `map`.
fn rewrite_sample_line(header: &str, map: &HashMap<String, String>) -> String {
    let mut out = String::with_capacity(header.len());
    for line in header.split_inclusive('\n') {
        let Some(rest) = line.strip_prefix("#CHROM") else {
            out.push_str(line);
            continue;
        };
        out.push_str("#CHROM");
        let (body, newline) = split_line_end(rest);
        // `rest` begins with the tab before POS; split past it so the first
        // field is POS (no leading empty field).
        let body = body.strip_prefix('\t').unwrap_or(body);
        // POS..FORMAT are eight fixed fields; the rest are sample names.
        let mut fields = body.split('\t');
        for _ in 0..8 {
            match fields.next() {
                Some(field) => {
                    out.push('\t');
                    out.push_str(field);
                }
                None => break,
            }
        }
        for name in fields {
            out.push('\t');
            out.push_str(map.get(name).map_or(name, String::as_str));
        }
        out.push_str(newline);
    }
    out
}

/// The `##contig` IDs of a header text, in file order. The index and the
/// record keys follow the *written* header (a `-h` replacement or a contig
/// rename changes the set and the order).
fn header_contig_names(header: &str) -> Vec<Vec<u8>> {
    let mut names = Vec::new();
    for line in header.lines() {
        let Some(rest) = line.strip_prefix("##contig=<ID=") else {
            continue;
        };
        let end = rest.find([',', '>']).unwrap_or(rest.len());
        names.push(rest.as_bytes()[..end].to_vec());
    }
    names
}

/// The per-chunk job: rewrite the CHROM/ID columns, pass the rest through.
/// Records are re-deflated in the worker (the filter's pattern), so the sink
/// only appends ordered blocks.
struct RenameJob {
    plan: Arc<CompiledPlan>,
    level: u8,
    compressed: bool,
    /// An index is being written, so the block writer consumes the record
    /// key. Without an index the key is never read, and computing it would
    /// reject perfectly legal input (no `##contig` lines).
    need_key: bool,
}

/// Worker scratch: the reusable block writer.
struct RenameScratch {
    writer: reed_io::write::block_writer::BlockWriter,
}

/// Worker → sink: the rewritten block plus the counters. With
/// `--rm-dup-ids` the raw text also rides along: the unique-ID suffix must
/// follow the global record order, which only the ordered sink knows.
#[derive(Default)]
struct RenamePartial {
    output: reed_io::write::block_writer::ChunkOutput,
    raw: Vec<u8>,
    n_records: u64,
    n_ids: u64,
}

impl ChunkJob for RenameJob {
    type Partial = RenamePartial;
    type Scratch = Option<RenameScratch>;

    fn run_chunk(
        &self,
        _chunk: &ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| RenameScratch {
            writer: reed_io::write::block_writer::BlockWriter::new(self.level, self.compressed)
                .expect("invariant: the level comes from the CLI clamp"),
        });
        let mut partial = RenamePartial {
            output: reed_io::write::block_writer::ChunkOutput::default(),
            raw: Vec::new(),
            n_records: 0,
            n_ids: 0,
        };
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            self.push_line(line, scratch, &mut partial)?;
        }
        partial.output = scratch.writer.finish_chunk().map_err(EngineError::from)?;
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: reed_core::ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| RenameScratch {
            writer: reed_io::write::block_writer::BlockWriter::new(self.level, self.compressed)
                .expect("invariant: the level comes from the CLI clamp"),
        });
        let mut partial = RenamePartial {
            output: reed_io::write::block_writer::ChunkOutput::default(),
            raw: Vec::new(),
            n_records: 0,
            n_ids: 0,
        };
        if !line.is_empty() && line[0] != b'#' {
            self.push_line(line, scratch, &mut partial)?;
        }
        partial.output = scratch.writer.finish_chunk().map_err(EngineError::from)?;
        Ok(partial)
    }
}

impl RenameJob {
    /// Rewrites one record and appends it to the worker's block.
    fn push_line(
        &self,
        line: &[u8],
        scratch: &mut RenameScratch,
        partial: &mut RenamePartial,
    ) -> Result<(), EngineError> {
        // A chunk body splits its records on `\n`, but a seam record arrives
        // with its terminator (pipeline.rs: "a complete line ending in \n").
        // Normalizing here keeps the rewrite and the dedup text on exactly
        // one terminator per record.
        let line = line.strip_suffix(b"\n").unwrap_or(line);
        let mut rewritten = self
            .rewrite_line(line, partial)
            .unwrap_or_else(|| line.to_vec());
        // The block writer keys each record for the index; parse the key
        // from the rewritten line (CHROM may have been renamed, so the key
        // must come after the edit). Without an index the key is unused.
        let key = if self.need_key {
            self.index_key(&rewritten)?
        } else {
            reed_io::write::block_writer::IndexKey {
                contig: reed_core::ContigId(0),
                start1: 0,
                ref_len: 0,
            }
        };
        rewritten.push(b'\n');
        scratch
            .writer
            .push_record(&rewritten, key)
            .map_err(EngineError::from)?;
        if self.plan.rm_dup_ids {
            partial.raw.extend_from_slice(&rewritten);
        }
        partial.n_records += 1;
        Ok(())
    }

    /// The block index key of a rewritten record (header contig order).
    fn index_key(&self, line: &[u8]) -> Result<reed_io::write::block_writer::IndexKey, EngineError> {
        let lr = LineRef::parse(line)
            .ok_or_else(|| EngineError::job("record line has fewer than 9 columns"))?;
        let site = reed_io::vcf::site::SiteView::parse(&lr).map_err(EngineError::from)?;
        let contig = self
            .plan
            .contig_id(lr.chrom())
            .ok_or_else(|| EngineError::job("record contig is not in the header"))?;
        Ok(reed_io::write::block_writer::IndexKey {
            contig,
            start1: site.pos1,
            ref_len: u32::try_from(lr.reference().len()).unwrap_or(u32::MAX),
        })
    }

    /// Rewrites one record line's CHROM and ID columns; `None` when the line
    /// is not a parsable VCF record (it passes through untouched — a byte
    /// editor never silently drops a record).
    fn rewrite_line(&self, line: &[u8], partial: &mut RenamePartial) -> Option<Vec<u8>> {
        let lr = LineRef::parse(line)?;
        let site = reed_io::vcf::site::SiteView::parse(&lr).ok()?;
        let mut chrom = lr.chrom().to_vec();
        if let Some(new) = self.plan.chrs.get(&chrom) {
            chrom = new.clone();
        } else if let Some(style) = self.plan.chr_style {
            let name = String::from_utf8_lossy(&chrom).into_owned();
            let styled = style.apply(&name);
            if styled.as_bytes() != chrom.as_slice() {
                chrom = styled.into_bytes();
            }
        }
        let id = self.new_id(&lr, &chrom, site.pos1, partial);
        // CHROM, POS, ID, REF, ALT are the first five columns; QUAL onwards
        // passes through byte-for-byte (including the fifth tab).
        let mut out = Vec::with_capacity(line.len() + 16);
        out.extend_from_slice(&chrom);
        out.push(b'\t');
        out.extend_from_slice(lr.pos_text());
        out.push(b'\t');
        out.extend_from_slice(&id);
        out.push(b'\t');
        out.extend_from_slice(lr.reference());
        out.push(b'\t');
        out.extend_from_slice(lr.alts());
        let mut col = 0_u32;
        let mut from = line.len();
        for (i, &b) in line.iter().enumerate() {
            if b == b'\t' {
                col += 1;
                if col == 5 {
                    // Start at the tab itself: it is the ALT|QUAL separator
                    // (the ALT text above carries no terminator).
                    from = i;
                    break;
                }
            }
        }
        out.extend_from_slice(&line[from..]);
        Some(out)
    }

    /// The record's new ID (no ID edit requested returns the original bytes).
    fn new_id(
        &self,
        lr: &LineRef<'_>,
        chrom: &[u8],
        pos1: u64,
        partial: &mut RenamePartial,
    ) -> Vec<u8> {
        let plan = &self.plan;
        let original = lr.id();
        if plan.clear_ids {
            partial.n_ids += 1;
            return b".".to_vec();
        }
        if let Some(template) = &plan.set_id {
            partial.n_ids += 1;
            return template.render(chrom, lr, pos1);
        }
        if let Some(template) = &plan.set_missing_var_ids
            && original == b"."
        {
            partial.n_ids += 1;
            return template.render(chrom, lr, pos1);
        }
        let mut id = original.to_vec();
        if let Some((new, slot)) = plan.update_name.get(&id) {
            id = new.clone();
            if let Some(hit) = plan.update_hits.get(*slot) {
                hit.store(true, Ordering::Relaxed);
            }
            partial.n_ids += 1;
        }
        if let Some(sep) = &plan.idspace_to {
            let mut replaced = Vec::with_capacity(id.len());
            for &b in &id {
                if b == b':' {
                    replaced.extend_from_slice(sep);
                } else {
                    replaced.push(b);
                }
            }
            if replaced != id {
                id = replaced;
                partial.n_ids += 1;
            }
        }
        id
    }
}

/// Renames one record's ID to a unique `id_1`, `id_2`, … form when it has
/// been seen before (B14-6: suffix, never drop). Returns `(bytes, dup)`.
/// `next_suffix` remembers where the last probe for a base ID stopped: a run
/// of N identical IDs is otherwise quadratic (each duplicate re-probed
/// `_2`, `_3`, … from the start).
fn dedup_line(
    line: &[u8],
    seen: &mut std::collections::HashSet<Vec<u8>>,
    next_suffix: &mut std::collections::HashMap<Vec<u8>, u32>,
) -> (Vec<u8>, bool) {
    let Some(lr) = LineRef::parse(line) else {
        return (line.to_vec(), false);
    };
    let id = lr.id().to_vec();
    if id == b"." || seen.insert(id.clone()) {
        return (line.to_vec(), false);
    }
    let mut n = next_suffix.get(&id).copied().unwrap_or(2);
    loop {
        let mut candidate = id.clone();
        candidate.extend_from_slice(format!("_{n}").as_bytes());
        if seen.insert(candidate.clone()) {
            next_suffix.insert(id.clone(), n.saturating_add(1));
            // Replace the ID column (the 3rd field): copy CHROM<TAB>POS<TAB>,
            // then the new ID, then the rest from the next tab on.
            let mut out = Vec::with_capacity(line.len() + 4);
            let mut col = 0_u32;
            for (i, &b) in line.iter().enumerate() {
                if b != b'\t' {
                    continue;
                }
                col += 1;
                if col == 2 {
                    // At the tab after POS: emit CHROM<TAB>POS<TAB>new ID.
                    out.extend_from_slice(&line[..i]);
                    out.push(b'\t');
                    out.extend_from_slice(&candidate);
                } else if col == 3 {
                    // At the tab after ID: the rest of the line follows.
                    out.extend_from_slice(&line[i..]);
                    return (out, true);
                }
            }
            return (line.to_vec(), false);
        }
        n += 1;
    }
}

/// The ordered sink: append the rewritten blocks and publish the report.
struct RenameSink {
    out: reed_io::write::ordered_output::OrderedOutput,
    report: RenameReport,
    slot: Arc<std::sync::Mutex<Option<RenameReport>>>,
    /// `--rm-dup-ids`: the IDs already emitted, in record order.
    seen_ids: std::collections::HashSet<Vec<u8>>,
    /// `--rm-dup-ids`: the next suffix to try per base ID (keeps a run of
    /// identical IDs linear instead of re-probing from `_2`).
    next_suffix: std::collections::HashMap<Vec<u8>, u32>,
    /// The contig table for the dedup path's block keys.
    plan: Arc<CompiledPlan>,
    /// `--rm-dup-ids`: the sink-side deflater for the deduplicated text.
    dedup_writer: Option<reed_io::write::block_writer::BlockWriter>,
}

impl OrderedSink<RenamePartial> for RenameSink {
    fn accept(
        &mut self,
        _id: reed_core::ChunkId,
        partial: &mut RenamePartial,
    ) -> Result<(), EngineError> {
        self.report.n_records += partial.n_records;
        self.report.n_ids += partial.n_ids;
        if let Some(writer) = &mut self.dedup_writer {
            // `--rm-dup-ids`: the suffix follows the global record order, so
            // the sink (the only ordered stage) rewrites the ID column here
            // and re-deflates the chunk.
            for line in partial.raw.split_inclusive(|&b| b == b'\n') {
                if line.is_empty() || line[0] == b'#' {
                    continue;
                }
                // The parser wants the record without its terminator.
                let bare = line.strip_suffix(b"\n").unwrap_or(line);
                let (bytes, dup) = dedup_line(bare, &mut self.seen_ids, &mut self.next_suffix);
                if dup {
                    self.report.n_dup_ids += 1;
                }
                let key = {
                    let lr = LineRef::parse(&bytes)
                        .ok_or_else(|| EngineError::job("dedup: unparsable record"))?;
                    let site = reed_io::vcf::site::SiteView::parse(&lr)
                        .map_err(EngineError::from)?;
                    let contig = self
                        .plan
                        .contig_id(lr.chrom())
                        .ok_or_else(|| EngineError::job("dedup: contig not in the header"))?;
                    reed_io::write::block_writer::IndexKey {
                        contig,
                        start1: site.pos1,
                        ref_len: u32::try_from(lr.reference().len()).unwrap_or(u32::MAX),
                    }
                };
                let mut with_newline = bytes;
                if !with_newline.ends_with(b"\n") {
                    with_newline.push(b'\n');
                }
                writer
                    .push_record(&with_newline, key)
                    .map_err(EngineError::from)?;
            }
            let output = writer.finish_chunk().map_err(EngineError::from)?;
            return self
                .out
                .accept_chunk(&output)
                .map_err(EngineError::from);
        }
        self.out
            .accept_chunk(&partial.output)
            .map_err(EngineError::from)
    }

    fn finish(self) -> Result<u64, EngineError> {
        let n = self.report.n_records;
        self.out.finish().map_err(EngineError::from)?;
        *self
            .slot
            .lock()
            .expect("invariant: rename slot not poisoned") = Some(self.report);
        Ok(n)
    }

    fn deflate_time(&self) -> std::time::Duration {
        self.out.deflate_time()
    }
}

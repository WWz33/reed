//! Text → BCF encoding through noodles (doc_v3 02-bcf.md B2).
//!
//! One `BcfEncoder` per worker: [`BcfEncoder::write_header`] primes noodles'
//! internal `StringMaps` once, and every surviving record is parsed into an
//! owned `RecordBuf` (whose missing FORMAT strings are patched to htslib's
//! literal `.` first) and encoded. The bytes then flow through the existing
//! `BlockWriter` like any other output.

use noodles::vcf::variant::io::Write as _;
use noodles::vcf::variant::record_buf::samples::{Keys, Samples, sample::Value};

use crate::error::{IoError, Locus};
use noodles::vcf;

/// The per-worker encoder: noodles' BCF writer over a byte buffer.
pub struct BcfEncoder {
    writer: noodles::bcf::io::Writer<std::io::Cursor<Vec<u8>>>,
    scratch: Vec<u8>,
    /// FORMAT keys the header types as String (noodles rejects an
    /// all-missing string field; htslib writes the literal `.`).
    string_keys: Vec<String>,
}

impl BcfEncoder {
    /// Creates the encoder; [`BcfEncoder::write_header`] must run before the
    /// first record (it builds the dictionary→id `StringMaps`).
    ///
    /// # Errors
    ///
    /// none today; the signature keeps room for setup failures.
    pub fn new() -> Result<Self, IoError> {
        let buffer = std::io::Cursor::new(Vec::new());
        Ok(Self {
            writer: noodles::bcf::io::Writer::from(buffer),
            scratch: Vec::with_capacity(1 << 16),
            string_keys: Vec::new(),
        })
    }

    /// Writes the header once, building the string maps every record is
    /// keyed against; the serialized bytes are dropped (the caller writes
    /// the BCF file header itself).
    ///
    /// # Errors
    ///
    /// noodles' header serialization.
    pub fn write_header(&mut self, header: &vcf::Header) -> Result<(), IoError> {
        self.string_keys = header
            .formats()
            .iter()
            .filter(|(key, map)| {
                key.as_str() != "GT"
                    && map.ty() == noodles::vcf::header::record::value::map::format::Type::String
            })
            .map(|(key, _)| key.to_string())
            .collect();
        self.writer
            .write_header(header)
            .map_err(|e| IoError::Io(format!("BCF header encode: {e}")))?;
        let buffer = self.writer.get_mut();
        buffer.set_position(0);
        buffer.get_mut().clear();
        Ok(())
    }

    /// Encodes one VCF text record (newline included) and appends the BCF
    /// bytes to `out`.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] when the text does not parse, and
    /// [`IoError::Io`] from the encoder.
    pub fn encode_record(
        &mut self,
        header: &vcf::Header,
        line: &[u8],
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        self.scratch.clear();
        self.scratch.extend_from_slice(line);
        // The borrowed record avoids per-value allocations; it only fails on
        // a String FORMAT field whose values are all missing (noodles builds
        // no width for it), which the owned fallback patches.
        let mut reader = vcf::io::Reader::new(&self.scratch[..]);
        let mut record = vcf::Record::default();
        reader
            .read_record(&mut record)
            .map_err(|e| IoError::MalformedRecord {
                locus: Locus {
                    contig_name: "<bcf-encode>".into(),
                    pos1: 0,
                },
                reason: leak(&e),
            })?;
        match self.writer.write_variant_record(header, &record) {
            Ok(()) => {}
            Err(e) if is_missing_strings(&e) => {
                let mut reader = vcf::io::Reader::new(&self.scratch[..]);
                let mut owned = vcf::variant::RecordBuf::default();
                reader.read_record_buf(header, &mut owned).map_err(|e| {
                    IoError::MalformedRecord {
                        locus: Locus {
                            contig_name: "<bcf-encode>".into(),
                            pos1: 0,
                        },
                        reason: leak(&e),
                    }
                })?;
                fix_missing_strings(&self.string_keys, &mut owned);
                self.writer
                    .write_variant_record(header, &owned)
                    .map_err(|e| IoError::Io(format!("BCF record encode: {e}")))?;
            }
            Err(e) => return Err(IoError::Io(format!("BCF record encode: {e}"))),
        }
        let buffer = self.writer.get_mut();
        out.extend_from_slice(buffer.get_ref().as_slice());
        buffer.set_position(0);
        buffer.get_mut().clear();
        Ok(())
    }

    /// A fresh encoder with no string maps; the worker primes it with
    /// [`BcfEncoder::write_header`].
    #[must_use]
    pub fn clone_for_worker(&self) -> Self {
        Self::new().expect("invariant: the encoder allocates")
    }
}

/// Whether the encoder refused a record because a String FORMAT field had
/// no values at all (noodles' wording).
fn is_missing_strings(e: &std::io::Error) -> bool {
    e.kind() == std::io::ErrorKind::InvalidInput && e.to_string().contains("missing String values")
}

/// htslib stores a missing FORMAT string as the literal `.`; noodles keeps
/// `None` (and rejects a field whose values are all missing). Patch the
/// values so both encoders agree.
fn fix_missing_strings(string_keys: &[String], record: &mut vcf::variant::RecordBuf) {
    if string_keys.is_empty() {
        return;
    }
    let samples = record.samples_mut();
    let (keys, mut values): (Keys, Vec<Vec<Option<Value>>>) = std::mem::take(samples).into();
    let key_names: Vec<&str> = keys.as_ref().iter().map(String::as_str).collect();
    for row in &mut values {
        for (i, value) in row.iter_mut().enumerate() {
            let is_string = key_names
                .get(i)
                .is_some_and(|k| string_keys.iter().any(|s| s == k));
            if value.is_none() && is_string {
                *value = Some(Value::String(".".into()));
            }
        }
    }
    *record.samples_mut() = Samples::new(keys, values);
}

/// An owned error string squeezed into a `&'static str` slot. The message is
/// leaked: encode errors are fatal per-record diagnostics, not a hot path.
fn leak(e: &std::io::Error) -> &'static str {
    Box::leak(e.to_string().into_boxed_str())
}

impl Default for BcfEncoder {
    fn default() -> Self {
        Self::new().expect("invariant: the encoder allocates")
    }
}

/// Renders the BCF file header text (NUL-terminated) with an explicit
/// `IDX=` on every FILTER/INFO/FORMAT line (doc_v3 02-bcf.md §3.3).
///
/// BCF records store dictionary *indices*, so the writer's encoder and
/// htslib's reader must agree on them. noodles builds the ids in
/// `StringMaps::try_from`; writing those ids out pins htslib to the same
/// numbering no matter how the text lines are ordered.
///
/// # Errors
///
/// Character rendering and noodles' header serialization.
pub fn bcf_header_text(header: &vcf::Header) -> Result<Vec<u8>, IoError> {
    let maps = vcf::header::StringMaps::try_from(header)
        .map_err(|e| IoError::Io(format!("BCF header maps: {e}")))?;
    let mut raw = Vec::new();
    vcf::io::Writer::new(&mut raw)
        .write_header(header)
        .map_err(|e| IoError::Io(format!("BCF header encode: {e}")))?;
    let mut out = Vec::with_capacity(raw.len() + 512);
    for line in raw.split_inclusive(|&b| b == b'\n') {
        const KINDS: [&[u8]; 3] = [b"##FILTER=<", b"##INFO=<", b"##FORMAT=<"];
        let kind = KINDS.iter().find(|k| line.starts_with(k)).copied();
        let Some(kind) = kind else {
            out.extend_from_slice(line);
            continue;
        };
        // The payload starts with `ID=value,` (or ends right after the value).
        let rest = &line[kind.len()..];
        let value = rest
            .strip_prefix(b"ID=")
            .ok_or_else(|| IoError::Io("BCF header line does not start with ID=".into()))?;
        let id_len = value
            .iter()
            .position(|&b| b == b',' || b == b'>')
            .unwrap_or(value.len());
        let id = std::str::from_utf8(&value[..id_len])
            .map_err(|_| IoError::Io("BCF header ID is not UTF-8".into()))?;
        let idx = maps.strings().get_index_of(id).ok_or_else(|| {
            IoError::Io(format!(
                "BCF header tag {id:?} is not in the string dictionary"
            ))
        })?;
        // Insert `,IDX=n` before the closing '>'.
        let close = line.iter().rposition(|&b| b == b'>').unwrap_or(line.len());
        out.extend_from_slice(&line[..close]);
        let _ = std::io::Write::write_fmt(&mut out, format_args!(",IDX={idx}"));
        out.extend_from_slice(&line[close..]);
    }
    out.push(b'\0');
    Ok(out)
}

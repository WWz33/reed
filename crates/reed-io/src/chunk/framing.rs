//! How a chunk's inflated text is delimited into records (doc_v3 02-bcf.md
//! §2.2). VCF text splits at newlines; BCF records are length-prefixed and
//! opaque to the engine, so the job owns framing and seam stitching.

/// The record delimiter a chunk source declares.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum RecordFraming {
    /// VCF text: the engine splits the chunk at newlines, hands the body to
    /// the job and stitches cross-chunk records from tail/head fragments.
    #[default]
    Newline,
    /// Length-prefixed records (BCF): the engine passes the whole inflated
    /// text; the job walks the length prefixes and stitches seams itself.
    Raw,
}

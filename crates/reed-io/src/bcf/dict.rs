//! The BCF dictionaries (hts-specs §6.2.1): FILTER/INFO/FORMAT share one
//! string dictionary in header order, contigs have their own in `##contig`
//! order. An explicit `IDX=` pins an entry's index; otherwise the next free
//! index is used — the same rule htslib applies when it reads a BCF header.

use crate::error::IoError;

/// Index → name lookups for one BCF header.
#[derive(Clone, Debug, Default)]
pub struct BcfDict {
    /// String dictionary: index → tag name (holes stay `None`).
    strings: Vec<Option<Vec<u8>>>,
    /// Contig names in `##contig` order.
    contigs: Vec<Vec<u8>>,
}

impl BcfDict {
    /// Builds the dictionaries from the header text embedded in the BCF.
    ///
    /// # Errors
    ///
    /// [`IoError::Header`] when a header line is malformed (no `ID=`).
    pub fn from_header(text: &[u8]) -> Result<Self, IoError> {
        let mut dict = Self::default();
        let mut next = 0_usize;
        for line in text.split(|&b| b == b'\n') {
            if let Some(rest) = line.strip_prefix(b"##contig=<") {
                let name = field(rest, b"ID");
                if let Some(name) = name {
                    dict.contigs.push(name.to_vec());
                }
                continue;
            }
            let kind = if line.starts_with(b"##FILTER=<") {
                Some(b"##FILTER=<".len())
            } else if line.starts_with(b"##INFO=<") {
                Some(b"##INFO=<".len())
            } else if line.starts_with(b"##FORMAT=<") {
                Some(b"##FORMAT=<".len())
            } else {
                None
            };
            let Some(prefix) = kind else { continue };
            let rest = &line[prefix..];
            let Some(name) = field(rest, b"ID") else {
                return Err(IoError::Header(format!(
                    "BCF header line has no ID=: {}",
                    String::from_utf8_lossy(&line[..line.len().min(60)])
                )));
            };
            let idx = field(rest, b"IDX")
                .and_then(|v| std::str::from_utf8(v).ok())
                .and_then(|v| v.parse::<usize>().ok())
                .unwrap_or_else(|| {
                    let i = next;
                    next += 1;
                    i
                });
            if dict.strings.len() <= idx {
                dict.strings.resize(idx + 1, None);
            }
            dict.strings[idx] = Some(name.to_vec());
            next = next.max(idx + 1);
        }
        Ok(dict)
    }

    /// The tag name at `index` (`None` outside the dictionary).
    #[must_use]
    pub fn tag(&self, index: u32) -> Option<&[u8]> {
        self.strings.get(index as usize).and_then(Option::as_deref)
    }

    /// The contig name at `index` (`None` outside the dictionary).
    #[must_use]
    pub fn contig(&self, index: i32) -> Option<&[u8]> {
        usize::try_from(index)
            .ok()
            .and_then(|i| self.contigs.get(i))
            .map(Vec::as_slice)
    }

    /// Number of contigs.
    #[must_use]
    pub fn n_contigs(&self) -> usize {
        self.contigs.len()
    }
}

/// The `KEY=value` inside a `##TAG=<...>` payload (value up to `,` or `>`).
fn field<'a>(rest: &'a [u8], key: &[u8]) -> Option<&'a [u8]> {
    let mut i = 0;
    while i + key.len() < rest.len() {
        if &rest[i..i + key.len()] == key && rest[i + key.len()] == b'=' {
            let start = i + key.len() + 1;
            let end = rest[start..]
                .iter()
                .position(|&b| b == b',' || b == b'>')
                .map_or(rest.len(), |p| start + p);
            return Some(&rest[start..end]);
        }
        // Skip to the next top-level comma (values may contain quoted commas).
        let mut quoted = false;
        while i < rest.len() {
            match rest[i] {
                b'"' => quoted = !quoted,
                b',' if !quoted => {
                    i += 1;
                    break;
                }
                b'>' if !quoted => return None,
                _ => {}
            }
            i += 1;
        }
        if i >= rest.len() {
            return None;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_the_dictionary_in_header_order() {
        let text = b"##fileformat=VCFv4.2\n\
                     ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
                     ##contig=<ID=chr1,length=1000>\n\
                     ##FILTER=<ID=LowQual,Description=\"Low quality\">\n\
                     ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"a\">\n\
                     ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
                     ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
                     ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q, with comma\">\n\
                     #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        assert_eq!(dict.tag(0).expect("0"), b"PASS");
        assert_eq!(dict.tag(1).expect("1"), b"LowQual");
        assert_eq!(dict.tag(2).expect("2"), b"AC");
        assert_eq!(dict.tag(3).expect("3"), b"DP");
        assert_eq!(dict.tag(4).expect("4"), b"GT");
        assert_eq!(dict.tag(5).expect("5"), b"GQ");
        assert_eq!(dict.contig(0).expect("c0"), b"chr1");
        assert!(dict.tag(6).is_none());
    }

    #[test]
    fn honors_explicit_idx() {
        let text = b"##FILTER=<ID=PASS,Description=\"p\">\n\
                     ##INFO=<ID=X,Number=1,Type=Integer,IDX=7,Description=\"x\">\n\
                     ##INFO=<ID=Y,Number=1,Type=Integer,Description=\"y\">\n";
        let dict = BcfDict::from_header(text).expect("invariant: parses");
        assert_eq!(dict.tag(7).expect("7"), b"X");
        assert_eq!(dict.tag(8).expect("8"), b"Y");
        assert!(dict.tag(1).is_none());
    }
}

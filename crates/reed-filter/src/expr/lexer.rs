//! Lexer for the bcftools expression subset (02-cli.md §4).

use std::fmt;

/// One token of an expression.
#[derive(Clone, Debug, PartialEq)]
pub enum Token {
    /// Bare identifier (`QUAL`, `DP`, `INFO/DP`, `FMT/DP`).
    Ident(String),
    /// Numeric literal.
    Number(f64),
    /// String literal (quotes stripped).
    Str(String),
    /// Comparison / match operator.
    Cmp(&'static str),
    /// Logical operator (`&&`, `||`, `&`, `|`).
    Logic(&'static str),
    /// `!`
    Not,
    /// `(`
    LParen,
    /// `)`
    RParen,
    /// `[`
    LBracket,
    /// `]`
    RBracket,
    /// `.` (the missing value).
    Missing,
}

impl fmt::Display for Token {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ident(s) => write!(f, "{s}"),
            Self::Number(n) => write!(f, "{n}"),
            Self::Str(s) => write!(f, "\"{s}\""),
            Self::Cmp(op) | Self::Logic(op) => write!(f, "{op}"),
            Self::Not => write!(f, "!"),
            Self::LParen => write!(f, "("),
            Self::RParen => write!(f, ")"),
            Self::LBracket => write!(f, "["),
            Self::RBracket => write!(f, "]"),
            Self::Missing => write!(f, "."),
        }
    }
}

/// Tokenizes `text`.
///
/// # Errors
///
/// The 1-based column of the first character that cannot start a token.
pub fn tokenize(text: &str) -> Result<Vec<Token>, String> {
    let bytes = text.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        match c {
            b' ' | b'\t' | b'\n' | b'\r' => i += 1,
            b'(' => {
                out.push(Token::LParen);
                i += 1;
            }
            b')' => {
                out.push(Token::RParen);
                i += 1;
            }
            b'[' => {
                out.push(Token::LBracket);
                i += 1;
            }
            b']' => {
                out.push(Token::RBracket);
                i += 1;
            }
            b'"' | b'\'' => {
                let quote = c;
                let start = i + 1;
                let mut end = start;
                while end < bytes.len() && bytes[end] != quote {
                    end += 1;
                }
                if end == bytes.len() {
                    return Err(format!("unterminated string literal at column {}", i + 1));
                }
                out.push(Token::Str(text[start..end].to_string()));
                i = end + 1;
            }
            b'.' if !bytes.get(i + 1).is_some_and(u8::is_ascii_digit) => {
                out.push(Token::Missing);
                i += 1;
            }
            b'=' => {
                if bytes.get(i + 1) == Some(&b'=') {
                    out.push(Token::Cmp("=="));
                    i += 2;
                } else {
                    out.push(Token::Cmp("="));
                    i += 1;
                }
            }
            b'!' => {
                if bytes.get(i + 1) == Some(&b'=') {
                    out.push(Token::Cmp("!="));
                    i += 2;
                } else if bytes.get(i + 1) == Some(&b'~') {
                    out.push(Token::Cmp("!~"));
                    i += 2;
                } else {
                    out.push(Token::Not);
                    i += 1;
                }
            }
            b'<' => {
                if bytes.get(i + 1) == Some(&b'=') {
                    out.push(Token::Cmp("<="));
                    i += 2;
                } else {
                    out.push(Token::Cmp("<"));
                    i += 1;
                }
            }
            b'>' => {
                if bytes.get(i + 1) == Some(&b'=') {
                    out.push(Token::Cmp(">="));
                    i += 2;
                } else {
                    out.push(Token::Cmp(">"));
                    i += 1;
                }
            }
            b'~' => {
                out.push(Token::Cmp("~"));
                i += 1;
            }
            b'&' => {
                if bytes.get(i + 1) == Some(&b'&') {
                    out.push(Token::Logic("&&"));
                    i += 2;
                } else {
                    out.push(Token::Logic("&"));
                    i += 1;
                }
            }
            b'|' => {
                if bytes.get(i + 1) == Some(&b'|') {
                    out.push(Token::Logic("||"));
                    i += 2;
                } else {
                    out.push(Token::Logic("|"));
                    i += 1;
                }
            }
            b'0'..=b'9' => {
                let start = i;
                while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                    i += 1;
                }
                // Exponent form (1e-3) is accepted too.
                if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
                    i += 1;
                    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
                        i += 1;
                    }
                    while i < bytes.len() && bytes[i].is_ascii_digit() {
                        i += 1;
                    }
                }
                let text_num = &text[start..i];
                let value: f64 = text_num
                    .parse()
                    .map_err(|_| format!("invalid number {text_num:?} at column {}", start + 1))?;
                out.push(Token::Number(value));
            }
            b'a'..=b'z' | b'A'..=b'Z' | b'_' => {
                let start = i;
                while i < bytes.len()
                    && (bytes[i].is_ascii_alphanumeric()
                        || matches!(bytes[i], b'_' | b'/' | b'.' | b'-' | b'*'))
                {
                    i += 1;
                }
                out.push(Token::Ident(text[start..i].to_string()));
            }
            _ => {
                return Err(format!(
                    "unexpected character {:?} at column {}",
                    c as char,
                    i + 1
                ))
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenizes_the_gold_expressions() {
        assert_eq!(
            tokenize("QUAL>=30").expect("invariant"),
            [Token::Ident("QUAL".into()), Token::Cmp(">="), Token::Number(30.0)]
        );
        assert_eq!(
            tokenize("QD<4 || FS>50").expect("invariant"),
            [
                Token::Ident("QD".into()),
                Token::Cmp("<"),
                Token::Number(4.0),
                Token::Logic("||"),
                Token::Ident("FS".into()),
                Token::Cmp(">"),
                Token::Number(50.0),
            ]
        );
        assert_eq!(
            tokenize("FILTER=\"PASS\" && INFO/DP[0] > 10").expect("invariant"),
            [
                Token::Ident("FILTER".into()),
                Token::Cmp("="),
                Token::Str("PASS".into()),
                Token::Logic("&&"),
                Token::Ident("INFO/DP".into()),
                Token::LBracket,
                Token::Number(0.0),
                Token::RBracket,
                Token::Cmp(">"),
                Token::Number(10.0),
            ]
        );
        assert_eq!(tokenize("AF != .").expect("invariant"), [
            Token::Ident("AF".into()),
            Token::Cmp("!="),
            Token::Missing,
        ]);
    }
}

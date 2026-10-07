//! AST and parser for the bcftools expression subset (02-cli.md §4).
//!
//! Precedence, lowest first: `|| |` → `&& &` → `!` → comparisons
//! (`== = != < <= > >= ~ !~`) → primaries (identifiers, literals, `.`,
//! parentheses, subscripts).

use super::lexer::{tokenize, Token};

/// A built-in site field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    /// `QUAL`
    Qual,
    /// `FILTER`
    Filter,
    /// `ID`
    Id,
    /// `CHROM`
    Chrom,
    /// `POS`
    Pos,
    /// `TYPE`
    Type,
    /// `N_ALT`
    NAlt,
    /// `N_SAMPLES`
    NSamples,
    /// `AC` (INFO when present, else genotypes)
    Ac,
    /// `AN`
    An,
    /// `AF`
    Af,
    /// `MAC`
    Mac,
    /// `MAF`
    Maf,
    /// `N_MISSING` (genotypes)
    NMissing,
    /// `F_MISSING` (genotypes)
    FMissing,
}

impl Field {
    /// The canonical name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Qual => "QUAL",
            Self::Filter => "FILTER",
            Self::Id => "ID",
            Self::Chrom => "CHROM",
            Self::Pos => "POS",
            Self::Type => "TYPE",
            Self::NAlt => "N_ALT",
            Self::NSamples => "N_SAMPLES",
            Self::Ac => "AC",
            Self::An => "AN",
            Self::Af => "AF",
            Self::Mac => "MAC",
            Self::Maf => "MAF",
            Self::NMissing => "N_MISSING",
            Self::FMissing => "F_MISSING",
        }
    }

    /// Whether the field is derived from genotypes (and therefore evaluated
    /// after the mask stage, 02-cli.md §4).
    #[must_use]
    pub const fn is_genotype_derived(self) -> bool {
        matches!(
            self,
            Self::Ac | Self::An | Self::Af | Self::Mac | Self::Maf | Self::NMissing | Self::FMissing
        )
    }

    /// Resolves a built-in field name.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Self> {
        Some(match name {
            "QUAL" => Self::Qual,
            "FILTER" => Self::Filter,
            "ID" => Self::Id,
            "CHROM" => Self::Chrom,
            "POS" => Self::Pos,
            "TYPE" => Self::Type,
            "N_ALT" => Self::NAlt,
            "N_SAMPLES" => Self::NSamples,
            "AC" => Self::Ac,
            "AN" => Self::An,
            "AF" => Self::Af,
            "MAC" => Self::Mac,
            "MAF" => Self::Maf,
            "N_MISSING" => Self::NMissing,
            "F_MISSING" => Self::FMissing,
            _ => return None,
        })
    }
}

/// A subscript on a value.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Subscript {
    /// `[n]`: the n-th element.
    Index(usize),
    /// `[*]`: every element (only meaningful per sample in bcftools).
    Any,
}

/// One expression node.
#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    /// Numeric literal.
    Number(f64),
    /// String literal.
    Str(String),
    /// `.` (missing).
    Missing,
    /// A built-in field.
    Field(Field),
    /// `INFO/TAG`.
    Info(String),
    /// `FMT/TAG` (parsed; `check` reports it as unsupported).
    Fmt(String),
    /// A bare identifier: INFO first, then the built-in names.
    Bare(String),
    /// Comparison or regex match.
    Cmp(&'static str, Box<Expr>, Box<Expr>),
    /// Logical and/or.
    Logic(&'static str, Box<Expr>, Box<Expr>),
    /// `!expr`
    Not(Box<Expr>),
    /// `value[subscript]`
    Subscript(Box<Expr>, Subscript),
}

impl Expr {
    /// Collects the `FORMAT` tags the expression reads, in first-seen order.
    pub fn fmt_tags(&self, out: &mut Vec<String>) {
        match self {
            Self::Fmt(tag) => {
                if !out.iter().any(|t| t == tag) {
                    out.push(tag.clone());
                }
            }
            Self::Cmp(_, a, b) | Self::Logic(_, a, b) => {
                a.fmt_tags(out);
                b.fmt_tags(out);
            }
            Self::Not(e) | Self::Subscript(e, _) => e.fmt_tags(out),
            _ => {}
        }
    }

    /// Whether the expression reads genotype-derived quantities.
    #[must_use]
    pub fn references_genotypes(&self) -> bool {
        match self {
            Self::Field(f) => f.is_genotype_derived(),
            Self::Bare(name) => Field::from_name(name).is_some_and(Field::is_genotype_derived),
            Self::Cmp(_, a, b) | Self::Logic(_, a, b) => {
                a.references_genotypes() || b.references_genotypes()
            }
            // A FORMAT read is per sample: the mask stage must evaluate the
            // expression (`--set-gts .`, 02-cli.md §3.3).
            Self::Fmt(_) => true,
            Self::Not(e) | Self::Subscript(e, _) => e.references_genotypes(),
            _ => false,
        }
    }
}

/// Parses an expression.
///
/// # Errors
///
/// A human-readable message with the offending column or token.
pub fn parse(text: &str) -> Result<Expr, String> {
    let tokens = tokenize(text)?;
    if tokens.is_empty() {
        return Err("empty expression".to_string());
    }
    let mut parser = Parser { tokens, pos: 0 };
    let expr = parser.logic_or()?;
    if parser.pos != parser.tokens.len() {
        return Err(format!(
            "unexpected trailing token {}",
            parser.tokens[parser.pos]
        ));
    }
    Ok(expr)
}

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<Token> {
        let token = self.tokens.get(self.pos).cloned();
        if token.is_some() {
            self.pos += 1;
        }
        token
    }

    fn logic_or(&mut self) -> Result<Expr, String> {
        let mut left = self.logic_and()?;
        while let Some(Token::Logic(op @ ("||" | "|"))) = self.peek() {
            let op = *op;
            self.next();
            let right = self.logic_and()?;
            left = Expr::Logic(op, Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn logic_and(&mut self) -> Result<Expr, String> {
        let mut left = self.unary()?;
        while let Some(Token::Logic(op @ ("&&" | "&"))) = self.peek() {
            let op = *op;
            self.next();
            let right = self.unary()?;
            left = Expr::Logic(op, Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn unary(&mut self) -> Result<Expr, String> {
        if let Some(Token::Not) = self.peek() {
            self.next();
            let inner = self.unary()?;
            return Ok(Expr::Not(Box::new(inner)));
        }
        self.comparison()
    }

    fn comparison(&mut self) -> Result<Expr, String> {
        let left = self.primary()?;
        let Some(Token::Cmp(op)) = self.peek() else {
            return Ok(left);
        };
        let op = *op;
        self.next();
        let right = self.primary()?;
        Ok(Expr::Cmp(op, Box::new(left), Box::new(right)))
    }

    fn primary(&mut self) -> Result<Expr, String> {
        let Some(token) = self.next() else {
            return Err("expression ends where a value was expected".to_string());
        };
        let mut expr = match token {
            Token::Number(n) => Expr::Number(n),
            Token::Str(s) => Expr::Str(s),
            Token::Missing => Expr::Missing,
            Token::LParen => {
                let inner = self.logic_or()?;
                match self.next() {
                    Some(Token::RParen) => {}
                    other => return Err(format!("expected `)`, found {other:?}")),
                }
                return Ok(inner);
            }
            Token::Ident(name) => {
                if let Some(rest) = name.strip_prefix("INFO/") {
                    Expr::Info(rest.to_string())
                } else if let Some(rest) = name.strip_prefix("FMT/") {
                    Expr::Fmt(rest.to_string())
                } else if let Some(field) = Field::from_name(&name) {
                    Expr::Field(field)
                } else {
                    Expr::Bare(name)
                }
            }
            other => return Err(format!("unexpected token {other}")),
        };
        if let Some(Token::LBracket) = self.peek() {
            self.next();
            let index = match self.next() {
                Some(Token::Number(n)) if n >= 0.0 && n.fract() == 0.0 => {
                    Subscript::Index(n as usize)
                }
                Some(Token::Ident(s)) if s == "*" => Subscript::Any,
                other => return Err(format!("invalid subscript {other:?}")),
            };
            match self.next() {
                Some(Token::RBracket) => {}
                other => return Err(format!("expected `]`, found {other:?}")),
            }
            expr = Expr::Subscript(Box::new(expr), index);
        }
        Ok(expr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_precedence_and_fields() {
        let e = parse("QD<4 || FS>50 && MQ>=40").expect("invariant");
        // && binds tighter than ||.
        let Expr::Logic("||", _, right) = &e else { panic!("{e:?}") };
        assert!(matches!(**right, Expr::Logic("&&", _, _)));

        let e = parse("QUAL>=30").expect("invariant");
        assert!(matches!(e, Expr::Cmp(">=", _, _)));
        assert!(!e.references_genotypes());

        let e = parse("F_MISSING<0.1 || INFO/DP>10").expect("invariant");
        assert!(e.references_genotypes());

        let e = parse("AC > 1").expect("invariant");
        assert!(e.references_genotypes(), "AC prefers INFO but can fall back to genotypes");

        assert!(parse("QUAL >= 30 &&").is_err());
        assert!(parse("(QUAL >= 30").is_err());
    }
}

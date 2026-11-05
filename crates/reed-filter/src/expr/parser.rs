//! AST and parser for the bcftools expression subset (02-cli.md §4).
//!
//! Precedence, lowest first: `|| |` → `&& &` → `!` → comparisons
//! (`== = != < <= > >= ~ !~`) → primaries (identifiers, literals, `.`,
//! parentheses, subscripts).

// Casts in this file are bounded by construction: Subscript indices are bounded by the expression text length.
#![allow(clippy::cast_possible_truncation)]

use super::lexer::{Token, tokenize};

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
            Self::Ac
                | Self::An
                | Self::Af
                | Self::Mac
                | Self::Maf
                | Self::NMissing
                | Self::FMissing
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

/// A comparison operator, resolved at parse time — the hot path never
/// re-matches the operator text (A6-18).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CmpOp {
    /// `==`/`=`
    Eq,
    /// `!=`
    Ne,
    /// `<`
    Lt,
    /// `<=`
    Le,
    /// `>`
    Gt,
    /// `>=`
    Ge,
    /// `~` (literal pattern match)
    Match,
    /// `!~` (literal pattern non-match)
    NotMatch,
}

impl CmpOp {
    /// The operator as written (error messages).
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Eq => "==",
            Self::Ne => "!=",
            Self::Lt => "<",
            Self::Le => "<=",
            Self::Gt => ">",
            Self::Ge => ">=",
            Self::Match => "~",
            Self::NotMatch => "!~",
        }
    }
}

/// A logical connective, resolved at parse time (A6-18).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LogicOp {
    /// `&&`/`&`
    And,
    /// `||`/`|`
    Or,
}

impl LogicOp {
    /// The operator as written (error messages).
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::And => "&&",
            Self::Or => "||",
        }
    }
}

/// Maps the lexer's operator token text to its [`CmpOp`].
fn cmp_op(text: &str) -> CmpOp {
    match text {
        "==" | "=" => CmpOp::Eq,
        "!=" => CmpOp::Ne,
        "<" => CmpOp::Lt,
        "<=" => CmpOp::Le,
        ">" => CmpOp::Gt,
        ">=" => CmpOp::Ge,
        "~" => CmpOp::Match,
        "!~" => CmpOp::NotMatch,
        _ => unreachable!("the lexer only produces the operators above"),
    }
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
    Cmp(CmpOp, Box<Expr>, Box<Expr>),
    /// Logical and/or.
    Logic(LogicOp, Box<Expr>, Box<Expr>),
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

    /// Whether the expression reads quantities that need the decoded
    /// genotypes, i.e. must run after the mask stage (`-i/-e` routing).
    ///
    /// The `INFO/` prefix is a plain INFO read: bcftools evaluates `-i/-e`
    /// before it subsets the samples (`vcfview.c:322` before `:341`), so
    /// `INFO/AC` names the record's own INFO value and never a recomputed
    /// one. A *bare* built-in is different — `AC` falls back to the genotype
    /// counts when the record carries no such INFO, so it needs the scanner
    /// and routes after the mask stage. V25-04 collapsed both spellings onto
    /// one path; that made `-s S1 -i 'INFO/AC>1'` read the post-subset count
    /// (0) and drop every record.
    #[must_use]
    pub fn references_genotypes(&self) -> bool {
        match self {
            Self::Field(f) => f.is_genotype_derived(),
            // An explicit INFO read stays in the head stage; only the bare
            // names fall back to the genotype counts.
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

    /// Whether the expression reads a genotype-derived **site** field (as
    /// opposed to a `FORMAT` column, which `--set-gts .` evaluates per
    /// sample). A `--set-gts` mask decision must not depend on the
    /// post-mask counts, so the compiler rejects these.
    ///
    /// This keeps both spellings here: `--set-gts` evaluates on the mask
    /// stage, where a bare `AC` would read the post-mask count, and an
    /// `INFO/AC` would read a stale cohort's INFO. Neither is a legal basis
    /// for a mask decision, so the distinction that matters for `-i/-e`
    /// routing ([`Self::references_genotypes`]) does not apply to this gate.
    #[must_use]
    pub fn reads_gt_fields(&self) -> bool {
        match self {
            Self::Field(f) => f.is_genotype_derived(),
            Self::Bare(name) | Self::Info(name) => {
                Field::from_name(name).is_some_and(Field::is_genotype_derived)
            }
            Self::Cmp(_, a, b) | Self::Logic(_, a, b) => {
                a.reads_gt_fields() || b.reads_gt_fields()
            }
            Self::Not(e) | Self::Subscript(e, _) => e.reads_gt_fields(),
            _ => false,
        }
    }
}

/// The deepest nesting the parser accepts, counting parentheses and `!`
/// runs. The grammar is recursive descent, so an unbounded expression is an
/// unbounded recursion: 10,000 nested parentheses overflowed the stack and
/// aborted the process (doc_v23 V23-05). The same limit bounds the AST's
/// recursive `Drop` and the evaluator's recursion.
///
/// The value comes from measurement, not taste: in a debug build the
/// `parse` + `Drop` pair overflows a 2 MiB stack just past 430 levels and a
/// 1 MiB stack past ~200 levels, so 128 leaves ~2x headroom on the smallest
/// stack a caller is likely to hand us. Real bcftools expressions nest a
/// handful of levels; the limit is a denial-of-service bound, not a
/// expressiveness one.
const MAX_NESTING: usize = 128;

/// The longest single `||`/`&&` chain. Chains are parsed iteratively but
/// build a left-deep tree, so the AST's `Drop` and the evaluator recurse once
/// per term; a 100,000-term expression would overflow there even though the
/// parser itself does not recurse. Measured: a 2,000-term chain parses and
/// drops on a 512 KiB stack, so 2,000 is safely inside the same envelope.
const MAX_CHAIN: usize = 2_000;

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
    let mut parser = Parser {
        tokens,
        pos: 0,
        depth: 0,
    };
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
    /// Current nesting depth (parentheses and `!` operators).
    depth: usize,
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

    /// Runs `f` one nesting level deeper, refusing to descend past
    /// [`MAX_NESTING`].
    fn nested<T>(
        &mut self,
        f: impl FnOnce(&mut Self) -> Result<T, String>,
    ) -> Result<T, String> {
        self.depth += 1;
        if self.depth > MAX_NESTING {
            self.depth -= 1;
            return Err(format!(
                "expression nesting exceeds {MAX_NESTING} levels"
            ));
        }
        let result = f(self);
        self.depth -= 1;
        result
    }

    fn logic_or(&mut self) -> Result<Expr, String> {
        let mut left = self.logic_and()?;
        let mut terms = 1_usize;
        while let Some(Token::Logic("||" | "|")) = self.peek() {
            terms += 1;
            if terms > MAX_CHAIN {
                return Err(format!("expression has more than {MAX_CHAIN} `||` terms"));
            }
            self.next();
            let right = self.logic_and()?;
            left = Expr::Logic(LogicOp::Or, Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn logic_and(&mut self) -> Result<Expr, String> {
        let mut left = self.unary()?;
        let mut terms = 1_usize;
        while let Some(Token::Logic("&&" | "&")) = self.peek() {
            terms += 1;
            if terms > MAX_CHAIN {
                return Err(format!("expression has more than {MAX_CHAIN} `&&` terms"));
            }
            self.next();
            let right = self.unary()?;
            left = Expr::Logic(LogicOp::And, Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn unary(&mut self) -> Result<Expr, String> {
        if let Some(Token::Not) = self.peek() {
            self.next();
            let inner = self.nested(Self::unary)?;
            return Ok(Expr::Not(Box::new(inner)));
        }
        self.comparison()
    }

    fn comparison(&mut self) -> Result<Expr, String> {
        let left = self.primary()?;
        let Some(Token::Cmp(op)) = self.peek() else {
            return Ok(left);
        };
        let op = cmp_op(op);
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
                let inner = self.nested(Self::logic_or)?;
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
        let Expr::Logic(LogicOp::Or, _, right) = &e else {
            panic!("{e:?}")
        };
        assert!(matches!(**right, Expr::Logic(LogicOp::And, _, _)));

        let e = parse("QUAL>=30").expect("invariant");
        assert!(matches!(e, Expr::Cmp(CmpOp::Ge, _, _)));
        assert!(!e.references_genotypes());

        let e = parse("F_MISSING<0.1 || INFO/DP>10").expect("invariant");
        assert!(e.references_genotypes());

        let e = parse("AC > 1").expect("invariant");
        assert!(
            e.references_genotypes(),
            "AC prefers INFO but can fall back to genotypes"
        );

        assert!(parse("QUAL >= 30 &&").is_err());
        assert!(parse("(QUAL >= 30").is_err());
    }

    /// doc_v23 V23-05: the grammar is recursive descent over parentheses, so
    /// an unbounded expression is an unbounded recursion. 10,000 nested
    /// parentheses used to abort the process with a stack overflow; the limit
    /// must sit at 500 levels (499 and 500 parse, 501 is an error).
    #[test]
    fn nesting_is_bounded_at_the_limit() {
        let at_limit = |n: usize| {
            let mut text = "(".repeat(n);
            text.push_str("QUAL>=30");
            text.push_str(&")".repeat(n));
            text
        };
        let e = parse(&at_limit(127)).expect("127 levels are within the limit");
        assert!(matches!(e, Expr::Cmp(CmpOp::Ge, _, _)));
        parse(&at_limit(128)).expect("128 levels are within the limit");
        let err = parse(&at_limit(129)).expect_err("129 levels exceed the limit");
        assert!(err.contains("nesting exceeds 128 levels"), "{err}");
        // `!` runs count as nesting too (the same recursion shape).
        let negated = format!("{}QUAL>=30", "!".repeat(128));
        parse(&negated).expect("128 `!` are within the limit");
        let err = parse(&format!("{}QUAL>=30", "!".repeat(129))).expect_err("129 `!` rejected");
        assert!(err.contains("nesting exceeds"), "{err}");
        // The original report's 10,000 levels is a clean error, not a crash.
        let err = parse(&at_limit(10_000)).expect_err("10,000 levels are refused");
        assert!(err.contains("nesting exceeds"), "{err}");
    }

    /// A long flat `||` chain builds a left-deep tree whose `Drop` recurses
    /// once per term. The parser stays iterative but the limit keeps the AST
    /// and the evaluator within the same bound.
    #[test]
    fn logic_chains_are_bounded() {
        let chain = |n: usize| vec!["QUAL>=30"; n].join(" || ");
        parse(&chain(2_000)).expect("2,000 terms are within the limit");
        let err = parse(&chain(2_001)).expect_err("2,001 terms exceed the limit");
        assert!(err.contains("more than 2000"), "{err}");
    }
}

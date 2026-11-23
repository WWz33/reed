//! The convert chunk job: per-record target dispatch and the
//! per-target fold helpers (split out of `job/mod.rs`, doc_v30 P4).

use super::msa::{cell_012, is_snp, raw_alt_count, resolve_msa_value, resolve_pick};
use super::{
    ChunkJob, ConvertJob, ConvertPartial, ConvertScratch, CountPartial, DrawKey, EngineError,
    GtSlots, HmpPartial, LineRef, MsaPartial, MsaRun, MsaSite, PhyloPartial, PhyloPlan,
    RecordClass, StatsJob, classify_line, iupac,
};

impl ChunkJob for ConvertJob {
    type Partial = ConvertPartial;
    type Scratch = ConvertScratch;

    fn run_chunk(
        &self,
        chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = ConvertPartial::default();
        // The composed SFS body is rebuilt only when a record must be
        // dropped: structural blanks need no dropping (the statistics job
        // skips them itself), but malformed lines (doc_v25 V25-03) and
        // records whose GT the strict scanner would reject — out-of-domain
        // alleles, non-digit tokens, polyploids (doc_v25 V25-01) — must not
        // reach it. Convert's policy for those is "count malformed, exit 0"
        // (B13-9), so handing them to the strict scanner used to turn the
        // same file into a fatal `--sfs` run. A clean chunk still hands the
        // original slice over (zero copy).
        let want_sfs = self.sfs.is_some();
        let sfs_samples = self.sfs.as_ref().map(StatsJob::active_samples);
        let mut dropped: Vec<usize> = Vec::new();
        let mut line_no: u32 = 0;
        for (index, line) in body.split(|&b| b == b'\n').enumerate() {
            let key = DrawKey {
                chunk: chunk.id.get(),
                line: line_no,
            };
            // Blank lines are not records: they do not consume a draw slot
            // (the key must stay aligned with the emitted rows).
            let outcome = self.push_line(line, key, &mut partial, scratch)?;
            if outcome != LineOutcome::Blank {
                line_no += 1;
            }
            if want_sfs {
                let keep = match outcome {
                    LineOutcome::Record => sfs_accepts_gt(scratch, sfs_samples),
                    // `Blank` passes through: the statistics job skips it.
                    LineOutcome::Blank => true,
                    LineOutcome::Malformed => false,
                };
                if !keep {
                    dropped.push(index);
                }
            }
        }
        if let Some(stats) = &self.sfs {
            let ConvertScratch {
                sanitized,
                stats: stats_scratch,
                ..
            } = scratch;
            let sfs_body = if dropped.is_empty() {
                body
            } else {
                sanitized.clear();
                let mut next_drop = dropped.iter().copied().peekable();
                for (index, line) in body.split(|&b| b == b'\n').enumerate() {
                    if next_drop.peek() == Some(&index) {
                        next_drop.next();
                        continue;
                    }
                    if !sanitized.is_empty() {
                        sanitized.push(b'\n');
                    }
                    sanitized.extend_from_slice(line);
                }
                sanitized.as_slice()
            };
            partial.sfs = Some(stats.run_chunk(chunk, sfs_body, stats_scratch)?);
        }
        Ok(partial)
    }

    fn run_seam(
        &self,
        id: reed_core::ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = ConvertPartial::default();
        // The seam record belongs to the previous chunk's stream.
        let key = DrawKey {
            chunk: id.get().saturating_sub(1),
            line: u32::MAX,
        };
        let outcome = self.push_line(line, key, &mut partial, scratch)?;
        // Feed the SFS target only what convert accepts (the same policy as
        // the body path: structural blanks/malformed and GT-domain rejects
        // never reach the strict statistics scanner).
        if let Some(stats) = &self.sfs
            && outcome == LineOutcome::Record
            && sfs_accepts_gt(scratch, Some(stats.active_samples()))
        {
            let ConvertScratch {
                stats: stats_scratch,
                ..
            } = scratch;
            partial.sfs = Some(stats.run_seam(id, line, stats_scratch)?);
        }
        Ok(partial)
    }
}

/// Whether the genotype just parsed by [`ConvertJob::parse_gts`] passes the
/// strict statistics scanner's grammar (doc_v25 V25-01): every *participating*
/// sample's GT is diploid-or-haploid, digit-only within the scanner's
/// five-digit token cap, and inside the record's allele domain.
///
/// `active` is the statistics sample view: excluded samples and samples
/// outside the SFS populations never reach the scanner, so neither may veto a
/// record here (v29 V29-02 / D29-03).
fn sfs_accepts_gt(scratch: &ConvertScratch, active: Option<&reed_core::sample::KeepMask>) -> bool {
    scratch.gts.iter().enumerate().all(|(sample, gt)| {
        active.is_some_and(|mask| !mask[sample])
            || !(gt.bad || gt.strict_bad || gt.overflow || gt.n > 2)
    })
}

impl ConvertJob {
    pub(super) fn has_row_targets(&self) -> bool {
        self.phylo.is_some() || self.count || self.hmp || self.msa.is_some()
    }

    /// Folds one line into every row target. Returns whether the line was a
    /// record slot (non-blank): malformed lines still consume one, exactly as
    /// vcf2phylip counts them; blank lines do not exist for it.
    fn push_line(
        &self,
        line: &[u8],
        key: DrawKey,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) -> Result<LineOutcome, EngineError> {
        scratch.record_malformed = false;
        let lr = match classify_line(line, self.n_samples) {
            RecordClass::Record(lr) => lr,
            // `line_is_sane` used to call a blank seam line sane and hand it
            // to the statistics job; a blank line is now simply not a record
            // (doc_v25 V25-03).
            RecordClass::Blank => return Ok(LineOutcome::Blank),
            RecordClass::Malformed(why) => {
                // vcf2phylip counts every non-empty, non-`#` line as
                // processed, malformed ones included (`snp_num` increments
                // before its checks).
                partial.counters.total_records += 1;
                self.bump_malformed(partial, scratch, why.label());
                return Ok(LineOutcome::Malformed);
            }
        };
        partial.counters.total_records += 1;
        let pos1 = reed_io::vcf::fixed::parse_pos1(lr.pos_text(), lr.chrom())
            .expect("invariant: classify_line validated POS");
        let pos0 = pos1
            .checked_sub(1)
            .expect("invariant: classify_line rejects POS 0");

        // The single per-sample GT parse (D13-7).
        self.parse_gts(&lr, scratch);
        // Polyploidy is outside the supported domain (v29 D29-06): when it
        // would reach `--sfs` or `--012`, count the record once and warn.
        let polyploid = scratch
            .gts
            .iter()
            .zip(self.keep.iter())
            .any(|(gt, &kept)| kept && gt.n > 2);
        if polyploid && (self.sfs.is_some() || self.count) {
            partial.counters.polyploid += 1;
            if !scratch.warned.polyploid {
                scratch.warned.polyploid = true;
                tracing::warn!(
                    "reed convert: records with polyploid genotypes are outside the supported \
                     domain; --sfs skips them and --012 keeps vcftools' -1, counted once \
                     (v29 D29-06)"
                );
            }
        }

        // Presence count over the kept samples (vcf2phylip counts samples
        // whose GT text does not start with `.`; v13 §3.1 step 2).
        let present = scratch
            .gts
            .iter()
            .zip(self.keep.iter())
            .filter(|(gt, kept)| **kept && gt.present)
            .count() as u32;
        // The phylo gate keeps vcf2phylip's stricter rule: `./1` counts as
        // missing there, unlike the count-family total above (doc_v30
        // §4.6(a); R01 had unified the two and changed phylo's site count).
        let phylo_present = scratch
            .gts
            .iter()
            .zip(self.keep.iter())
            .filter(|(gt, kept)| **kept && gt.phylo_present)
            .count() as u32;

        // --- The phylogenetic stream (nucleotide + bin, vcf2phylip order) --
        // Its gates skip *this stream only*: a phylo early exit must not
        // starve the count/hapmap/MSA folds below (doc_v15 D15-4).
        'phylo: {
            let Some(plan) = &self.phylo else {
                break 'phylo;
            };
            if phylo_present < plan.gate {
                partial.counters.excess_missing += 1;
                break 'phylo;
            }
            self.build_alleles(&lr, true, scratch);
            if !is_snp(&lr, &scratch.allele_ranges) {
                partial.counters.mnps += 1;
                break 'phylo;
            }
            if plan.formats.any() {
                match self.phylo_column(scratch, key, plan.options.resolve_iupac) {
                    Ok(()) => {
                        let write_used = plan.options.write_used_sites;
                        let phy = Self::phylo_partial(partial, plan, self.n_kept);
                        for (i, byte) in scratch.column.iter().copied().enumerate() {
                            phy.cols[i].push(byte);
                        }
                        phy.n_sites += 1;
                        if write_used {
                            let used = &mut phy.used;
                            used.extend_from_slice(lr.chrom());
                            used.push(b'\t');
                            used.extend_from_slice(lr.pos_text());
                            used.push(b'\t');
                            used.extend_from_slice(present.to_string().as_bytes());
                            used.push(b'\n');
                        }
                        partial.counters.accepted += 1;
                    }
                    Err(Malformed) => {
                        // vcf2phylip's malformed skip also drops the record
                        // from the bin stream (the nucleotide check runs
                        // first and `continue`s).
                        self.bump_malformed(
                            partial,
                            scratch,
                            "a genotype with an unknown allele letter or an index beyond the record's alleles",
                        );
                        break 'phylo;
                    }
                }
            }
            // --- Binary NEXUS (nested under the SNP gate, v13 §3.4) -------
            if plan.bin && lr.alts().len() == 1 {
                let samples = lr.samples_text();
                let mut cells = Vec::with_capacity(self.n_kept);
                for (gt, &kept) in scratch.gts.iter().zip(self.keep.iter()) {
                    if !kept {
                        continue;
                    }
                    let raw = gt
                        .raw
                        .as_ref()
                        .and_then(|span| samples.get(span.start as usize..span.end as usize))
                        .unwrap_or(&b"."[..]);
                    cells.push(iupac::gen_bin_lookup(raw).unwrap_or(b'?'));
                }
                {
                    let phy = Self::phylo_partial(partial, plan, self.n_kept);
                    for (i, byte) in cells.into_iter().enumerate() {
                        phy.bin_cols[i].push(byte);
                    }
                    phy.n_bin_sites += 1;
                }
                partial.counters.biallelic_for_bin += 1;
            }
        }

        self.push_count_matrix(&lr, present, partial, scratch);
        self.push_hapmap(&lr, present, partial, scratch);
        self.push_msa(&lr, pos0, partial, scratch)?;
        Ok(LineOutcome::Record)
    }

    /// The phylogenetic partial, created with the shape both branches agree
    /// on (nucleotide columns and/or binary columns, per the plan).
    fn phylo_partial<'a>(
        partial: &'a mut ConvertPartial,
        plan: &PhyloPlan,
        n_kept: usize,
    ) -> &'a mut PhyloPartial {
        partial.phylo.get_or_insert_with(|| PhyloPartial {
            n_sites: 0,
            n_bin_sites: 0,
            cols: if plan.formats.any() {
                vec![Vec::new(); n_kept]
            } else {
                Vec::new()
            },
            bin_cols: if plan.bin {
                vec![Vec::new(); n_kept]
            } else {
                Vec::new()
            },
            used: Vec::new(),
        })
    }

    fn bump_malformed(
        &self,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
        why: &str,
    ) {
        if scratch.record_malformed {
            return;
        }
        scratch.record_malformed = true;
        partial.counters.malformed += 1;
        if !scratch.warned.malformed {
            scratch.warned.malformed = true;
            tracing::warn!(
                "reed convert: skipping a malformed record ({why}); further ones are counted in PREFIX.log"
            );
        }
    }

    /// Parses every sample's GT tokens into `scratch.gts` (the D13-7 view).
    /// An allele index above the record's ALT count marks the genotype bad
    /// (doc_v25 V25-01): every target then counts the record as malformed —
    /// the `--012` path used to code a haploid out-of-domain call as `-1`
    /// silently.
    fn parse_gts(&self, lr: &LineRef<'_>, scratch: &mut ConvertScratch) {
        // The record's ALT count (0 for `ALT=.`): legal indices are 0..=n_alt.
        let n_alt: u32 = if lr.alts() == b"." {
            0
        } else {
            u32::try_from(memchr::memchr_iter(b',', lr.alts()).count() + 1).unwrap_or(u32::MAX)
        };
        // One GT-position lookup per record (doc_v25 V25-02): the old
        // `memchr(b':')` took the first subfield, so `FORMAT=DP:GT` handed
        // the DP value to the genotype parser.
        let gt_index = reed_io::vcf::format_layout::gt_index(lr.format_text());
        scratch.gts.clear();
        scratch.gts.resize(self.n_samples, GtSlots::default());
        // doc_v28 V28-24: one borrowed-span pass. The GT range is stored
        // relative to the sample region (what the binary-NEXUS path slices);
        // it points at the GT subfield, not the whole sample column.
        let base = lr.samples_start();
        for (gt, span) in scratch.gts.iter_mut().zip(lr.sample_ranges()) {
            let column = &lr.line[span.clone()];
            let begin = span.start - base;
            let raw_span =
                gt_index.and_then(|i| reed_io::vcf::format_layout::subfield_span(column, i));
            // A column shorter than the GT position has no GT text; treat it
            // as missing, never as the first subfield.
            let raw = raw_span
                .as_ref()
                .map_or(&b"."[..], |range| &column[range.clone()]);
            gt.raw = raw_span.map(|range| (begin + range.start) as u32..(begin + range.end) as u32);
            let mut tokens = 0_usize;
            for token in raw.split(|&b| b == b'/' || b == b'|') {
                if tokens >= 4 {
                    gt.overflow = true;
                    break;
                }
                if token == b"." {
                    gt.slots[tokens] = None;
                } else if !token.is_empty() && token.iter().all(u8::is_ascii_digit) {
                    // The row targets' grammar: any digit run (bcftools
                    // accepts zero-padded alleles), but the strict scanner
                    // caps a token at five digits — record that separately
                    // instead of failing this record for every target.
                    if token.len() > 5 {
                        gt.strict_bad = true;
                    }
                    match std::str::from_utf8(token)
                        .ok()
                        .and_then(|t| t.parse::<u32>().ok())
                    {
                        Some(idx) if idx <= n_alt => gt.slots[tokens] = Some(idx),
                        _ => gt.bad = true,
                    }
                } else {
                    gt.bad = true;
                }
                tokens += 1;
            }
            gt.n = u8::try_from(tokens.min(4)).unwrap_or(4);
            gt.present = reed_core::genotype::has_called_allele(raw);
            // vcf2phylip's own presence rule for the `--min-samples-locus`
            // gate: any GT text starting with `.` is missing, so `./1`
            // does not count (doc_v30 §4.6(a); R01 had unified this with
            // stats' half-call rule and silently changed the phylo gate).
            gt.phylo_present = raw.first().is_some_and(|&b| b != b'.');
        }
    }

    /// Builds the processed allele texts: range 0 = REF, `i + 1` = `ALT[i]`.
    /// `uppercase` selects the phylogenetic alphabet (uppercase, `-` → `*`,
    /// `<NON_REF>` → REF) versus MSA's literal alleles (`*` → `-`).
    fn build_alleles(&self, lr: &LineRef<'_>, uppercase: bool, scratch: &mut ConvertScratch) {
        scratch.alleles.clear();
        scratch.allele_ranges.clear();
        let push = |scratch: &mut ConvertScratch, bytes: &[u8]| {
            let start = scratch.alleles.len() as u32;
            scratch.alleles.extend_from_slice(bytes);
            scratch.allele_ranges.push((start, bytes.len() as u32));
        };
        if uppercase {
            let mut reference: Vec<u8> = lr
                .reference()
                .iter()
                .map(|&b| {
                    if b == b'-' {
                        b'*'
                    } else {
                        b.to_ascii_uppercase()
                    }
                })
                .collect();
            push(scratch, &reference);
            let ref_processed = std::mem::take(&mut reference);
            for alt in lr.alts().split(|&b| b == b',') {
                let mut processed: Vec<u8> = Vec::with_capacity(alt.len());
                let mut i = 0;
                while i < alt.len() {
                    if alt[i] == b'-' {
                        processed.push(b'*');
                        i += 1;
                    } else if alt[i..].starts_with(b"<NON_REF>") {
                        processed.extend_from_slice(&ref_processed);
                        i += b"<NON_REF>".len();
                    } else {
                        processed.push(alt[i].to_ascii_uppercase());
                        i += 1;
                    }
                }
                push(scratch, &processed);
            }
        } else {
            push(scratch, lr.reference());
            for alt in lr.alts().split(|&b| b == b',') {
                if alt == b"*" {
                    push(scratch, b"-");
                } else {
                    push(scratch, alt);
                }
            }
        }
    }

    /// vcf2phylip's SNP gate on the *processed* allele table: REF is one
    /// character and every ALT allele is one character.
    fn phylo_column(
        &self,
        scratch: &mut ConvertScratch,
        key: DrawKey,
        resolve: bool,
    ) -> Result<(), Malformed> {
        let ConvertScratch {
            gts,
            allele_ranges,
            alleles,
            column,
            ..
        } = scratch;
        column.clear();
        for (sample, gt) in gts.iter().enumerate() {
            if !self.keep[sample] {
                continue;
            }
            if gt.bad || gt.overflow {
                return Err(Malformed);
            }
            let letter_at = |index: u32| -> Result<&[u8], Malformed> {
                let &(start, len) = allele_ranges.get(index as usize).ok_or(Malformed)?;
                Ok(&alleles[start as usize..(start + len) as usize])
            };
            if resolve {
                if gt.n == 0 {
                    return Err(Malformed);
                }
                let n = usize::from(gt.n);
                let token = resolve_pick(key, sample as u32, n, self.seed);
                let letter = match gt.slots[token] {
                    None => b'N',
                    Some(index) => {
                        let allele = letter_at(index)?;
                        if allele.len() != 1 {
                            return Err(Malformed);
                        }
                        if allele[0] == b'*' { b'-' } else { allele[0] }
                    }
                };
                column.push(letter);
                continue;
            }
            let mut set = [0_u8; 8];
            let n = usize::from(gt.n).min(8);
            for (t, slot) in set.iter_mut().enumerate().take(n) {
                *slot = match gt.slots[t] {
                    None => b'N',
                    Some(index) => {
                        let allele = letter_at(index)?;
                        if allele.len() != 1 {
                            return Err(Malformed);
                        }
                        allele[0]
                    }
                };
            }
            set[..n].sort_unstable();
            let mut deduped = 0;
            for j in 0..n {
                if j == 0 || set[j] != set[j - 1] {
                    set[deduped] = set[j];
                    deduped += 1;
                }
            }
            match iupac::ambig_lookup(&set[..deduped]) {
                Some(letter) => column.push(letter),
                // The reference implementation raises here (its AMBIG lookup
                // is outside the malformed guard); reed skips + counts
                // (B13-9).
                None => return Err(Malformed),
            }
        }
        Ok(())
    }

    fn push_count_matrix(
        &self,
        lr: &LineRef<'_>,
        present: u32,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) {
        if !self.count {
            return;
        }
        // vcftools' biallelic gate: `>2` alleles is skipped; `ALT='.'`
        // counts one allele and is kept.
        if raw_alt_count(lr) > 1 {
            if !scratch.warned.multiallelic_012 {
                scratch.warned.multiallelic_012 = true;
                tracing::warn!(
                    "reed convert --012: skipping multiallelic loci (vcftools behavior)"
                );
            }
            return;
        }
        if self.count_gate.is_some_and(|gate| present < gate) {
            return;
        }
        let mut cells = Vec::with_capacity(self.n_kept);
        for (gt, &kept) in scratch.gts.iter().zip(self.keep.iter()) {
            if !kept {
                continue;
            }
            let Some(cell) = cell_012(gt) else {
                self.bump_malformed(partial, scratch, "a genotype index beyond the ALT alleles");
                return;
            };
            cells.push(cell);
        }
        let out = partial.count.get_or_insert_with(|| CountPartial {
            n_sites: 0,
            cells: vec![Vec::new(); self.n_kept],
            pos: Vec::new(),
        });
        for (i, cell) in cells.into_iter().enumerate() {
            out.cells[i].push(cell);
        }
        out.n_sites += 1;
        out.pos.extend_from_slice(lr.chrom());
        out.pos.push(b'\t');
        out.pos.extend_from_slice(lr.pos_text());
        out.pos.push(b'\n');
    }

    fn push_hapmap(
        &self,
        lr: &LineRef<'_>,
        present: u32,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) {
        if !self.hmp {
            return;
        }
        if self.hmp_gate.is_some_and(|gate| present < gate) {
            return;
        }
        // The SNP gate (shared with the phylogenetic stream, v13 §3.9).
        self.build_alleles(lr, true, scratch);
        if !is_snp(lr, &scratch.allele_ranges) {
            return;
        }
        match self.hapmap_row(lr, scratch) {
            Ok(row) => {
                let out = partial.hmp.get_or_insert_with(HmpPartial::default);
                out.rows.extend_from_slice(&row);
                out.n_sites += 1;
            }
            Err(Malformed) => {
                self.bump_malformed(partial, scratch, "a hapmap row without an IUPAC letter");
            }
        }
    }

    /// Renders one hapmap site row. Cells use TASSEL's diploid alphabet
    /// (v13 §3.9): the IUPAC letter when TASSEL reads it (A/C/G/T/R/Y/S/W/
    /// K/M/N and `-`), else `N` (its reader has no lowercase letters and no
    /// three-base codes V/H/D/B).
    fn hapmap_row(
        &self,
        lr: &LineRef<'_>,
        scratch: &mut ConvertScratch,
    ) -> Result<Vec<u8>, Malformed> {
        self.phylo_column(scratch, DrawKey { chunk: 0, line: 0 }, false)?;
        for letter in &mut scratch.column {
            *letter = match *letter {
                b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                | b'-' => *letter,
                _ => b'N',
            };
        }
        let mut row = Vec::with_capacity(64 + self.n_kept);
        // rs#: the VCF ID, or `chrom:pos` when absent.
        if lr.id() == b"." {
            row.extend_from_slice(lr.chrom());
            row.push(b':');
            row.extend_from_slice(lr.pos_text());
        } else {
            row.extend_from_slice(lr.id());
        }
        row.push(b'\t');
        // alleles: the declared REF/ALT set (B13-5).
        let mut first = true;
        let mut seen: Vec<&[u8]> = Vec::new();
        for (i, &(start, len)) in scratch.allele_ranges.iter().enumerate() {
            let mut text = &scratch.alleles[start as usize..(start + len) as usize];
            // Skip an ALT `.` (no alternative allele).
            if i > 0 && text == b"." {
                continue;
            }
            // The GATK deletion symbol is TASSEL's gap.
            if text == b"*" {
                text = b"-";
            }
            if seen.contains(&text) {
                continue;
            }
            seen.push(text);
            if !first {
                row.push(b'/');
            }
            row.extend_from_slice(text);
            first = false;
        }
        row.push(b'\t');
        row.extend_from_slice(lr.chrom());
        row.push(b'\t');
        row.extend_from_slice(lr.pos_text());
        row.extend_from_slice(b"\t+\tNA\tNA\tNA\tNA\tNA\tNA\t");
        for (i, &letter) in scratch.column.iter().enumerate() {
            row.push(letter);
            if i + 1 < self.n_kept {
                row.push(b'\t');
            }
        }
        row.push(b'\n');
        Ok(row)
    }

    fn push_msa(
        &self,
        lr: &LineRef<'_>,
        pos0: u64,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) -> Result<(), EngineError> {
        let Some(plan) = &self.msa else {
            return Ok(());
        };
        if !plan.allows(lr.chrom(), pos0) {
            return Ok(());
        }
        // D29-01 rule 5: the record's whole REF span must exist in the
        // reference and match it (case-insensitive).
        let reference = lr.reference();
        let Some(sequence) = plan.reference.get(lr.chrom()) else {
            return Err(EngineError::job(format!(
                "MSA record {}:{} names a contig missing from the reference FASTA",
                String::from_utf8_lossy(lr.chrom()),
                pos0 + 1
            )));
        };
        let ref_len = u64::try_from(reference.len()).unwrap_or(u64::MAX);
        let end = pos0.checked_add(ref_len).ok_or_else(|| {
            EngineError::job(format!(
                "MSA record {}:{}: REF span overflows the coordinate space",
                String::from_utf8_lossy(lr.chrom()),
                pos0 + 1
            ))
        })?;
        let contig_len = sequence.len() as u64;
        if end > contig_len {
            return Err(EngineError::job(format!(
                "MSA record {}:{}: REF span of {ref_len} bp exceeds the reference contig length {contig_len}",
                String::from_utf8_lossy(lr.chrom()),
                pos0 + 1
            )));
        }
        if !sequence[pos0 as usize..end as usize].eq_ignore_ascii_case(reference) {
            return Err(EngineError::job(format!(
                "MSA record {}:{}: REF {} does not match the reference FASTA",
                String::from_utf8_lossy(lr.chrom()),
                pos0 + 1,
                String::from_utf8_lossy(reference)
            )));
        }
        let ref_len = u32::try_from(ref_len).map_err(|_| {
            EngineError::job(format!(
                "MSA record {}:{}: REF span is longer than u32",
                String::from_utf8_lossy(lr.chrom()),
                pos0 + 1
            ))
        })?;
        // Literal allele texts (vcf2msa's genotype_resolve domain).
        self.build_alleles(lr, false, scratch);
        match self.msa_values(scratch) {
            Ok((vals, width)) => {
                let out = partial.msa.get_or_insert_with(MsaPartial::default);
                let site = MsaSite {
                    pos0,
                    ref_len,
                    vals,
                    width,
                };
                match out.runs.last_mut() {
                    Some(run) if run.contig == lr.chrom() => run.sites.push(site),
                    _ => out.runs.push(MsaRun {
                        contig: lr.chrom().to_vec(),
                        sites: vec![site],
                    }),
                }
            }
            Err(Malformed) => {
                self.bump_malformed(partial, scratch, "an MSA genotype with an unknown allele");
            }
        }
        Ok(())
    }

    /// Resolves every kept sample's MSA value for one record (vcf2msa's
    /// `genotype_resolve` minus the external aligner). The block width is
    /// `max(len(REF), len(ALTs))`; every value is right-padded with `-` to
    /// that width, while the reference cursor advances by `len(REF)` only
    /// (v29 D29-01).
    fn msa_values(&self, scratch: &mut ConvertScratch) -> Result<(Vec<u8>, u32), Malformed> {
        let ConvertScratch {
            gts,
            allele_ranges,
            alleles,
            msa_vals,
            msa_starts,
            msa_lens,
            msa_padded,
            warned,
            ..
        } = scratch;
        msa_vals.clear();
        msa_starts.clear();
        msa_lens.clear();
        for (sample, gt) in gts.iter().enumerate() {
            if !self.keep[sample] {
                continue;
            }
            if gt.bad || gt.overflow {
                return Err(Malformed);
            }
            msa_starts.push(msa_vals.len() as u32);
            resolve_msa_value(gt, allele_ranges, alleles, msa_vals)?;
            let start = *msa_starts.last().expect("invariant: just pushed");
            msa_lens.push(msa_vals.len() as u32 - start);
        }
        let width = allele_ranges
            .iter()
            .map(|&(_, len)| len)
            .max()
            .unwrap_or(1)
            .max(1);
        if width == 1 {
            return Ok((std::mem::take(msa_vals), 1));
        }
        if !warned.padded_msa {
            warned.padded_msa = true;
            tracing::warn!(
                "reed convert --msa: indel sites are gap-padded deterministically \
                 (no external aligner; B13-6, v29 D29-01)"
            );
        }
        msa_padded.clear();
        for (i, &len) in msa_lens.iter().enumerate() {
            let start = msa_starts[i] as usize;
            msa_padded.extend_from_slice(&msa_vals[start..start + len as usize]);
            for _ in len..width {
                msa_padded.push(b'-');
            }
        }
        Ok((std::mem::take(msa_padded), width))
    }
}

/// What one input line was for the row targets (`push_line`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LineOutcome {
    /// Only line terminators: not a record, consumes no draw slot.
    Blank,
    /// A non-blank line that fails the structural checks: counted as
    /// malformed (vcf2phylip's `snp_num` slot) and skipped.
    Malformed,
    /// A structurally valid record; GT-domain rejection is decided per
    /// target from the parsed `gts`.
    Record,
}

/// A malformed record for one target (skipped with a counter).
pub(super) struct Malformed;

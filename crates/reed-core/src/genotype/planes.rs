//! Genotype bit planes (`docs/architecture-spec.md` §5.4).
//!
//! One diploid genotype is 4 bits of state, stored as four per-sample planes:
//! `slot0_alt`, `slot0_miss`, `slot1_alt`, `slot1_miss` (plus an optional
//! `phased` plane). Invariants: `slotX_alt & slotX_miss == 0` for every
//! sample; trailing padding bits are zero.
//!
//! Block layout is site-major with the planes of one site contiguous:
//! `words[(site * n_planes + plane) * W + w]`, where `W = ceil(n_samples/64)`.

/// Planes per site, excluding the optional phased plane.
pub const N_BASE_PLANES: usize = 4;
/// Plane index: slot 0 carries ALT.
pub const PLANE_SLOT0_ALT: usize = 0;
/// Plane index: slot 0 is missing.
pub const PLANE_SLOT0_MISS: usize = 1;
/// Plane index: slot 1 carries ALT.
pub const PLANE_SLOT1_ALT: usize = 2;
/// Plane index: slot 1 is missing.
pub const PLANE_SLOT1_MISS: usize = 3;
/// Plane index: genotype is phased (only present when requested).
pub const PLANE_PHASED: usize = 4;

/// Borrowed view of one site's planes (word slices of length `W`).
#[derive(Clone, Copy, Debug)]
pub struct SitePlanes<'a> {
    /// Slot 0 ALT bits.
    pub slot0_alt: &'a [u64],
    /// Slot 0 missing bits.
    pub slot0_miss: &'a [u64],
    /// Slot 1 ALT bits.
    pub slot1_alt: &'a [u64],
    /// Slot 1 missing bits.
    pub slot1_miss: &'a [u64],
    /// Phased bits, present only in phased blocks.
    pub phased: Option<&'a [u64]>,
}

impl SitePlanes<'_> {
    /// Number of `u64` words per plane.
    #[must_use]
    pub const fn n_words(&self) -> usize {
        self.slot0_alt.len()
    }
}

/// Genotype planes for a run of consecutive sites (one chunk).
#[derive(Clone, Debug, Default)]
pub struct PlaneBlock {
    words: Box<[u64]>,
    n_sites: usize,
    n_samples: usize,
    n_words_per_plane: usize,
    n_planes: usize,
}

impl PlaneBlock {
    /// Allocates a zeroed block for `n_sites` sites over `n_samples` samples.
    #[must_use]
    pub fn new(n_sites: usize, n_samples: usize, with_phase: bool) -> Self {
        let n_words_per_plane = n_samples.div_ceil(64).max(1);
        let n_planes = if with_phase {
            N_BASE_PLANES + 1
        } else {
            N_BASE_PLANES
        };
        let words = vec![0; n_sites * n_planes * n_words_per_plane].into_boxed_slice();
        Self {
            words,
            n_sites,
            n_samples,
            n_words_per_plane,
            n_planes,
        }
    }

    /// Number of sites the block holds.
    #[must_use]
    pub const fn n_sites(&self) -> usize {
        self.n_sites
    }

    /// Number of samples per site.
    #[must_use]
    pub const fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Whether the block carries the phased plane.
    #[must_use]
    pub const fn with_phase(&self) -> bool {
        self.n_planes == N_BASE_PLANES + 1
    }

    /// `u64` words per plane.
    #[must_use]
    pub const fn n_words_per_plane(&self) -> usize {
        self.n_words_per_plane
    }

    /// The plane slice for one site and plane index.
    ///
    /// # Panics
    ///
    /// Panics if `site >= n_sites` or `plane >= n_planes` (caller invariant).
    #[must_use]
    pub fn plane(&self, site: usize, plane: usize) -> &[u64] {
        assert!(site < self.n_sites && plane < self.n_planes, "invariant: plane index out of range");
        let start = (site * self.n_planes + plane) * self.n_words_per_plane;
        &self.words[start..start + self.n_words_per_plane]
    }

    /// Mutable plane slice for one site and plane index.
    ///
    /// # Panics
    ///
    /// Panics if the indices are out of range.
    pub fn plane_mut(&mut self, site: usize, plane: usize) -> &mut [u64] {
        assert!(site < self.n_sites && plane < self.n_planes, "invariant: plane index out of range");
        let start = (site * self.n_planes + plane) * self.n_words_per_plane;
        let w = self.n_words_per_plane;
        &mut self.words[start..start + w]
    }

    /// Sets one sample's slot state (`alt`, `miss`) for a site.
    ///
    /// # Panics
    ///
    /// Panics if `site >= n_sites` or `sample >= n_samples`.
    pub fn set_slot(&mut self, site: usize, sample: usize, slot: usize, alt: bool, miss: bool) {
        assert!(site < self.n_sites && sample < self.n_samples, "invariant: genotype index out of range");
        let (alt_plane, miss_plane) = match slot {
            0 => (PLANE_SLOT0_ALT, PLANE_SLOT0_MISS),
            _ => (PLANE_SLOT1_ALT, PLANE_SLOT1_MISS),
        };
        let word = sample / 64;
        let bit = 1u64 << (sample % 64);
        let write = |block: &mut Self, plane: usize, set: bool| {
            let w = &mut block.plane_mut(site, plane)[word];
            if set {
                *w |= bit;
            } else {
                *w &= !bit;
            }
        };
        write(self, alt_plane, alt);
        write(self, miss_plane, miss);
    }

    /// Sets one sample's phased bit for a site; a no-op on non-phased blocks.
    ///
    /// # Panics
    ///
    /// Panics if `site >= n_sites` or `sample >= n_samples`.
    pub fn set_phased(&mut self, site: usize, sample: usize, phased: bool) {
        assert!(site < self.n_sites && sample < self.n_samples, "invariant: genotype index out of range");
        if !self.with_phase() {
            return;
        }
        let word = sample / 64;
        let bit = 1u64 << (sample % 64);
        let w = &mut self.plane_mut(site, PLANE_PHASED)[word];
        if phased {
            *w |= bit;
        } else {
            *w &= !bit;
        }
    }

    /// Grows the block to `n_sites` sites, preserving existing content
    /// (worker scratch reuse: one allocation warm-up per chunk).
    pub fn resize_sites(&mut self, n_sites: usize) {
        if n_sites <= self.n_sites {
            return;
        }
        let site_words = self.n_planes * self.n_words_per_plane;
        let mut words = vec![0_u64; n_sites * site_words].into_boxed_slice();
        words[..self.n_sites * site_words].copy_from_slice(&self.words);
        self.words = words;
        self.n_sites = n_sites;
    }

    /// Borrowed view of one site's planes.
    ///
    /// # Panics
    ///
    /// Panics if `site >= n_sites`.
    #[must_use]
    pub fn site_planes(&self, site: usize) -> SitePlanes<'_> {
        assert!(site < self.n_sites, "invariant: site index out of range");
        SitePlanes {
            slot0_alt: self.plane(site, PLANE_SLOT0_ALT),
            slot0_miss: self.plane(site, PLANE_SLOT0_MISS),
            slot1_alt: self.plane(site, PLANE_SLOT1_ALT),
            slot1_miss: self.plane(site, PLANE_SLOT1_MISS),
            phased: if self.with_phase() {
                Some(self.plane(site, PLANE_PHASED))
            } else {
                None
            },
        }
    }

    /// Checks the per-sample plane invariants over every site (tests and
    /// debug paths; the decoders maintain them by construction).
    #[must_use]
    pub fn invariants_hold(&self) -> bool {
        (0..self.n_sites).all(|site| {
            let planes = self.site_planes(site);
            planes
                .slot0_alt
                .iter()
                .zip(planes.slot0_miss.iter())
                .chain(
                    planes
                        .slot1_alt
                        .iter()
                        .zip(planes.slot1_miss.iter()),
                )
                .all(|(&a, &m)| a & m == 0)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_and_read_slots() {
        let mut block = PlaneBlock::new(3, 130, true);
        block.set_slot(1, 0, 0, true, false); // s0 slot0 = ALT
        block.set_slot(1, 0, 1, false, false); // s0 slot1 = REF
        block.set_slot(1, 129, 0, false, true); // s129 slot0 missing
        block.set_phased(1, 0, true);

        let planes = block.site_planes(1);
        assert!(planes.slot0_alt[0] & 1 == 1);
        assert!(planes.slot0_miss[0] & 1 == 0);
        assert!(planes.slot1_alt[0] & 1 == 0);
        assert!(planes.phased.expect("invariant: phased plane present")[0] & 1 == 1);
        // Last sample lives in the second word.
        assert!(planes.slot0_miss[2] >> (129 - 128) & 1 == 1);
        assert!(block.invariants_hold());
    }

    #[test]
    fn non_phased_block_has_no_phase_plane() {
        let block = PlaneBlock::new(2, 10, false);
        assert!(!block.with_phase());
        assert!(block.site_planes(0).phased.is_none());
    }
}

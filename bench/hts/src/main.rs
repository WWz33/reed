use std::env;
use std::time::Instant;

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Read as BcfRead, Reader as BcfReader};

fn main() {
    let path = env::args().nth(1).expect("usage: bench-hts <file>");
    let t0 = Instant::now();
    let mut n_rec = 0u64;
    let mut n_alt = 0u64;
    let mut n_miss = 0u64;
    let mut n_phased = 0u64;

    let mut reader = BcfReader::from_path(&path).unwrap();
    for rec in reader.records() {
        let rec = rec.unwrap();
        n_rec += 1;
        let gts = rec.genotypes().unwrap();
        for i in 0..rec.sample_count() as usize {
            let g = gts.get(i);
            let mut phased = false;
            for a in g.iter() {
                match a {
                    GenotypeAllele::Phased(idx) => {
                        n_alt += (*idx > 0) as u64;
                        phased = true;
                    }
                    GenotypeAllele::Unphased(idx) => n_alt += (*idx > 0) as u64,
                    GenotypeAllele::PhasedMissing => {
                        n_miss += 1;
                        phased = true;
                    }
                    GenotypeAllele::UnphasedMissing => n_miss += 1,
                }
            }
            n_phased += phased as u64;
        }
    }

    let dt = t0.elapsed();
    println!(
        "n_rec={n_rec}\nn_alt={n_alt}\nn_miss={n_miss}\nn_phased={n_phased}\nsecs={:.3}",
        dt.as_secs_f64()
    );
}

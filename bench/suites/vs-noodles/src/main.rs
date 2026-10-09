use std::env;
use std::fs::File;
use std::io::BufReader;
use std::time::Instant;

use noodles::bcf as bcf;
use noodles::bgzf as bgzf;
use noodles::vcf as vcf;
use vcf::variant::record::samples::series::Value;
use vcf::variant::record::samples::series::value::genotype::Phasing;
use vcf::variant::record::samples::Series as _;

fn tally(
    gt: &dyn vcf::variant::record::samples::series::value::Genotype,
    n_alt: &mut u64,
    n_miss: &mut u64,
    n_phased: &mut u64,
) {
    let mut phased = false;
    for pair in gt.iter() {
        let (pos, phasing) = pair.unwrap();
        match pos {
            Some(i) => {
                *n_alt += (i > 0) as u64;
                if phasing == Phasing::Phased {
                    phased = true;
                }
            }
            None => *n_miss += 1,
        }
    }
    *n_phased += phased as u64;
}

fn main() {
    let path = env::args().nth(1).expect("usage: bench-noo <file>");
    let t0 = Instant::now();
    let mut n_rec = 0u64;
    let mut n_alt = 0u64;
    let mut n_miss = 0u64;
    let mut n_phased = 0u64;

    if path.ends_with(".bcf") {
        let mut reader = bcf::io::Reader::new(BufReader::new(File::open(&path).unwrap()));
        let header = reader.read_header().unwrap();
        for result in reader.records() {
            let rec = result.unwrap();
            n_rec += 1;
            let samples = rec.samples().unwrap();
            for s in samples.series() {
                let series = s.unwrap();
                if series.name(&header).unwrap() == "GT" {
                    for v in series.iter(&header) {
                        match v.unwrap() {
                            Some(Value::Genotype(gt)) => {
                                tally(gt.as_ref(), &mut n_alt, &mut n_miss, &mut n_phased)
                            }
                            _ => {}
                        }
                    }
                }
            }
        }
    } else {
        let mut reader =
            vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(File::open(&path).unwrap())));
        let header = reader.read_header().unwrap();
        for result in reader.records() {
            let rec = result.unwrap();
            n_rec += 1;
            if let Some(series) = rec.samples().select("GT") {
                for v in series.iter(&header) {
                    match v.unwrap() {
                        Some(Value::Genotype(gt)) => {
                            tally(gt.as_ref(), &mut n_alt, &mut n_miss, &mut n_phased)
                        }
                        _ => {}
                    }
                }
            }
        }
    }

    let dt = t0.elapsed();
    println!(
        "n_rec={n_rec}\nn_alt={n_alt}\nn_miss={n_miss}\nn_phased={n_phased}\nsecs={:.3}",
        dt.as_secs_f64()
    );
}

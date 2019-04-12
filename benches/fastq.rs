#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate criterion;

use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use rand::{Rng, SeedableRng};
use rand_distr::Normal;
use rand_isaac::isaac64::Isaac64Rng;
use seq_io::prelude::*;

use std::iter::repeat;

/// number of records for all benchmarks
const N: usize = 25_000;
/// standard deviation of sequence lengths relative to mean sequence length
const SEQLEN_SD_FRAC: f64 = 0.2;

/// generates 'nrecords' FASTQ records with given properties
fn gen_fastq(
    nrecords: usize,
    id_len: usize,
    desc_len: usize,
    seq_len: usize,
    sep_ids: bool,
    cr: bool,
) -> Vec<u8> {
    let newline = if cr { b"\r\n".to_vec() } else { b"\n".to_vec() };
    let mut rec: Vec<u8> = vec![];
    rec.push(b'@');
    let id: Vec<_> = repeat(b'i').take(id_len).collect();
    rec.extend(&id);
    rec.push(b' ');
    rec.extend(repeat(b'd').take(desc_len));
    rec.extend(&newline);

    let norm = Normal::new(seq_len as f64, seq_len as f64 * SEQLEN_SD_FRAC).unwrap();
    let rng = Isaac64Rng::from_seed([5; 32]);

    rng.sample_iter(&norm)
        .map(|slen| {
            let slen = slen.round() as usize;
            let mut r = rec.clone();
            r.extend(repeat(b'A').take(slen));
            r.extend(&newline);
            r.extend(b"+");
            if sep_ids {
                r.extend(Some(b' ').into_iter().chain(id.iter().cloned()));
            }
            r.extend(&newline);
            r.extend(repeat(66).take(slen));
            r.extend(&newline);
            r
        })
        .take(nrecords)
        .flat_map(|r| r)
        .collect()
}

/// generates 'nrecords' FASTQ records with fixed ID / description lengths (20 and 50), but configurable otherwise
fn with_seqlen(nrecords: usize, seq_len: usize, sep_ids: bool, cr: bool) -> Vec<u8> {
    gen_fastq(nrecords, 20, 50, seq_len, sep_ids, cr)
}

// data to be used with parallel readers that require 'static
lazy_static! {
    static ref L300: Vec<u8> = with_seqlen(N, 300, false, false);
}

macro_rules! bench {
    ($c:expr, $name:expr, $data:ident, $code:block) => {
        let name = format!("fastq {}", $name);
        let id = BenchmarkId::new(name, 0);
        $c.bench_with_input(id, $data, move |b, $data| b.iter(|| $code));
    };
}

macro_rules! fastq {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let mut reader = seq_io::fastq::Reader::new($data); //.set_store::<seq_io::fastq::multiline::MultiRangeStore>();
            while let Some(r) = reader.next() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

macro_rules! fastq_multi {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let mut reader = seq_io::fastq::multiline::Reader::new($data);
            while let Some(r) = reader.next() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

macro_rules! fastx {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let mut reader = seq_io::fastx::Reader::new($data);
            while let Some(r) = reader.next() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

macro_rules! fastx_dynamic {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let mut reader = seq_io::fastx::dynamic::reader($data, false)
                .unwrap()
                .unwrap();
            while let Some(r) = reader.next_fastx() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

macro_rules! bio {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let reader = bio::io::fastq::Reader::new($data);
            for r in reader.records() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

macro_rules! fastq_rs {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let reader = fastq::Parser::new($data);
            reader.each(|$rec| {
                $code
                true
            }).unwrap();
        });
     };
}

fn readers(c: &mut Criterion) {
    bench_readers(c, &L300);
}

fn bench_readers(c: &mut Criterion, data: &'static [u8]) {
    let mut group = c.benchmark_group("fastq");
    group.throughput(Throughput::Bytes(data.len() as u64));

    // warm up (first measurement seems not always stable)
    fastq!(group, "seqio discard", data, r, {
        black_box(r);
    });
    fastq_multi!(group, "seqio_multi discard", data, r, {
        black_box(r);
    });

    // simple parsing
    fastq!(group, "seqio borrow", data, r, {
        black_box(r);
    });
    fastq_multi!(group, "seqio_multi borrow", data, r, {
        black_box(r);
    });
    fastx!(group, "seqio_fastx borrow", data, r, {
        black_box(r);
    });
    fastx_dynamic!(group, "seqio_fastx_dynamic borrow", data, r, {
        black_box(r);
    });
    fastq_rs!(group, "fastq_rs borrow", data, r, {
        black_box(r);
    });

    // owned
    fastq_rs!(group, "fastq_rs owned", data, r, {
        let r = r.to_owned_record();
        black_box(r);
    });
    fastq!(group, "seqio owned", data, r, {
        let r = r.to_owned_record();
        black_box(r);
    });
    bench!(group, "seqio_records owned", data, {
        for res in seq_io::fastq::Reader::new(data).into_records() {
            let rec = res.unwrap();
            black_box(rec);
        }
    });
    bench!(group, "seqio_clone_into owned", data, {
        let mut reader = seq_io::fastq::Reader::new(data);
        let mut owned = seq_io::fastq::OwnedRecord::default();
        while let Some(res) = reader.next() {
            let rec = res.unwrap();
            rec.clone_into_owned(&mut owned);
        }
    });
    bio!(group, "bio owned", data, r, {
        black_box(r);
    });

    // access sequence
    fastq!(group, "seqio seq", data, r, {
        let s: &[u8] = &r.full_seq();
        black_box(s);
    });
    bench!(group, "seqio_seq_given seq", data, {
        let mut reader = seq_io::fastq::Reader::new(data);
        let mut seq = vec![];
        while let Some(res) = reader.next() {
            let rec = res.unwrap();
            rec.full_seq_given(|| &mut seq);
        }
    });
    fastx!(group, "seqio_fastx seq", data, r, {
        let s: &[u8] = &r.full_seq();
        black_box(s);
    });
    bio!(group, "bio seq", data, r, {
        let s = r.seq();
        black_box(s);
    });

    // iterate sequence
    fastq!(group, "seqio seq_iter", data, r, {
        for s in r.seq_lines() {
            black_box(s);
        }
    });
    fastx!(group, "seqio_fastx seq_iter", data, r, {
        for s in r.seq_lines() {
            black_box(s);
        }
    });

    // get sequence ID
    fastq!(group, "seqio_bytes id", data, r, {
        let id = r.id_bytes();
        black_box(id);
    });
    fastq!(group, "seqio_str id", data, r, {
        let id = r.id().unwrap();
        black_box(id);
    });

    // read into record sets without parallelism
    bench!(group, "seqio recordset", data, {
        let mut reader = seq_io::fastq::Reader::new(data);
        let mut rset = seq_io::fastq::RecordSet::default();
        while reader.read_record_set(&mut rset).unwrap() {
            for r in &rset {
                black_box(r);
            }
        }
    });

    bench!(group, "fastq_rs parallel", data, {
        let parser = fastq::Parser::new(data);
        let res: Vec<_> = parser
            .parallel_each(2, |rsets| {
                for rset in rsets {
                    for rec in rset.iter() {
                        black_box(rec);
                    }
                }
            })
            .unwrap();
        black_box(res);
    });

    // parallel
    bench!(group, "seqio_rset parallel", data, {
        let reader = seq_io::fastq::Reader::new(data);
        seq_io::parallel::read_process_recordsets(
            reader,
            2,
            2,
            |rset, _| {
                for r in &*rset {
                    black_box(r);
                }
            },
            |mut rsets| {
                while let Some(result) = rsets.next() {
                    let (rset, _): (_, &mut ()) = result.unwrap();
                    for r in &*rset {
                        black_box(r);
                    }
                }
            },
        );
    });

    bench!(group, "seqio_record parallel", data, {
        let reader = seq_io::fastq::Reader::new(data);
        seq_io::parallel::read_process_fastq_records::<_, _, _, _, (), ()>(
            reader,
            2,
            2,
            |r, _| {
                black_box(r);
            },
            |r, _| {
                black_box(r);
                None
            },
        )
        .unwrap();
    });
}

// compare different buffer capacities

fn readers_cap(c: &mut Criterion) {
    let k = 1024;
    let caps = [8 * k, 16 * k, 32 * k, 64 * k, 128 * k, 256 * k];
    let lengths = [100, 250, 500, 1000, 10000];

    for &seqlen in &lengths {
        let n = N * 250 / seqlen;
        let input = with_seqlen(n, seqlen, false, false);
        let input = input.as_slice();
        let mut group = c.benchmark_group(format!("fastq_cap_{}", seqlen));
        group.throughput(Throughput::Bytes(input.len() as u64));

        for &cap in &caps {
            let name = format!("{} {}", seqlen, cap);
            bench!(group, name, input, {
                let mut reader = seq_io::fastq::Reader::with_capacity(input, cap);
                while let Some(r) = reader.next() {
                    let _ = r.unwrap();
                }
            });
        }
    }
}

criterion_group!(benches, readers, readers_cap);
criterion_main!(benches);

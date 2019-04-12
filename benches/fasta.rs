#[macro_use]
extern crate criterion;

use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use rand::{Rng, SeedableRng};
use rand_distr::Normal;
use rand_isaac::isaac64::Isaac64Rng;
use std::iter::repeat;

use seq_io::prelude::*;
use seq_io::{fasta, fastx};

/// number of records for all benchmarks
const N: usize = 30_000;
const SEQLEN_SD_FRAC: f64 = 0.2;

/// generates 'nrecords' FASTA records with given properties
fn gen_fasta(
    nrecords: usize,
    id_len: usize,
    desc_len: usize,
    seq_len: usize,
    break_seq: Option<usize>,
    cr: bool,
) -> Vec<u8> {
    let newline = if cr { b"\r\n".to_vec() } else { b"\n".to_vec() };
    let mut rec: Vec<u8> = vec![];
    rec.push(b'>');
    rec.extend(repeat(b'i').take(id_len));
    rec.push(b' ');
    rec.extend(repeat(b'd').take(desc_len));
    rec.extend(&newline);

    let norm = Normal::new(seq_len as f64, seq_len as f64 * SEQLEN_SD_FRAC).unwrap();
    let rng = Isaac64Rng::from_seed([5; 32]);

    rng.sample_iter(&norm)
        .map(|slen| {
            let slen = slen.round() as usize;
            let mut r = rec.clone();
            let seq_line = break_seq.unwrap_or(slen);
            let rest = slen % seq_line;
            for n in repeat(seq_line).take(slen / seq_line).chain(if rest > 0 {
                Some(rest)
            } else {
                None
            }) {
                r.extend(repeat(b'A').take(n));
                r.extend(&newline);
            }
            assert!(
                r.len() - rec.len()
                    == slen
                        + slen / seq_line * newline.len()
                        + if slen % seq_line > 0 {
                            newline.len()
                        } else {
                            0
                        }
            );
            r
        })
        .take(nrecords)
        .flat_map(|r| r)
        .collect()
}

/// generates 'nrecords' FASTA with fixed ID / description lengths (20 and 50), but configurable otherwise
fn with_seqlen(nrecords: usize, seq_len: usize, break_seq: Option<usize>, cr: bool) -> Vec<u8> {
    gen_fasta(nrecords, 20, 50, seq_len, break_seq, cr)
}

macro_rules! bench {
    ($c:expr, $name:expr, $data:ident, $code:block) => {
        let name = format!("fasta {}", $name);
        let id = BenchmarkId::new(name, 0);
        $c.bench_with_input(id, $data, move |b, $data| b.iter(|| $code));
    };
}

macro_rules! fasta {
    ($c:expr, $name:expr, $data:ident, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let mut reader = fasta::Reader::new($data);
            while let Some(r) = reader.next() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

macro_rules! fasta_single {
    ($c:expr, $name:expr, $data:ident, $Store:ty, $rec:ident, $code:block) => {
        bench!($c, $name, $data, {
            let mut reader: fasta::single_line::Reader<_, _, $Store> =
                fasta::single_line::Reader::new($data).set_store();
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
            let mut reader = fastx::Reader::new($data);
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
            let mut reader = fastx::dynamic::reader($data, false).unwrap().unwrap();
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
            let reader = bio::io::fasta::Reader::new($data);
            for r in reader.records() {
                let $rec = r.unwrap();
                $code
            }
        });
    };
}

fn readers(c: &mut Criterion) {
    let data = with_seqlen(N, 300, None, false);
    reader_bench(c, "fasta", &data, true);

    let data = with_seqlen(N, 300, Some(80), false);
    reader_bench(c, "fasta_multiline", &data, false);
}

fn reader_bench(c: &mut Criterion, name: &str, data: &[u8], single_line: bool) {
    let mut group = c.benchmark_group(name);
    group.throughput(Throughput::Bytes(data.len() as u64));

    // simple parsing
    fasta!(group, "seqio borrow", data, r, {
        black_box(r);
    });
    fastx!(group, "seqio_fastx borrow", data, r, {
        black_box(r);
    });
    fastx_dynamic!(group, "seqio_fastx_dynamic borrow", data, r, {
        black_box(r);
    });

    if single_line {
        fasta_single!(
            group,
            "seqio_single_linestore borrow",
            data,
            fasta::LineStore,
            r,
            {
                black_box(r);
            }
        );
        fasta_single!(
            group,
            "seqio_single borrow",
            data,
            fasta::single_line::RangeStore,
            r,
            {
                black_box(r);
            }
        );
    }

    // owned
    bench!(group, "seqio_clone_into owned", data, {
        let mut reader = fasta::Reader::new(data);
        let mut owned = fasta::OwnedRecord::default();
        while let Some(res) = reader.next() {
            let rec = res.unwrap();
            rec.clone_into_owned(&mut owned);
        }
    });
    bio!(group, "bio owned", data, r, {
        black_box(r);
    });
    fasta!(group, "seqio owned", data, r, {
        let r = r.to_owned_record();
        black_box(r);
    });
    bench!(group, "seqio_records owned", data, {
        for res in fasta::Reader::new(data).into_records() {
            let rec = res.unwrap();
            black_box(rec);
        }
    });

    // access sequence
    fasta!(group, "seqio seq", data, r, {
        let s: &[u8] = &r.full_seq();
        black_box(s);
    });
    bench!(group, "seqio_seq_given seq", data, {
        let mut reader = seq_io::fasta::Reader::new(data);
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
    fasta!(group, "seqio seq_iter", data, r, {
        for s in r.seq_lines() {
            black_box(s);
        }
    });
    fastx!(group, "seqio_fastx seq_iter", data, r, {
        for s in r.seq_lines() {
            black_box(s);
        }
    });
    fastx_dynamic!(group, "seqio_fastx_dynamic seq_iter", data, r, {
        for s in r.seq_lines() {
            black_box(s);
        }
    });

    // get sequence ID
    fasta!(group, "seqio_bytes id", data, r, {
        let id = r.id_bytes();
        black_box(id);
    });
    fasta!(group, "seqio_str id", data, r, {
        let id = r.id().unwrap();
        black_box(id);
    });

    // read into record sets without parallelism
    bench!(group, "seqio recordset", data, {
        let mut reader = fasta::Reader::new(data);
        let mut rset = fasta::RecordSet::default();
        while reader.read_record_set(&mut rset).unwrap() {
            for r in &rset {
                black_box(r);
            }
        }
    });

    // parallel
    bench!(group, "seqio_rset parallel", data, {
        let reader = fasta::Reader::new(data);
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
        let reader = fasta::Reader::new(data);
        seq_io::parallel::read_process_fasta_records::<_, _, _, _, (), ()>(
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
        let input = with_seqlen(n, seqlen, None, false);
        let input = input.as_slice();
        let mut group = c.benchmark_group(format!("fasta_cap_{}", seqlen));
        group.throughput(Throughput::Bytes(input.len() as u64));

        for &cap in &caps {
            let name = format!("{} {}k", seqlen, cap);
            bench!(group, name, input, {
                let mut reader = seq_io::fasta::Reader::with_capacity(input, cap);
                while let Some(r) = reader.next() {
                    let _ = r.unwrap();
                }
            });
        }
    }
}

criterion_group!(benches, readers, readers_cap);
criterion_main!(benches);

#![feature(test)]
#![allow(non_snake_case, unused_variables)]

extern crate test;
extern crate seq_io;
extern crate fastq;
extern crate bio;
#[macro_use]
extern crate lazy_static;

use std::io::Cursor;
use std::iter::repeat;
use test::Bencher;


/// generates 'nrecords' FASTQ records with given properties
fn gen_fastq(nrecords: usize, id_len: usize, desc_len: usize,
             seq_len: usize, sep_ids: bool, cr: bool) -> Vec<u8> {

    let newline = if cr { b"\r\n".to_vec() } else { b"\n".to_vec() };
    let mut rec: Vec<u8> = vec![];
    rec.push(b'@');
    let id: Vec<_> = repeat(b'i').take(id_len).collect();
    rec.extend(&id);
    rec.push(b' ');
    rec.extend(repeat(b'd').take(desc_len));
    rec.extend(&newline);
    rec.extend(repeat(b'A').take(seq_len));
    rec.extend(&newline);
    rec.push(b'+');
    if sep_ids {
        rec.extend(Some(b' ').into_iter().chain(id.into_iter()));
    }
    rec.extend(&newline);
    rec.extend(repeat(66).take(seq_len));
    rec.extend(&newline);

    (0..nrecords).flat_map(|_| rec.clone()).collect()
}

/// generates 'nrecords' FASTQ records with fixed ID / description lengths (20 and 50), but configurable otherwise
fn with_seqlen(nrecords: usize, seq_len: usize, sep_ids: bool, cr: bool) -> Vec<u8> {
    gen_fastq(nrecords, 20, 50, seq_len, sep_ids, cr)
}


/// number of records for all benchmarks
const N: usize = 100_000;

/// data to be used with parallel readers (require 'static)
lazy_static! {
    static ref L500: Vec<u8> = { with_seqlen(N, 500, false, false) };
}

macro_rules! bench {
    ($name:ident, $seqlen:expr, $lbreak:expr, $cursor:ident, $code:block) => {
        #[bench]
        fn $name(b: &mut Bencher) {
            let data = with_seqlen(N, $seqlen, $lbreak, false);
            b.bytes = data.len() as u64;
            b.iter(|| {
                let $cursor = Cursor::new(&data);
                $code
            });
        }
     };
}

macro_rules! bench_static500 {
    ($name:ident, $cursor:ident, $code:block) => {
        #[bench]
        fn $name(b: &mut Bencher) {
            let data: &'static Vec<u8> = &L500 as &Vec<u8>;
            b.bytes = data.len() as u64;
            b.iter(move || {
                let $cursor = Cursor::new(data);
                $code
            });
        }
     };
}

macro_rules! fastq {
    ($name:ident, $seqlen:expr, $lbreak:expr, $rec:ident, $code:block) => {
        bench!($name, $seqlen, $lbreak, cursor, {
            let mut reader = seq_io::fastq::Reader::new(cursor);
            while let Some($rec) = reader.next() {
                let $rec = $rec.unwrap();
                $code
            }
        });
     };
}

macro_rules! fastq_rs {
    ($name:ident, $seqlen:expr, $lbreak:expr, $rec:ident, $code:block) => {
        bench!($name, $seqlen, $lbreak, cursor, {
            let reader = fastq::Parser::new(cursor);
            reader.each(|$rec| {
                $code
                true
            }).unwrap();
        });
     };
}

macro_rules! bio {
    ($name:ident, $seqlen:expr, $lbreak:expr, $rec:ident, $code:block) => {
        bench!($name, $seqlen, $lbreak, cursor, {
            let reader = bio::io::fastq::Reader::new(cursor);
            for r in reader.records() {
                let $rec = r.unwrap();
                $code
            }
        });
     };
}


fastq!(fastq_iter_200____seqio, 200, false, r, {});
fastq!(fastq_iter_500____seqio, 500, false, r, {});
fastq!(fastq_iter_1000____seqio, 1000, false, r, {});

fastq_rs!(fastq_iter_200____fastqrs, 200, false, r, {});
fastq_rs!(fastq_iter_500____fastqrs, 500, false, r, {});
fastq_rs!(fastq_iter_1000____fastqrs, 1000, false, r, {});

bio!(fastq_iter_200__owned__bio,  200,  false, r, {});
bio!(fastq_iter_500__owned__bio,  500,  false, r, {});
bio!(fastq_iter_1000__owned__bio, 1000, false, r, {});


// bench!(fastq_iter_500__buf5k__seqio, 500, false, cursor, {
//     let mut reader = seq_io::fastq::Reader::with_capacity(cursor, 5*1024);
//     while let Some(rec) = reader.next() {
//         rec.unwrap();
//     }
// });
//
// bench!(fastq_iter_500__buf5k__fastqrs, 500, false, cursor, {
//     let reader = fastq::Parser::with_capacity(cursor, 5*1024);
//     reader.each(|_| {true}).unwrap();
// });


fastq!(fastq_iter_500__owned__seqio,  500,  false, r, {
    let _ = r.to_owned_record();
});

fastq_rs!(fastq_iter_500__owned__fastqrs,  500,  false, r, {
    let _ = r.to_owned_record();
});


// parallel processing (just iterate)

bench!(fastq_iter_500_recset__parallel_seqio, 500, false, cursor, {
    let reader = seq_io::fastq::Reader::new(cursor);
    seq_io::parallel::read_parallel(reader, 2, 2, |rset| {
        for _ in &*rset {}
    }, |rsets| {
        while let Some(result) = rsets.next() {
            let (rset, _) = result.unwrap();
            for _ in &*rset {}
        }
    });
});

bench!(fastq_iter_500_records__parallel_seqio, 500, false, cursor, {
    let reader = seq_io::fastq::Reader::new(cursor);
    seq_io::parallel::parallel_fastq::<_, (), _, _, ()>(reader, 2, 2, |_, _| {}, |_, _| {None}).unwrap();
});

// bench_static500!(fastq_iter_500_static__parallel_seqio, cursor, {
//     let reader = seq_io::fastq::Reader::new(cursor);
//     seq_io::parallel::read_parallel(reader, 2, 2, |rset| { for _ in &*rset {} }, |rsets| {
//         while let Some(result) = rsets.next() {
//             let (rset, _) = result.unwrap();
//             for _ in &*rset {}
//         }
//     });
// });

// input needs static lifetime, but impact of using lazy_static is neglectable
bench_static500!(fastq_iter_500_recset__parallel_fastqrs, cursor, {
    let parser = fastq::Parser::new(cursor);
    let res: Vec<_> = parser.parallel_each(2, |rsets| {
        for rset in rsets {
            for _ in rset.iter() {}
        }
    }).unwrap();
});


// iterate over record sets without parallelism

bench!(fastq_iter_500_recset___seqio, 500, false, cursor, {
    let mut reader = seq_io::fastq::Reader::new(cursor);
    let mut rset = seq_io::fastq::RecordSet::default();
    while let Some(result) = reader.read_record_set(&mut rset) {
        result.unwrap();
        for _ in &rset {}
    }
});

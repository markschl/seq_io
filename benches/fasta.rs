#![feature(test)]
#![allow(non_snake_case, unused_variables)]

extern crate test;
extern crate seq_io;
extern crate bio;

use std::io::Cursor;
use std::iter::repeat;
use test::Bencher;

use seq_io::fasta;


/// generates 'nrecords' FASTA records with given properties
fn gen_fasta(nrecords: usize, id_len: usize, desc_len: usize,
             seq_len: usize, break_seq: Option<usize>, cr: bool) -> Vec<u8> {

    let newline = if cr { b"\r\n".to_vec() } else { b"\n".to_vec() };
    let mut rec: Vec<u8> = vec![];
    rec.push(b'>');
    rec.extend(repeat(b'i').take(id_len));
    rec.push(b' ');
    rec.extend(repeat(b'd').take(desc_len));
    rec.extend(&newline);

    let seq: Vec<_> = repeat(b'A').take(seq_len).collect();
    for s in seq.chunks(break_seq.unwrap_or(seq_len)) {
        rec.extend(s);
        rec.extend(&newline);
    }

    (0..nrecords).flat_map(|_| rec.clone()).collect()
}

/// generates 'nrecords' FASTA with fixed ID / description lengths (20 and 50), but configurable otherwise
fn with_seqlen(nrecords: usize, seq_len: usize, break_seq: Option<usize>, cr: bool) -> Vec<u8> {
    gen_fasta(nrecords, 20, 50, seq_len, break_seq, cr)
}

/// number of records for all benchmarks
const N: usize = 100_000;



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

macro_rules! fasta {
    ($name:ident, $seqlen:expr, $lbreak:expr, $rec:ident, $code:block) => {
        bench!($name, $seqlen, $lbreak, cursor, {
            let mut reader = fasta::Reader::new(cursor);
            while let Some(Ok($rec)) = reader.next() $code
        });
     };
}

macro_rules! bio {
    ($name:ident, $seqlen:expr, $lbreak:expr, $rec:ident, $code:block) => {
        bench!($name, $seqlen, $lbreak, cursor, {
            let reader = bio::io::fasta::Reader::new(cursor);
            for r in reader.records() {
                let $rec = r.unwrap();
                $code
            }

        });
     };
}




fasta!(fasta_iter_200____seqio,  200,  None, r, {});
fasta!(fasta_iter_500____seqio,  500,  None, r, {});
fasta!(fasta_iter_1000____seqio, 1000, None, r, {});

bio!(fasta_iter_200__owned__bio,  200,  None, r, {});
bio!(fasta_iter_500__owned__bio,  500,  None, r, {});
bio!(fasta_iter_1000__owned__bio, 1000, None, r, {});

fasta!(fasta_iter_500_multiline___seqio,  500,  Some(100), r, {});
bio!(fasta_iter_500_multiline_owned__bio, 500,  Some(100), r, {});

fasta!(fasta_iter_500__owned__seqio,  500,  None, r, {
    let _ = r.to_owned_record();
});

fasta!(fasta_iter_500_multiline_owned__seqio,  500,  Some(100), r, {
    let _ = r.to_owned_record();
});


// parallel

bench!(fasta_iter_500_recset__parallel_seqio, 500, None, cursor, {
    let reader = fasta::Reader::new(cursor);
    seq_io::parallel::read_parallel(reader, 2, 2, |rset| { for _ in &*rset {} }, |rsets| {
        while let Some(result) = rsets.next() {
            let (rset, _) = result.unwrap();
            for _ in &*rset {}
        }
    });
});

bench!(fasta_iter_500_records__parallel_seqio, 500, None, cursor, {
    let reader = fasta::Reader::new(cursor);
    seq_io::parallel::parallel_fasta::<_, (), _, _>(reader, 2, 2, |_, _| {}, |_, _| {true}).unwrap();
});

// read into record sets without parallelism

bench!(fasta_iter_500_recset___seqio, 500, None, cursor, {
    let mut reader = fasta::Reader::new(cursor);
    let mut rset = fasta::RecordSet::default();
    while let Some(result) = reader.read_record_set(&mut rset) {
        result.unwrap();
        for _ in &rset {}
    }
});

fasta!(fasta_seq_500____seqio,  500,  None, r, {
    for _ in r.seq_lines() {}
});

bio!(fasta_seq_500____bio,  500,  None, r, {
    let _ = r.seq();
});

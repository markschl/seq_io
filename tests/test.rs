#![allow(dead_code)]

extern crate seq_io;

use seq_io::fasta::{self,Record};
use seq_io::fastq::{self,Record as FastqRecord};



const FASTA: &'static [&'static [u8]; 11] = &[
   b">id desc",
   b"ACCGTAGGCT",
   b"CCGTAGGCTG",
   b"CGTAGGCTGA",
   b"GTAGGCTGAA",
   b"CCCC",
   b">id2",
   b"ATTGTTGTTT",
   b"ATTGTTGTTT",
   b"ATTGTTGTTT",
   b"GGGG"
];


fn concat_lines(lines: &[&[u8]], terminator: &[u8], last: bool) -> Vec<u8> {
   let mut out: Vec<_> = lines.iter()
        .flat_map(|s| s.iter().chain(terminator))
        .cloned()
        .collect();
   if !last {
      let l = out.len();
      out.truncate(l - terminator.len());
   }
   out
}

#[test]
fn test_fasta_reader() {
   let expected = [
      (Ok("id"), Some(Ok("desc")), (1, 6)),
      (Ok("id2"), None, (7, 11))
   ];
   let lterms: [&[u8]; 2] = [b"\n", b"\r\n"];

   // try different line endings
   for t in lterms.into_iter() {

      let fasta = concat_lines(FASTA, *t, true);
      let exp_seqs: Vec<_> = expected.iter().map(|&(_, _, (start, end))| (
         // raw sequence
         concat_lines(&FASTA[start..end], *t, false),
         // concatenated sequence
         FASTA[start..end].concat().to_vec()
      )).collect();

      // try different initial capacities to test
      // buffer growing feature
      for cap in 3 .. 100 {
         let mut exp_iter = expected.iter().zip(&exp_seqs);
         let mut reader = fasta::Reader::with_capacity(fasta.as_slice(), cap);
          while let Some((&(ref id, ref desc, _), &(ref raw_seq, ref seq))) = exp_iter.next() {
              let record = reader.next().unwrap().expect(&format!("Error reading record at cap. {}", cap));

              assert_eq!(record.id(), *id, "ID mismatch at cap. {}", cap);
              assert_eq!(record.desc(), *desc, "desc mismatch at cap. {}", cap);
              assert_eq!(record.seq(), raw_seq.as_slice(), "raw seq mismatch at cap. {}", cap);
              assert_eq!(record.owned_seq().as_slice(), seq.as_slice(), "seq mismatch at cap. {}", cap);

              let owned = record.to_owned_record();
              assert_eq!(owned.id(), *id, "ID mismatch at cap. {}", cap);
              assert_eq!(owned.desc(), *desc, "desc mismatch at cap. {}", cap);
              assert_eq!(owned.seq(), seq.as_slice(), "seq mismatch at cap. {}", cap);

          }
      }
   }
}

#[test]
fn test_fasta_invalid_start() {
    let mut reader = fasta::Reader::new(&b"id\nATGC"[..]);
    let rec = reader.next().unwrap();
    assert!(rec.is_err() && format!("{}", rec.err().unwrap()).contains("Expected > at file start"));
}

#[test]
fn test_fasta_truncated() {
    let mut reader = fasta::Reader::new(&b">id\n"[..]);
    let rec = reader.next().unwrap();
    assert!(rec.is_err() && format!("{}", rec.err().unwrap()).contains("Unexpected end of input"));
}


#[test]
fn test_fasta_none_after_err() {
    let mut reader = fasta::Reader::new(&b"id\nATGC"[..]);
    assert!(reader.next().unwrap().is_err());
    assert!(reader.next().is_none());
}


/// FASTQ

const FASTQ: &'static [u8] = b"@id desc
ATGC
+
~~~~
@id2
ATGC
+
~~~~";


#[test]
fn test_fastq_reader() {
   let expected = [
      (Ok("id"), Some(Ok("desc")), b"ATGC", b"~~~~"),
      (Ok("id2"), None, b"ATGC", b"~~~~")
   ];

   // try different initial capacities to test
   // buffer growing feature
   for cap in 3 .. 100 {
      let mut exp_iter = expected.iter();
      let mut reader = fastq::Reader::with_capacity(FASTQ, cap);
      while let Some(&(ref id, ref desc, ref seq, ref qual)) = exp_iter.next() {
           let record = reader.next().unwrap().expect(&format!("Error reading record at cap. {}", cap));

           assert_eq!(record.id(), *id, "ID mismatch at cap. {}", cap);
           assert_eq!(record.desc(), *desc, "desc mismatch at cap. {}", cap);
           assert_eq!(record.seq(), *seq, "seq at cap. {}", cap);
           assert_eq!(record.qual(), *qual, "qual mismatch at cap. {}", cap);

           let owned = record.to_owned_record();
           assert_eq!(owned.id(), *id, "ID mismatch at cap. {}", cap);
           assert_eq!(owned.desc(), *desc, "desc mismatch at cap. {}", cap);
           assert_eq!(owned.seq(), *seq, "seq at cap. {}", cap);
           assert_eq!(owned.qual(), *qual, "qual mismatch at cap. {}", cap);
       }
   }
}


#[test]
fn test_fastq_invalid_start() {
    let mut reader = fastq::Reader::new(&b"id\nATGC"[..]);
    let rec = reader.next().unwrap();
    assert!(rec.is_err() && format!("{}", rec.err().unwrap()).contains("Expected '@' but found 'i'"));
}

#[test]
fn test_fastq_truncated() {
    let mut reader = fastq::Reader::new(&b"@id\nATGC\n+"[..]);
    let rec = reader.next().unwrap();
    assert!(rec.is_err() && format!("{}", rec.err().unwrap()).contains("Unexpected end of input"));
}

#[test]
fn test_fastq_unequal() {
    let mut reader = fastq::Reader::new(&b"@id\nATGC\n+\n~~"[..]);
    let rec = reader.next().unwrap();
    assert!(rec.is_err() && format!("{}", rec.err().unwrap()).contains("Unequal lengths"));
}

#[test]
fn test_fastq_no_sep() {
    let mut reader = fastq::Reader::new(&b"@id\nATGC\n~~~~\n@id2"[..]);
    let rec = reader.next().unwrap();
    assert!(rec.is_err() && format!("{}", rec.err().unwrap()).contains("Expected '+' but found '~'"));
}

#[test]
fn test_fastq_none_after_err() {
    let mut reader = fastq::Reader::new(&b"@id\nATGC"[..]);
    assert!(reader.next().unwrap().is_err());
    assert!(reader.next().is_none());
}

// FASTA writing

#[test]
fn test_fasta_write_head() {
    let mut out = vec![];
    fasta::write_head(&mut out, b"id desc").unwrap();
    assert_eq!(&out, b">id desc\n");
}

#[test]
fn test_fasta_write_seq() {
    let mut out = vec![];
    fasta::write_seq(&mut out, b"ATGC").unwrap();
    assert_eq!(&out, b"ATGC\n");
}

#[test]
fn test_fasta_write_seq_wrap() {
    let mut out = vec![];
    fasta::write_wrap_seq(&mut out, b"ATGCA", 2).unwrap();
    assert_eq!(&out, b"AT\nGC\nA\n");
}

#[test]
fn test_fasta_write_seq_iter() {
    let mut out = vec![];
    fasta::write_seq_iter(&mut out, b"ATGCA".chunks(2)).unwrap();
    assert_eq!(&out, b"ATGCA\n");
}

#[test]
fn test_fasta_write_seq_iter_wrap() {
    for size in 1..11 {
        let mut out = vec![];
        fasta::write_wrap_seq_iter(&mut out, b"AAAATTTTGGG".chunks(size), 3).unwrap();
        assert_eq!(&out, b"AAA\nATT\nTTG\nGG\n");

        let mut out = vec![];
        fasta::write_wrap_seq_iter(&mut out, b"AAAATTTTGGG".chunks(size), 4).unwrap();
        assert_eq!(&out, b"AAAA\nTTTT\nGGG\n");
    }
}


// FASTQ writing

#[test]
fn test_fastq_write_record() {
    let mut out = vec![];
    let mut rdr = fastq::Reader::new(FASTQ);
    while let Some(Ok(r)) = rdr.next() {
        r.write(&mut out).unwrap();
    }
    assert_eq!(&out, &FASTQ);
}

#[test]
fn test_fastq_write_record_unchanged() {
    let mut out = vec![];
    let mut rdr = fastq::Reader::new(FASTQ);
    while let Some(Ok(r)) = rdr.next() {
        r.write_unchanged(&mut out).unwrap();
    }
    assert_eq!(&out, &FASTQ);
}


extern crate seq_io;

use seq_io::fastq::{self, Record};


const FASTQ: &'static [u8] = b"@id desc
ATGC
+
~~~~
@id2
ATGC
+
~~~~
";


#[test]
fn test_fastq_reader() {
    let expected = [
        (
            Ok("id"),
            Some(Ok("desc")),
            &b"id desc"[..],
            b"ATGC",
            b"~~~~",
        ),
        (Ok("id2"), None, &b"id2"[..], b"ATGC", b"~~~~"),
    ];

    // try different initial capacities to test
    // buffer growing feature
    for cap in 3..100 {
        let mut exp_iter = expected.iter();
        let mut reader = fastq::Reader::with_capacity(FASTQ, cap);
        while let Some(&(id, desc, head, seq, qual)) = exp_iter.next() {
            let record = reader
                .next()
                .unwrap()
                .expect(&format!("Error reading record at cap. {}", cap));

            assert_eq!(record.id(), id, "ID mismatch at cap. {}", cap);
            assert_eq!(record.desc(), desc, "desc mismatch at cap. {}", cap);
            assert_eq!(record.head(), head, "head mismatch at cap. {}", cap);
            assert_eq!(record.seq(), seq, "seq at cap. {}", cap);
            assert_eq!(record.qual(), qual, "qual mismatch at cap. {}", cap);

            let owned = record.to_owned_record();
            assert_eq!(owned.id(), id, "ID mismatch at cap. {}", cap);
            assert_eq!(owned.desc(), desc, "desc mismatch at cap. {}", cap);
            assert_eq!(owned.head(), head, "head mismatch at cap. {}", cap);
            assert_eq!(owned.seq(), seq, "seq at cap. {}", cap);
            assert_eq!(owned.qual(), qual, "qual mismatch at cap. {}", cap);
        }
    }
}


#[test]
fn test_fastq_invalid_start() {
    let mut reader = fastq::Reader::new(&b"id\nATGC\n+\n~~~~"[..]);
    let rec = reader.next().unwrap();
    assert!(
        rec.is_err() && format!("{}", rec.err().unwrap()).contains("Expected '@' but found 'i'")
    );
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
    assert!(
        rec.is_err() && format!("{}", rec.err().unwrap()).contains("Expected '+' but found '~'")
    );
}

#[test]
fn test_fastq_none_after_err() {
    let mut reader = fastq::Reader::new(&b"@id\nATGC"[..]);
    assert!(reader.next().unwrap().is_err());
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_empty() {
    let mut reader = fastq::Reader::new(&b""[..]);
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_empty_lines_end() {
    // up to 3 newlines at end possible
    let mut reader = fastq::Reader::new(&b"@id\nATGC\n+\n~~~~\n\n\n"[..]);
    assert!(reader.next().unwrap().unwrap().id_bytes() == b"id");
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_no_newline_end() {
    let mut reader = fastq::Reader::new(&b"@id\nATGC\n+\n~~~~"[..]);
    assert!(reader.next().unwrap().unwrap().id_bytes() == b"id");
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_recset() {
    let expected = [
        (
            Ok("id"),
            Some(Ok("desc")),
            &b"id desc"[..],
            b"ATGC",
            b"~~~~",
        ),
        (Ok("id2"), None, &b"id2"[..], b"ATGC", b"~~~~"),
    ];
    let mut rset = fastq::RecordSet::default();
    let mut reader = fastq::Reader::new(FASTQ);
    reader.read_record_set(&mut rset).unwrap().unwrap();
    let mut rset_iter = rset.into_iter();
    for &(id, desc, head, seq, qual) in expected.into_iter() {
        let rec = rset_iter.next().unwrap();
        assert_eq!(rec.id(), id);
        assert_eq!(rec.desc(), desc);
        assert_eq!(rec.head(), head);
        assert_eq!(rec.seq(), seq);
        assert_eq!(rec.qual(), qual);
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

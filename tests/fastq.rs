extern crate seq_io;
#[macro_use]
extern crate matches;
#[macro_use]
extern crate lazy_static;

use seq_io::fastq::*;
use std::io;
use std::str::Utf8Error;

const FASTQ: &[u8] = b"@id desc
ATGC
+
~~~~
@id2
ATGC
+
~~~~
";

type ExpectedRecord = (
    Result<&'static str, Utf8Error>,
    Option<Result<&'static str, Utf8Error>>,
    &'static [u8],
    &'static [u8; 4],
    &'static [u8; 4],
);

lazy_static! {
    static ref EXPECTED: [ExpectedRecord; 2] = [
        (
            Ok("id"),
            Some(Ok("desc")),
            &b"id desc"[..],
            b"ATGC",
            b"~~~~"
        ),
        (Ok("id2"), None, &b"id2"[..], b"ATGC", b"~~~~"),
    ];
}

#[test]
fn test_fastq_reader() {
    // try different initial capacities to test
    // buffer growing feature
    for cap in 3..100 {
        let mut exp_iter = EXPECTED.iter();
        let mut reader = Reader::with_capacity(FASTQ, cap);
        while let Some(&(id, desc, head, seq, qual)) = exp_iter.next() {
            let record = reader
                .next()
                .unwrap()
                .unwrap_or_else(|_| panic!("Error reading record at cap. {}", cap));

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
        assert!(reader.next().is_none());
    }
}

#[test]
fn test_fastq_policy() {
    struct Limited(usize);

    impl seq_io::BufPolicy for Limited {
        fn grow_to(&mut self, current_size: usize) -> Option<usize> {
            if current_size > self.0 {
                return None
            }
            Some(current_size * 2)
        }
    }

    let mut reader = Reader::with_capacity(&b">id\nAT\nGC\n"[..], 3)
        .set_policy(Limited(5));
    let res = reader.next().unwrap();
    assert_matches!(res, Err(Error::BufferLimit));
}

#[test]
fn test_fastq_invalid_start() {
    let mut reader = Reader::new(&b"@id1\nA\n+\n~\nid\nATGC\n+\n~~~~"[..]);
    reader.next().unwrap().unwrap();
    let rec = reader.next().unwrap();
    assert_matches!(
        rec,
        Err(Error::InvalidStart {
            found: b'i',
            pos: ErrorPosition { line: 5, id: None }
        })
    );
}

#[test]
#[allow(unused_variables)]
fn test_fastq_truncated() {
    let mut reader = Reader::new(&b"@id\nATGC\n+"[..]);
    let rec = reader.next().unwrap();
    let p = ErrorPosition {
        line: 3,
        id: Some("id".to_string()),
    };
    assert_matches!(rec, Err(Error::UnexpectedEnd { pos: p }));
}

#[test]
#[allow(unused_variables)]
fn test_fastq_unequal() {
    let mut reader = Reader::new(&b"@id\nATGC\n+\n~~"[..]);
    let rec = reader.next().unwrap();
    let p = ErrorPosition {
        line: 1,
        id: Some("id".to_string()),
    };
    assert_matches!(
        rec,
        Err(Error::UnequalLengths {
            seq: 4,
            qual: 2,
            pos: p
        })
    );
}

#[test]
#[allow(unused_variables)]
fn test_fastq_no_sep() {
    let mut reader = Reader::new(&b"@id\nATGC\n~~~~\n"[..]);
    let rec = reader.next().unwrap();
    let p = ErrorPosition {
        line: 3,
        id: Some("id".to_string()),
    };
    assert_matches!(
        rec,
        Err(Error::InvalidSep {
            found: b'~',
            pos: p
        })
    );
}

#[test]
fn test_fastq_none_after_err() {
    let mut reader = Reader::new(&b"@id\nATGC"[..]);
    assert!(reader.next().unwrap().is_err());
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_empty() {
    let mut reader = Reader::new(&b""[..]);
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_empty_lines_end() {
    // up to 3 newlines at end possible
    let mut reader = Reader::new(&b"@id\nATGC\n+\n~~~~\n\n\n"[..]);
    assert_eq!(reader.next().unwrap().unwrap().id_bytes(), b"id");
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_no_newline_end() {
    let mut reader = Reader::new(&b"@id\nATGC\n+\n~~~~"[..]);
    assert_eq!(reader.next().unwrap().unwrap().id_bytes(), b"id");
    assert!(reader.next().is_none());
}

#[test]
fn test_fastq_recset() {
    for cap in 3..400 {
        let mut reader = Reader::with_capacity(FASTQ, cap);
        let mut rsets = vec![];
        loop {
            rsets.push(RecordSet::default());
            let rset = rsets.last_mut().unwrap();
            if let Some(r) = reader.read_record_set(rset) {
                r.unwrap();
            } else {
                break;
            }
        }
        let mut rset_iter = rsets.iter().flat_map(|r| r.into_iter());

        let mut reader = Reader::new(FASTQ);
        while let Some(r0) = reader.next() {
            let rec = rset_iter.next().unwrap();
            let r0 = r0.unwrap();
            assert_eq!(rec.id(), r0.id());
            assert_eq!(rec.desc(), r0.desc());
            assert_eq!(rec.head(), r0.head());
            assert_eq!(rec.seq(), r0.seq());
            assert_eq!(rec.qual(), r0.qual());
        }
    }
}

#[test]
fn test_fastq_parallel() {
    for cap in 3..400 {
        let par_reader = Reader::with_capacity(FASTQ, cap);
        let mut reader = Reader::new(FASTQ);

        seq_io::parallel::parallel_fastq(
            par_reader,
            1,
            2,
            |_, out| {
                *out = ();
            },
            |rec, _| {
                // runs in main thread
                let r0 = reader.next().unwrap().unwrap();
                assert_eq!(rec.id(), r0.id());
                assert_eq!(rec.desc(), r0.desc());
                assert_eq!(rec.head(), r0.head());
                assert_eq!(rec.seq(), r0.seq());
                assert_eq!(rec.qual(), r0.qual());
                None::<()>
            },
        ).unwrap();
    }
}

#[test]
fn test_fastq_seek() {
    for cap in 3..31 {
        let mut reader = Reader::with_capacity(
            io::Cursor::new(b"@s1\nA\n+\n~\n@s2\nA\n+\n~\n@s3\nA\n+\n~\n"),
            cap,
        );
        let pos0 = reader.position().clone();
        let rec1 = reader.next().unwrap().unwrap().to_owned_record();
        let pos1 = reader.position().clone();
        reader.next().unwrap().unwrap();
        let rec3 = reader.next().unwrap().unwrap().to_owned_record();
        let pos3 = reader.position().clone();

        assert_eq!(pos0.line(), 1);
        assert_eq!(pos1.line(), 1);
        assert_eq!(pos3.line(), 9);
        assert_eq!(pos0.byte(), 0);
        assert_eq!(pos1.byte(), 0);
        assert_eq!(pos3.byte(), 20);

        reader.seek(&pos0).unwrap();
        assert_eq!(reader.next().unwrap().unwrap().id_bytes(), rec1.id_bytes());
        reader.next().unwrap().unwrap();
        assert_eq!(reader.next().unwrap().unwrap().id_bytes(), rec3.id_bytes());

        reader.seek(&pos1).unwrap();
        assert_eq!(reader.next().unwrap().unwrap().id_bytes(), rec1.id_bytes());
        reader.next().unwrap().unwrap();
        assert_eq!(reader.next().unwrap().unwrap().id_bytes(), rec3.id_bytes());

        reader.seek(&pos3).unwrap();
        assert_eq!(reader.next().unwrap().unwrap().id_bytes(), rec3.id_bytes());
        assert!(reader.next().is_none());

        reader.seek(&pos3).unwrap();
        assert_eq!(reader.next().unwrap().unwrap().id_bytes(), rec3.id_bytes());
        assert!(reader.next().is_none());
    }
}

#[test]
fn test_fastq_seek_err() {
    for cap in 3..20 {
        let mut reader =
            Reader::with_capacity(io::Cursor::new(&b"@s1\nA\n+\n~\n@s2\nA\n~\n"[..]), cap);

        let pos0 = reader.position().clone();
        reader.next().unwrap().unwrap();
        assert_matches!(
            reader.next().unwrap(),
            Err(Error::InvalidSep {
                found: b'~',
                pos:
                    ErrorPosition {
                        line: 7,
                        id: Some(_),
                    },
            })
        );

        let err_pos = reader.position().clone();
        assert!(reader.next().is_none());

        reader.seek(&pos0).unwrap();
        assert_matches!(reader.next(), Some(Ok(_)));

        reader.seek(&err_pos).unwrap();
        assert_matches!(
            reader.next().unwrap(),
            Err(Error::InvalidSep {
                found: b'~',
                pos:
                    ErrorPosition {
                        line: 7,
                        id: Some(_),
                    },
            })
        );
    }
}

#[test]
fn test_fastq_seek_none() {
    for cap in 3..20 {
        let mut reader = Reader::with_capacity(io::Cursor::new(&b"@s1\nA\n+\n~\n"[..]), cap);

        let pos0 = reader.position().clone();
        reader.next().unwrap().unwrap();
        assert!(reader.next().is_none());

        let end_pos = reader.position().clone();
        reader.seek(&pos0).unwrap();
        reader.seek(&end_pos).unwrap();
        assert!(reader.next().is_none());
    }
}

// FASTQ writing

#[test]
fn test_fastq_write_record() {
    let mut out = vec![];
    let mut rdr = Reader::new(FASTQ);
    while let Some(Ok(r)) = rdr.next() {
        r.write(&mut out).unwrap();
    }
    assert_eq!(&out, &FASTQ);
}

#[test]
fn test_fastq_write_record_unchanged() {
    let mut out = vec![];
    let mut rdr = Reader::new(FASTQ);
    while let Some(Ok(r)) = rdr.next() {
        r.write_unchanged(&mut out).unwrap();
    }
    assert_eq!(&out, &FASTQ);
}

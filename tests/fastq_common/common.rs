#[macro_export]
macro_rules! impl_common_fastq_tests {
    ($input:expr, $expected:expr, $Reader:ident, $PositionStore:ty, $RecordSet:ty, $ErrorKind:ident, $validate_ref:path,
        $next_fn:ident, $read_record_set_fn:ident, $seek_fn:ident, $parallel_reader_func:path) => {


use seq_io::Position as _Position;

macro_rules! test_reader {
    ($fastq:expr, $reader:ident, $block:block) => {
        for cap in 5..80 {
            if let Err(_) = std::panic::catch_unwind(|| {
                #[allow(unused_mut)]
                {
                    let mut $reader: $Reader<_, _, $PositionStore> = $Reader::with_capacity($fastq, cap).set_store();
                    $block
                }
            }) {
                panic!("Reader failed at capacity {}", cap);
            }
        }
    }
}

impl_common_tests!($input, $expected, RecordSet, test_reader, $validate_ref,
                   $next_fn, $read_record_set_fn, $seek_fn, $parallel_reader_func);


#[test]
#[should_panic(expected = "capacity smaller than")]
fn policy() {
    let fq = &b"@id\nATGC\n+\nIIII\n"[..];
    for cap in 5..80 {
        let policy = seq_io::policy::DoubleUntilLimited::new(2, 5);
        let mut reader: $Reader<_, _, $PositionStore> = $Reader::with_capacity(fq, cap).set_store()
            .set_policy(policy);
        let res = reader.next().unwrap();
        let err = res.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::BufferLimit);
        assert!(err.position().is_none());
    }
}


#[test]
#[allow(unused_variables)] // TODO: remove?
fn truncated1() {
    let fq = &b"\n@id"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnexpectedEnd { pos: _ });
        let pos = err.position().unwrap();
        let exp_pos = _Position::new().set_line(1).set_byte(1).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        assert_eq!(pos.error_offset().unwrap().line(), 0);
        assert_eq!(pos.error_offset().unwrap().byte(), 2);
        let exp_pos = _Position::new().set_line(1).set_byte(3).set_record(0).clone();
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.record_id().is_none());
    });
}

#[test]
#[allow(unused_variables)]
fn truncated2() {
    let fq = &b"\n@id\nATGC"[..];
    test_reader!(fq, reader, {
        let res = reader.next().unwrap();
        let err = res.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnexpectedEnd { pos: _ });
        let pos = err.position().unwrap();
        let exp_pos = _Position::new().set_line(1).set_byte(1).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        assert_eq!(pos.error_offset().unwrap().line(), 1);
        assert_eq!(pos.error_offset().unwrap().byte(), 7);
        let exp_pos = _Position::new().set_line(2).set_byte(8).set_record(0).clone();
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.record_id().unwrap() == "id");
    });
}

#[test]
#[allow(unused_variables)]
fn truncated3() {
    let fq = &b"\n\n@id\nATGC\n+"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnexpectedEnd { pos: _ });
        let pos = err.position().unwrap();
        let exp_pos = _Position::new().set_line(2).set_byte(2).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        assert_eq!(pos.error_offset().unwrap().line(), 2);
        assert_eq!(pos.error_offset().unwrap().byte(), 9);
        let exp_pos = _Position::new().set_line(4).set_byte(11).set_record(0).clone();
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.record_id().unwrap() == "id");
    });
}

#[test]
#[allow(unused_variables)]
fn truncated4() {
    let fq = &b"\n\n@id\nATGC\n+\n"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnexpectedEnd { pos: _ });
        let pos = err.position().unwrap();
        let exp_pos = _Position::new().set_line(2).set_byte(2).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        assert_eq!(pos.error_offset().unwrap().line(), 2);
        assert_eq!(pos.error_offset().unwrap().byte(), 10);
        let exp_pos = _Position::new().set_line(4).set_byte(12).set_record(0).clone();
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.record_id().unwrap() == "id");
    });
}

#[test]
#[allow(unused_variables)]
fn unequal() {
    let fq = &b"@id\nATGC\n+\nIII"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnequalLengths { pos: _, seq: 4, qual: 3 });
        let pos = err.position().unwrap();
        let exp_pos = _Position::new().set_line(0).set_byte(0).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.error_offset().is_none());
        assert!(pos.record_id().expect("should have an ID") == "id");
    });
}


#[test]
fn none_after_err() {
    let fq = &b"@id\nATGC"[..];
    test_reader!(fq, reader, {
        assert!(reader.next().unwrap().is_err());
        assert!(reader.next().is_none());
    });
}

#[test]
fn no_newline_end() {
    let fq = &b"@id\nATGC\n+\nIIII"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id_bytes(), b"id");
        assert_eq!(rec.seq(), b"ATGC");
        assert_eq!(rec.opt_qual(), Some(&b"IIII"[..]));
        assert!(reader.next().is_none());
    });
}


// TODO: check b"@id\nSS\n+\nQ\nQ"

#[test]
fn seek_err() {
    let fq = &b"@s1\nA\n+\nI\n@s2\nA\n"[..];
    test_reader!(std::io::Cursor::new(fq), reader, {
        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.id(), Ok("s1"));
        let pos1 = reader.position().clone();

        // advance to errorneous record
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnexpectedEnd { pos: _ });
        let error_pos = err.position().unwrap().to_owned();
        let exp_pos = _Position::new().set_line(4).set_byte(10).set_record(1).clone();
        validate_position!(error_pos.record_position().unwrap(), exp_pos);
        assert_eq!(error_pos.error_offset().unwrap().line(), 1);
        assert_eq!(error_pos.error_offset().unwrap().byte(), 5);
        let exp_pos = _Position::new().set_line(5).set_byte(15).set_record(1).clone();
        validate_position!(error_pos.position().unwrap(), exp_pos);
        assert_eq!(error_pos.record_id(), Some("s2"));

        // no next record after error
        assert!(reader.next().is_none());
        // current position should be the same as the one reported in error
        assert_eq!(&reader.position(), error_pos.record_position().unwrap());
        // seek back to first record
        reader.seek(&pos1).unwrap();
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id(), Ok("s1"));
        // seek to error position -> same error should happen again
        reader.seek(error_pos.record_position().unwrap()).unwrap();
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::UnexpectedEnd { pos: _ });
        assert_eq!(err.position().unwrap(), &error_pos);
    });
}

}
}

#[macro_export]
macro_rules! impl_common_fasta_tests {
    ($input:expr, $expected:expr, $ReaderBuilder:ident, $PositionStore:path, $RecordSet:ty, $ErrorKind:ident,
        $next_fn:ident, $read_record_set_fn:ident, $seek_fn:ident, $parallel_reader_func:path) => {

use seq_io::Position as _Position;
use seq_io::fasta;

macro_rules! test_reader {
    ($fastq:expr, $reader:ident, $block:block) => {
        for cap in 3..200 {
            if let Err(_) = std::panic::catch_unwind(|| {
                #[allow(unused_mut)]
                {
                    let mut $reader = make_reader!($ReaderBuilder, $fastq, $PositionStore, cap);
                    $block
                }
            }) {
                panic!("Reader failed at capacity {}", cap);
            }
        }
    }
}

macro_rules! validate_ref {
    ($record:expr, $exp:expr) => {
        let seq_lines: Vec<_> = $record.seq_lines().collect();
        assert_eq!(seq_lines, $exp.seq_lines, "sequence line mismatch");
    }
}


impl_common_tests!($input, $expected, RecordSet, test_reader, validate_ref, $next_fn, $read_record_set_fn, $seek_fn, $parallel_reader_func);


#[test]
fn invalid_start() {
    let mut reader = make_reader!($ReaderBuilder, &b"\r\nid\nATGC\n"[..], $PositionStore);
    let rec = reader.next().unwrap();
    let err = rec.err().expect("Should be an error");
    assert_matches!(err.kind(), $ErrorKind::InvalidStart { pos: _, found: b'i' });
    let pos = err.position().unwrap();
    let exp_pos = _Position::new().set_line(1).set_byte(2).set_record(0).clone();
    validate_position!(pos.record_position().unwrap(), exp_pos);
    assert!(pos.error_offset().is_none());
    validate_position!(pos.position().unwrap(), exp_pos);
    assert!(pos.record_id().is_none());
}


#[test]
fn policy() {
    let p = seq_io::policy::DoubleUntilLimited::new(2, 5);
    let mut reader = make_reader!($ReaderBuilder, &b">id\nAT\nGC\n"[..], $PositionStore, 3, p);
    let err = reader.next().unwrap().unwrap_err();
    assert_matches!(err.kind(), $ErrorKind::BufferLimit);
}

#[test]
fn none_after_err() {
    let p = seq_io::policy::DoubleUntilLimited::new(2, 3);
    let mut reader = make_reader!($ReaderBuilder, &b">id\nATGC\n"[..], $PositionStore, 3, p);
    assert!(reader.next().unwrap().is_err());
    assert!(reader.next().is_none());
}

#[test]
fn empty_lines_end() {
    let mut reader = make_reader!($ReaderBuilder, &b">id\nATGC\n\n\n\n\n\n\n\n\n"[..], $PositionStore);
    assert_eq!(reader.next().unwrap().unwrap().id_bytes(), b"id");
    assert!(reader.next().is_none());
}

#[test]
fn no_newline_end() {
    let mut reader = make_reader!($ReaderBuilder, &b">id\nATGC"[..], $PositionStore);
    assert_eq!(reader.next().unwrap().unwrap().id_bytes(), b"id");
    assert!(reader.next().is_none());
}

#[test]
fn write_head() {
    let mut out = vec![];
    fasta::write_head(&mut out, b"id desc").unwrap();
    assert_eq!(&out, b">id desc\n");
}

#[test]
fn write_seq() {
    let mut out = vec![];
    fasta::write_seq(&mut out, b"ATGC").unwrap();
    assert_eq!(&out, b"ATGC\n");
}

#[test]
fn write_seq_wrap() {
    let mut out = vec![];
    fasta::write_wrap_seq(&mut out, b"ATGCA", 2).unwrap();
    assert_eq!(&out, b"AT\nGC\nA\n");
}

#[test]
fn write_seq_iter() {
    let mut out = vec![];
    fasta::write_seq_iter(&mut out, b"ATGCA".chunks(2)).unwrap();
    assert_eq!(&out, b"ATGCA\n");
}

#[test]
fn write_seq_iter_wrap() {
    for size in 1..11 {
        let mut out = vec![];
        fasta::write_wrap_seq_iter(&mut out, b"AAAATTTTGGG".chunks(size), 3).unwrap();
        assert_eq!(&out, b"AAA\nATT\nTTG\nGG\n");

        let mut out = vec![];
        fasta::write_wrap_seq_iter(&mut out, b"AAAATTTTGGG".chunks(size), 4).unwrap();
        assert_eq!(&out, b"AAAA\nTTTT\nGGG\n");
    }
}

}
}

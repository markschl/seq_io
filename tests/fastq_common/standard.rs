#[macro_export]
macro_rules! impl_fastq_standard_tests {
    ($Reader:ident, $PositionStore:ty, $RecordSet:ident, $ErrorKind:ident, $validate_ref:ident, $parallel_reader_func:path) => {

use seq_io::Position;
//use seq_io::parallel::parallel_fastq;
use crate::fastq_common::{FASTQ, FASTQ_EXPECTED};

impl_common_fastq_tests!(FASTQ, FASTQ_EXPECTED, $Reader, $PositionStore, $RecordSet, $ErrorKind, $validate_ref,
    next, read_record_set, seek, $parallel_reader_func);

#[test]
fn invalid_start() {
    let fq = &b"@id1\nA\n+\nI\nid\nATGC\n+\nIIII"[..];
    test_reader!(fq, reader, {
        reader.next().unwrap().unwrap();
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), $ErrorKind::InvalidStart { pos: _, found: b'i' });
        let pos = err.position().unwrap();
        validate_position!(reader.position(), Position::new().set_line(4).set_byte(11).set_record(1));
        assert!(pos.error_offset().is_none());
        assert!(pos.record_id().is_none());
    });
}


#[test]
fn no_sep() {
    let fq = &b"\n@id\nATGC\nIIII\n"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), ErrorKind::InvalidSep { pos: _, found: Some(b'I') });
        let pos = err.position().unwrap();
        let exp_pos = Position::new().set_line(1).set_byte(1).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        assert_eq!(pos.error_offset().unwrap().line(), 2);
        assert_eq!(pos.error_offset().unwrap().byte(), 9);
        let exp_pos = Position::new().set_line(3).set_byte(10).set_record(0).clone();
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.record_id().unwrap() == "id");
    });
}


// same as single-line FASTQ version, but with two newlines at end
#[test]
fn write_record() {
    let fastq_in = &b"\n\n@id\nSEQQ\n+ id\r\nQUAL\n\n@id2\r\nSEQ\r\n+\nQUA"[..];
    // lines joined and CR removed
    let fastq_out = &b"@id\nSEQQ\n+\nQUAL\n@id2\nSEQ\n+\nQUA\n"[..];
    // not changed apart from newlines before and after record
    let fastq_out_unchanged = &b"@id\nSEQQ\n+ id\r\nQUAL\n@id2\r\nSEQ\r\n+\nQUA\n"[..];
    test_reader!(fastq_in, reader, {
        let mut out = vec![];
        let mut out_unchanged = vec![];
        while let Some(res) = reader.next() {
            let rec = res.unwrap();
            rec.write(&mut out).unwrap();
            rec.write_unchanged(&mut out_unchanged).unwrap();
        }
        assert_eq!(out.as_slice(), fastq_out);
        assert_eq!(out_unchanged.as_slice(), fastq_out_unchanged);
    });
}

}
}

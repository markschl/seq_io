#[macro_export]
macro_rules! impl_fastq_multi_tests {
    ($Reader:ident, $PositionStore:path, $RecordSet:ty, $ErrorKind:ident, $validate_ref:path, $parallel_reader_func:path) => {

use seq_io::Position;
//use seq_io::parallel::parallel_fastq;
use crate::fastq_common::{FASTQ_MULTILINE, FASTQ_MULTI_EXPECTED};

impl_common_fastq_tests!(FASTQ_MULTILINE, FASTQ_MULTI_EXPECTED, $Reader, $PositionStore, $RecordSet, $ErrorKind, $validate_ref,
                        next, read_record_set, seek, $parallel_reader_func);

// TODO: explain
#[test]
#[allow(unused_variables)]
fn no_sep() {
    let fq = &b"\n@id\nATGC\nIII\nI"[..];
    test_reader!(fq, reader, {
        let rec = reader.next().unwrap();
        let err = rec.err().expect("Should be an error");
        assert_matches!(err.kind(), ErrorKind::UnexpectedEnd { pos: _ });
        let pos = err.position().unwrap();
        let exp_pos = Position::new().set_line(1).set_byte(1).set_record(0).clone();
        validate_position!(pos.record_position().unwrap(), exp_pos);
        assert_eq!(pos.error_offset().unwrap().line(), 3);
        assert_eq!(pos.error_offset().unwrap().byte(), 13);
        let exp_pos = Position::new().set_line(4).set_byte(14).set_record(0).clone();
        validate_position!(pos.position().unwrap(), exp_pos);
        assert!(pos.record_id().unwrap() == "id");
    });
}

// same as single-line FASTQ version, but with two newlines at end
#[test]
fn write_record() {
    let fastq_in = &b"\n\n@id\nSEQQ\n+ id\r\nQUAL\n\n@id2\r\nS\nE\nQ\r\n+\nQU\nA"[..];
    // lines joined and CR removed
    let fastq_out = &b"@id\nSEQQ\n+\nQUAL\n@id2\nSEQ\n+\nQUA\n"[..];
    // not changed apart from newlines before and after record
    let fastq_out_unchanged = &b"@id\nSEQQ\n+ id\r\nQUAL\n\n@id2\r\nS\nE\nQ\r\n+\nQU\nA\n"[..];
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

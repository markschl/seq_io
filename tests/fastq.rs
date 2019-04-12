#[macro_use]
extern crate matches;
#[macro_use]
extern crate lazy_static;

#[macro_use]
mod common;
#[macro_use]
mod fastq_common;
mod seq;

macro_rules! validate_ref {
    ($record:expr, $exp:expr) => {
        validate_record!($record, $exp);
        let seq_lines: Vec<_> = $record.seq_lines().collect();
        assert_eq!(seq_lines, $exp.seq_lines, "sequence lines mismatch");
        let qual_lines: Vec<_> = $record.qual_lines().collect();
        assert_eq!(Some(qual_lines), $exp.qual_lines, "quality lines mismatch");
    };
}

mod standard {
    use seq_io::fastq::{ErrorKind, RangeStore, Reader, RecordSet};
    use seq_io::prelude::*;
    impl_fastq_standard_tests!(
        Reader,
        RangeStore,
        RecordSet,
        ErrorKind,
        validate_ref,
        seq_io::parallel::read_process_fastq_records
    );
}

mod multi_line {
    use seq_io::fastq::multiline::{MultiRangeStore, Reader};
    use seq_io::fastq::{ErrorKind, RecordSet};
    use seq_io::prelude::*;
    impl_fastq_multi_tests!(
        Reader,
        MultiRangeStore,
        RecordSet,
        ErrorKind,
        validate_ref,
        seq_io::parallel::read_process_fastq_records
    );
}

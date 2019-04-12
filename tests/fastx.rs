#[macro_use]
extern crate matches;
#[macro_use]
extern crate lazy_static;

//use seq_io::parallel::parallel_fastq;

#[macro_use]
mod common;
#[macro_use]
mod fasta_common;
#[macro_use]
mod fastq_common;
mod seq;

mod fasta {
    use seq_io::fastx::{ErrorKind, Reader, RecordSet};
    use seq_io::prelude::*;
    impl_fasta_standard_tests!(
        Reader,
        seq_io::fasta::LineStore,
        RecordSet,
        ErrorKind,
        seq_io::parallel::read_process_fastx_records
    );
}

mod fasta_multi_qual {
    use seq_io::fastx::multiline_qual::Reader;
    use seq_io::fastx::{ErrorKind, LineStore, RecordSet};
    use seq_io::prelude::*;
    impl_fasta_standard_tests!(
        Reader,
        LineStore,
        RecordSet,
        ErrorKind,
        seq_io::parallel::read_process_fastx_records
    );
}

macro_rules! validate_fastx_ref {
    ($record:expr, $exp:expr) => {
        validate_record!($record, $exp);
        let seq_lines: Vec<_> = $record.seq_lines().collect();
        assert_eq!(seq_lines, $exp.seq_lines, "sequence lines mismatch");
        let qual_lines: Vec<_> = $record.opt_qual_lines().unwrap().collect();
        assert_eq!(Some(qual_lines), $exp.qual_lines, "quality lines mismatch");
    };
}

mod fastq {
    use seq_io::fastx::{ErrorKind, Reader, RecordSet};
    use seq_io::prelude::*;
    impl_fastq_standard_tests!(
        Reader,
        seq_io::fastq::RangeStore,
        RecordSet,
        ErrorKind,
        validate_fastx_ref,
        seq_io::parallel::read_process_fastx_records
    );
}

mod fastq_dyn {
    use seq_io::fastx::{ErrorKind, Reader, RecordSet};
    use seq_io::prelude::*;
    impl_fastq_standard_tests!(
        Reader,
        seq_io::fastq::RangeStore,
        RecordSet,
        ErrorKind,
        validate_fastx_ref,
        seq_io::parallel::read_process_fastx_records
    );
}

mod fastq_linestore {
    use seq_io::fastx::{ErrorKind, Reader, RecordSet};
    use seq_io::prelude::*;
    impl_fastq_standard_tests!(
        Reader,
        seq_io::fastx::LineStore,
        RecordSet,
        ErrorKind,
        validate_fastx_ref,
        seq_io::parallel::read_process_fastx_records
    );
}

mod fastq_multiline {
    use seq_io::fastx::multiline_qual::Reader;
    use seq_io::fastx::{ErrorKind, LineStore, RecordSet};
    use seq_io::prelude::*;
    impl_fastq_multi_tests!(
        Reader,
        LineStore,
        RecordSet,
        ErrorKind,
        validate_fastx_ref,
        seq_io::parallel::read_process_fastx_records
    );
}

mod fastx_dynamic_fastq {
    use seq_io::fastx::dynamic::ReaderBuilder;
    use seq_io::fastx::SeqFormat;
    use seq_io::prelude::*;

    macro_rules! test_dyn_fastx_reader {
        ($content:expr, $reader:ident, $block:block) => {
            for cap in 5..80 {
                if let Err(_) = std::panic::catch_unwind(|| {
                    let mut $reader = ReaderBuilder::new()
                        .set_capacity(cap)
                        .set_store::<seq_io::fastx::LineStore>()
                        .from_reader($content)
                        .unwrap()
                        .expect("Input empty");
                    $block
                }) {
                    panic!("Reader failed at capacity {}", cap);
                }
            }
        };
    }

    #[test]
    fn reader() {
        use crate::fastq_common::{FASTQ, FASTQ_EXPECTED};
        test_dyn_fastx_reader!(FASTQ, reader, {
            let mut exp_iter = FASTQ_EXPECTED.iter();
            while let Some(exp) = exp_iter.next() {
                let record = reader.next_fastx().unwrap().unwrap();
                validate_fastx_ref!(record, exp);
                validate_record!(record.to_owned_record(), exp);
                dbg!(reader.position_fastx());
                dbg!(&exp.pos);
                validate_position!(reader.position_fastx(), exp.pos);
            }
            assert!(reader.next_fastx().is_none());
        });
    }

    #[test]
    fn fasta_recognition() {
        let fasta = b"\r\n\n\n\n\n\n>id\n";
        test_dyn_fastx_reader!(&fasta[..], reader, {
            assert_eq!(reader.fastx_records().count(), 1);
            assert_eq!(reader.format(), Some(SeqFormat::FASTA));
        });
    }

    #[test]
    fn fastq_recognition() {
        let fastq = b"\n\n\r\n\n\n\n@id\n\n+\n\n";
        test_dyn_fastx_reader!(&fastq[..], reader, {
            assert_eq!(reader.fastx_records().count(), 1);
            assert_eq!(reader.format(), Some(SeqFormat::FASTQ));
        });
    }

    #[test]
    fn empty_recognition() {
        let empty = b"\n\n\r\n\n\n\n\n\n\n\n";
        let reader = ReaderBuilder::new()
            .set_capacity(3)
            .set_store::<seq_io::fastx::LineStore>()
            .from_reader(&empty[..])
            .unwrap();
        assert!(reader.is_none());
    }
}

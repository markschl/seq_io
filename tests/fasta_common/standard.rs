#[macro_export]
macro_rules! impl_fasta_standard_tests {
    ($ReaderBuilder:ident, $PositionStore:path, $RecordSet:ty, $ErrorKind:ident, $parallel_reader_func:path) => {
        //use seq_io::parallel::parallel_fastq;
        use crate::fasta_common::{FASTA, FASTA_EXPECTED};

        impl_common_fasta_tests!(
            FASTA,
            FASTA_EXPECTED,
            $ReaderBuilder,
            $PositionStore,
            $RecordSet,
            $ErrorKind,
            next,
            read_record_set,
            seek,
            $parallel_reader_func
        );

        #[test]
        fn no_seq() {
            let mut reader = make_reader!($ReaderBuilder, &b">id1\n>id2\n"[..], $PositionStore);
            {
                let r = reader.next().unwrap().unwrap();
                assert_eq!(r.id_bytes(), b"id1");
                assert!(r.seq().is_empty());
                assert!(r.num_seq_lines() == 0);
                assert!(r.seq_lines().count() == 0);
                assert!(r.seq_lines().rev().count() == 0);
            }
            {
                let r = reader.next().unwrap().unwrap();
                assert_eq!(r.id_bytes(), b"id2");
                assert!(r.seq().is_empty());
                assert!(r.num_seq_lines() == 0);
                assert!(r.seq_lines().count() == 0);
                assert!(r.seq_lines().rev().count() == 0);
            }
        }

        #[test]
        fn test_write_fasta() {
            let fasta_in = &b"\n\n>id desc\r\nACGT\nACGT\r\nACGT\n\r\n\n\r\n"[..];
            let fasta_out = &b">id desc\nACG\nTAC\nGTA\nCGT\n"[..];
            test_reader!(fasta_in, reader, {
                let mut out = vec![];
                while let Some(res) = reader.next() {
                    res.unwrap().write_wrap(&mut out, 3).unwrap();
                }
                assert_eq!(out.as_slice(), fasta_out);
            });
        }
    };
}

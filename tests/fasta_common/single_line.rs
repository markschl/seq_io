#[macro_export]
macro_rules! impl_fasta_single_tests {
    ($Reader:ident, $PositionStore:path, $RecordSet:ty, $ErrorKind:ident, $parallel_reader_func:path) => {
        //use seq_io::parallel::parallel_fastq;
        use crate::fasta_common::{FASTA_SINGLE, FASTA_SINGLE_EXPECTED};

        impl_common_fasta_tests!(
            FASTA_SINGLE,
            FASTA_SINGLE_EXPECTED,
            $Reader,
            $PositionStore,
            $RecordSet,
            $ErrorKind,
            next,
            read_record_set,
            seek,
            $parallel_reader_func
        );

        #[test]
        fn test_write_fasta() {
            let fasta_in = &b"\n\n>id desc\r\nACGTACGTACGT\n\r\n\n\r\n"[..];
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

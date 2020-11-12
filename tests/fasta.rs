#[macro_use]
extern crate matches;
#[macro_use]
extern crate lazy_static;

#[macro_use]
mod common;
mod seq;
#[macro_use]
mod fasta_common;

mod standard {
    use seq_io::fasta::{ErrorKind, LineStore, ReaderBuilder, RecordSet};
    use seq_io::prelude::*;
    impl_fasta_standard_tests!(
        ReaderBuilder,
        LineStore,
        RecordSet,
        ErrorKind,
        seq_io::parallel::read_process_fasta_records
    );
}

mod single_line {
    use seq_io::fasta::single_line::ReaderBuilder;
    use seq_io::fasta::{ErrorKind, LineStore, RecordSet};
    use seq_io::prelude::*;
    impl_fasta_single_tests!(
        ReaderBuilder,
        LineStore,
        RecordSet,
        ErrorKind,
        seq_io::parallel::read_process_fasta_records
    );
}

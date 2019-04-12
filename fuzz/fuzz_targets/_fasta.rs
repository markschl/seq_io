use seq_io::{fasta, fastq, fastx};
use simple_reader::{compare_simple, compare_readers, compare_recset};

mod simple_reader;

pub fn evaluate(data: &[u8]) {
    // FASTA
    // normal
    let reader = fasta::Reader::with_capacity(data, 3);
    let simple_rdr = simple_reader::Reader::new_fasta(data, true);
    compare_simple(reader, simple_rdr, false);
    compare_readers(
        fasta::Reader::with_capacity(data, 3).set_store::<fasta::LineStore>(),
        fasta::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
    );
    compare_readers(
        // TODO: line numbers not correct
        fasta::Reader::with_capacity(data, 3).set_store::<fasta::LineStore>(),
        fasta::Reader::with_capacity(data, 3).set_store::<fastq::RangeStore>(),
    );
    compare_readers(
        fasta::Reader::with_capacity(data, 3).set_store::<fasta::LineStore>(),
        fasta::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
    );
    compare_recset(
        fasta::Reader::with_capacity(data, 3), 
        fasta::Reader::with_capacity(data, 3)
    );

    // single-line
    let reader = fasta::single_line::Reader::with_capacity(data, 3);
    let simple_rdr = simple_reader::Reader::new_fasta(data, false);
    compare_simple(reader, simple_rdr, false);
    compare_readers(
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fasta::single_line::RangeStore>(),
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fasta::LineStore>(),
    );
    compare_readers(
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fasta::single_line::RangeStore>(),
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
    );
    compare_readers(
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fasta::single_line::RangeStore>(),
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fastq::RangeStore>(),
    );
    compare_readers(
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fasta::single_line::RangeStore>(),
        fasta::single_line::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
    );
    compare_recset(
        fasta::single_line::Reader::with_capacity(data, 3), 
        fasta::single_line::Reader::with_capacity(data, 3)
    );
    
    // FASTX <-> FASTA
    // normal
    let reader = fasta::Reader::with_capacity(data, 3);
    let simple_rdr = simple_reader::Reader::new_fastx(data, true, false);
    compare_simple(reader, simple_rdr, true);

    // FASTX <-> FASTX
    // normal
    let reader = fastx::Reader::with_capacity(data, 3);
    let simple_rdr = simple_reader::Reader::new_fastx(data, true, false);
    compare_simple(reader, simple_rdr, false);
    compare_readers(
        fastx::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
        fastx::Reader::with_capacity(data, 3).set_store::<fastq::RangeStore>(),
    );
    compare_readers(
        fastx::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
        fastx::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
    );
    compare_recset(
        fastx::Reader::with_capacity(data, 3), 
        fastx::Reader::with_capacity(data, 3)
    );

    // multi-line quality
    let reader = fastx::multiline_qual::Reader::with_capacity(data, 3);
    let simple_rdr = simple_reader::Reader::new_fastx(data, true, true);
    compare_simple(reader, simple_rdr, false);
    compare_readers(
        fastx::multiline_qual::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
        fastx::multiline_qual::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
    );
    compare_recset(
        fastx::multiline_qual::Reader::with_capacity(data, 3), 
        fastx::multiline_qual::Reader::with_capacity(data, 3)
    );

}

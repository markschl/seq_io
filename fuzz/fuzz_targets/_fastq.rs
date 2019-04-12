use seq_io::{fastq, fastx};
use simple_reader::{compare_simple, compare_readers, compare_recset};

mod simple_reader;

pub fn evaluate(data: &[u8]) {
    // FASTQ
    // normal
    let reader = fastq::Reader::with_capacity(data, 3);
    let simple_reader = simple_reader::Reader::new_fastq(data, false);
    compare_simple(reader, simple_reader, false);
    compare_readers(
        fastq::Reader::with_capacity(data, 3).set_store::<fastq::RangeStore>(),
        fastq::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
    );
    compare_readers(
        fastq::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
        fastq::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
    );
    compare_recset(
        fastq::Reader::with_capacity(data, 3), 
        fastq::Reader::with_capacity(data, 3)
    );

    // multi-line
    let reader = fastq::multiline::Reader::with_capacity(data, 3);
    let simple_reader = simple_reader::Reader::new_fastq(data, true);
    compare_simple(reader, simple_reader, false);
    compare_readers(
        // num_lines not correct
        fastq::Reader::with_capacity(data, 3).set_store::<fastq::RangeStore>(),
        fastq::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
    );
    compare_readers(
        fastq::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
        fastq::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
    );
    compare_recset(
        fastq::multiline::Reader::with_capacity(data, 3), 
        fastq::multiline::Reader::with_capacity(data, 3)
    );

    // FASTQ <-> FASTX
    // normal
    let reader = fastq::Reader::with_capacity(data, 3);
    let simple_reader = simple_reader::Reader::new_fastx(data, true, false);
    compare_simple(reader, simple_reader, true);

    // multi-line qualities
    let reader = fastq::multiline::Reader::with_capacity(data, 3);
    let simple_reader = simple_reader::Reader::new_fastx(data, true, true);
    compare_simple(reader, simple_reader, true);

    // FASTX <-> FASTX
    // normal
    let reader = fastx::Reader::with_capacity(data, 3);
    let simple_reader = simple_reader::Reader::new_fastx(data, true, false);
    compare_simple(reader, simple_reader, false);
    compare_readers(
        // TODO: num_lines...
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

    // multi-line qualities
    let reader = fastx::multiline_qual::Reader::with_capacity(data, 3);
    let simple_reader = simple_reader::Reader::new_fastx(data, true, true);
    compare_simple(reader, simple_reader, false);
    compare_readers(
        fastx::multiline_qual::Reader::with_capacity(data, 3).set_store::<fastx::LineStore>(),
        fastx::multiline_qual::Reader::with_capacity(data, 3).set_store::<fastq::multiline::MultiRangeStore>(),
    );
    // TODO: could fastq::RangeStore work?
    compare_recset(
        fastx::multiline_qual::Reader::with_capacity(data, 3), 
        fastx::multiline_qual::Reader::with_capacity(data, 3)
    );
}


# FASTA and FASTQ parsing and writing in Rust.

![docs.rs](https://docs.rs/seq_io/badge.svg)
[![crates.io](https://img.shields.io/crates/v/seq_io.svg)](https://crates.io/crates/seq_io)
[![Build status](https://github.com/markschl/seq_io/workflows/ci/badge.svg)](https://github.com/markschl/seq_io/actions)

This library provides an(other) attempt at parsing of the sequence formats FASTA and FASTQ, as well as writing.

**Features:**

* Fast readers that minimize the use of allocations and copying of memory
* Flexible methods for writing FASTA and FASTQ
* Informative errors
* Support for seeking
* Serde support (for owned data structures)
* Functions for parallel processing
* Tested using fuzzing techniques [see here](fuzz/README.md)

The FASTA parser can read and write multi-line files and allows
iterating over the sequence lines without doing any allocation or
copying. The FASTQ parser does not support multiple sequence / quality lines.

### Documentation

[**Documentation for the stable version (0.3.x)**](https://docs.rs/seq_io)

The [v0.4 branch](https://github.com/markschl/seq_io/tree/v0.4) contains code for a new
version, which includes a FASTX reader. Although it works and has been tested to some
extent, there will be further large changes, which are not quite ready yet.

[**Documentation for development version (0.4.0-alpha.x)**](https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/index.html)

### Example
Reads FASTA sequences from STDIN and writes them to STDOUT
if long enough. Otherwise it prints a message. This should
be very fast because the sequence is not allocated (`seq_lines()`).
```rust
use seq_io::fasta::{Reader,Record};
use std::io;

let mut reader = Reader::new(io::stdin());
let mut stdout = io::stdout();

while let Some(result) = reader.next() {
    let record = result.unwrap();
    // determine sequence length
    let seqlen = record.seq_lines()
                       .fold(0, |l, seq| l + seq.len());
    if seqlen > 100 {
        record.write_wrap(&mut stdout, 80).unwrap();
    } else {
        eprintln!("{} is only {} long", record.id().unwrap(), seqlen);
    }
}
```

Records are directly borrowing data from the internal buffered reader,
therefore the `while let` is required. By default, the buffer will automatically
grow if a record is too large to fit in. How it grows can be configured, it is
also possible to set a size limit. Iterators over owned records are also provided.

**Note:** [LTO](https://doc.rust-lang.org/cargo/reference/profiles.html#lto)
might be evaluated to see, whether it improves performance. But generally, library should
work reasonably fast without LTO, too.

### Multi-threaded processing

The `parallel` module contains functions for sending FASTQ/FASTA
records to a thread pool where expensive calculations are done.
Sequences are processed in batches (`RecordSet`) because sending across
channels has a performance impact. FASTA/FASTQ records can be accessed in
both the 'worker' function and (after processing) a function running in the
main thread.

### Similar projects in Rust

* *[Rust-Bio](https://rust-bio.github.io)*: Binformatics library that provides
  simple FASTA and FASTQ readers.
* *[fastq-rs](https://github.com/aseyboldt/fastq-rs)*: FASTQ parser with
  comparable performance (see below). `seq_io` was inspired by `fastq_rs`.
* *[Needletail](https://github.com/onecodex/needletail)* has a FASTA parser.
* *[fasten](https://github.com/lskatz/fasten)* (FASTQ parser)

### Performance comparisons

The FASTQ reader from this crate performs similar to the
[fastq-rs](https://github.com/aseyboldt/fastq-rs) reader.
The [rust-bio](http://rust-bio.github.io/) readers are slower due
to allocations, copying, and UTF-8 validity checks.

All comparisons were run on a set of 100,000 auto-generated, synthetic sequences
with lengths normally distributed around 500 bp and loaded into memory.
The parsers from this crate (*seq_io*) are compared with [fastq-rs](https://github.com/aseyboldt/fastq-rs) (*fastq_rs*)
and [Rust-Bio](https://rust-bio.github.io/) (*bio*).
The bars represent the throughput in GB/s (+/- standard error of the mean).
Run on a Thinkpad X1 Carbon (i7-5500U) with a fixed frequency of 2.3 GHz
using Rust 1.31 nightly

![benchmark results](bench_results/reader_comparison_simple.png)

**Explanation of labels**:

* *Top bars*: Iteration over all records without further action.
* *owned*: An owned copy of each record is created for comparison with *Rust-Bio*,
  which does not provide zero copy parsing.
* *multiline*: The FASTA sequence is split into 5 x 100 bp lines.
* *recordset*: Records are parsed into record sets using `read_record_set()` (involves some copying).
* *parallel*: Record sets are are sent to worker threads for parallel processing
  where they are being iterated over and then sent back to the main thread
  where there is another iteration over the records (the latter only in seq_io)

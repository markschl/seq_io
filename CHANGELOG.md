# seq_io change log

## v0.3.4 (Mar 14, 2025)

This release includes a new method called `read_record_set_exact()`, which allows
specifying the number of records that should be present in a record set.
This feature is required for parallel processing of paired-end reads
([see also here](https://github.com/markschl/seq_io/issues/24)).
Thanks @nh13 for the [contribution](https://github.com/markschl/seq_io/pull/23)
to this feature.

More changes:

* Added `buf_capacity` and `shrink_buffer_to_fit` to `fasta::RecordSet` and
  `fastq::RecordSet`, allowing more control over the buffer size.


## v0.3.3 (Jan 28, 2025)

* Some internal refactoring with the aim to always correctly handle mixed calls
  to `next()`, `read_record_set()` and `seek()`.
* Added `len()` and `is_empty()` methods to `fasta::RecordSet` and `fastq::RecordSet`
* Added `from_path_with_capacity()` to the FASTA and FASTQ readers
* Thanks to [changes in buf_redux 1.0.1](https://github.com/dignifiedquire/buffer-redux/pull/2)
  (now enforced in Cargo.toml), LTO is not necessarily required anymore for full performance.

**Bug fixes**
* Fix issue where calling `read_record_set()` after `next()` lead to incorrect
  reading ([#20](https://github.com/markschl/seq_io/issues/20))

## v0.3.2 (Aug 07, 2023)

This version mostly updated all the dependencies and switched to using
Rust edition 2021. The minimum supported Rust version is still 1.56.1,
as it was before.

The unmaintained `buf_redux` was replaced by the modernized fork `buffer_redux`
(fixes #15).

In addition, a4c50de adds some `#[inline]` statements, which may improve performance
in non-LTO builds (but there [remains an issue](https://github.com/dignifiedquire/buffer-redux/pull/2)
still requiring LTO for maximum performance).

## v0.3.1 (Oct 17, 2020)

This version got many improvements and three bugfixes (related to FASTA parsing).
An update is recommended since certain invalid FASTA files were parsed without an error
before.
In addition, the FASTA parser is now less restrictive (7e713ab, 70d98da) and 
better tested thanks to fuzzing (3950cd8).

**Note**: in hindsight, the changes d6de345 and 61bda2b are potentially breaking
API changes, it would have been better to release 0.4.0 instead.

**Bug fixes** (finished Apr 11, 2019, but unfortunately no release made):
* f3c647a/ac8cc8c FASTA reader did not correctly validate the first
    record, resulting in malformed FASTA files passing without error.
* 07274b7 *bugfix*: Fix fasta::RefRecord::seq() panic with empty sequences.

Other important commits:
* bd99053 Updated dependencies (especially crossbeam, see #4)
* d6de345 Moved buffer policy types into separate module.
* Updated dependencies (fixing #4) and silenced warnings
* 61bda2b Improve writing API: take io::Write by value, not by mutable reference 
* 2b4a674 Implement ExactSizeIterator for seq_io::fasta::SeqLines 
* 187cd9e Add 'DoubleUntilLimited' buffer policy and rename default policy 

## v0.3.0 (Aug 21, 2018)

Most important change:

ae5f71c Renamed buffer growth "strategy" to "policy" module and refactored the API.

## v0.2.6 (Aug 14, 2018)

Most important change:

022872c Fixed edge case that could occur in parallel module.

## v0.2.5 (Aug 13, 2018)

* Added cargo-fuzz infrastructure (thanks @aseyboldt) and fixed a panic detected
  by fuzzing

## v0.2.4 (Apr 16, 2018)

* Updated dependencies and fixed warning

## v0.2.2 + v0.2.3 (Jan 13, 2018)

* Make proceed() private
* Added iterators of owned records and a few other methods: full_seq(),
  write_wrap_seq_iter(), full_seq() and num_seq_lines()
* 

## v0.2.1 (Dec 16, 2017)

* Remove Record impl from readers
* Add more flexible methods for parallel processing. 

## v0.2.0 (Nov 21, 2017)

* 77b6627 Added a simple test crate allowing to compare the parsed sequences with
    the parsers of the rust-bio crate and thus validating seq_io on selected
    FASTA/FASTQ files.
* Better errors and error messages
* added position() and seek() reader methods
* Improved `parallel` module.
* More options for FASTA writing
* Fixes FASTQ writer and fastq::Reader::qual()
* Correctly deal with missing newlines at EOF

## v0.1.0 (Jul 31, 2017)

Initial release

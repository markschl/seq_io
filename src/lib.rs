//! This library provides an(other) attempt at high performance FASTA and FASTQ parsing and writing.
//! The FASTA parser can read and write multi-line files. The FASTQ parser supports only single
//! lines. The sequence length of records in the FASTA/FASTQ files
//! is not limited by the size of the buffer. Instead, the buffer will grow until
//! the record fits, allowing parsers with a minimum amount of copying required.
//! How it grows can be configured (see [`BufStrategy`](trait.BufStrategy.html)).
//!
//! See also the documentation for the [FASTA Reader](fasta/struct.Reader.html) and the
//! [FASTQ Reader](fastq/struct.Reader.html). The methods for writing are documented
//! [here](fasta/index.html#functions) for FASTA and [here](fastq/index.html#functions)
//! for FASTQ.
//!
//! # Example FASTQ parser:
//! This code prints the ID string from each FASTQ record.
//!
//! ```no_run
//! use seq_io::fastq::{Reader,Record};
//!
//! let mut reader = Reader::from_path("seqs.fasta").unwrap();
//!
//! while let Some(record) = reader.next() {
//!     let record = record.expect("Error reading record");
//!     println!("{}", record.id().unwrap());
//! }
//! ```
//!
//! # Example FASTA parser calculating mean sequence length:
//! The FASTA reader works just the same. One challenge with the FASTA
//! format is that the sequence can be broken into multiple lines.
//! Therefore, it is not always possible to get a slice to the whole sequence
//! without copying the data. But it is possible to use `seq_lines()`
//! for efficiently iterating over each sequence line:
//!
//! ```no_run
//! use seq_io::fasta::{Reader,Record};
//!
//! let mut reader = Reader::from_path("seqs.fasta").unwrap();
//!
//! let mut n = 0;
//! let mut sum = 0;
//! while let Some(record) = reader.next() {
//!     let record = record.expect("Error reading record");
//!     for s in record.seq_lines() {
//!         sum += s.len();
//!     }
//!     n += 1;
//! }
//! println!("mean sequence length of {} records: {:.1} bp", n, sum as f32 / n as f32);
//! ```
//! If the whole sequence is required at once, there is the
//! [`full_seq`](fasta/struct.RefRecord.html#method.full_seq),
//! which will only allocate the sequence if there are multiple lines.
//! use seq_io::fasta::{Reader,OwnedRecord};
//!
//! # Owned records
//! Both readers also provide iterators similar to *Rust-Bio*, which return owned data. This
//! is slower, but make sense, e.g. if the records are collected in to a vector:
//!
//! ```no_run
//! use seq_io::fasta::Reader;
//!
//! let mut reader = Reader::from_path("input.fasta").unwrap();
//!
//! let records: Result<Vec<_>, _> = reader.records().collect();
//! ```
//!
//! # Parallel processing
//! Functions for parallel processing can be found in the [`parallel`](parallel/index.html) module

extern crate buf_redux;
extern crate memchr;

#[macro_use]
extern crate serde_derive;
extern crate serde;

use std::error;
use std::fmt;
use std::io;

pub use strategy::*;

mod strategy;

macro_rules! try_opt {
    ($expr: expr) => {
        match $expr {
            Ok(item) => item,
            Err(e) => return Some(Err(::std::convert::From::from(e))),
        }
    };
}

macro_rules! unwrap_or {
    ($expr:expr, $or:block) => {
        match $expr {
            Some(item) => item,
            None => $or,
        }
    };
}

pub mod fasta;
pub mod fastq;
pub mod parallel;

/// used by more than one module

#[derive(Default, Debug)]
struct ReadAlways;

impl buf_redux::strategy::ReadStrategy for ReadAlways {
    fn should_read(&self, _: &buf_redux::Buffer) -> bool {
        true
    }
}

/// Remove a final '\r' from a byte slice
#[inline]
fn trim_cr(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', remaining)) = line.split_last() {
        remaining
    } else {
        line
    }
}

/// Makes sure the buffer is full after this call (unless EOF reached)
/// code adapted from `io::Read::read_exact`
fn fill_buf<R, Rs, Ms>(reader: &mut buf_redux::BufReader<R, Rs, Ms>) -> io::Result<usize>
where
    R: io::Read,
    Rs: buf_redux::strategy::ReadStrategy,
    Ms: buf_redux::strategy::MoveStrategy,
{
    let mut num_read = reader.get_buf().len();
    while num_read < reader.capacity() {
        match reader.read_into_buf() {
            Ok(0) => break,
            Ok(n) => num_read += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(num_read)
}

//! This library provides an(other) attempt at high performance FASTA and FASTQ parsing.
//! There are many similarities to the excellent [fastq-rs](https://github.com/aseyboldt/fastq-rs).
//! However, the API that provides streaming iterators where possible.
//! Additionally, the sequence length of records in the FASTA/FASTQ files
//! is not limited by the size of the buffer. Instead, the buffer will grow until
//! the record fits, allowing parsers with a minimum amount of copying required.
//! How it grows can be configured (see [`BufGrowStrategy`](trait.BufGrowStrategy.html)).
//!
//! # Example FASTQ parser:
//! This code prints the ID string from each FASTQ record.
//!
//! ```no_run
//! use seq_io::fastq::{Reader,Record};
//!
//! let mut reader = Reader::from_path("seqs.fastq").unwrap();
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
//! Therefore, it is not possible to get a slice to the whole sequence
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
//!
//! # Parallel processing
//! Functions for parallel processing can be found in the [`parallel`](parallel/index.html) module



//use std::io;
//use std::str::{self,Utf8Error};
extern crate buf_redux;
extern crate memchr;

use std::error;
use std::fmt;
use std::io;

pub use strategy::*;

mod strategy;


macro_rules! try_opt {
    ($expr: expr) => {
        match $expr {
            Ok(item) => item,
            Err(e) => return Some(Err(::std::convert::From::from(e)))
        }
     };
}


macro_rules! unwrap_or {
    ($expr:expr, $or:block) => {
        match $expr {
            Some(item) => item,
            None => $or
        }
     };
}


pub mod parallel;
pub mod fasta;
pub mod fastq;



/// used by more than one module

#[derive(Default, Debug)]
struct ReadAlways;

impl buf_redux::strategy::ReadStrategy for ReadAlways {
    fn should_read(&self, _: &buf_redux::Buffer) -> bool { true }
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
    where R: io::Read,
          Rs: buf_redux::strategy::ReadStrategy,
          Ms: buf_redux::strategy::MoveStrategy
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

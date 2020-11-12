//! This module provides types and functions for FASTX reading via dynamic
//! dispatch.
//!
//! In contrast to [`fastx::Reader`](crate::fastx::Reader), which has an API
//! equivalent to [`fasta::Reader`](crate::fasta::Reader) and
//! [`fastq::Reader`](crate::fastq::Reader), this module provides the traits
//! [`FastxReader`](FastxReader) and [`FastxSeekReader`](FastxSeekReader).
//! These are implemented for all readers in the `fasta` and `fastq` modules.
//!
//! # Example
//!
//! This example uses the [`reader`](reader) function to obtain a
//! [`Box<dyn FastxReader<...>`](FastxReader). The format is recognized based
//! on the first byte of the first non-empty line, and an appropriate reader
//! for the format is chosen and returned:
//!
//! ```rust
//! use seq_io::prelude::*;
//! use seq_io::fastx::dynamic::reader;
//! use seq_io::fastx::SeqFormat;
//!
//! # fn main() {
//! let fasta = b">id
//! SEQUENCE
//! ";
//!
//! let mut fasta_rdr = reader(&fasta[..], false).unwrap().expect("input empty");
//! // The format should be correctly recognized
//! assert_eq!(fasta_rdr.format(), Some(SeqFormat::FASTA));
//!
//! // The 'next' method is called 'next_fastx' here to prevent name clashes:
//! let record = fasta_rdr.next_fastx().unwrap().unwrap();
//! assert_eq!(record.id(), Ok("id"));
//! # }
//! ```
//!
//! An example with FASTQ:
//!
//! ```rust
//! # use seq_io::prelude::*;
//! # use seq_io::fastx::dynamic::reader;
//! # use seq_io::fastx::SeqFormat;
//! # fn main() {
//! let fastq = b"
//! @id
//! SEQUENCE
//! +
//! QUALITYY
//! ";
//!
//! let mut fastq_rdr = reader(&fastq[..], false).unwrap().expect("input empty");
//! // The format should be correctly recognized
//! assert_eq!(fastq_rdr.format(), Some(SeqFormat::FASTQ));
//!
//! // The 'next' method is called 'next_fastx' here to prevent name clashes:
//! let record = fastq_rdr.next_fastx().unwrap().unwrap();
//! assert_eq!(record.id(), Ok("id"));
//! # }
//! ```
//!
//! Empty input will lead to `None`:
//!
//! ```rust
//! # use seq_io::prelude::*;
//! # use seq_io::fastx::dynamic::reader;
//! # use seq_io::fastx::SeqFormat;
//! # fn main() {
//! let empty_lines = b"
//!
//!
//! ";
//!
//! assert!(reader(&empty_lines[..], false).unwrap().is_none());
//! # }
//! ```
//!
//! # Using ReaderBuilder
//!
//! [`ReaderBuilder`](ReaderBuilder) allows for applying different configuration
//! options before the actual reader is created.
//!
//! This example applies several non-default settings at once:
//!
//! ```rust
//! use seq_io::fastx::dynamic::ReaderBuilder;
//! use seq_io::fastx::SeqFormat;
//! use seq_io::policy::DoubleUntil;
//!
//! # fn main() {
//! let fastq = b"
//! @id
//! SEQUENCE
//! +
//! QUALITYY
//! ";
//!
//! let policy = DoubleUntil(16 * 1024 * 1024);
//! let capacity = 32 * 1024 * 1024;
//!
//! let mut reader = ReaderBuilder::new()
//!     .buf_policy(policy)
//!     .capacity(capacity)
//!     .multiline_fastq(true)
//!     .from_reader(&fastq[..])
//!     .unwrap()
//!     .expect("empty input");
//!
//! assert_eq!(reader.format(), Some(SeqFormat::FASTQ));
//! # }
//! ```
//!
//! Another use of `ReaderBuilder` is to construct a
//! [`FastxSeekReader`](FastxSeekReader) from a source other than [`File`](std::fs::File)
//! (in which case you could use [`from_path`](from_path)):
//!
//! ```rust
//! use seq_io::prelude::*;
//! use seq_io::fastx::dynamic::ReaderBuilder;
//! use std::io::Cursor;
//!
//! # fn main() {
//! let fastq = b"
//! >id1
//! SEQUENCE1
//! >id2
//! SEQUENCE2
//! ";
//! let input = Cursor::new(&fastq[..]);
//!
//! let mut reader = ReaderBuilder::new()
//!     .from_seekable(input)
//!     .unwrap()
//!     .expect("empty input");
//!
//! // first record
//! let rec1 = reader.next_fastx().unwrap().expect("no record");
//! let id1 = rec1.id().unwrap().to_owned();
//! let pos1 = reader.position_fastx();
//!
//! // proceed to second record
//! let _ = reader.next_fastx().unwrap().expect("no record");
//! assert!(reader.next_fastx().is_none());
//!
//! // now seek back
//! reader.seek_fastx(&pos1);
//! let rec1_again = reader.next_fastx().unwrap().expect("no record");
//! assert_eq!(rec1_again.id().unwrap(), &id1);
//! # }
//! ```

use std::convert::AsRef;
use std::fs::File;
use std::io;
use std::marker::PhantomData;
use std::path::Path;

use crate::core::{BufReader, QualRecordPosition, BUFSIZE};
use crate::fastx::{recognize_format, LineStore, SeqFormat};
use crate::policy::{BufPolicy, StdPolicy};
use crate::{fasta, fastq};

pub use super::error::{Error, Result};
pub use super::record::{OwnedRecord, Record, RecordSet, RefRecord};

/// Recognizes the sequence format of the input and returns an appropriate FASTA
/// or FASTQ reader as [`Box<dyn FastxReader<...>`](FastxReader) or `None` if
/// the input has only empty lines.
///
/// `reader` is a convenience function for the most frequent case. If options
/// need to be changed, use [`ReaderBuilder`](ReaderBuilder). The following two
/// calls are equivalent:
///
/// ```rust ignore
/// // using the reader() function:
/// let reader = seq_io::fastx::dynamic::reader(rdr, multiline_fastq).unwrap();
///
/// // using the builder API:
/// let reader = ReaderBuilder::new()
///     .set_multiline_fastq(multiline_fastq)
///     .from_reader(rdr)
///     .unwrap();
/// ```
///
/// If the first non-empty line starts with `>`, a
/// [`fasta::Reader`](crate::fasta::Reader) is returned. If it starts with
/// `@`, [`fastq::Reader`](crate::fastq::Reader) is returned if
/// `multiline_fastq` is `false`, otherwise
/// [`fastq::multiline::Reader`](crate::fastq::multiline::Reader) is returned,
/// If the first non-empty line contains an invalid start byte, an error with
/// `ErrorKind::InvalidStart` will be returned.
pub fn reader<'s, R>(
    reader: R,
    multiline_fastq: bool,
) -> Result<Option<Box<dyn FastxReader<R, StdPolicy, LineStore> + 's>>>
where
    R: io::Read + 's,
{
    ReaderBuilder::new()
        .multiline_fastq(multiline_fastq)
        .from_reader(reader)
}

/// Recognizes the sequence format of the file and returns an appropriate FASTA
/// or FASTQ reader as [`Box<dyn FastxSeekReader<...>`](FastxSeekReader) or `None` if
/// the input has only empty lines. See [`reader`](reader) for more information
/// about the behaviour.
pub fn from_path<'s, P>(
    path: P,
    multiline_fastq: bool,
) -> Result<Option<Box<dyn FastxSeekReader<File, StdPolicy, LineStore> + 's>>>
where
    P: AsRef<Path> + 's,
{
    ReaderBuilder::new()
        .multiline_fastq(multiline_fastq)
        .from_path(path)
}

macro_rules! get_reader {
    ($builder:ident, $reader:expr, $ReaderType:ident) => {{
        let multiline_fastq = $builder.multiline_fastq;
        let mut buf_reader =
            BufReader::with_capacity($reader, $builder.capacity).set_policy($builder.buf_policy);
        recognize_format(&mut buf_reader)?
            .map(|(fmt, (byte, line))| {
                let out: Box<dyn $ReaderType<_, _, _>> = match fmt {
                    SeqFormat::FASTA => {
                        Box::new(fasta::Reader::from_buf_reader(buf_reader, byte, line))
                    }
                    SeqFormat::FASTQ => {
                        if multiline_fastq {
                            Box::new(fastq::multiline::Reader::from_buf_reader(
                                buf_reader, byte, line,
                            ))
                        } else {
                            //  println!("new fastq {} {}", byte, line);
                            Box::new(fastq::Reader::from_buf_reader(buf_reader, byte, line))
                        }
                    }
                };
                Ok(out)
            })
            .transpose()
    }};
}

/// Allows building a [`FastxReader`](FastxReader)
/// or [`FastxSeekReader`](FastxReader)
/// with various configuration options.
#[derive(Debug)]
pub struct ReaderBuilder<P: BufPolicy = StdPolicy, S: QualRecordPosition = LineStore> {
    buf_policy: P,
    _store: PhantomData<S>,
    capacity: usize,
    multiline_fastq: bool,
}

impl ReaderBuilder {
    /// Creates a new `ReaderBuilder`
    #[inline]
    pub fn new() -> Self {
        ReaderBuilder {
            buf_policy: StdPolicy,
            _store: PhantomData,
            capacity: BUFSIZE,
            multiline_fastq: false,
        }
    }
}

impl<P, S> ReaderBuilder<P, S>
where
    P: BufPolicy,
    S: QualRecordPosition,
{
    /// Creates a new [`FastxReader`](FastxReader)
    /// with the current configuration of the builder.
    /// See [`reader`](reader) for more information about the behaviour.
    #[inline]
    pub fn from_reader<'s, R>(self, reader: R) -> Result<Option<Box<dyn FastxReader<R, P, S> + 's>>>
    where
        R: io::Read + 's,
        P: 's,
        S: 's,
    {
        get_reader!(self, reader, FastxReader)
    }

    /// Creates a new [`FastxSeekReader`](FastxSeekReader) from the given path
    /// with the current configuration of the builder.
    #[inline]
    pub fn from_path<'s, Q>(
        self,
        path: Q,
    ) -> Result<Option<Box<dyn FastxSeekReader<File, P, S> + 's>>>
    where
        Q: AsRef<Path> + 's,
        P: 's,
        S: 's,
    {
        let f = File::open(path)?;
        self.from_seekable(f)
    }

    /// Creates a new [`FastxSeekReader`](FastxSeekReader)
    /// with the current configuration of the builder.
    #[inline]
    pub fn from_seekable<'s, R>(
        self,
        reader: R,
    ) -> Result<Option<Box<dyn FastxSeekReader<R, P, S> + 's>>>
    where
        R: io::Read + io::Seek + 's,
        P: 's,
        S: 's,
    {
        get_reader!(self, reader, FastxSeekReader)
    }

    /// Creates a new reader with a given
    /// [`core::QualRecordPosition`](crate::core::QualRecordPosition)
    /// as defined in the type argument of the method. The method consumes
    /// the reader and returns a new `Reader` instance.
    #[inline]
    pub fn pos_store<T: QualRecordPosition>(self) -> ReaderBuilder<P, T> {
        ReaderBuilder {
            buf_policy: self.buf_policy,
            _store: PhantomData,
            capacity: self.capacity,
            multiline_fastq: self.multiline_fastq,
        }
    }

    /// Applies a [`BufPolicy`](crate::policy::BufPolicy) to the
    /// current reader. The method consumes the reader and returns a new
    /// `Reader` instance.
    #[inline]
    pub fn buf_policy<T: crate::policy::BufPolicy>(self, buf_policy: T) -> ReaderBuilder<T, S> {
        ReaderBuilder {
            buf_policy,
            _store: self._store,
            capacity: self.capacity,
            multiline_fastq: self.multiline_fastq,
        }
    }

    /// Sets a reader capacity.
    #[inline]
    pub fn capacity(mut self, cap: usize) -> Self {
        self.capacity = cap;
        self
    }

    /// Sets whether multi-line FASTQ should be recognized by the reader.
    /// If `true`, the final reader will be
    /// [`fastq::multiline::Reader`](crate::fastq::multiline::Reader)
    /// in the case of FASTQ input.
    #[inline]
    pub fn multiline_fastq(mut self, multiline: bool) -> Self {
        self.multiline_fastq = multiline;
        self
    }
}

/// Trait for FASTX reading.
///
/// Provides the same methods as individual `Reader` types, but in contrast of
/// returning / operating on different `RefRecord` / `RecordSet` /
/// `RecordsIter` ... types, types from the `fastx` module are used, allowing
/// for FASTX functionality with dynamic dispatch.
pub trait FastxReader<R, P, S>
where
    R: io::Read,
    P: BufPolicy,
    S: QualRecordPosition,
{
    /// Returns the next [`fastx::RefRecord`](crate::fastx::RefRecord), if any.
    fn next_fastx(&mut self) -> Option<Result<RefRecord<S>>>;

    /// Updates a [`fastx::RecordSet`](crate::fastx::RecordSet) with new data.
    fn read_record_set_fastx(&mut self, record_set: &mut RecordSet<S>) -> Result<bool>;

    fn read_record_set_exact_fastx(
        &mut self,
        record_set: &mut crate::fastx::RecordSet<S>,
        n_records: usize,
    ) -> crate::fastx::Result<bool>;

    /// Returns the sequence format (`SeqFormat::FASTA` or `SeqFormat::FASTQ`)
    /// if known. For FASTA / FASTQ readers, the format is known in the
    /// beginning already, for FASTX readers, a first record has to be returned
    /// successfully (without error) before the format can be known.
    /// If not able to call this method in the middle of a parsing process
    /// because of borrowing issues, there is also
    /// [`has_quality()`](crate::BaseRecord::has_quality)
    /// method, which can be called on individual sequence records and
    /// essentially gives the same information.
    fn format(&self) -> Option<SeqFormat>;

    /// Returns a borrowed iterator over all records.
    fn fastx_records(&mut self) -> RecordsIter<R, P, S>;

    /// Returns an iterator over all records like `FastxReader::records()`,
    /// but with the difference that it owns the underlying reader.
    fn into_fastx_records<'a>(self) -> RecordsIntoIter<'a, R, P, S>
    where
        R: 'a,
        P: 'a,
        S: 'a;

    /// Returns the position of the current record (found by the previous call
    /// to [`next_fastx()`](#next_fastx)). Equivalent to the `position()` method
    /// of individual readers. Useful with [`seek_fastx()`](FastxSeekReader).
    fn position_fastx(&self) -> crate::Position;
}

// TODO: implementation for box not needed?
// impl<'a, R, P, S> FastxReader<R, P, S> for Box<dyn FastxReader<R, P, S> + 'a>
// where
//     R: io::Read,
//     P: BufPolicy,
//     S: QualRecordPosition,
// {
//     fn next_fastx(&mut self) -> Option<Result<RefRecord<S>>> {
//         (**self).next_fastx()
//     }

//     fn read_record_set_fastx(&mut self, record_set: &mut RecordSet<S>) -> Result<bool> {
//         (**self).read_record_set_fastx(record_set)
//     }

//     fn format(&self) -> Option<SeqFormat> {
//         (**self).format()
//     }

//     fn fastx_records(&mut self) -> RecordsIter<R, P, S> {
//         (**self).fastx_records()
//     }

//     fn into_fastx_records<'s>(self) -> RecordsIntoIter<'s, R, P, S>
//     where
//         R: 's,
//         P: 's,
//         S: 's
//     {
//         <Self as FastxReader<R, P, S>>::into_fastx_records(self)
//     }

//     fn position_fastx(&self) -> crate::Position {
//         (**self).position_fastx()
//     }
// }

pub trait FastxSeekReader<R, P, S>: FastxReader<R, P, S>
where
    R: io::Read + io::Seek,
    P: BufPolicy,
    S: QualRecordPosition,
{
    /// Seeks to a specified position. Equivalent to the `seek()` method of the
    /// individual readers.
    fn seek_fastx(&mut self, pos: &crate::Position) -> io::Result<()>;
}

/// Borrowed iterator of `OwnedRecord`
pub struct RecordsIter<'a, R, P, S>
where
    P: crate::policy::BufPolicy + 'a,
    R: std::io::Read + 'a,
    S: QualRecordPosition + 'a,
{
    rdr: &'a mut dyn FastxReader<R, P, S>,
}

impl<'a, R, P, S> RecordsIter<'a, R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: QualRecordPosition,
{
    #[inline]
    pub(crate) fn new(rdr: &'a mut dyn FastxReader<R, P, S>) -> Self {
        RecordsIter { rdr }
    }
}

impl<'a, R, P, S> Iterator for RecordsIter<'a, R, P, S>
where
    P: crate::policy::BufPolicy + 'a,
    R: std::io::Read + 'a,
    S: QualRecordPosition,
{
    type Item = Result<OwnedRecord>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.rdr
            .next_fastx()
            .map(|rec| rec.map(|r| r.to_owned_record()))
    }
}

/// Iterator of `OwnedRecord` that owns the underlying reader
pub struct RecordsIntoIter<'a, R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: QualRecordPosition,
{
    rdr: Box<dyn FastxReader<R, P, S> + 'a>,
    _s: std::marker::PhantomData<S>,
}

impl<'a, R, P, S> RecordsIntoIter<'a, R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: QualRecordPosition,
{
    #[inline]
    pub(crate) fn new(rdr: Box<dyn FastxReader<R, P, S> + 'a>) -> Self {
        RecordsIntoIter {
            rdr,
            _s: std::marker::PhantomData,
        }
    }
}

impl<'a, R, P, S> Iterator for RecordsIntoIter<'a, R, P, S>
where
    P: crate::policy::BufPolicy,
    R: std::io::Read,
    S: QualRecordPosition,
{
    type Item = Result<OwnedRecord>;
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        self.rdr
            .next_fastx()
            .map(|rec| rec.map(|r| r.to_owned_record()))
    }
}

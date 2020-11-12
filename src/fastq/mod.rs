//! FASTQ reading and writing
//!
//! # Flavours
//!
//! There are two flavours of this parser:
//!
//! * [`fastq::Reader`](crate::fastq::Reader) in this module parses standard
//!   single-line FASTQ.
//! * [`Reader`](crate::fastq::multiline::Reader) in [`fastq::multiline`](crate::fastq::multiline)
//!   parses multi-line FASTQ as well. This parser runs slightly slower on
//!   single-line FASTQ than `fastq::Reader`.
//!
//! # Example
//!
//! The following example shows how to use [`Reader`](Reader).
//!
//! ```rust
//! use seq_io::prelude::*;  // needed to import necessary traits
//! use seq_io::fastq::Reader;
//!
//! # fn main() {
//! let seq = b"@id1 some description
//! SEQUENCE
//! +
//! IIIIIIII
//! @id2
//! SEQUENCE
//! +
//! IIIIIIII
//! ";
//!
//! // Construct the reader
//! let mut reader = Reader::new(&seq[..]);
//!
//! // We'll write the records back to this vector
//! let mut output = vec![];
//!
//! while let Some(result) = reader.next() {
//!     let rec = result.unwrap();
//!
//!     // Access the ID and the description parts of the header (separated by a space)
//!     let id = rec.id().unwrap();
//!     let desc = rec.desc().transpose().unwrap();
//!     println!("ID: {}, description: {:?}", id, desc);
//!
//!     // Print the sequence and quality scores
//!     println!("seq:  {}", std::str::from_utf8(rec.seq()).unwrap());
//!     println!("qual: {}", std::str::from_utf8(rec.qual()).unwrap());
//!
//!     // Write the record to 'output'
//!     rec.write(&mut output).unwrap();
//! }
//!
//! // The output is identical
//! assert_eq!(&seq[..], output.as_slice());
//! # }
//! ```
//!
//! The output will be:
//!
//! ```text
//! ID: id1, description: Some("some description")
//! seq:  SEQUENCE
//! qual: IIIIIIII
//! ID: id2, description: None
//! seq:  SEQUENCE
//! qual: IIIIIIII
//! ```
//!
//! As the record returned by the [`next()`](Reader::next)
//! method borrows its data from the underlying buffer, it is not possible to
//! use a `for` loop for iterating. Therefore, we use the `while let ...`
//! construct.
//!
//! # Sequence record types
//!
//! Similarly to [`fasta`](crate::fasta)`::Reader`, there are two record types,
//! which both implement the common [`BaseRecord`](crate::BaseRecord)
//! trait, [and `fastq::Record`](Record) providing additional FASTQ
//! specific methods:
//!
//! * [`RefRecord`](RefRecord), the type returned by
//!   `Reader::next()`, only remembers the position of the record in the buffer
//!   without copying any data.
//! * [`OwnedRecord`](OwnedRecord) owns its data.
//!
//! # Writing FASTQ
//!
//! Records can be written to output using
//! [`BaseRecord::write()`](crate::BaseRecord::write).
//! `RefRecord` additionally has the method
//! [`write_unchanged`](RefRecord::write_unchanged), which may be
//! faster.
//!
//! It is also possible to write data not part of a FASTQ record directly using
//! a set of different functions [listed here](#functions).
//!
//! # Details on parsing and writing
//!
//! * Like all parsers in this crate, `fasta::Reader` handles UNIX (LF) and
//!   Windows (CRLF) line endings, but not old Mac-style (CR) endings. LF and
//!   CRLF may be mixed within the same file.
//! * FASTQ writing currently always uses UNIX line endings.
//! * The first non-empty line should start with `@`, indicating the first
//!   header. If not, an error with `ErrorKind::InvalidStart` is returned.
//! * Empty lines are allowed before and after records, but *not within*
//!   records.
//! * Whitespace at the end of header and sequence lines is never removed.
//! * Empty input will result in `None` being returned immediately by
//!   `Reader::next()` and in empty iterators for `RecordsIter` /
//!   `RecordsIntoIter`.
//! * `Reader::next()` compares sequence and quality line lengths and returns
//!   an error of `ErrorKind::UnequalLengths` if different. It is possible to
//!   omit this check by calling
//!   [`Reader::next_unchecked_len()`](struct.Reader::next_unchecked_len).
//!   The lengths can be checked later with [`Record::check_lengths()`](Record::check_lengths)
//!   or [`RefRecord::check_lengths_strict()`](RefRecord::check_lengths_strict)
//! * The quality line of the last record should either terminated by a line ending.
//!   If not, an error of `ErrorKind::UnexpectedEnd` is returned.
//!
//! ## Error priority
//!
//! Validity checks are done in the following order:
//!
//! * Is the start byte correct (`@`)? If not: [`InvalidStart`](ErrorKind::InvalidStart).
//! * Do the the header, sequence and separator lines have a line terminator?
//!   If not: [`UnexpectedEnd`](ErrorKind::UnexpectedEnd).
//! * Is the separator byte correct (`+`)? If not: [`InvalidSep`](ErrorKind::InvalidSep).
//! * Do the quality scores have a line terminator *or* are they
//!   non-empty? If not: [`UnexpectedEnd`](ErrorKind::UnexpectedEnd).
#[macro_use]
mod error;
#[macro_use]
mod reader;
pub mod multiline;
mod position;
mod record;
mod write;

pub use self::error::*;
pub use self::position::*;
pub use self::reader::*;
pub use self::record::*;
pub use self::write::*;

//! FASTA reading and writing
//!
//! # Flavours
//!
//! There are two flavours of this parser:
//!
//! * [`fasta::Reader`](crate::fasta::Reader) in this module parses standard
//!   multi-line FASTA.
//! * [`Reader`](crate::fasta::single_line::Reader) in [`fasta::single_line`](crate::fasta::single_line)
//!   accepts *only* single-line FASTA. This may not be of general use,
//!   the only advantage of this parser is that it runs faster than the multi-
//!   line parser.
//!
//! # Example
//!
//! The following example shows how to use [`Reader`](Reader).
//!
//! ```rust
//! use seq_io::prelude::*;  // needed to import necessary traits
//! use seq_io::fasta::Reader;
//!
//! # fn main() {
//! let seq = b">id1 some description
//! SEQUENCE
//! .ANOTHER.LINE
//! >id2
//! SEQUENCE
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
//!     println!("ID: '{}', description: {:?}", id, desc);
//!
//!     // Print the sequence
//!     println!("{}", std::str::from_utf8(rec.seq()).unwrap());
//!
//!     // Write the record to 'output'
//!     rec.write(&mut output).unwrap();
//! }
//!
//! println!("\nWritten:\n{}", std::str::from_utf8(&output).unwrap());
//! # }
//! ```
//!
//! The output will be:
//!
//! ```text
//! ID: 'id1', description: Some("some description")
//! SEQUENCE
//! .ANOTHER.LINE
//! ID: 'id2', description: None
//! SEQUENCE
//!
//! Written:
//! >id1 some description
//! SEQUENCE.ANOTHER.LINE
//! >id2
//! SEQUENCE
//! ```
//!
//! As the record returned by the [`next()`](Reader::next)
//! method borrows its data from the underlying buffer, it is not possible to
//! use a `for` loop for iterating. Therefore, we use the `while let ...`
//! construct.
//!
//! The ID of the record can be accessed using `record.id()`, which first has
//! to check the validity of the UTF-8 first.
//! [`record.description()`](crate::BaseRecord::id) returns
//! everything after the space in the header, if present.
//! The whole header line can be accessed using
//! [`record.head()`](crate::BaseRecord::head). This returns a byte
//! slice, which needs to be manually converted to `String`.
//!
//! # Sequence record types
//!
//! * The type of record obtained by calling `Reader::next()` is
//! [`RefRecord`](RefRecord), which is a simple wrapper
//! borrowing its data from the reader and its buffer. Essentially, the position
//! of the record in the buffer is stored without any further copying necessary.
//!
//! * [`OwnedRecord`](OwnedRecord) owns its data
//!  ([see more here](#owned-records)).
//!
//! Both record types implement [`BaseRecord`](crate::BaseRecord), the
//! common trait for all records in this crate. FASTA specific methods are
//! additionally defined in the [`fasta::Record`](Record) trait.
//! It is easiest to simply import these traits (along others) with
//! `use seq_io::prelude::*`.
//!
//! # Accessing the sequence
//!
//! There are different possibibilities to access the sequence in `RefRecord`.
//! The fastest are presented first:
//!
//! ## Taking a slice of everything
//!
//! In the [above example](#simple-parsing), printing the output of `rec.seq()`
//! results in two sequence lines, like in the input. This is on purpose, not a
//! bug. The [`record.seq()`](crate::BaseRecord::seq) method simply
//! takes a slice from the buffer without transforming the sequence in any way.
//! This may make sense in some cases, but if multiple lines are expected in the
//! FASTA input, the sequence will be interrupted by line breaks.
//!
//! ## Iterating over sequence lines
//!
//! Since the position of each line is remembered by the record, iterating over
//! sequence lines is also very fast:
//!
//! ```rust ignore
//! for line in rec.seq_lines() {
//!     println!("line: {}", std::str::from_utf8(&line)?);
//! }
//! ```
//!
//! This would be the output for the first record from the above example:
//!
//! ```text
//! line: SEQUENCE
//! line: .ANOTHER.LINE
//! ```
//!
//! ## Obtaining the full sequence
//!
//! If it is required to have the whole, contiguous sequence in one slice,
//! there are a two possibilities:
//!
//! * [`full_seq()`](crate::BaseRecord::full_seq) returns
//!   [`Cow<[u8]>`](std::borrow::Cow), meaning that if there is only one
//!   sequence line, no copying will be done. Multiple lines are instead copied
//!   into a newly allocated `Vec<u8>`.
//! * [`full_seq_given()`](crate::BaseRecord::full_seq_given)
//!   allows reusing allocations in different calls by letting the user supply
//!   a `Vec<u8>`:
//!
//! ```rust ignore
//! let mut seq = vec![];
//! rec.full_seq_given(|| &mut seq);
//! println!("{}", std::str::from_utf8(&seq)?);
//! ```
//!
//! This would print `SEQUENCE.ANOTHER.LINE` given the first sequence record.
//! The reasoning behind the closure was to also allow it to be used with arena
//! allocators. The closure is only called if needed due to multiple sequence
//! lines.
//!
//! # Owned records
//!
//! With the method shown above it is not possible to iterate using a `for` loop
//! since `Reader` is a streaming iterator, which does not implement `Iterator`.
//! Records cannot be stored in a vector for later reuse. In this case, it is
//! necessary to create an owned copy, which stores the header and sequence in
//! allocated vectors:
//!
//! ```rust ignore
//! let records = vec![];
//! while let Some(result) = reader.next() {
//!     let rec = result?;
//!     records.push(rec.to_owned_record());
//! }
//! ```
//!
//! An even easier way is to use the [`records()`](Reader::records)
//! or [`into_records`](Reader::into_records) iterators:
//!
//! ```rust ignore
//! let records: Vec<_> = reader.records().collect()?;
//! ```
//!
//! Of course, this slows down everything because allocations take time.
//! However, it is also possible to use
//! [`clone_into_owned`](RefRecord::clone_into_owned)
//! in order to reuse the `OwnedRecord`s:
//!
//! ```rust ignore
//! let records = vec![];
//! while let Some(result) = reader.next() {
//!     let rec = result?;
//!     // Obtain the record from somewhere or create a new one using
//!     // OwnedRecord::default()
//!     let mut owned_record = ...
//!     // Update it with new data
//!     rec.clone_into_owned(&mut owned_record);
//!     records.push(owned_record);
//! }
//! ```
//!  
//! # Writing FASTA
//!
//! Records can be written to output using
//! [`BaseRecord::write()`](crate::BaseRecord::write).
//! [`fasta::Record::write_wrap()`](Record::write_wrap)
//! Writes FASTA, wrapping lines at a given width.
//! [`RefRecord`](RefRecord) additionally has the method
//! [`write_unchanged`](RefRecord::write_unchanged), which should be
//! faster, but does not remove line wrap if present or re-wrap lines.
//!
//! It is also possible to write data not part of a FASTA record directly using
//! a set of different functions [listed here](#functions).
//!
//! # Details on parsing and writing
//!
//! Valid FASTA records require the header line to start with `>` and terminated
//! with a line terminator (UNIX-style `\n` or Windows-style `\r\n`).
//! Sequence lines are optional.
//! Therefore, the following is also accepted and interpreted as two records
//! with an empty header and no sequence:
//!
//! ```text
//! >
//! >
//! ```
//!
//! Accordingly, [`BaseRecord::num_seq_lines()`](crate::BaseRecord::num_seq_lines)
//! returns `0` and [`RefRecord::seq_lines()`](RefRecord::seq_lines) an
//! empty iterator.
//!
//! More details:
//!
//! * Like all parsers in this crate, `fasta::Reader` handles UNIX (LF) and
//!   Windows (CRLF) line endings, but not old Mac-style (CR) endings. LF and
//!   CRLF may be mixed within the same file.
//! * FASTA writing currently always uses UNIX line endings.
//! * The first non-empty line should start with `>`, indicating the first
//!   header. If not, an error with `ErrorKind::InvalidStart` is returned.
//! * Whitespace at the end of header and sequence lines is never removed.
//! * Empty input will result in `None` being returned immediately by
//!   `fasta::Reader::next()` and in empty iterators for `RecordsIter` /
//!   `RecordsIntoIter`.
//! * Comment lines starting with `;` are not recognized.
//!   If at the start of a file, there will be an error, since `>` is expected.
//!   Intermediate comments are interpreted as belonging to the sequence.
//! * The last record header requires to be terminated by a line ending.
//!   If not, an error of `ErrorKind::UnexpectedEnd` is returned.
#[macro_use]
mod error;
#[macro_use]
mod reader;
mod position;
mod record;
pub mod single_line;
mod write;

pub use self::error::*;
pub use self::position::*;
pub use self::reader::*;
pub use self::record::*;
pub use self::write::*;

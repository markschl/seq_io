//! FASTX reading and writing
//!
//! The `fastx` reader guesses the sequence format based on examination of the
//! first non-empty line (should start with `>` for FASTA and `@` for FASTQ).
//! Once the format is known, the file is parsed according to the behaviour
//! described in the [`fasta`](crate::fasta) and [`fastq`](crate::fastq).
//!
//! # Flavours
//!
//! There are two flavours of this parser:
//!
//! * [`fastx::Reader`](crate::fastx::Reader) in this module parses FASTA and
//!   and single-line FASTQ
//! * [`Reader`](crate::fastx::multiline_qual::Reader) in
//!   [`fastx::multiline_qual`](crate::fastx::multiline_qual)
//!   parses FASTA and and multi-line FASTQ.
//!
//! # Trait object approach
//!
//! In addition to the above readers, an approach to FASTX reading using trait
//! objects was implemented. It can be found in the [`fastx::dynamic`](crate::fastx::dynamic)
//! module. The API is different, but has the advantage that it theoretically
//! allows integrating sequence parsers not defined in this crate.
//!
//! # Example
//!
//! This example parses FASTQ and transcribes it to FASTA.
//!
//! ```rust
//! use seq_io::prelude::*;
//! use seq_io::fastx::Reader;
//! use seq_io::fastx::SeqFormat;
//!
//! # fn main() {
//! let fastq = b"@id1
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
//! let mut reader = Reader::new(&fastq[..]);
//!
//! // Records are written here
//! let mut output = vec![];
//!
//! while let Some(result) = reader.next() {
//!     let rec = result.unwrap();
//!     rec.write_as(&mut output, SeqFormat::FASTA, None).unwrap();
//! }
//!
//! let fasta_equivalent = b">id1
//! SEQUENCE
//! >id2
//! SEQUENCE
//! ";
//!
//! assert_eq!(output, fasta_equivalent);
//! # }
//! ```

#[macro_use]
mod error;
#[macro_use]
mod record;
#[macro_use]
pub mod dynamic;
#[macro_use]
mod reader;
pub mod multiline_qual;
mod position;
mod recognition;

pub use self::error::*;
pub use self::position::*;
pub use self::reader::*;
pub use self::recognition::*;
pub use self::record::*;

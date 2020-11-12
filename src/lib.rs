//! This library provides parsing and writing of FASTA, FASTQ and FASTX at a
//! high performance.
//!
//! For a detailed documentation of the components, please refer to the
//! module pages listed here:
//!
//! * **FASTA**: See [`fasta`](crate::fasta) module for an introduction.
//! * **FASTQ**: See [`fastq`](crate::fastq) module. A separate parser supporting
//!   multi-line FASTQ is found in [`fastq::multiline`](crate::fastq::multiline).
//! * **FASTX**: There are two approaches:
//!     - The parsers from the [`fastx`](crate::fastx) and
//!       [`fastx::multiline_qual`](crate::fastx::multiline_qual) module
//!     - An approach based on trait objects and dynamic dispatch found in the
//!       [`fastx::dynamic`](crate::fastx::dynamic) module.
//!
//! # Special features
//!
//! The following features are not covered in the module docs for the individual
//! formats:
//!
//! ## Parallel processing
//!
//! All readers allow reading sequence records into chunks called record sets.
//! Effectively, they are just copies of the internal buffer with associated
//! positional information. These record sets can be sent around in channels
//! without the overhead that sending single records would have. The idea for
//! this was borrowed from [fastq-rs](https://github.com/aseyboldt/fastq-rs).
//!
//! The [`parallel`](crate::parallel) module offers a few functions for
//! processing sequence records and record sets in a worker pool and then
//! sending them along with the processing results to the main thread.
//! The functions work; the API design may not be optimal yet.
//!
//! ## Position tracking and seeking
//!
//! All readers keep track of the byte offset, line number and record number
//! while parsing. The current position can be stored and used later for seeking
//! back to the same position. See [`here`](crate::fasta::Reader::seek) for
//! an example.
//!
//! It is not yet possible to restore a record completely given positional
//! information (such as from a `.fai` file). All that is done currently is
//! to set the position, so `next()` will return the correct record.
//!
//! # Design notes
//!
//! Apart from `R: io::Read`, all readers have two additional generic
//! parameters. It is normally not necessary to change the defaults, but in
//! some cases this may be relevant.
//!
//! ## Buffer growth policy
//!
//! The parsers avoid allocations and copying as much as possible.
//! To achieve this, each sequence record must fit into the underlying
//! buffer as a whole. This may not be possible if dealing with large sequences.
//! Therefore, the internal buffer of the reader will grow automatically to fit
//! the whole sequence record again. The buffer may grow until it reaches 1 GiB;
//! larger records will cause an error.
//!
//! The behaviour of buffer growth can be further configured by applying
//! a different policy. This is documented in the [`policy`](policy) module.
//!
//! ## Position stores
//!
//! At the core of the different parsers is the same code, which is called
//! with different parameters. While searching the buffer for sequence records,
//! the position of the different features is stored. This allows to later
//! access the header, sequence and quality features directly as slices taken
//! from the internal buffer. Which positional information needs to be stored
//! depends on the format. For example, the [`fasta`](crate::fasta) reader
//! stores the position of every sequence line in order to allow fast iteration
//! lines later. The [`fastq`](crate::fastq) reader needs to remember the
//! position of the quality scores, but doesn't need to store information about
//! multiple lines, which allows for a simpler data structure. In turn,
//! the [`fastx`](crate::fastx) reader needs to store FASTA lines *and* quality
//! scores.
//!
//! Therefore, all readers have a third generic parameter, which allows
//! assigning a specific "storage backend" implementing the traits
//! [`SeqRecordPosition`](crate::SeqRecordPosition) and
//! [`QualRecordPosition`](crate::QualRecordPosition). Usually, it is
//! not necessary to deal with this parameter since each parser has a reasonable
//! default. The only case where it is changed in this crate is with the trait
//! object approach implemented in [`fastx::dynamic`](crate::fastx::dynamic).

pub use self::error::*;
pub use self::helpers::*;
pub use self::record::*;

mod helpers;
mod record;
#[macro_use]
mod error;

#[macro_use]
pub mod core;
#[macro_use]
pub mod fastx;
pub mod fasta;
pub mod fastq;
pub mod parallel;
pub mod policy;
pub mod prelude;
// pub mod paired;

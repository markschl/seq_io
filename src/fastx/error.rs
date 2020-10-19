use crate::ErrorPosition;
use std::convert::From;
use std::fmt;
use std::io;

/// The kind of error. Currently it has the same options as `seq_io::fastq::ErrorKind`.
#[derive(Debug)]
pub enum ErrorKind {
    /// `std::io::Error`
    Io(io::Error),
    /// Invalid start byte encountered (expected `>` or `@`, depending on the format).
    /// This error occurs either at the beginning of the file if the first byte of the
    /// first non-empty line does matches neither of the two possible characters,
    /// required for either FASTA or FASTQ.
    /// In later records, this error will be returned if the start byte of a record
    /// does not match the expected byte chosen at the beginning.
    InvalidStart {
        /// Position, where the error occurred. `ErrorPosition::position()`
        /// returns the record start (same as `ErrorPosition::record_position()`),
        /// `ErrorPosition::error_offset()` will return `None`.
        pos: ErrorPosition,
        /// Byte found instead.
        found: u8,
    },
    /// Invalid separator byte encountered (expected `+`). This error cannot
    /// occur with multi-line FASTQ.
    InvalidSep {
        /// Position, where the error occurred. `ErrorPosition::position()`
        /// returns the position of the byte that should be `+`.
        pos: ErrorPosition,
        /// Byte found instead (`Some` if any, otherwise `None`)
        found: Option<u8>,
    },
    /// Truncated record found at the end of the input.
    UnexpectedEnd {
        /// Position, where the error occurred. `ErrorPosition::position()`
        /// returns a `Position` referring to the last byte of the file.
        /// `ErrorPosition::record_id()` only returns an ID if the
        /// unexpected end did not occur in within the header line. In this
        /// case, the ID may be truncated.
        pos: ErrorPosition,
    },
    /// Sequence and qualitiy lengths found to be different.
    UnequalLengths {
        /// Position, where the error occurred. `ErrorPosition::position()`
        /// returns the record start (same as `ErrorPosition::record_position()`),
        /// `ErrorPosition::error_offset()` will return `None`.
        /// The position will be `None` if the error occurred in a call to
        /// `BaseRecord::check_lengths`.
        pos: ErrorPosition,
        /// Length of sequence
        seq: usize,
        /// Length of quality information
        qual: usize,
    },
    /// Size limit of buffer was reached, which happens if `policy::BufPolicy::grow_to()` returned
    /// `None`. This does not happen with the default `DoubleUntil` policy,
    /// since this policy does not impose a memory limit.
    BufferLimit,
    /// Hints that destructuring should not be exhaustive,
    /// makes sure that adding new variants will not break the code.
    #[doc(hidden)]
    __Nonexhaustive,
}

impl ErrorKind {
    /// Returns the position for this error, if one exists.
    pub fn position(&self) -> Option<&ErrorPosition> {
        match self {
            ErrorKind::InvalidStart { pos, .. } => Some(pos),
            ErrorKind::InvalidSep { pos, .. } => Some(pos),
            ErrorKind::UnexpectedEnd { pos } => Some(pos),
            ErrorKind::UnequalLengths { pos, .. } => Some(pos),
            _ => None,
        }
    }
}

impl_error!(ErrorKind);

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.kind() {
            ErrorKind::Io(ref e) => e.fmt(f),
            ErrorKind::InvalidStart { pos, found } => write!(
                f,
                "FASTQ parse error: expected '@' at record start but found '{}' ({})",
                (*found as char).escape_default(),
                pos
            ),
            ErrorKind::InvalidSep { pos, found } => {
                write!(f, "FASTQ parse error: Expected '+' separator but found ")?;
                if let Some(b) = found {
                    write!(f, "'{}'", (*b as char).escape_default())?;
                } else {
                    write!(f, "empty line")?;
                }
                write!(f, " ({}).", pos)
            }
            ErrorKind::UnexpectedEnd { pos } => {
                write!(f, "FASTQ parse error: unexpected end of input ({})", pos)
            }
            ErrorKind::UnequalLengths { pos, seq, qual } => write!(
                f,
                "FASTQ parse error: sequence length is {}, but quality length is {} ({})",
                seq, qual, pos
            ),
            ErrorKind::BufferLimit => write!(f, "FASTQ parse error: Buffer limit reached"),
            _ => Ok(()),
        }
    }
}

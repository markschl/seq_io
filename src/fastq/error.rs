use crate::fastx;
use crate::ErrorPosition;
use std::fmt;
use std::io;

#[derive(Debug)]
pub enum ErrorKind {
    /// `std::io::Error`
    Io(io::Error),
    /// Invalid start byte encountered (expected `@`). With multi-line FASTQ,
    /// this error may only occur at the beginning of the file.
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
    /// // TODO: copy-paste
    /// This is a convenience function that permits callers to easily access
    /// the position on an error without doing case analysis on `ErrorKind`.
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
            ErrorKind::UnequalLengths { pos, seq, qual } => write!(
                f,
                "FASTQ parse error: sequence length is {}, but quality length is {} ({})",
                seq, qual, pos
            ),
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
            ErrorKind::BufferLimit => write!(f, "FASTQ parse error: Buffer limit reached"),
            _ => Ok(()),
        }
    }
}

impl From<Error> for fastx::Error {
    fn from(e: Error) -> Self {
        let kind = match e.into_kind() {
            ErrorKind::Io(e) => fastx::ErrorKind::Io(e),
            ErrorKind::InvalidStart { pos, found } => fastx::ErrorKind::InvalidStart { pos, found },
            ErrorKind::InvalidSep { pos, found } => fastx::ErrorKind::InvalidSep { pos, found },
            ErrorKind::UnexpectedEnd { pos } => fastx::ErrorKind::UnexpectedEnd { pos },
            ErrorKind::UnequalLengths { pos, seq, qual } => {
                fastx::ErrorKind::UnequalLengths { pos, seq, qual }
            }
            ErrorKind::BufferLimit => fastx::ErrorKind::BufferLimit,
            ErrorKind::__Nonexhaustive => fastx::ErrorKind::__Nonexhaustive,
        };
        fastx::Error::new(kind)
    }
}

impl From<fastx::Error> for Error {
    fn from(e: fastx::Error) -> Self {
        let kind = match e.into_kind() {
            fastx::ErrorKind::Io(e) => ErrorKind::Io(e),
            fastx::ErrorKind::InvalidStart { pos, found } => ErrorKind::InvalidStart { pos, found },
            fastx::ErrorKind::InvalidSep { pos, found } => ErrorKind::InvalidSep { pos, found },
            fastx::ErrorKind::UnexpectedEnd { pos } => ErrorKind::UnexpectedEnd { pos },
            fastx::ErrorKind::UnequalLengths { pos, seq, qual } => {
                ErrorKind::UnequalLengths { pos, seq, qual }
            }
            fastx::ErrorKind::BufferLimit => ErrorKind::BufferLimit,
            fastx::ErrorKind::__Nonexhaustive => ErrorKind::__Nonexhaustive,
        };
        Error::new(kind)
    }
}

use crate::fastx;
use crate::ErrorPosition;
use std::fmt;
use std::io;

#[derive(Debug)]
pub enum ErrorKind {
    /// `std::io::Error`
    Io(io::Error),
    /// Invalid start byte encountered (expected `>`)
    InvalidStart {
        /// Position, where the error occurred. `ErrorPosition::position()`
        /// returns the record start (same as `ErrorPosition::record_position()`),
        /// `ErrorPosition::error_offset()` will return `None`.
        pos: ErrorPosition,
        /// Byte found instead.
        found: u8,
    },
    /// Truncated record found at the end of the input. This occurs if
    /// there is no newline character `\n` after a sequence header.
    UnexpectedEnd {
        /// Position, where the error occurred. `ErrorPosition::position()`
        /// returns a `Position` referring to the last byte of the file.
        /// `ErrorPosition::record_id()` only returns an ID if the
        /// unexpected end did not occur in within the header line. In this
        /// case, the ID may be truncated.
        pos: ErrorPosition,
    },
    /// Size limit of buffer was reached, which happens if `policy::BufPolicy::grow_to()` returned
    /// `None`. This does not happen with the default `DoubleUntil` policy.
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
                "FASTA parse error: expected '>' at record start but found '{}' ({}).",
                (*found as char).escape_default(),
                pos
            ),
            ErrorKind::BufferLimit => write!(f, "FASTA parse error: Buffer limit reached."),
            _ => Ok(()),
        }
    }
}

impl From<Error> for fastx::Error {
    fn from(e: Error) -> Self {
        let kind = match e.into_kind() {
            ErrorKind::Io(e) => fastx::ErrorKind::Io(e),
            ErrorKind::InvalidStart { pos, found } => fastx::ErrorKind::InvalidStart { pos, found },
            ErrorKind::UnexpectedEnd { pos } => fastx::ErrorKind::UnexpectedEnd { pos },
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
            fastx::ErrorKind::InvalidSep { .. } => unimplemented!(),
            fastx::ErrorKind::UnexpectedEnd { pos } => ErrorKind::UnexpectedEnd { pos },
            fastx::ErrorKind::UnequalLengths { .. } => unimplemented!(),
            fastx::ErrorKind::BufferLimit => ErrorKind::BufferLimit,
            fastx::ErrorKind::__Nonexhaustive => ErrorKind::__Nonexhaustive,
        };
        Error::new(kind)
    }
}

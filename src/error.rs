use crate::Position;
use std::fmt;

/// Position of a parsing error within the sequence record
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ErrorOffset {
    /// Line, at which the error occurred, specified
    /// as offset (0-based) from the start of the record.
    pub(crate) line: u64,
    /// Byte offset from the start of the record
    pub(crate) byte: u64,
}

impl ErrorOffset {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn line(&self) -> u64 {
        self.line
    }

    pub fn byte(&self) -> u64 {
        self.byte
    }

    pub fn set_line(&mut self, line: u64) -> &mut Self {
        self.line = line;
        self
    }

    pub fn set_byte(&mut self, byte: u64) -> &mut Self {
        self.byte = byte;
        self
    }
}

/// Position of a parsing error within the file
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ErrorPosition {
    record_pos: Option<Position>,
    err_offset: Option<ErrorOffset>,
    id: Option<String>,
}

impl ErrorPosition {
    pub fn new(
        record_pos: Option<Position>,
        err_offset: Option<ErrorOffset>,
        id: Option<String>,
    ) -> Self {
        ErrorPosition {
            record_pos,
            err_offset,
            id,
        }
    }

    /// Returns the exact position, where the error occurred or the
    /// record start if there is no specific position.
    /// `None` can be returned if the position is not known because the error
    /// occurred in a call to `BaseRecord::check_lengths()`. Errors returned
    /// by any of the readers directly (not record instances) are guaranteed
    /// to always have a defined position.
    ///
    /// In contrast to the position returned by the function
    /// `ErrorPosition::record_position()`, an optional `ErrorOffset`
    /// is applied to the sequence record position
    /// (see `ErrorPosition::error_offset()`) to obtain
    /// the exact line and byte position of the error.
    #[inline]
    pub fn position(&self) -> Option<Position> {
        self.record_pos.as_ref().map(|pos| {
            if let Some(o) = self.error_offset() {
                return Position {
                    line: pos.line() + o.line(),
                    byte: pos.byte() + o.byte(),
                    record: pos.record(),
                };
            }
            pos.clone()
        })
    }

    /// Returns the ID of the record where the error occurred,
    /// if available. If parsing already failed within the sequence header,
    /// `None` is returned. Invalid UTF-8 bytes are replaced
    /// (see `String::from_utf8_lossy`).
    #[inline]
    pub fn record_id(&self) -> Option<&str> {
        self.id.as_ref().map(|s| s.as_str())
    }

    /// Returns a reference to the `Position` describing the
    /// line and byte offset of the sequence record, where the
    /// error occurred (if known, see `position()`).
    #[inline]
    pub fn record_position(&self) -> Option<&Position> {
        self.record_pos.as_ref()
    }

    /// Returns the offset (line / byte) from the record start,
    /// at which the error occurred, if available.
    #[inline]
    pub fn error_offset(&self) -> Option<&ErrorOffset> {
        self.err_offset.as_ref()
    }
}

impl fmt::Display for ErrorPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(id) = self.id.as_ref() {
            write!(f, "record '{}'", id)?;
        }
        if let Some(pos) = self.record_pos.as_ref() {
            write!(f, " at line {}", pos.line() + 1)?;
        }
        Ok(())
    }
}

macro_rules! impl_error {
    ($ErrorKind:ident) => {
        pub type Result<T> = std::result::Result<T, Error>;

        /// Parsing error
        #[derive(Debug)]
        pub struct Error {
            kind: Box<$ErrorKind>,
        }

        impl Error {
            #[inline]
            pub fn new(kind: $ErrorKind) -> Self {
                Error {
                    kind: Box::new(kind),
                }
            }

            /// Returns a reference to the [`ErrorKind`](ErrorKind)
            /// associated with the error.
            #[inline]
            pub fn kind(&self) -> &$ErrorKind {
                &self.kind
            }

            /// Returns the [`ErrorKind`](ErrorKind) associated with
            /// the error, thereby consuming the error.
            #[inline]
            pub fn into_kind(self) -> $ErrorKind {
                *self.kind
            }

            /// Returns the [`ErrorPosition`](ErrorPosition) of
            /// the error within file.
            #[inline]
            pub fn position(&self) -> Option<&ErrorPosition> {
                self.kind().position()
            }
        }

        impl From<std::io::Error> for Error {
            fn from(e: std::io::Error) -> Error {
                Error {
                    kind: Box::new($ErrorKind::Io(e)),
                }
            }
        }

        impl From<Error> for std::io::Error {
            fn from(err: Error) -> std::io::Error {
                std::io::Error::new(std::io::ErrorKind::Other, err)
            }
        }

        impl std::error::Error for Error {
            fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
                match self.kind() {
                    $ErrorKind::Io(ref err) => Some(err),
                    _ => None,
                }
            }
        }
    };
}

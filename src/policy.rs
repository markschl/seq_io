//! This module defines the [`BufPolicy`](BufPolicy) trait,
//! which configures how the internal buffer of the parsers should grow upon
//! encountering large sequences that don't fit into the buffer.
//!
//! The standard policy ([`DoubleUntil`](DoubleUntil))
//! causes the initial buffer to double its size until it reaches 32 MiB, then
//! further grow linearly until it reaches a size of 1 GiB. Larger sequences
//! will cause an error of `ErrorKind::BufferLimit` in readers.
//!
//! These numbers can be changed by using [`DoubleUntilLimited`](DoubleUntilLimited),
//! or the limit can be completely removed by using [`DoubleUntil`](DoubleUntil).
//!
//! # Example
//!
//! The following code applies a policy that lets the buffer double until it
//! reaches 64 MiB, and then further grow linearly without any limit:
//!
//! ```no_run
//! use seq_io::fasta::Reader;
//! use seq_io::policy::DoubleUntil;
//! use std::io::stdin;
//!
//! let policy = DoubleUntil(64 * 1024 * 1024);
//! let mut reader = Reader::with_buf_policy(stdin(), policy);
//! // (...)
//! ```
//!
//! # Custom policy
//!
//! This example shows how to implement a custom buffer policy:
//!
//! ```no_run
//! # extern crate seq_io;
//! # fn main() {
//! use seq_io::prelude::*;
//! use seq_io::policy::BufPolicy;
//! use seq_io::fasta::Reader;
//! use std::io::stdin;
//!
//! struct Max1G;
//!
//! // This policy lets the buffer double each time, but
//! // limits the buffer size to 1 GiB. Note that this is similar to how
//! // `DoubleUntilLimited` works.
//! impl BufPolicy for Max1G {
//!     fn grow(&mut self, current_size: usize) -> usize {
//!         current_size * 2
//!     }
//!
//!     fn limit(&self) -> Option<usize> {
//!         Some(1 << 30)
//!     }
//! }
//!
//! let mut reader = Reader::with_buf_policy(stdin(), Max1G);
//!
//! while let Some(record) = reader.next() {
//!     println!("{}", record.unwrap().id().unwrap());
//! }
//! # }
//! ```

/// Policy that configures how the internal buffer grows upon
/// encountering large sequences that don't fit into the current buffer.
pub trait BufPolicy: Send + Sync {
    /// Takes the current buffer size in bytes and returns the new
    /// size the the buffer should grow to. This function is called every time
    /// the buffer has to be enlarged.
    fn grow(&mut self, current_size: usize) -> usize;

    /// Returns a buffer limit, if any. Called every time the buffer has to be
    /// enlarged. If the new buffer size (as calculated based on the call to
    /// `grow()`) exceeds the given limit, the readers will return an error
    /// of `ErrorKind::BufferLimit`.
    fn limit(&self) -> Option<usize> {
        None
    }

    /// Combines `grow()` and `limit()` into one call. Takes the current buffer
    /// size and returns the new size, unless it is larger than the limit.
    fn grow_limited(&mut self, current_size: usize) -> Option<usize> {
        let new_size = self.grow(current_size);
        if let Some(l) = self.limit() {
            if new_size > l {
                return None;
            }
        }
        Some(new_size)
    }
}

/// Standard buffer policy: This policy corresponds to
/// `DoubleUntilLimited(1 << 25, 1 << 30)`, meaning that buffer size
/// doubles until it reaches 32 MiB, then increases in steps of 32 MiB until
/// it reaches 1 GiB, which is the limit.
pub struct StdPolicy;

impl BufPolicy for StdPolicy {
    fn grow(&mut self, current_size: usize) -> usize {
        if current_size < 1 << 25 {
            current_size * 2
        } else {
            current_size + (1 << 25)
        }
    }

    fn limit(&self) -> Option<usize> {
        Some(1 << 30)
    }
}

/// Buffer size doubles until it reaches a given limit
/// (in bytes). Above, it will increase linearly in
/// steps of 'limit'. Buffer size is not limited,
/// it could theoretically grow indefinitely.
pub struct DoubleUntil(pub usize);

impl BufPolicy for DoubleUntil {
    fn grow(&mut self, current_size: usize) -> usize {
        if current_size < self.0 {
            current_size * 2
        } else {
            current_size + self.0
        }
    }
}

/// Buffer size doubles until it reaches a given limit
/// (in bytes). Above, it will increase linearly in
/// steps of 'double_until'. Buffer size is additionally
/// limited to `limit` bytes. Readers will return an error
/// if this limit is .
pub struct DoubleUntilLimited {
    double_until: usize,
    limit: usize,
}

impl DoubleUntilLimited {
    pub fn new(double_until: usize, limit: usize) -> Self {
        DoubleUntilLimited {
            double_until,
            limit,
        }
    }
}

impl BufPolicy for DoubleUntilLimited {
    fn grow(&mut self, current_size: usize) -> usize {
        if current_size < self.double_until {
            current_size * 2
        } else {
            current_size + self.double_until
        }
    }

    fn limit(&self) -> Option<usize> {
        Some(self.limit)
    }
}

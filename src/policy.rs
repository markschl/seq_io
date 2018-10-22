/// Policy that decides how a buffer should grow
///
/// Takes the current buffer size in bytes and returns the new
/// size the the buffer should grow to. Returning `None` instead will indicate
/// that the buffer has grown too big. In this case, the FASTA and FASTQ readers
/// will return `Error::BufferLimit`.
///
/// # Example
///
/// ```no_run
/// # extern crate seq_io;
/// # fn main() {
/// use seq_io::BufPolicy;
/// use seq_io::fasta::{Reader,Record};
/// use std::io::stdin;
///
/// struct Max1G;
///
/// // This policy lets the buffer double each time, but
/// // limits the buffer size to 1 GiB. Note that this is similar to how
/// // `DoubleUntilLimited` works.
/// impl BufPolicy for Max1G {
///     fn grow_to(&mut self, current_size: usize) -> Option<usize> {
///         if current_size > 1 << 30 {
///             return None
///         }
///         Some(current_size * 2)
///     }
/// }
///
/// let mut reader = Reader::new(stdin()).set_policy(Max1G);
///
/// while let Some(record) = reader.next() {
///     println!("{}", record.unwrap().id().unwrap());
/// }
/// # }
/// ```
pub trait BufPolicy {
    fn grow_to(&mut self, current_size: usize) -> Option<usize>;
}

/// Standard buffer policy: This policy corresponds to
/// `DoubleUntil(8 * 1024 * 1024)`, meaning that buffer size
/// doubles until it reaches 8 MiB. Above, it will
/// increase in steps of 8 MiB. Buffer size is not limited,
/// it could theoretically grow indefinitely.
pub struct StdPolicy;

impl BufPolicy for StdPolicy {
    fn grow_to(&mut self, current_size: usize) -> Option<usize> {
        Some(if current_size < 1 << 23 {
            current_size * 2
        } else {
            current_size + (1 << 23)
        })
    }
}

/// Buffer size doubles until it reaches a given limit
/// (in bytes). Above, it will increase linearly in
/// steps of 'limit'. Buffer size is not limited,
/// it could theoretically grow indefinitely.
pub struct DoubleUntil(pub usize);

impl BufPolicy for DoubleUntil {
    fn grow_to(&mut self, current_size: usize) -> Option<usize> {
        Some(if current_size < self.0 {
            current_size * 2
        } else {
            current_size + self.0
        })
    }
}

/// Buffer size doubles until it reaches a given limit
/// (in bytes). Above, it will increase linearly in
/// steps of 'double_until'. Buffer size is additionally
/// limited to `limit` bytes. Readers will return an error
/// if this limit is .
pub struct DoubleUntilLimited {
    double_until: usize,
    limit: usize
}

impl DoubleUntilLimited {
    pub fn new(double_until: usize, limit: usize) -> Self {
        DoubleUntilLimited { double_until, limit }
    }
}

impl BufPolicy for DoubleUntilLimited {
    fn grow_to(&mut self, current_size: usize) -> Option<usize> {
        let new_size = if current_size < self.double_until {
            current_size * 2
        } else {
            current_size + self.double_until
        };
        if new_size <= self.limit {
            Some(new_size)
        } else {
            None
        }
    }
}

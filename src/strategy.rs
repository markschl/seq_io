/// Strategy that decides how a buffer should grow
///
/// Returns the number of additional bytes given the
/// current size. Returning None instead will indicate
/// that the buffer has grown too big.
/// Creates a new reader with a given buffer capacity and growth strategy
///
/// # Example
///
/// ```no_run
/// # extern crate seq_io;
/// # fn main() {
/// use seq_io::BufStrategy;
/// use seq_io::fasta::{Reader,Record};
/// use std::io::stdin;
///
/// pub struct Max1G;
///
/// // This BufStrategy limits the buffer size to 1 GB
/// impl BufStrategy for Max1G {
///     fn grow_to(&mut self, current_size: usize) -> Option<usize> {
///         if current_size > 1 << 30 {
///             return None
///         }
///         Some(current_size * 2)
///     }
/// }
///
/// let mut reader = Reader::with_cap_and_strategy(stdin(), 68 * 1024, Max1G);
/// while let Some(record) = reader.next() {
///     println!("{}", record.unwrap().id().unwrap());
/// }
/// # }
/// ```
pub trait BufStrategy {
    fn grow_to(&mut self, current_size: usize) -> Option<usize>;
}

/// Buffer size doubles until it
/// reaches 8 MB. Above, it will
/// increase in steps of 8 MB
pub struct DoubleUntil8M;

impl BufStrategy for DoubleUntil8M {
    fn grow_to(&mut self, current_size: usize) -> Option<usize> {
        Some(if current_size < 1 << 23 {
            current_size * 2
        } else {
            current_size + 1 << 23
        })
    }
}


/// Buffer size doubles until it reaches
/// `double_size_limit` (in bytes). Above,
/// it increases in steps of `double_size_limit`
pub struct DoubleUntil(pub usize);

impl BufStrategy for DoubleUntil {
    fn grow_to(&mut self, current_size: usize) -> Option<usize> {
        Some(if current_size < self.0 {
            current_size * 2
        } else {
            current_size + self.0
        })
    }
}

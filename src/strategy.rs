/// Strategy that decides how a buffer should grow
///
/// Returns the number of additional bytes given the
/// current size. Returning None instead will indicate
/// that the buffer has grown too big.
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

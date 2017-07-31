/// Strategy that decides how a buffer should grow
///
/// Returns the number of additional bytes given the
/// current size. It's implementors can also choose to
/// panic at some point.
/// **TODO:** should Result be used instead?
pub trait BufGrowStrategy {
    fn new_size(&self, current_size: usize) -> usize;
}

/// Buffer size doubles until it
/// reaches 8 MB. Above, it will
/// increase in steps of 8 MB
pub struct DoubleUntil8M;

impl BufGrowStrategy for DoubleUntil8M {
    fn new_size(&self, current_size: usize) -> usize {
        if current_size < 1 << 23 {
            return current_size * 2;
        }
        current_size + 1 << 23
    }
}


/// Buffer size doubles until it reaches
/// `double_size_limit` (in bytes). Above,
/// it increases in steps of `double_size_limit`
pub struct DoubleUntil(usize);

impl DoubleUntil {
    pub fn new(double_size_limit: usize) -> DoubleUntil {
        DoubleUntil(double_size_limit)
    }
}

impl BufGrowStrategy for DoubleUntil {
    fn new_size(&self, current_size: usize) -> usize {
        if current_size < self.0 {
            return current_size * 2;
        }
        current_size + self.0
    }
}

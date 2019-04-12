use crate::policy::{BufPolicy, StdPolicy};
use std::fs::File;
use std::io::{self, BufRead, Seek};
use std::path::Path;

/// The default initial buffer size for readers.
pub const BUFSIZE: usize = 64 * 1024;

/// Wraps buf_redux::BufReader, managing buffer growth
/// based on BufPolicy and the relocation of buffer contents
/// with make_room(). Does not implement `std::io::BufRead`
pub struct BufReader<R, P = StdPolicy>
where
    R: std::io::Read,
{
    buf_reader: buf_redux::BufReader<R>,
    // Buffer resizing policy
    buf_policy: P,
    // Search position within file in total (line and byte offset)
    file_offset: u64,
}

impl<P> BufReader<File, P>
where
    P: BufPolicy,
{
    /// Creates a buffered reader from a file path.
    #[inline]
    pub fn from_path<F: AsRef<Path>>(path: F) -> io::Result<BufReader<File>> {
        File::open(path).map(BufReader::new)
    }
}

impl<R> BufReader<R>
where
    R: io::Read,
{
    #[inline]
    pub fn new(reader: R) -> Self {
        Self::with_capacity(reader, BUFSIZE)
    }

    #[inline]
    pub fn with_capacity(reader: R, capacity: usize) -> Self {
        BufReader {
            buf_reader: buf_redux::BufReader::with_capacity(capacity, reader),
            buf_policy: StdPolicy,
            file_offset: 0,
        }
    }
}

impl<R, P> BufReader<R, P>
where
    R: io::Read,
    P: BufPolicy,
{
    #[inline]
    pub fn policy(&self) -> &P {
        &self.buf_policy
    }

    #[inline]
    pub fn set_policy<T: BufPolicy>(self, buf_policy: T) -> BufReader<R, T> {
        BufReader {
            buf_reader: self.buf_reader,
            buf_policy,
            file_offset: self.file_offset,
        }
    }

    #[inline]
    pub fn buffer(&self) -> &[u8] {
        self.buf_reader.buffer()
    }

    #[inline]
    pub fn capacity(&self) -> usize {
        self.buf_reader.capacity()
    }

    #[inline]
    pub fn file_offset(&self) -> u64 {
        self.file_offset
    }

    // grow buffer based on policy
    #[inline]
    pub fn grow(&mut self) -> bool {
        let cap = self.buf_reader.capacity();
        if let Some(new_size) = self.buf_policy.grow_limited(cap) {
            //println!("grow {} -> {}", cap, new_size);
            let additional = new_size - cap;
            assert!(additional > 0);
            self.buf_reader.reserve(additional);
            return true;
        }
        false
    }

    // move incomplete bytes to start of buffer and retry
    #[inline]
    pub fn make_room(&mut self, offset: usize) {
        self.file_offset += offset as u64;
        //println!("make room -> {} at bufsize {} / {}", self.file_offset, self.buf_reader.buffer().len(), self.buf_reader.capacity());
        self.buf_reader.consume(offset as usize);
        self.buf_reader.make_room();
    }

    /// Makes sure the buffer is full after this call (unless EOF reached)
    /// code adapted from `io::Read::read_exact`
    #[inline]
    pub fn fill_buf(&mut self) -> io::Result<usize> {
        let initial_size = self.buf_reader.buffer().len();
        let mut num_read = 0;
        while initial_size + num_read < self.buf_reader.capacity() {
            match self.buf_reader.read_into_buf() {
                Ok(0) => break,
                Ok(n) => num_read += n,
                Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(num_read)
    }
}

impl<R, P> BufReader<R, P>
where
    R: io::Read + std::io::Seek,
    P: BufPolicy,
{
    #[inline]
    pub fn seek_to(&mut self, pos: u64) -> io::Result<()> {
        self.file_offset = pos;
        self.buf_reader.seek(io::SeekFrom::Start(pos))?;
        Ok(())
    }
}

//! Efficient FASTA reading and writing
//!

use std::io::{self,BufRead};
use std::fs::File;
use std::path::Path;
use std::slice;
use std::iter;
use std::str::{self,Utf8Error};

use memchr::memchr;
use buf_redux;

use super::*;


type DefaultBufGrowStrategy = DoubleUntil8M;


const BUFSIZE: usize = 68 * 1024;


/// Parser for FASTA files.
pub struct Reader<R: io::Read, S = DefaultBufGrowStrategy> {
    buffer: buf_redux::BufReader<R, ReadAlways, buf_redux::strategy::NeverMove>,
    position: RecordPosition,
    n_searched: usize,
    finished: bool,
    grow_strategy: S,
}


impl<R> Reader<R, DefaultBufGrowStrategy>
    where R: io::Read
{
    pub fn new(reader: R) -> Reader<R, DoubleUntil8M> {
        Reader::with_cap_and_strategy(reader, BUFSIZE, DoubleUntil8M)
    }

    pub fn with_capacity(reader: R, capacity: usize) -> Reader<R, DoubleUntil8M> {
        Reader::with_cap_and_strategy(reader, capacity, DoubleUntil8M)
    }
}

impl Reader<File, DefaultBufGrowStrategy> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Reader<File>> {
        File::open(path).map(Reader::new)
    }
}

impl<R, S> Reader<R, S>
    where R: io::Read,
          S: BufGrowStrategy
{

    #[inline]
    pub fn with_cap_and_strategy(reader: R, cap: usize, grow_strategy: S) -> Reader<R, S> {
        assert!(cap >= 3);
        Reader {
            buffer: buf_redux::BufReader::with_cap_and_strategies(
                reader, cap,
                ReadAlways, buf_redux::strategy::NeverMove
            ),
            position: RecordPosition {
                start: 0,
                seq_starts: Vec::with_capacity(2),
            },
            n_searched: 0,
            finished: false,
            grow_strategy: grow_strategy,
        }
    }

    #[inline]
    pub fn proceed(&mut self) -> Option<Result<(), ParseError>> {

        if self.finished {
            return None;
        }

        if ! try_opt!(self.read_next()) {
            try_opt!(self.next_complete());
        }

        Some(Ok(()))
    }

    /// Search the next FASTA record and return a `RefRecord` that
    /// borrows it's data from the underlying buffer of this reader
    pub fn next<'a>(&'a mut self) -> Option<Result<RefRecord<'a>, ParseError>> {
        self.proceed().map(|r| r.map(
            move |_| RefRecord { buffer: self.get_buf(), position: &self.position }
        ))
    }

    /// Updates a `RecordSet` with a new buffer and searches for records. Old data will be erased.
    /// Returns `None` if the input reached its end
    pub fn read_record_set(&mut self, rset: &mut RecordSet) -> Option<Result<(), ParseError>> {

        // note: remember to always check for 'self.finished', or subsequent calls to
        // next_complete will have an incorrect 'next_start' prop

        if self.finished {
            return None;
        }

        try_opt!(self.next_complete());

        // copy buffer AFTER call to next_complete (initialization of buffer is done there)
        rset.buffer.clear();
        rset.buffer.extend(self.get_buf());

        // Update records that are already in the positions vector
        let mut n = 0;
        for pos in rset.positions.iter_mut() {
            n += 1;
            pos.update(&self.position);

            if self.finished || ! try_opt!(self.read_next()) {
                rset.npos = n;
                return Some(Ok(()));
            }
        }

        // Add more positions if necessary
        loop {
            n += 1;
            rset.positions.push(self.position.clone());

            if self.finished || ! try_opt!(self.read_next()) {
                rset.npos = n;
                return Some(Ok(()));
            }
        }
    }

    #[inline]
    fn read_next(&mut self) -> Result<bool, ParseError> {
        // initialize (clear)
        self.position.start = self.n_searched;
        self.position.seq_starts.clear();
        self.search()
    }

    #[inline]
    fn get_buf(&self) -> &[u8] {
        self.buffer.get_buf()
    }

    /// Finds the position of the next record
    /// and returns true if found; false if end of buffer reached.
    #[inline]
    fn search(&mut self) -> Result<bool, ParseError> {

        let bufsize = self.get_buf().len();

        while let Some(pos) = memchr(b'\n', &self.get_buf()[self.n_searched .. ]) {
            // don't search the last byte, since we need to look forward one byte, looking for the next record
            let pos = self.n_searched + pos;
            let line_start = pos + 1;

            // start of next record is not actually a sequence start, but indicates end of Record
            self.n_searched = line_start;

            if line_start == bufsize {
                // cannot check next byte -> treat as incomplete
                self.n_searched -= 1;  // make sure last byte is re-searched next time
                return self.check_end();
            }

            self.position.seq_starts.push(line_start);
            if self.position.seq_starts.len() > 1 && self.get_buf()[line_start] == b'>' {
                // next found
                return Ok(true);
            }
        }

        // record end not found
        self.n_searched = bufsize;

        self.check_end()
    }

    // checks if the buffer size is smaller than expected, indicating EOF
    #[inline]
    fn check_end(&mut self) -> Result<bool, ParseError> {
        let bufsize = self.get_buf().len();
        if bufsize < self.buffer.capacity() && bufsize > 0 { // bufsize == 0 means that init() has not yet been executed
            // EOF reached, there will be no next record
            self.finished = true;
            self.position.seq_starts.push(bufsize);
            if self.position.seq_starts.len() == 1 {
                return Err(ParseError::UnexpectedEnd);
            }
            return Ok(true);
        }
        Ok(false)
    }

    // moves to the first record posiion, ignoring newline / carriage return characters
    #[inline]
    fn init(&mut self) -> Result<(), ParseError> {
        loop {
            let n = fill_buf(&mut self.buffer)?;
            if n == 0 {
                // TODO: should this really be an ParseError?
                self.finished = true;
                return Err(ParseError::EmptyInput);
            }
            if let Some((i, c)) = self.buffer.get_buf()
                                      .iter()
                                      .enumerate()
                                      .skip_while(|&(_, c)| *c == b'\r' || *c == b'\n')
                                      .next() {
                if *c != b'>' {
                    self.finished = true;
                    return Err(ParseError::InvalidStart(i));
                }
                self.position.start = i;
                return Ok(());
            }
            // whole buffer consists of newlines (unlikely)
            let cap = self.buffer.capacity();
            self.buffer.consume(cap);
        }
    }

    /// To be called when the end of the buffer is reached and `next_pos` does not find
    /// the next record. Incomplete bytes will be moved to the start of the buffer.
    /// If the record still doesn't fit in, the buffer will be enlarged.
    /// After calling this function, the position will therefore always be 'complete'.
    /// this function assumes that the buffer was fully searched
    #[inline]
    fn next_complete(&mut self) -> Result<(), ParseError> {

        loop {

            if self.get_buf().len() == 0 {
                // not yet initialized
                self.init()?;

            } else if self.position.start == 0 {
                // first record -> buffer too small, grow until big enough
                let cap = self.buffer.capacity();
                let new_size = self.grow_strategy.new_size(cap).ok_or(ParseError::BufferOverflow)?;
                let additional = new_size - cap;
                self.buffer.grow(additional);

            } else {
                // not the first record -> buffer may be big enough
                // move incomplete bytes to start of buffer and retry
                let consumed = self.position.start;
                self.buffer.consume(consumed);
                self.buffer.make_room();
                self.position.start = 0;
                self.n_searched -= consumed;
                for s in &mut self.position.seq_starts {
                    *s -= consumed;
                }
            }

            // fill up remaining buffer
            fill_buf(&mut self.buffer)?;

            if self.search()? {
                return Ok(());
            }
        }
    }

    // TODO: these methods are shared with RefRecord -> another trait?

    pub fn seq_lines(&self) -> SeqLines {
        self.position.seq_lines(self.get_buf())
    }

    pub fn owned_seq(&self) -> Vec<u8> {
        self.position.owned_seq(self.get_buf())
    }

    pub fn to_owned_record(&self) -> OwnedRecord {
        self.position.get_owned_record(self.get_buf())
    }

    pub fn write_unchanged<W: io::Write>(&self, writer: &mut W) -> io::Result<usize> {
        self.position.write(writer, self.get_buf())
    }
}


#[derive(Debug)]
pub enum ParseError {
    Io(io::Error),
    EmptyInput,
    InvalidStart(usize),
    UnexpectedEnd,
    BufferOverflow,
}


impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ParseError::Io(ref e) => write!(f, "{}", e),
            ParseError::EmptyInput => write!(f, "Input empty"),
            ParseError::InvalidStart(pos) => write!(f, "Expected > at file start (position: {})", pos),
            ParseError::UnexpectedEnd => write!(f, "Unexpected end of input"),
            ParseError::BufferOverflow => write!(f, "Buffer overflow"),
        }
    }
}

impl From<io::Error> for ParseError {
    fn from(e: io::Error) -> ParseError {
        ParseError::Io(e)
    }
}

impl error::Error for ParseError {
    fn description(&self) -> &str { "FASTA parsing ParseError" }
}




#[derive(Clone, Debug)]
struct RecordPosition {
    start: usize,
    seq_starts: Vec<usize>,
}

impl RecordPosition {
    #[inline]
    fn head<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.start + 1 .. *self.seq_starts.first().unwrap() - 1])
    }

    #[inline]
    fn seq<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[*self.seq_starts.first().unwrap() .. *self.seq_starts.last().unwrap() - 1])
    }

    #[inline]
    fn write<W: io::Write>(&self, writer: &mut W, buffer: &[u8]) -> io::Result<usize> {
        let data = &buffer[self.start .. *self.seq_starts.last().unwrap()];
        writer.write_all(data)?;
        Ok(data.len())
    }

    #[inline]
    fn update(&mut self, other: &Self) {
        self.start = other.start;
        self.seq_starts.clear();
        self.seq_starts.extend(&other.seq_starts);
    }

    #[inline]
    fn seq_lines<'a>(&'a self, buffer: &'a [u8]) -> SeqLines<'a> {
        SeqLines {
            data: buffer,
            pos_iter: self.seq_starts.iter().zip(self.seq_starts.iter().skip(1))
        }
    }

    #[inline]
    fn owned_seq(&self, buffer: &[u8]) -> Vec<u8> {
        let mut seq = Vec::new();
        for segment in self.seq_lines(buffer) {
            seq.extend(segment);
        }
        return seq;
    }

    #[inline]
    pub fn get_owned_record(&self, buffer: &[u8]) -> OwnedRecord {
        OwnedRecord {
            head: self.head(buffer).to_vec(),
            seq: self.owned_seq(buffer),
        }
    }
}


pub trait Record {
    /// Return the header line of the record as byte slice
    fn head(&self) -> &[u8];
    /// Return the FASTA sequence as byte slice
    fn seq(&self) -> &[u8];

    fn id_bytes(&self) -> &[u8] {
        self.head().split(|b| *b == b' ').next().unwrap()
    }

    /// Return the ID of the record (everything before an optional space) as string slice
    fn id(&self) -> Result<&str, Utf8Error> {
        str::from_utf8(self.id_bytes())
    }

    fn desc_bytes(&self) -> Option<&[u8]> {
        self.head().splitn(2, |b| *b == b' ').nth(1)
    }

    /// Return the description of the record as string slice, if present. Otherwise, `None` is returned.
    fn desc(&self) -> Option<Result<&str, Utf8Error>> {
        self.desc_bytes().map(str::from_utf8)
    }

    /// Return both the ID and the description of the record (if present)
    /// This should be faster than calling `id()` and `desc()` separately.
    fn id_desc_bytes(&self) -> (&[u8], Option<&[u8]>) {
        let mut h = self.head().splitn(2, |c| *c == b' ');
        (h.next().unwrap(), h.next())
    }

    /// Return both the ID and the description of the record (if present)
    /// This should be faster than calling `id()` and `desc()` separately.
    fn id_desc(&self) -> Result<(&str, Option<&str>), Utf8Error> {
        let mut h = str::from_utf8(self.head())?.splitn(2, ' ');
        Ok((h.next().unwrap(), h.next()))
    }

    fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<usize> {
        write_to(writer, self.head(), self.seq())
    }
}


impl<R, S> Record for Reader<R, S>
    where R: io::Read,
          S: BufGrowStrategy
{
    #[inline]
    fn head(&self) -> &[u8] {
        self.position.head(self.get_buf())
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        self.position.seq(self.get_buf())
    }
}



/// A FASTA record that borrows data from a buffer.
#[derive(Debug, Clone)]
pub struct RefRecord<'a> {
    buffer: &'a [u8],
    position: &'a RecordPosition,
}


impl<'a> Record for RefRecord<'a> {

    #[inline]
    fn head(&self) -> &[u8] {
        self.position.head(self.buffer)
    }

    /// Return the FASTA sequence as byte slice.
    /// Note that this method of `RefRecord` returns
    /// the **raw** sequence, which may contain line breaks.
    /// Use `seq_lines()` to iterate over all lines witout
    /// breaks, or collect into an owned sequence using `owned_seq()`
    #[inline]
    fn seq(&self) -> &[u8] {
        self.position.seq(self.buffer)
    }
}



impl<'a> RefRecord<'a> {
    /// Return an iterator over all sequence lines in the data
    pub fn seq_lines(&self) -> SeqLines {
        self.position.seq_lines(self.buffer)
    }

    /// Returns the sequence as owned `Vec`. **Note**: This function
    /// must be called in order to obtain a sequence that does not contain
    /// line endings (as returned by `seq()`)
    pub fn owned_seq(&self) -> Vec<u8> {
        self.position.owned_seq(self.buffer)
    }

    /// Creates an owned copy of the record.
    pub fn to_owned_record(&self) -> OwnedRecord {
        self.position.get_owned_record(self.buffer)
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `RefRecord::write`
    pub fn write_unchanged<W: io::Write>(&self, writer: &mut W) -> io::Result<usize> {
        self.position.write(writer, self.buffer)
    }
}


pub struct SeqLines<'a> {
    data: &'a [u8],
    pos_iter: iter::Zip<slice::Iter<'a, usize>, iter::Skip<slice::Iter<'a, usize>>>,
}

impl<'a> Iterator for SeqLines<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<&'a [u8]> {
        self.pos_iter.next().map(
            |(start, next_start)| trim_cr(&self.data[*start .. *next_start - 1])
        )
    }
}


/// A FASTA record that ownes its data (requiring two allocations)
#[derive(Debug, Clone)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
}


impl Record for OwnedRecord {
    #[inline]
    fn head(&self) -> &[u8] {
        &self.head
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        &self.seq
    }
}



/// Set of FASTA records that owns it's buffer
/// and knows the positions of each record.
#[derive(Clone, Debug)]
pub struct RecordSet {
    buffer: Vec<u8>,
    positions: Vec<RecordPosition>,
    npos: usize,
}

impl Default for RecordSet {
    fn default() -> RecordSet {
        RecordSet {
            buffer: vec![],
            positions: vec![],
            npos: 0,
        }
    }
}



impl<'a> iter::IntoIterator for &'a RecordSet {
    type Item = RefRecord<'a>;
    type IntoIter = RecordSetIter<'a>;
     fn into_iter(self) -> Self::IntoIter {
         RecordSetIter {
             buffer: &self.buffer,
             pos: self.positions.iter().take(self.npos),
         }
     }
}


pub struct RecordSetIter<'a> {
    buffer: &'a [u8],
    pos: iter::Take<slice::Iter<'a, RecordPosition>>
}

impl<'a> Iterator for RecordSetIter<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<RefRecord<'a>> {
        self.pos.next().map(|p| {
            RefRecord {
                buffer: self.buffer,
                position: p,
            }
        })
    }
}


/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTA format
pub fn write_to<W: io::Write>(writer: &mut W, head: &[u8], seq: &[u8]) -> io::Result<usize> {
    let mut n = 0;
    n += writer.write(b">")?;
    n += writer.write(head)?;
    n += writer.write(b"\n")?;
    n += writer.write(seq)?;
    n += writer.write(b"\n")?;
    Ok(n)
}

/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTA format. The ID and description parts of the header are supplied separately
/// instead of a whole header line
pub fn write_parts<W: io::Write>(writer: &mut W, id: &[u8], desc: Option<&[u8]>, seq: &[u8]) -> io::Result<usize> {
    let mut n = 0;
    n += writer.write(b">")?;
    n += writer.write(id)?;
    if let Some(d) = desc {
        n += writer.write(b" ")?;
        n += writer.write(d)?;
    }
    n += writer.write(b"\n")?;
    n += writer.write(seq)?;
    n += writer.write(b"\n")?;
    Ok(n)
}

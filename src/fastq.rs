//! Efficient FASTQ reading and writing
//!

use std::io::{self,BufRead};
use std::fs::File;
use std::path::Path;
use memchr::memchr;
use std::slice;
use std::iter;
use std::char;
use std::str::{self,Utf8Error};

use buf_redux;

use super::*;


type DefaultBufGrowStrategy = DoubleUntil8M;

const BUFSIZE: usize = 68 * 1024;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
enum CursorPos {
    HEAD,
    SEQ,
    SEP,
    QUAL
}


/// FASTQ parser.
pub struct Reader<R: io::Read, S = DefaultBufGrowStrategy> {
    buffer: buf_redux::BufReader<R, ReadAlways, buf_redux::strategy::NeverMove>,
    position: RecordPosition,
    cursor_pos: CursorPos,
    finished: bool,
    grow_strategy: S,
}


impl<R> Reader<R, DefaultBufGrowStrategy>
    where R: io::Read
{
    /// Creates a new reader with the default buffer size of 68 KB
    pub fn new(reader: R) -> Reader<R, DoubleUntil8M> {
        Reader::with_cap_and_strategy(reader, BUFSIZE, DoubleUntil8M)
    }

    /// Creates a reader with the given buffer size
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

    /// Creates a reader given buffer size (cap) and growth strategy
    #[inline]
    pub fn with_cap_and_strategy(reader: R, cap: usize, grow_strategy: S) -> Reader<R, S> {
        assert!(cap >= 3);
        Reader {
            buffer: buf_redux::BufReader::with_cap_and_strategies(
                reader, //unsafe {buf_redux::AssertTrustRead::new(reader)},
                cap,
                ReadAlways, buf_redux::strategy::NeverMove
            ),
            position: RecordPosition::default(),
            cursor_pos: CursorPos::HEAD,
            finished: false,
            grow_strategy: grow_strategy,
        }
    }

    #[inline]
    pub fn proceed(&mut self) -> Option<Result<(), ParseError>> {
        if ! try_opt!(self.find_next()) {
            if try_opt!(self.next_complete()).is_none() {
                return None;
            }
        }
        Some(Ok(()))
    }

    /// Search the next FASTQ record and return a `RefRecord` that
    /// borrows it's data from the underlying buffer of this reader
    #[inline]
    pub fn next<'a>(&'a mut self) -> Option<Result<RefRecord<'a>, ParseError>> {
        self.proceed().map(|r| r.map(
            move |_| RefRecord { buffer: self.get_buf(), position: &self.position }
        ))
    }

    #[inline]
    fn init(&mut self) -> Result<(), ParseError> {
        let n = fill_buf(&mut self.buffer)?;
        if n == 0 {
            // TODO: should this really be an ParseError?
            self.finished = true;
            return Err(ParseError::EmptyInput);
        }
        self.require_byte(0, b'@')
    }

    /// Updates a `RecordSet` with a new buffer and new positions.
    /// Present data will be erased. Returns `None` if the input reached its end
    pub fn read_record_set(&mut self, rset: &mut RecordSet) -> Option<Result<(), ParseError>> {

        if try_opt!(self.next_complete()).is_none() {
            return None;
        }

        // fill buffer AFTER call to next_complete(); initialization of buffer is done there
        rset.buffer.clear();
        rset.buffer.extend(self.get_buf());

        rset.positions.clear();
        rset.positions.push(self.position.clone());

        while try_opt!(self.find_next()) {
            rset.positions.push(self.position.clone());
        }

        Some(Ok(()))
    }


    #[inline(always)]
    fn get_buf(&self) -> &[u8] {
        self.buffer.get_buf()
    }

    // Reads the current record  and returns true if found
    // (false if incomplete because end of buffer reached)
    // meaning that the last record may be incomplete
    // search_start >= self.position.start
    // updates the position
    #[inline(always)]
    fn find_next(&mut self) -> Result<bool, ParseError> {

        self.position.pos.0 = self.position.pos.1;

        self.position.seq = unwrap_or!(self.find_line(self.position.pos.0), {
            self.cursor_pos = CursorPos::HEAD;
            return Ok(false);
        });

        self.position.sep = unwrap_or!(self.find_line(self.position.seq), {
            self.cursor_pos = CursorPos::SEQ;
            return Ok(false);
        });

        self.position.qual = unwrap_or!(self.find_line(self.position.sep), {
            self.cursor_pos = CursorPos::SEP;
            return Ok(false);
        });

        self.position.pos.1 = unwrap_or!(self.find_line(self.position.qual), {
            self.cursor_pos = CursorPos::QUAL;
            return Ok(false);
        });

        self.validate()?;

        // caution: this fuction does not reset `self.cursor_pos` to `CursorPos::HEAD` =>
        // once it returned false, find_incomplete has to be called next, otherwise
        // there will be problems
        Ok(true)
    }

    // Continues reading an incomplete record without
    // re-searching positions that were already found.
    // The resulting position may still be incomplete.
    #[inline]
    fn find_incomplete(&mut self) -> Result<bool, ParseError> {

        if self.cursor_pos == CursorPos::HEAD {
            self.position.seq = unwrap_or!(self.find_line(self.position.pos.0), {
                self.cursor_pos = CursorPos::HEAD;
                return Ok(false);
            });
        }

        if self.cursor_pos <= CursorPos::SEQ {
            self.position.sep = unwrap_or!(self.find_line(self.position.seq), {
                self.cursor_pos = CursorPos::SEQ;
                return Ok(false);
            });
        }

        if self.cursor_pos <= CursorPos::SEP {
            self.position.qual = unwrap_or!(self.find_line(self.position.sep), {
                self.cursor_pos = CursorPos::SEP;
                return Ok(false);
            });
        }

        if self.cursor_pos <= CursorPos::QUAL {
            self.position.pos.1 = unwrap_or!(self.find_line(self.position.qual), {
                self.cursor_pos = CursorPos::QUAL;
                return Ok(false);
            });
        }

        self.cursor_pos = CursorPos::HEAD;

        self.validate()?;

        Ok(true)
    }

    // should only be called on a complete RecordPosition
    #[inline(always)]
    fn validate(&mut self) -> Result<(), ParseError> {

        let p0 = self.position.pos.0;
        self.require_byte(p0, b'@')?;
        let ps = self.position.sep;
        self.require_byte(ps, b'+')?;

        let qual_len = self.position.pos.1 - self.position.qual;
        let seq_len = self.position.sep - self.position.seq;
        if seq_len != qual_len {
            self.finished = true;
            return Err(ParseError::UnequalLengths(seq_len, qual_len));
        }
        Ok(())
    }

    #[inline(always)]
    fn find_line(&self, search_start: usize) -> Option<usize> {
        memchr(b'\n', &self.get_buf()[search_start.. ]).map(|pos| search_start + pos + 1)
    }

    #[inline(always)]
    fn require_byte(&mut self, pos: usize, expected: u8) -> Result<(), ParseError> {
        let found = self.get_buf()[pos];
        if found != expected {
            self.finished = true;
            return Err(ParseError::Unexpected(expected, found));
        }
        Ok(())
    }

    // To be called when the end of the buffer is reached and `next_pos` does not find
    // the next record. Incomplete bytes will be moved to the start of the buffer.
    // If the record still doesn't fit in, the buffer will be enlarged.
    // After calling this function, the position will therefore always be 'complete'.
    #[inline]
    fn next_complete(&mut self) -> Result<Option<()>, ParseError> {

        if self.get_buf().len() == 0 {
            // not yet initialized
            self.init()?;
            if self.find_next()? {
                return Ok(Some(()));
            }
        }

        loop {

            // this function assumes that the buffer was fully searched
            let bufsize = self.get_buf().len();
            if bufsize < self.buffer.capacity() {
                if self.finished {
                    return Ok(None);
                }
                // EOF reached, there will be no next record
                self.finished = true;
                if self.cursor_pos == CursorPos::QUAL {
                    // no line ending at end of last record
                    self.position.pos.1 = bufsize + 1;
                    self.validate()?;
                    return Ok(Some(()));
                }

                let rest = &self.get_buf()[self.position.pos.0..];
                if rest.split(|c| *c == b'\n').all(|l| trim_cr(l).len() == 0)  {
                    // allow up to 4 newlines after last record (more will cause an Unexpected error)
                    return Ok(None);
                }

                return Err(ParseError::UnexpectedEnd);

            } else if self.position.pos.0 == 0 {
                // first record already incomplete -> buffer too small, grow until big enough
                let cap = self.buffer.capacity();
                let new_size = self.grow_strategy.new_size(cap).ok_or(ParseError::BufferOverflow)?;
                let additional = new_size - cap;
                self.buffer.grow(additional);
                fill_buf(&mut self.buffer)?;

            } else {
                // not the first record -> buffer may be big enough
                // move incomplete bytes to start of buffer and retry
                let consumed = self.position.pos.0;
                self.buffer.consume(consumed);
                self.buffer.make_room();

                self.position.pos.0 = 0;

                if self.cursor_pos >= CursorPos::SEQ {
                    self.position.seq -= consumed;
                }
                if self.cursor_pos >= CursorPos::SEP {
                    self.position.sep -= consumed;
                }
                if self.cursor_pos >= CursorPos::QUAL {
                    self.position.qual -= consumed;
                }

                fill_buf(&mut self.buffer)?;
            }

            // self.position.pos.1 = 0;
            // self.cursor_pos = CursorPos::HEAD;
            if self.find_incomplete()? {
                return Ok(Some(()));
            }
        }
    }

    // TODO: duplicate code
    #[inline]
    pub fn get_owned_record(&self) -> OwnedRecord {
        self.position.get_owned_record(self.get_buf())
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `RefRecord::write`
    #[inline]
    pub fn write_unchanged<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        self.position.write_unchanged(writer, self.get_buf())
    }
}


#[derive(Debug)]
pub enum ParseError {
    Io(io::Error),
    EmptyInput,
    UnequalLengths(usize, usize),
    Unexpected(u8, u8),
    UnexpectedEnd,
    BufferOverflow,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            ParseError::Io(ref e) => write!(f, "{}", e),
            ParseError::EmptyInput => write!(f, "Input empty"),
            ParseError::UnequalLengths(seq, qual) => write!(f, "Unequal lengths: sequence length is {}, but quality length is {}", seq, qual),
            ParseError::Unexpected(exp, found) => write!(f, "Expected '{}' but found '{}'", exp as char, found as char),
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
    fn description(&self) -> &str { "FASTQ parsing ParseError" }
}



#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct RecordPosition {
    // (start, stop), but might include \r at the end
    pos: (usize, usize),
    seq: usize,
    sep: usize,
    qual: usize,
}

impl RecordPosition {

    #[inline]
    fn head<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.pos.0 + 1 .. self.seq - 1])
    }

    #[inline]
    fn seq<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.seq .. self.sep - 1])
    }

    #[inline]
    fn qual<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.qual .. self.pos.1 - 1])
    }

    #[inline]
    fn write_unchanged<W: io::Write>(&self, writer: &mut W, buffer: &[u8]) -> io::Result<()> {
        let data = &buffer[self.pos.0 .. self.pos.1];
        writer.write_all(data)?;
        Ok(())
    }

    #[inline]
    pub fn get_owned_record(&self, buffer: &[u8]) -> OwnedRecord {
        OwnedRecord {
            head: self.head(buffer).to_vec(),
            seq: self.seq(buffer).to_vec(),
            qual: self.qual(buffer).to_vec(),
        }
    }
}


pub trait Record {
    /// Return the header line of the record as byte slice
    fn head(&self) -> &[u8];
    /// Return the FASTA sequence as byte slice
    fn seq(&self) -> &[u8];
    /// Return the FASTA qualities as byte slice
    fn qual(&self) -> &[u8];

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

    /// Writes a record to the given `io::Write` instance
    #[inline]
    fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        write_to(writer, self.head(), self.seq(), self.qual())
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

    #[inline]
    fn qual(&self) -> &[u8] {
        self.position.qual(self.get_buf())
    }
}


/// A FASTQ record that borrows data from a buffer
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

    #[inline]
    fn seq(&self) -> &[u8] {
        self.position.seq(self.buffer)
    }

    #[inline]
    fn qual(&self) -> &[u8] {
        self.position.qual(self.buffer)
    }
}



impl<'a> RefRecord<'a> {
    #[inline]
    pub fn to_owned_record(&self) -> OwnedRecord {
        self.position.get_owned_record(self.buffer)
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `RefRecord::write`
    #[inline]
    pub fn write_unchanged<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        self.position.write_unchanged(writer, self.buffer)
    }
}


/// A FASTQ record that ownes its data (requires allocations)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}


impl Record for OwnedRecord {
    #[inline]
    fn head(&self) -> &[u8] { &self.head }
    #[inline]
    fn seq(&self) -> &[u8]  { &self.seq  }
    #[inline]
    fn qual(&self) -> &[u8] { &self.qual }
}


/// Set of FASTQ records that owns it's buffer
/// and knows the positions of each record.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RecordSet {
    buffer: Vec<u8>,
    positions: Vec<RecordPosition>,
}

impl Default for RecordSet {
    fn default() -> RecordSet {
        RecordSet {
            buffer: vec![],
            positions: vec![],
        }
    }
}

impl<'a> iter::IntoIterator for &'a RecordSet {
    type Item = RefRecord<'a>;
    type IntoIter = RecordSetIter<'a>;
     fn into_iter(self) -> Self::IntoIter {
         RecordSetIter {
             buffer: &self.buffer,
             pos: self.positions.iter(),
         }
     }
}


pub struct RecordSetIter<'a> {
    buffer: &'a [u8],
    pos: slice::Iter<'a, RecordPosition>
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
/// to the FASTQ format
pub fn write_to<W: io::Write>(
    writer: &mut W,
    head: &[u8], seq: &[u8], qual: &[u8]) -> io::Result<()>
{
    writer.write_all(b"@")?;
    writer.write_all(head)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}


/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTQ format. The ID and description parts of the header are supplied separately
/// instead of a whole header line
pub fn write_parts<W: io::Write>(
    writer: &mut W,
    id: &[u8], desc: Option<&[u8]>,
    seq: &[u8], qual: &[u8]) -> io::Result<()>
{
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    if let Some(d) = desc {
        writer.write_all(b" ")?;
        writer.write_all(d)?;
    }
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

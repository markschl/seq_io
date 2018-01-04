//! Efficient FASTQ reading and writing
//!

use std::io::{self,BufRead,Seek};
use std::fs::File;
use std::path::Path;
use memchr::memchr;
use std::slice;
use std::iter;
use std::char;
use std::str::{self,Utf8Error};

use buf_redux;

use std::error::Error as StdError;
use super::*;


type DefaultBufStrategy = DoubleUntil8M;

const BUFSIZE: usize = 68 * 1024;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
enum SearchPos {
    HEAD,
    SEQ,
    SEP,
    QUAL
}


/// FASTQ parser.
pub struct Reader<R: io::Read, S = DefaultBufStrategy> {
    buffer: buf_redux::BufReader<R, ReadAlways, buf_redux::strategy::NeverMove>,
    buf_pos: BufferPosition,
    search_pos: SearchPos,
    position: Position,
    finished: bool,
    buf_strategy: S,
}


impl<R> Reader<R, DefaultBufStrategy>
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

impl Reader<File, DefaultBufStrategy> {
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Reader<File>> {
        File::open(path).map(Reader::new)
    }
}


impl<R, S> Reader<R, S>
    where R: io::Read,
          S: BufStrategy
{
    /// Creates a new reader with a given buffer capacity and growth strategy. See
    /// [See here](trait.BufStrategy.html) for an example using the FASTA reader, but otherwise
    /// equivalent.
    pub fn with_cap_and_strategy(reader: R, cap: usize, buf_strategy: S) -> Reader<R, S> {
        assert!(cap >= 3);
        Reader {
            buffer: buf_redux::BufReader::with_cap_and_strategies(
                reader,
                cap,
                ReadAlways, buf_redux::strategy::NeverMove
            ),
            buf_pos: BufferPosition::default(),
            search_pos: SearchPos::HEAD,
            position: Position::new(1, 0),
            finished: false,
            buf_strategy: buf_strategy,
        }
    }

    fn proceed(&mut self) -> Option<Result<(), Error>> {

        if self.finished || ! self.initialized() && ! try_opt!(self.init()) {
            return None;
        }

        if ! self.buf_pos.is_new() {
            self.next_pos();
        }

        if ! try_opt!(self.find()) && ! try_opt!(self.next_complete()) {
            return None;
        }

        Some(Ok(()))
    }

    /// Search the next FASTQ record and return a `RefRecord` that
    /// borrows it's data from the underlying buffer of this reader
    pub fn next<'a>(&'a mut self) -> Option<Result<RefRecord<'a>, Error>> {
        self.proceed().map(|r| r.map(
            move |_| RefRecord { buffer: self.get_buf(), buf_pos: &self.buf_pos }
        ))
    }

    #[inline(never)]
    fn init(&mut self) -> Result<bool, Error> {
        let n = fill_buf(&mut self.buffer)?;
        if n == 0 {
            self.finished = true;
            return Ok(false);
        }
        Ok(true)
    }

    #[inline]
    fn initialized(&self) -> bool {
        self.get_buf().len() != 0
    }

    /// Updates a `RecordSet` with a new buffer and searches for records. Old data will be erased.
    /// Returns `None` if the input reached its end
    pub fn read_record_set(&mut self, rset: &mut RecordSet) -> Option<Result<(), Error>> {

        if self.finished {
            return None;
        }

        if  ! self.initialized() {
            if ! try_opt!(self.init()) {
                return None;
            }
            if ! try_opt!(self.find()) {
                return Some(Ok(()));
            }
        } else {
            if ! try_opt!(self.next_complete()) {
                return None;
            }
        };

        rset.buffer.clear();
        rset.buffer.extend(self.get_buf());

        rset.buf_positions.clear();
        rset.buf_positions.push(self.buf_pos.clone());

        loop {
            self.next_pos();
            if ! try_opt!(self.find()) {
                return Some(Ok(()))
            }
            rset.buf_positions.push(self.buf_pos.clone());
        }
    }

    #[inline]
    fn get_buf(&self) -> &[u8] {
        self.buffer.get_buf()
    }

    // Sets starting points for next position
    #[inline]
    fn next_pos(&mut self) {
        self.position.byte += (self.buf_pos.pos.1 + 1 - self.buf_pos.pos.0) as u64;
        self.position.line += 4;
        self.buf_pos.pos.0 = self.buf_pos.pos.1 + 1;
    }

    // Reads the current record  and returns true if found
    // (false if incomplete because end of buffer reached)
    // meaning that the last record may be incomplete
    // search_start >= self.buf_pos.start
    // updates the position.
    fn find(&mut self) -> Result<bool, Error> {

        self.buf_pos.seq = unwrap_or!(self.find_line(self.buf_pos.pos.0), {
            self.search_pos = SearchPos::HEAD;
            return Ok(false);
        });

        self.buf_pos.sep = unwrap_or!(self.find_line(self.buf_pos.seq), {
            self.search_pos = SearchPos::SEQ;
            return Ok(false);
        });

        self.buf_pos.qual = unwrap_or!(self.find_line(self.buf_pos.sep), {
            self.search_pos = SearchPos::SEP;
            return Ok(false);
        });

        self.buf_pos.pos.1 = unwrap_or!(self.find_line(self.buf_pos.qual), {
            self.search_pos = SearchPos::QUAL;
            return Ok(false);
        }) - 1;

        self.validate()?;

        Ok(true)
    }

    // Resumes reading an incomplete record without
    // re-searching positions that were already found.
    // The resulting position may still be incomplete (-> false).
    fn find_incomplete(&mut self) -> Result<bool, Error> {

        if self.search_pos == SearchPos::HEAD {
            self.buf_pos.seq = unwrap_or!(self.find_line(self.buf_pos.pos.0), {
                self.search_pos = SearchPos::HEAD;
                return Ok(false);
            });
        }

        if self.search_pos <= SearchPos::SEQ {
            self.buf_pos.sep = unwrap_or!(self.find_line(self.buf_pos.seq), {
                self.search_pos = SearchPos::SEQ;
                return Ok(false);
            });
        }

        if self.search_pos <= SearchPos::SEP {
            self.buf_pos.qual = unwrap_or!(self.find_line(self.buf_pos.sep), {
                self.search_pos = SearchPos::SEP;
                return Ok(false);
            });
        }

        if self.search_pos <= SearchPos::QUAL {
            self.buf_pos.pos.1 = unwrap_or!(self.find_line(self.buf_pos.qual), {
                self.search_pos = SearchPos::QUAL;
                return Ok(false);
            }) - 1;
        }

        self.search_pos = SearchPos::HEAD;

        self.validate()?;

        Ok(true)
    }

    // should only be called on a complete BufferPosition
    #[inline(always)]  // has performance impact and would not be inlined otherwise
    fn validate(&mut self) -> Result<(), Error> {
        let start_byte = self.get_buf()[self.buf_pos.pos.0];
        if start_byte != b'@' {
            self.finished = true;
            return Err(Error::InvalidStart {
                found: start_byte,
                pos: self.get_error_pos(0, false)
            });
        }

        let sep_byte = self.get_buf()[self.buf_pos.sep];
        if sep_byte != b'+' {
            self.finished = true;
            return Err(Error::InvalidSep {
                found: sep_byte,
                pos: self.get_error_pos(2, true)
            });
        }

        let qual_len = self.buf_pos.pos.1 - self.buf_pos.qual + 1;
        let seq_len = self.buf_pos.sep - self.buf_pos.seq;
        if seq_len != qual_len {
            self.finished = true;
            return Err(Error::UnequalLengths {
                seq: self.buf_pos.seq(self.get_buf()).len(),
                qual: self.buf_pos.qual(self.get_buf()).len(),
                pos: self.get_error_pos(0, true)
            });
        }
        Ok(())
    }

    #[inline(never)]
    fn get_error_pos(&self, offset: u64, parse_id: bool) -> ErrorPosition {
        let id =
            if parse_id {
                let id = self.buf_pos
                    .head(self.get_buf())
                    .split(|b| *b == b' ')
                    .next()
                    .unwrap();
                Some(String::from_utf8_lossy(id).into())
            } else {
                None
            };
        ErrorPosition {
            line: self.position.line + offset,
            id: id
        }
    }

    #[inline]
    fn find_line(&self, search_start: usize) -> Option<usize> {
        memchr(b'\n', &self.get_buf()[search_start.. ]).map(|pos| search_start + pos + 1)
    }

    // To be called when the end of the buffer is reached and `next_pos` does not find
    // the next record. Incomplete bytes will be moved to the start of the buffer.
    // If the record still doesn't fit in, the buffer will be enlarged.
    // After calling this function, the position will therefore always be 'complete'.
    #[inline(never)]
    fn next_complete(&mut self) -> Result<bool, Error> {

        loop {
            let bufsize = self.get_buf().len();
            if bufsize < self.buffer.capacity() {
                // EOF reached, there will be no next record
                self.finished = true;
                if self.search_pos == SearchPos::QUAL {
                    // no line ending at end of last record
                    self.buf_pos.pos.1 = bufsize;
                    self.validate()?;
                    return Ok(true);
                }

                let rest = &self.get_buf()[self.buf_pos.pos.0..];
                if rest.split(|c| *c == b'\n').all(|l| trim_cr(l).len() == 0)  {
                    // allow up to 3 newlines after last record (more will cause an Unexpected error)
                    return Ok(false);
                }

                return Err(Error::UnexpectedEnd {
                    pos: self.get_error_pos(self.search_pos as u64, self.search_pos > SearchPos::HEAD)
                });

            } else if self.buf_pos.pos.0 == 0 {
                // first record already incomplete -> buffer too small
                self.grow()?;
            } else {
                // not the first record -> buffer may be big enough
                self.make_room();
            }

            fill_buf(&mut self.buffer)?;

            // self.buf_pos.pos.1 = 0;
            // self.search_pos = SearchPos::HEAD;
            if self.find_incomplete()? {
                return Ok(true);
            }
        }
    }

    // grow buffer based on strategy
    fn grow(&mut self) -> Result<(), Error> {
        let cap = self.buffer.capacity();
        let new_size = self.buf_strategy.grow_to(cap).ok_or(Error::BufferLimit)?;
        let additional = new_size - cap;
        self.buffer.grow(additional);
        Ok(())
    }

    // move incomplete bytes to start of buffer and retry
    fn make_room(&mut self) {
        let consumed = self.buf_pos.pos.0;
        self.buffer.consume(consumed);
        self.buffer.make_room();

        self.buf_pos.pos.0 = 0;

        if self.search_pos >= SearchPos::SEQ {
            self.buf_pos.seq -= consumed;
        }
        if self.search_pos >= SearchPos::SEP {
            self.buf_pos.sep -= consumed;
        }
        if self.search_pos >= SearchPos::QUAL {
            self.buf_pos.qual -= consumed;
        }
    }

    /// Returns the current position (useful with `seek()`)
    #[inline]
    pub fn position(&self) -> &Position {
        &self.position
    }


    /// Returns a borrowed iterator over all FASTA records. The records
    /// are owned (`OwnedRecord`), this is therefore slower than using
    /// the `Reader::next()`, but makes sense if an owned copy is required.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate seq_io;
    /// # fn main() {
    /// use seq_io::fastq::{Reader,OwnedRecord};
    ///
    /// let fastq = b"@id1
    /// ACGT
    /// +
    /// IIII
    /// @id2
    /// TGCA
    /// +
    /// IIII";
    ///
    /// let mut reader = Reader::new(&fastq[..]);
    ///
    /// let records: Result<Vec<_>, _> = reader
    ///     .records()
    ///     .collect();
    ///
    /// assert_eq!(records.unwrap(),
    ///     vec![
    ///         OwnedRecord {head: b"id1".to_vec(), seq: b"ACGT".to_vec(), qual: b"IIII".to_vec()},
    ///         OwnedRecord {head: b"id2".to_vec(), seq: b"TGCA".to_vec(), qual: b"IIII".to_vec()}
    ///     ]
    /// );
    /// # }
    /// ```
    pub fn records(&mut self) -> RecordsIter<R, S> {
        RecordsIter { rdr: self }
    }

    /// Returns an iterator over all FASTA records like `Reader::records()`,
    /// but with the difference that it owns the underlying reader.
    /// Returns a borrowed iterator over all FASTA records. The records
    /// are owned (`OwnedRecord`), this is therefore slower than using
    /// the `Reader::next()`, but makes sense if an owned copy is required.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate seq_io;
    /// # fn main() {
    /// use seq_io::fastq::{Reader,OwnedRecord};
    ///
    /// let fastq = b"@id1
    /// ACGT
    /// +
    /// IIII
    /// @id2
    /// TGCA
    /// +
    /// IIII";
    ///
    /// let mut reader = Reader::new(&fastq[..]);
    ///
    /// let records: Result<Vec<_>, _> = reader
    ///     .into_records()
    ///     .collect();
    ///
    /// assert_eq!(records.unwrap(),
    ///     vec![
    ///         OwnedRecord {head: b"id1".to_vec(), seq: b"ACGT".to_vec(), qual: b"IIII".to_vec()},
    ///         OwnedRecord {head: b"id2".to_vec(), seq: b"TGCA".to_vec(), qual: b"IIII".to_vec()}
    ///     ]
    /// );
    /// # }
    /// ```
    pub fn into_records(self) -> RecordsIntoIter<R, S> {
        RecordsIntoIter { rdr: self }
    }
}


impl<R, S> Reader<R, S>
    where R: io::Read + Seek,
          S: BufStrategy
{
    /// Seeks to a specified position.
    /// Keep the underyling buffer if the seek position is found within it, otherwise it has to be
    /// discarded.
    /// If an error was returned before, seeking to that position will return the same error.
    /// The same is not always true with `None`. If there is no newline character at the end of the
    /// file, the last record will be returned instead of `None`.
    pub fn seek(&mut self, pos: &Position) -> Result<(), Error> {
        self.finished = false;
        let diff = pos.byte as i64 - self.position.byte as i64;
        let endpos = self.buf_pos.pos.0 as i64 + diff;
        self.position = pos.clone();

        if endpos >= 0 && endpos < (self.get_buf().len() as i64) {
            // position reachable within buffer -> no actual seeking necessary
            self.buf_pos.reset(endpos as usize); // is_new() will return true
            return Ok(());
        }

        self.buffer.seek(io::SeekFrom::Start(pos.byte))?;
        self.buf_pos.reset(0);
        Ok(())
    }
}


/// Borrowed iterator of `OwnedRecord`
pub struct RecordsIter<'a, R, S = DefaultBufStrategy>
    where S: 'a,
          R: io::Read + 'a
{
    rdr: &'a mut Reader<R, S>
}

impl<'a, R, S> Iterator for RecordsIter<'a, R, S>
    where S: BufStrategy + 'a,
          R: io::Read + 'a
{
    type Item = Result<OwnedRecord, Error>;
    fn next(&mut self) -> Option<Self::Item> {
        self.rdr.next().map(|rec| rec.map(|r| r.to_owned_record()))
    }
}


/// Iterator of `OwnedRecord` that owns the underlying reader
pub struct RecordsIntoIter<R: io::Read, S = DefaultBufStrategy> {
    rdr: Reader<R, S>
}

impl<R, S> Iterator for RecordsIntoIter<R, S>
    where S: BufStrategy,
          R: io::Read
{
    type Item = Result<OwnedRecord, Error>;
    fn next(&mut self) -> Option<Self::Item> {
        self.rdr.next().map(|rec| rec.map(|r| r.to_owned_record()))
    }
}


#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Position {
    line: u64,
    byte: u64,
}

impl Position {
    pub fn new(line: u64, byte: u64) -> Position {
        Position { line: line, byte: byte }
    }

    /// Line number (starting with 1)
    pub fn line(&self) -> u64 {
        self.line
    }

    /// Byte offset within the file
    pub fn byte(&self) -> u64 {
        self.byte
    }
}


#[derive(Debug)]
pub enum Error {
    Io(io::Error),
    /// sequence and qualitiy lengths are not equal
    UnequalLengths {
        /// Length of sequence
        seq: usize,
        /// Length of qualities
        qual: usize,
        /// Position within file.
        /// `ErrorPosition::line` has the position of the header, not sequence/qualities
        pos: ErrorPosition,
    },
    /// Invalid start byte encountered (expected '@')
    InvalidStart {
        /// Byte found instead.
        found: u8,
        /// Position within file. `ErrorPosition::id` will be `None`.
        pos: ErrorPosition,
    },
    /// Invalid separator byte encountered (expected '+')
    InvalidSep {
        /// Byte found instead.
        found: u8,
        /// Position within file
        pos: ErrorPosition,
    },
    /// Truncated record found
    UnexpectedEnd {
        /// Position within file.
        pos: ErrorPosition,
    },
    /// Size limit of buffer was reached, which happens if `BufStrategy::new_size()` returned
    /// `None` (not the case by default).
    BufferLimit,
}


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ErrorPosition {
    /// Line number where the error occurred (starting with 1)
    pub line: u64,
    /// ID of record if available
    pub id: Option<String>,
}


impl fmt::Display for ErrorPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(id) = self.id.as_ref() {
            write!(f, "record '{}' at ", id)?;
        }
        write!(f, "line {}", self.line)
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Io(ref e) => e.fmt(f),
            Error::UnequalLengths { seq, qual, ref pos } => write!(f,
                "FASTQ parse error: sequence length is {}, but quality length is {} ({}).",
                seq, qual, pos
            ),
            Error::InvalidStart { found, ref pos } => write!(f,
                "FASTQ parse error: expected '@' at record start but found '{}' ({}).",
                (found as char).escape_default(), pos
            ),
            Error::InvalidSep { found, ref pos } => write!(f,
                "FASTQ parse error: Expected '+' separator but found '{}' ({}).",
                (found as char).escape_default(), pos
            ),
            Error::UnexpectedEnd { ref pos } => write!(f,
                "FASTQ parse error: unexpected end of input ({}).", pos
            ),
            Error::BufferLimit => write!(f,
                "FASTQ parse error: Buffer limit reached."
            ),
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Error {
        Error::Io(e)
    }
}

impl StdError for Error {
    fn description(&self) -> &str {
        match *self {
            Error::Io(ref e) => e.description(),
            Error::UnequalLengths {..} => "sequence and quality lengths are different",
            Error::InvalidStart {..} => "invalid record start",
            Error::InvalidSep {..} => "invalid record separator",
            Error::UnexpectedEnd {..} => "unexpected end of input",
            Error::BufferLimit => "buffer limit reached",
        }
    }

    fn cause(&self) -> Option<&StdError> {
        match *self {
            Error::Io(ref err) => Some(err),
            _ => None
        }
    }
}


/// Represents the position of a record within a buffer
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct BufferPosition {
    // (start, stop), but might include \r at the end
    pos: (usize, usize),
    seq: usize,
    sep: usize,
    qual: usize,
}

impl BufferPosition {

    #[inline]
    fn is_new(&self) -> bool {
        self.pos.1 == 0
    }

    #[inline]
    fn reset(&mut self, start: usize) {
        self.pos.0 = start;
        self.pos.1 = 0;
    }

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
        trim_cr(&buffer[self.qual .. self.pos.1])
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


/// A FASTQ record that borrows data from a buffer
#[derive(Debug, Clone)]
pub struct RefRecord<'a> {
    buffer: &'a [u8],
    buf_pos: &'a BufferPosition,
}


impl<'a> Record for RefRecord<'a> {

    #[inline]
    fn head(&self) -> &[u8] {
        self.buf_pos.head(self.buffer)
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        self.buf_pos.seq(self.buffer)
    }

    #[inline]
    fn qual(&self) -> &[u8] {
        self.buf_pos.qual(self.buffer)
    }
}


impl<'a> RefRecord<'a> {
    #[inline]
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            head: self.head().to_vec(),
            seq: self.seq().to_vec(),
            qual: self.qual().to_vec(),
        }
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `RefRecord::write`
    #[inline]
    pub fn write_unchanged<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        #[inline]
        let data = &self.buffer[self.buf_pos.pos.0 .. self.buf_pos.pos.1];
        writer.write_all(data)?;
        writer.write_all(b"\n")
    }
}


/// A FASTQ record that ownes its data (requires allocations)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
    buf_positions: Vec<BufferPosition>,
}

impl Default for RecordSet {
    fn default() -> RecordSet {
        RecordSet {
            buffer: vec![],
            buf_positions: vec![],
        }
    }
}

impl<'a> iter::IntoIterator for &'a RecordSet {
    type Item = RefRecord<'a>;
    type IntoIter = RecordSetIter<'a>;
     fn into_iter(self) -> Self::IntoIter {
         RecordSetIter {
             buffer: &self.buffer,
             pos: self.buf_positions.iter(),
         }
     }
}


pub struct RecordSetIter<'a> {
    buffer: &'a [u8],
    pos: slice::Iter<'a, BufferPosition>
}

impl<'a> Iterator for RecordSetIter<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<RefRecord<'a>> {
        self.pos.next().map(|p| {
            RefRecord {
                buffer: self.buffer,
                buf_pos: p,
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

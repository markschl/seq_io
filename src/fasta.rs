//! Efficient FASTA reading and writing
//!
//! # Example
//!
//! This example reads some content, writes it back and compares the output
//! (should be the same):
//!
//! ```
//! use seq_io::fasta::{Reader, Record};
//!
//! let input = b">id1
//! ACGT
//! ACGT
//! >id2
//! TGCA
//! TGCA
//! ";
//!
//! let mut reader = Reader::new(&input[..]);
//! let mut output = vec![];
//!
//! while let Some(record) = reader.next() {
//!     let record = record.expect("Error reading record");
//!     println!("id: {}", record.id().unwrap());
//!     record.write_wrap(&mut output, 4);
//! }
//!
//! assert_eq!(input, output.as_slice());
//! ```
//!
//! # Details on parsing behaviour
//!
//! * The parser handles UNIX (LF) and Windows (CRLF) line endings, but not old
//!   Mac-style (CR) endings. However, FASTA writing currently always uses UNIX
//!   line endings.
//! * Empty lines are allowed anywhere in the file, they will just be ignored.
//!   The first non-empty line must start with `>`, indicating the first header.
//! * Whitespace at the end of header and sequence lines is never removed.
//! * If two consecutive FASTA header lines (starting with `>`) are encountered
//!   without intermediate sequence line, the first record will have an empty
//!   sequence. The same is true if the input ends with a header line.
//! * Empty input will result in `None` being returned immediately by
//!   `fasta::Reader::next()` and in empty iterators for `RecordsIter` /
//!   `RecordsIntoIter`.
//! * Comment lines starting with `;` are not supported.
//!   If at the start of a file, there will be an error, since `>` is expected.
//!   Intermediate comments are appended to the sequence.


use std::borrow::Cow;
use std::fs::File;
use std::io::{self, BufRead, Seek};
use std::iter;
use std::path::Path;
use std::slice;
use std::str::{self, Utf8Error};

use buf_redux;
use memchr::Memchr;

use super::*;
use super::policy::{BufPolicy, StdPolicy};

type DefaultPolicy = StdPolicy;

const BUFSIZE: usize = 64 * 1024;

/// Parser for FASTA files.
pub struct Reader<R: io::Read, P = DefaultPolicy> {
    buffer: buf_redux::BufReader<R>,
    buf_pos: BufferPosition,
    position: Position,
    search_pos: usize,
    finished: bool,
    buf_policy: P,
}

impl<R> Reader<R, DefaultPolicy>
where
    R: io::Read,
{
    /// Creates a new reader with the default buffer size of 64 KiB
    ///
    /// # Example:
    ///
    /// ```
    /// use seq_io::fasta::{Reader,Record};
    /// let fasta = b">id\nSEQUENCE";
    ///
    /// let mut reader = Reader::new(&fasta[..]);
    /// let record = reader.next().unwrap().unwrap();
    /// assert_eq!(record.id(), Ok("id"))
    /// ```
    #[inline]
    pub fn new(reader: R) -> Reader<R, StdPolicy> {
        Reader::with_capacity(reader, BUFSIZE)
    }

    /// Creates a new reader with a given buffer capacity. The minimum allowed
    /// capacity is 3.
    #[inline]
    pub fn with_capacity(reader: R, capacity: usize) -> Reader<R, DefaultPolicy> {
        assert!(capacity >= 3);
        Reader {
            buffer: buf_redux::BufReader::with_capacity(capacity, reader),
            buf_pos: BufferPosition {
                start: 0,
                seq_pos: Vec::with_capacity(2),
            },
            position: Position::new(0, 0),
            search_pos: 0,
            finished: false,
            buf_policy: StdPolicy,
        }
    }
}

impl Reader<File, DefaultPolicy> {
    /// Creates a reader from a file path.
    ///
    /// # Example:
    ///
    /// ```no_run
    /// use seq_io::fasta::Reader;
    ///
    /// let mut reader = Reader::from_path("seqs.fasta").unwrap();
    ///
    /// // (... do something with the reader)
    /// ```
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Reader<File>> {
        File::open(path).map(Reader::new)
    }
}

impl<R, P> Reader<R, P>
where
    R: io::Read,
    P: BufPolicy,
{
    /// Returns a reader with the given buffer policy applied
    #[inline]
    pub fn set_policy<T: BufPolicy>(self, policy: T) -> Reader<R, T> {
        Reader {
            buffer: self.buffer,
            buf_pos: self.buf_pos,
            position: self.position,
            search_pos: self.search_pos,
            finished: self.finished,
            buf_policy: policy,
        }
    }

    /// Returns the `BufPolicy` of the reader
    #[inline]
    pub fn policy(&self) -> &P {
        &self.buf_policy
    }

    /// Searches the next FASTA record and returns a [RefRecord](struct.RefRecord.html) that
    /// borrows its data from the underlying buffer of this reader.
    ///
    /// # Example:
    ///
    /// ```no_run
    /// use seq_io::fasta::{Reader,Record};
    ///
    /// let mut reader = Reader::from_path("seqs.fasta").unwrap();
    ///
    /// while let Some(record) = reader.next() {
    ///     let record = record.unwrap();
    ///     println!("{}", record.id().unwrap());
    /// }
    /// ```
    pub fn next(&mut self) -> Option<Result<RefRecord, Error>> {
        if self.finished || !self.initialized() && !try_opt!(self.init()) {
            return None;
        }

        if !self.buf_pos.is_new() {
            self.next_pos();
        }

        if !try_opt!(self.search()) && !try_opt!(self.next_complete()) {
            return None;
        }

        Some(Ok(
            RefRecord {
                buffer: self.get_buf(),
                buf_pos: &self.buf_pos,
            }
        ))
    }

    /// Updates a [RecordSet](struct.RecordSet.html) with new data. The contents of the internal
    /// buffer are just copied over to the record set and the positions of all records are found.
    /// Old data will be erased. Returns `None` if the input reached its end.
    pub fn read_record_set(&mut self, rset: &mut RecordSet) -> Option<Result<(), Error>> {
        if self.finished {
            return None;
        }

        if !self.initialized() {
            if !try_opt!(self.init()) {
                return None;
            }
            if !try_opt!(self.search()) {
                return Some(Ok(()));
            }
        } else if !try_opt!(self.next_complete()) {
            return None;
        };

        // copy buffer AFTER call to next_complete (initialization of buffer is done there)
        rset.buffer.clear();
        rset.buffer.extend(self.get_buf());

        // Update records that are already in the positions vector
        let mut n = 0;
        for pos in &mut rset.positions {
            n += 1;
            pos.update(&self.buf_pos);

            self.next_pos();
            if self.finished || !try_opt!(self.search()) {
                rset.npos = n;
                return Some(Ok(()));
            }
        }

        // Add more positions if necessary
        loop {
            n += 1;
            rset.positions.push(self.buf_pos.clone());

            self.next_pos();
            if self.finished || !try_opt!(self.search()) {
                rset.npos = n;
                return Some(Ok(()));
            }
        }
    }

    // Sets starting points for next position
    #[inline]
    fn next_pos(&mut self) {
        self.position.line += self.buf_pos.seq_pos.len() as u64;
        self.position.byte += (self.search_pos - self.buf_pos.start) as u64;
        self.buf_pos.start = self.search_pos;
        self.buf_pos.seq_pos.clear();
    }

    #[inline(always)]
    fn get_buf(&self) -> &[u8] {
        self.buffer.buffer()
    }

    #[inline(always)]
    fn initialized(&self) -> bool {
        self.position.line != 0
    }

    // moves to the first record positon, ignoring newline characters
    fn init(&mut self) -> Result<bool, Error> {
        if let Some((line_num, pos, byte)) = self.first_byte()? {
            if byte == b'>' {
                self.buf_pos.start = pos;
                self.position.byte = pos as u64;
                self.position.line = line_num as u64;
                self.search_pos = pos + 1;
                return Ok(true);
            } else {
                self.finished = true;
                return Err(Error::InvalidStart {
                    line: line_num,
                    found: byte,
                });
            }
        }
        self.finished = true;
        Ok(false)
    }

    fn first_byte(&mut self) -> Result<Option<(usize, usize, u8)>, Error> {

        let mut line_num = 0;

        while fill_buf(&mut self.buffer)? > 0 {
            let mut pos = 0;

            for line in self.get_buf().split(|b| *b == b'\n') {
                line_num += 1;
                if !line.is_empty() && line != b"\r" {
                    return Ok(Some((line_num, pos, line[0])));
                }
                pos += line.len() + 1;
            }
            self.buffer.consume(pos - 1);
        }
        Ok(None)
    }

    /// Finds the position of the next record
    /// and returns true if found; false if end of buffer reached.
    #[inline]
    fn search(&mut self) -> Result<bool, Error> {
        if self._search() {
            return Ok(true);
        }

        // nothing found
        if self.get_buf().len() < self.buffer.capacity() {
            // EOF reached, there will be no next record
            self.finished = true;
            self.buf_pos.seq_pos.push(self.search_pos);
            return Ok(true);
        }

        Ok(false)
    }

    // returns true if complete position found, false if end of buffer reached.
    #[inline]
    fn _search(&mut self) -> bool {
        let bufsize = self.get_buf().len();

        for pos in Memchr::new(b'\n', &self.buffer.buffer()[self.search_pos..]) {
            let pos = self.search_pos + pos;
            let next_line_start = pos + 1;

            if next_line_start == bufsize {
                // cannot check next byte -> treat as incomplete
                self.search_pos = pos; // make sure last byte is re-searched next time
                return false;
            }

            self.buf_pos.seq_pos.push(pos);
            if self.get_buf()[next_line_start] == b'>' {
                // complete record was found
                self.search_pos = next_line_start;
                return true;
            }
        }

        // record end not found
        self.search_pos = bufsize;

        false
    }

    /// To be called when the end of the buffer is reached and `next_pos` does not find
    /// the next record. Incomplete bytes will be moved to the start of the buffer.
    /// If the record still doesn't fit in, the buffer will be enlarged.
    /// After calling this function, the position will therefore always be 'complete'.
    /// this function assumes that the buffer was fully searched
    fn next_complete(&mut self) -> Result<bool, Error> {
        loop {
            if self.buf_pos.start == 0 {
                // first record -> buffer too small
                self.grow()?;
            } else {
                // not the first record -> buffer may be big enough
                self.make_room();
            }

            // fill up remaining buffer
            fill_buf(&mut self.buffer)?;

            if self.search()? {
                return Ok(true);
            }
        }
    }

    // grow buffer
    fn grow(&mut self) -> Result<(), Error> {
        let cap = self.buffer.capacity();
        let new_size = self.buf_policy.grow_to(cap).ok_or(Error::BufferLimit)?;
        let additional = new_size - cap;
        self.buffer.reserve(additional);
        Ok(())
    }

    // move incomplete bytes to start of buffer
    fn make_room(&mut self) {
        let consumed = self.buf_pos.start;
        self.buffer.consume(consumed);
        self.buffer.make_room();
        self.buf_pos.start = 0;
        self.search_pos -= consumed;
        for s in &mut self.buf_pos.seq_pos {
            *s -= consumed;
        }
    }

    /// Returns the current position (useful with `seek()`).
    /// If `next()` has not yet been called, `None` will be returned.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate seq_io;
    /// # fn main() {
    /// use seq_io::fasta::{Reader,Position};
    ///
    /// let fasta = b">id1
    /// ACGT
    /// >id2
    /// TGCA";
    ///
    /// let mut reader = Reader::new(&fasta[..]);
    ///
    /// // skip one record
    /// reader.next().unwrap();
    /// // second position
    /// reader.next().unwrap();
    ///
    /// assert_eq!(reader.position(), Some(&Position::new(3, 10)));
    /// # }
    /// ```
    #[inline]
    pub fn position(&self) -> Option<&Position> {
        if self.buf_pos.is_new() {
            return None;
        }
        Some(&self.position)
    }

    /// Returns a borrowed iterator over all FASTA records. The records
    /// are owned (`OwnedRecord`), this is therefore slower than using
    /// `Reader::next()`.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate seq_io;
    /// # fn main() {
    /// use seq_io::fasta::{Reader,OwnedRecord};
    ///
    /// let fasta = b">id1
    /// ACGT
    /// >id2
    /// TGCA";
    ///
    /// let mut reader = Reader::new(&fasta[..]);
    ///
    /// let records: Result<Vec<_>, _> = reader
    ///     .records()
    ///     .collect();
    ///
    /// assert_eq!(records.unwrap(),
    ///     vec![
    ///         OwnedRecord {head: b"id1".to_vec(), seq: b"ACGT".to_vec()},
    ///         OwnedRecord {head: b"id2".to_vec(), seq: b"TGCA".to_vec()}
    ///     ]
    /// );
    /// # }
    /// ```
    pub fn records(&mut self) -> RecordsIter<R, P> {
        RecordsIter { rdr: self }
    }

    /// Returns an iterator over all FASTA records like `Reader::records()`,
    /// but with the difference that it owns the underlying reader.
    pub fn into_records(self) -> RecordsIntoIter<R, P> {
        RecordsIntoIter { rdr: self }
    }
}

impl<R, P> Reader<R, P>
where
    R: io::Read + Seek,
    P: BufPolicy,
{
    /// Seeks to a specified position.  Keeps the underyling buffer if the seek position is
    /// found within it, otherwise it has to be discarded.
    /// If an error was returned before, seeking to that position will return the same error.
    /// The same is not always true with `None`. If there is no newline character at the end of the
    /// file, the last record will be returned instead of `None`.
    ///
    /// # Example
    ///
    /// ```
    /// # extern crate seq_io;
    /// # fn main() {
    /// use seq_io::fasta::{Reader,Position,OwnedRecord};
    /// use std::io::Cursor;
    ///
    /// let fasta = b">id1
    /// ACGT
    /// >id2
    /// TGCA";
    ///
    /// let mut cursor = Cursor::new(&fasta[..]);
    /// let mut reader = Reader::new(cursor);
    ///
    /// // read the first record and get its position
    /// let record1 = reader.next().unwrap().unwrap().to_owned_record();
    /// let pos1 = reader.position().unwrap().to_owned();
    ///
    /// // read the second record
    /// reader.next().unwrap().unwrap();
    ///
    /// // now seek to position of first record
    /// reader.seek(&pos1);
    /// assert_eq!(reader.next().unwrap().unwrap().to_owned_record(), record1);
    /// # }
    /// ```
    pub fn seek(&mut self, pos: &Position) -> Result<(), Error> {
        self.finished = false;
        let diff = pos.byte as i64 - self.position.byte as i64;
        let rel_pos = self.buf_pos.start as i64 + diff;
        if rel_pos >= 0 && rel_pos < (self.get_buf().len() as i64) {
            // position reachable within buffer -> no actual seeking necessary
            self.search_pos = rel_pos as usize;
            self.buf_pos.reset(rel_pos as usize);
            self.position = pos.clone();
            return Ok(());
        }
        self.position = pos.clone();
        self.search_pos = 0;
        self.buffer.seek(io::SeekFrom::Start(pos.byte))?;
        fill_buf(&mut self.buffer)?;
        self.buf_pos.reset(0);
        Ok(())
    }
}

/// Borrowed iterator of `OwnedRecord`
pub struct RecordsIter<'a, R, P = DefaultPolicy>
where
    P: 'a,
    R: io::Read + 'a,
{
    rdr: &'a mut Reader<R, P>,
}

impl<'a, R, P> Iterator for RecordsIter<'a, R, P>
where
    P: BufPolicy + 'a,
    R: io::Read + 'a,
{
    type Item = Result<OwnedRecord, Error>;
    fn next(&mut self) -> Option<Self::Item> {
        self.rdr.next().map(|rec| rec.map(|r| r.to_owned_record()))
    }
}

/// Iterator of `OwnedRecord` that owns the underlying reader
pub struct RecordsIntoIter<R: io::Read, P = DefaultPolicy> {
    rdr: Reader<R, P>,
}

impl<R, P> Iterator for RecordsIntoIter<R, P>
where
    P: BufPolicy,
    R: io::Read,
{
    type Item = Result<OwnedRecord, Error>;
    fn next(&mut self) -> Option<Self::Item> {
        self.rdr.next().map(|rec| rec.map(|r| r.to_owned_record()))
    }
}

/// Holds line number and byte offset of a FASTA record
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Position {
    line: u64,
    byte: u64,
}

impl Position {
    pub fn new(line: u64, byte: u64) -> Position {
        Position { line, byte }
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

/// FASTA parsing error
#[derive(Debug)]
pub enum Error {
    /// io::Error
    Io(io::Error),
    /// First non-empty line does not start with `>`
    InvalidStart {
        /// line number (1-based)
        line: usize,
        /// byte that was found instead
        found: u8,
    },
    /// Size limit of buffer was reached, which happens if `policy::BufPolicy::grow_to()` returned
    /// `None`. This does not happen with the default `struct.DoubleUntil.html` policy.
    BufferLimit,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Io(ref e) => e.fmt(f),
            Error::InvalidStart { line, found } => write!(
                f,
                "FASTA parse error: expected '>' but found '{}' at file start, line {}.",
                (found as char).escape_default(),
                line
            ),
            Error::BufferLimit => write!(f, "FASTA parse error: buffer limit reached."),
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Error {
        Error::Io(e)
    }
}

impl error::Error for Error {
    fn description(&self) -> &str {
        match *self {
            Error::Io(ref e) => e.description(),
            Error::InvalidStart { .. } => "invalid record start",
            Error::BufferLimit => "buffer limit reached",
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct BufferPosition {
    /// index of '>'
    start: usize,
    /// Indicate line start, but actually it is one byte before (start - 1), which is usually
    /// the line terminator of the header (if there is one). The last index in the Vec is always
    /// the last byte of the last sequence line (including line terminator if present).
    /// Therefore, the length of this Vec should never be 0.
    seq_pos: Vec<usize>,
}

impl BufferPosition {
    #[inline]
    fn is_new(&self) -> bool {
        self.seq_pos.is_empty()
    }

    #[inline]
    fn reset(&mut self, start: usize) {
        self.seq_pos.clear();
        self.start = start;
    }

    #[inline]
    fn update(&mut self, other: &Self) {
        self.start = other.start;
        self.seq_pos.clear();
        self.seq_pos.extend(&other.seq_pos);
    }
}

/// FASTA record trait implemented by both `RefRecord` and `OwnedRecord`
pub trait Record {
    /// Return the header line of the record as byte slice
    fn head(&self) -> &[u8];
    /// Return the FASTA sequence as byte slice
    fn seq(&self) -> &[u8];
    /// Write the record to the given `io::Write` instance. The sequence will occupy one line only.
    fn write<W: io::Write>(&self, writer: W) -> io::Result<()>;
    /// Write the record to the given `io::Write` instance. The sequence is wrapped to produce
    ///  multi-line FASTA with a maximum width specified by `wrap`.
    fn write_wrap<W: io::Write>(&self, writer: W, wrap: usize) -> io::Result<()>;

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
}

/// A FASTA record that borrows data from a buffer.
#[derive(Debug, Clone)]
pub struct RefRecord<'a> {
    buffer: &'a [u8],
    buf_pos: &'a BufferPosition,
}

impl<'a> Record for RefRecord<'a> {
    #[inline]
    fn head(&self) -> &[u8] {
        trim_cr(&self.buffer[self.buf_pos.start + 1..*self.buf_pos.seq_pos.first().unwrap()])
    }

    /// Return the FASTA sequence as byte slice.
    /// Note that this method of `RefRecord` returns
    /// the **raw** sequence, which may contain line breaks.
    /// Use `seq_lines()` to iterate over all lines without
    /// breaks, or use [`full_seq()`](struct.RefRecord.html#method.full_seq)
    /// to access the whole sequence at once.
    #[inline]
    fn seq(&self) -> &[u8] {
        if self.buf_pos.seq_pos.len() > 1 {
            let start = *self.buf_pos.seq_pos.first().unwrap() + 1;
            let end = *self.buf_pos.seq_pos.last().unwrap();
            trim_cr(&self.buffer[start..end])
        } else {
            b""
        }
    }

    #[inline]
    fn write<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        write_head(&mut writer, self.head())?;
        write_seq_iter(&mut writer, self.seq_lines())
    }

    #[inline]
    fn write_wrap<W: io::Write>(&self, mut writer: W, wrap: usize) -> io::Result<()> {
        write_head(&mut writer, self.head())?;
        write_wrap_seq_iter(&mut writer, self.seq_lines(), wrap)
    }
}

impl<'a> RefRecord<'a> {
    /// Return an iterator over all sequence lines in the data
    #[inline]
    pub fn seq_lines(&self) -> SeqLines {
        SeqLines {
            data: &self.buffer,
            len: self.buf_pos.seq_pos.len() - 1,
            pos_iter: self
                .buf_pos
                .seq_pos
                .iter()
                .zip(self.buf_pos.seq_pos.iter().skip(1)),
        }
    }

    /// Returns the number of sequence lines.
    /// Equivalent to `self.seq_lines().len()`
    #[inline]
    pub fn num_seq_lines(&self) -> usize {
        self.seq_lines().len()
    }

    /// Returns the full sequence. If the sequence consists of a single line,
    /// then the sequence will be borrowed from the underlying buffer
    /// (equivalent to calling `RefRecord::seq()`). If there are multiple
    /// lines, an owned copy will be created (equivalent to `RefRecord::owned_seq()`).
    pub fn full_seq(&self) -> Cow<[u8]> {
        if self.num_seq_lines() == 1 {
            // only one line
            self.seq().into()
        } else {
            self.owned_seq().into()
        }
    }

    /// Returns the sequence as owned `Vec`. **Note**: This function
    /// must be called in order to obtain a sequence that does not contain
    /// line endings (as returned by `seq()`)
    pub fn owned_seq(&self) -> Vec<u8> {
        let mut seq = Vec::new();
        for segment in self.seq_lines() {
            seq.extend(segment);
        }
        seq
    }

    /// Creates an owned copy of the record.
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            head: self.head().to_vec(),
            seq: self.owned_seq(),
        }
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `RefRecord::write`
    pub fn write_unchanged<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        let data = &self.buffer[self.buf_pos.start..*self.buf_pos.seq_pos.last().unwrap()];
        writer.write_all(data)?;
        if *data.last().unwrap() != b'\n' {
            writer.write_all(&[b'\n'])?;
        }
        Ok(())
    }
}

/// Iterator over sequence the lines of a FASTA record.
pub struct SeqLines<'a> {
    data: &'a [u8],
    len: usize,
    pos_iter: iter::Zip<slice::Iter<'a, usize>, iter::Skip<slice::Iter<'a, usize>>>,
}

impl<'a> Iterator for SeqLines<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<&'a [u8]> {
        self.pos_iter
            .next()
            .map(|(start, next_start)| trim_cr(&self.data[*start + 1..*next_start]))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        let l = self.len();
        (l, Some(l))
    }
}

impl<'a> DoubleEndedIterator for SeqLines<'a> {
    #[inline]
    fn next_back(&mut self) -> Option<&'a [u8]> {
        self.pos_iter
            .next_back()
            .map(|(start, next_start)| trim_cr(&self.data[*start + 1..*next_start]))
    }
}

impl<'a> ExactSizeIterator for SeqLines<'a> {
    #[inline]
    fn len(&self) -> usize {
        self.len
    }
}

/// A FASTA record that ownes its data (requiring two allocations)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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

    #[inline]
    fn write<W: io::Write>(&self, writer: W) -> io::Result<()> {
        write_to(writer, &self.head, &self.seq)
    }

    #[inline]
    fn write_wrap<W: io::Write>(&self, mut writer: W, wrap: usize) -> io::Result<()> {
        write_head(&mut writer, &self.head)?;
        write_wrap_seq(&mut writer, &self.seq, wrap)
    }
}

/// Set of FASTA records that owns it'P buffer
/// and knows the positions of each record.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RecordSet {
    buffer: Vec<u8>,
    positions: Vec<BufferPosition>,
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

/// Iterator over record sets
pub struct RecordSetIter<'a> {
    buffer: &'a [u8],
    pos: iter::Take<slice::Iter<'a, BufferPosition>>,
}

impl<'a> Iterator for RecordSetIter<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<RefRecord<'a>> {
        self.pos.next().map(|p| RefRecord {
            buffer: self.buffer,
            buf_pos: p,
        })
    }
}

/// Writes data (not necessarily stored in a `Record` instance) to the FASTA format.
#[inline]
pub fn write_to<W>(mut writer: W, head: &[u8], seq: &[u8]) -> io::Result<()>
where
    W: io::Write,
{
    write_head(&mut writer, head)?;
    write_seq(&mut writer, seq)
}

/// Writes data to the FASTA format. ID and description parts of the header are supplied
/// separately instead of a whole header line.
#[inline]
pub fn write_parts<W>(mut writer: W, id: &[u8], desc: Option<&[u8]>, seq: &[u8]) -> io::Result<()>
where
    W: io::Write,
{
    write_id_desc(&mut writer, id, desc)?;
    write_seq(&mut writer, seq)
}

/// Writes data to the FASTA format. Wraps the sequence to produce multi-line FASTA
/// with a maximum width specified by the `wrap` parameter.
#[inline]
pub fn write_wrap<W>(
    mut writer: W,
    id: &[u8],
    desc: Option<&[u8]>,
    seq: &[u8],
    wrap: usize,
) -> io::Result<()>
where
    W: io::Write,
{
    write_id_desc(&mut writer, id, desc)?;
    write_wrap_seq(&mut writer, seq, wrap)
}

/// Writes only the sequence header.
#[inline]
pub fn write_head<W>(mut writer: W, head: &[u8]) -> io::Result<()>
where
    W: io::Write,
{
    writer.write_all(b">")?;
    writer.write_all(head)?;
    writer.write_all(b"\n")
}

/// Writes only the sequence header given ID and description parts.
#[inline]
pub fn write_id_desc<W>(mut writer: W, id: &[u8], desc: Option<&[u8]>) -> io::Result<()>
where
    W: io::Write,
{
    writer.write_all(b">")?;
    writer.write_all(id)?;
    if let Some(d) = desc {
        writer.write_all(b" ")?;
        writer.write_all(d)?;
    }
    writer.write_all(b"\n")
}

/// Writes only the sequence line.
#[inline]
pub fn write_seq<W>(mut writer: W, seq: &[u8]) -> io::Result<()>
where
    W: io::Write,
{
    writer.write_all(seq)?;
    writer.write_all(b"\n")
}

/// Writes the sequence line, and wraps the output to a maximum width specified by `wrap`.
#[inline]
pub fn write_wrap_seq<W>(mut writer: W, seq: &[u8], wrap: usize) -> io::Result<()>
where
    W: io::Write,
{
    assert!(wrap > 0);
    for chunk in seq.chunks(wrap) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// Writes the sequence line from an iterator (such as `SeqLines`)
#[inline]
pub fn write_seq_iter<'a, W, P>(mut writer: W, seq: P) -> io::Result<()>
where
    W: io::Write,
    P: Iterator<Item = &'a [u8]>,
{
    for subseq in seq {
        writer.write_all(subseq)?;
    }
    writer.write_all(b"\n")
}

/// Writes the sequence line from an iterator (such as `SeqLines`) and wraps the output
/// to a maximum width specified by `wrap`.
#[inline]
pub fn write_wrap_seq_iter<'a, W, P>(mut writer: W, seq: P, wrap: usize) -> io::Result<()>
where
    W: io::Write,
    P: IntoIterator<Item = &'a [u8]>,
{
    assert!(wrap > 0);
    let mut n_line = 0;
    for subseq in seq {
        let mut chunk = subseq;
        loop {
            let remaining = wrap - n_line;
            if chunk.len() <= remaining {
                writer.write_all(chunk)?;
                n_line += chunk.len();
                break;
            }
            // chunk longer than line -> break
            let (line, rest) = chunk.split_at(remaining);
            chunk = rest;
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
            n_line = 0;
        }
    }
    writer.write_all(b"\n")?;
    Ok(())
}

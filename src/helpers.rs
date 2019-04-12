use crate::core::{trim_cr, trim_end, SimpleLines};
use std::clone::Clone;
use std::default::Default;
use std::fmt::Debug;
use std::io;
use std::slice;

/// Helper trait  used to allow supplying either the whole head or separate ID
/// and description parts to the  `write_...()` functions in the `fasta` and
/// `fastq` modules.
pub trait HeadWriter {
    /// Writes the header line to output.
    fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write;
}

impl<'a> HeadWriter for &'a [u8] {
    fn write_head<W>(&self, mut writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        write!(writer, "{}", start_byte as char)?;
        writer.write_all(self)?;
        writer.write_all(b"\n")
    }
}

macro_rules! impl_write_head {
    ($t:ty) => {
        impl HeadWriter for $t {
            fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
            where
                W: io::Write,
            {
                (&self[..]).write_head(writer, start_byte)
            }
        }
    };
}

impl_write_head!(&Vec<u8>);
impl_write_head!(&[u8; 1]);
impl_write_head!(&[u8; 2]);
impl_write_head!(&[u8; 3]);
impl_write_head!(&[u8; 4]);
impl_write_head!(&[u8; 5]);
impl_write_head!(&[u8; 6]);
impl_write_head!(&[u8; 7]);
impl_write_head!(&[u8; 8]);
impl_write_head!(&[u8; 9]);
impl_write_head!(&[u8; 10]);
impl_write_head!(&[u8; 11]);
impl_write_head!(&[u8; 12]);
impl_write_head!(&[u8; 13]);
impl_write_head!(&[u8; 14]);
impl_write_head!(&[u8; 15]);
impl_write_head!(&[u8; 16]);
impl_write_head!(&[u8; 17]);
impl_write_head!(&[u8; 18]);
impl_write_head!(&[u8; 19]);
impl_write_head!(&[u8; 20]);

impl<'a> HeadWriter for &'a str {
    fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        self.as_bytes().write_head(writer, start_byte)
    }
}

impl<'a> HeadWriter for &String {
    fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        self.as_str().write_head(writer, start_byte)
    }
}

impl<'a> HeadWriter for (&'a [u8], Option<&'a [u8]>) {
    fn write_head<W>(&self, mut writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        write!(writer, "{}", start_byte)?;
        writer.write_all(self.0)?;
        writer.write_all(b" ")?;
        if let Some(desc) = self.1 {
            writer.write_all(desc)?;
        }
        Ok(())
    }
}

impl<'a> HeadWriter for (&Vec<u8>, Option<&Vec<u8>>) {
    fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        (self.0.as_slice(), self.1.map(|d| d.as_slice())).write_head(writer, start_byte)
    }
}

impl<'a> HeadWriter for (&'a str, Option<&str>) {
    fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        (self.0.as_bytes(), self.1.map(|d| d.as_bytes())).write_head(writer, start_byte)
    }
}

impl<'a> HeadWriter for (&String, Option<&String>) {
    fn write_head<W>(&self, writer: W, start_byte: u8) -> io::Result<()>
    where
        W: io::Write,
    {
        (self.0.as_str(), self.1.map(|d| d.as_str())).write_head(writer, start_byte)
    }
}

/// Holds line number and byte offset of a FASTQ record
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Position {
    pub(crate) line: u64,
    pub(crate) byte: u64,
    pub(crate) record: u64,
}

impl Position {
    #[inline]
    pub fn new() -> Position {
        Position::default()
    }

    /// Line index (0-based)
    #[inline]
    pub fn line(&self) -> u64 {
        self.line
    }

    /// Byte offset within the input
    #[inline]
    pub fn byte(&self) -> u64 {
        self.byte
    }

    /// Record index (0-based) in the input
    #[inline]
    pub fn record(&self) -> u64 {
        self.record
    }

    /// Sets the line index (0-based)
    #[inline]
    pub fn set_line(&mut self, line: u64) -> &mut Self {
        self.line = line;
        self
    }

    /// Sets the byte offset
    #[inline]
    pub fn set_byte(&mut self, byte: u64) -> &mut Self {
        self.byte = byte;
        self
    }

    /// Sets the record index (0-based)
    #[inline]
    pub fn set_record(&mut self, idx: u64) -> &mut Self {
        self.record = idx;
        self
    }
}

// /// Iterator over the lines of a text slice, whose positions
// /// have been searched before. Therefore, iterating is very fast.
// pub struct LinePositionIter<'a> {
//     data: &'a [u8],
//     len: usize,
//     pos_iter: iter::Zip<slice::Iter<'a, usize>, iter::Skip<slice::Iter<'a, usize>>>,
// }

// impl<'a> LinePositionIter<'a> {
//     /// Assumes that line_ends are sorted by position, and no duplicate
//     /// values are present.
//     #[inline]
//     pub fn new(buffer: &'a [u8], line_ends: &'a [usize]) -> Self {
//         let start_iter = line_ends.iter();
//         let end_iter = line_ends.iter().skip(1);
//         LinePositionIter {
//             data: buffer,
//             len: line_ends.len() - 1,
//             pos_iter: start_iter.zip(end_iter),
//         }
//     }
// }

// impl<'a> Iterator for LinePositionIter<'a> {
//     type Item = &'a [u8];

//     #[inline]
//     fn next(&mut self) -> Option<&'a [u8]> {
//         self.pos_iter
//             .next()
//             .map(|(&start, &next_start)| trim_cr(&self.data[start..next_start - 1]))
//     }

//     #[inline]
//     fn size_hint(&self) -> (usize, Option<usize>) {
//         (self.len, Some(self.len))
//     }
// }

// impl<'a> DoubleEndedIterator for LinePositionIter<'a> {
//     #[inline]
//     fn next_back(&mut self) -> Option<&'a [u8]> {
//         self.pos_iter
//             .next_back()
//             .map(|(start, next_start)| trim_cr(&self.data[*start + 1..*next_start]))
//     }
// }

// impl<'a> ExactSizeIterator for LinePositionIter<'a> {}

// TODO: compare performance
/// Iterator over the lines of a text slice, whose positions have been searched
/// before and stored. Iterating is very fast.
pub struct LinePositionIter<'a> {
    data: &'a [u8],
    pos_iter: slice::Windows<'a, usize>,
}

impl<'a> LinePositionIter<'a> {
    /// Assumes that line_ends are sorted by position, and no duplicate
    /// values are present.
    #[inline]
    pub fn new(buffer: &'a [u8], positions: &'a [usize]) -> Self {
        LinePositionIter {
            data: buffer,
            pos_iter: positions.windows(2),
        }
    }
}

impl<'a> Iterator for LinePositionIter<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<&'a [u8]> {
        self.pos_iter
            .next()
            .map(|pos| trim_cr(&self.data[pos[0]..pos[1] - 1]))
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        self.pos_iter.size_hint()
    }
}

impl<'a> DoubleEndedIterator for LinePositionIter<'a> {
    #[inline]
    fn next_back(&mut self) -> Option<&'a [u8]> {
        self.pos_iter
            .next_back()
            .map(|pos| trim_cr(&self.data[pos[0]..pos[1] - 1]))
    }
}

impl<'a> ExactSizeIterator for LinePositionIter<'a> {}

#[derive(Debug)]
enum ParseType<'a> {
    One(Option<&'a [u8]>),
    Many(SimpleLines<'a>),
}

/// Iterator over the lines of a text slice. The line endings are searched in
/// the text while iterating, except for the case where it is known that there
/// is only one line.
///
/// Iterating `> 1` line is thus slow compared to `LinePositionIter`. In case,
/// of a single line, the two iterators should be equally fast.
// TODO: implement ExactSizeIterator?
#[derive(Debug)]
pub struct LineSearchIter<'a> {
    inner: ParseType<'a>,
}

impl<'a> LineSearchIter<'a> {
    #[inline]
    pub fn new(text: &'a [u8], one_line: bool) -> Self {
        let inner = if one_line {
            // trim_end also removes single \r, which should only be fed to this
            // constructor if at the end of the input.
            ParseType::One(Some(trim_end(text)))
        } else {
            // SimpleLines also trims orphan \r at the end of the input.
            // The reason for all this is consistency within the crate
            ParseType::Many(SimpleLines::new(text))
        };
        LineSearchIter { inner }
    }
}

impl<'a> Iterator for LineSearchIter<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self.inner {
            ParseType::One(ref mut l) => l.take(),
            ParseType::Many(ref mut parser) => parser.next(),
        }
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        match self.inner {
            ParseType::One(_) => (1, Some(1)),
            ParseType::Many(ref parser) => parser.size_hint(),
        }
    }
}

impl<'a> DoubleEndedIterator for LineSearchIter<'a> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        match self.inner {
            ParseType::One(ref mut l) => l.take(),
            ParseType::Many(ref mut parser) => parser.next_back(),
        }
    }
}

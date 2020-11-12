use std::borrow::Cow;
use std::io;
use std::iter;
use std::slice;
use std::str;

use super::{write, write_iter, write_wrap, write_wrap_iter, LineStore};
use crate::core::{join_lines, join_lines_given, record_head, record_seq, SeqRecordPosition};
use crate::BaseRecord;
use serde::{Deserialize, Serialize};

/// FASTA record trait implemented by both `RefRecord` and `OwnedRecord`,
/// which adds more methods to [`BaseRecord`](crate::BaseRecord).
pub trait Record: BaseRecord {
    /// Writes the record to the given `io::Write` instance.
    /// The sequence is wrapped to produce multi-line FASTA with a maximum width
    /// specified by `wrap`.
    fn write_wrap<W>(&self, writer: W, wrap: usize) -> io::Result<()>
    where
        W: io::Write,
        Self: Sized;
}

impl<'a, R: Record> Record for &'a R {
    fn write_wrap<W>(&self, writer: W, wrap: usize) -> io::Result<()>
    where
        W: io::Write,
        Self: Sized,
    {
        (**self).write_wrap(writer, wrap)
    }
}

/// A FASTA record that borrows data from a buffer.
/// It implements the traits [`BaseRecord`](crate::BaseRecord) and
/// [`Record`](Record).
#[derive(Debug, Clone)]
pub struct RefRecord<'a, S = LineStore>
where
    S: SeqRecordPosition,
{
    pub(crate) buffer: &'a [u8],
    pub(crate) buf_pos: &'a S,
}

impl<'a, S> BaseRecord for RefRecord<'a, S>
where
    S: SeqRecordPosition,
{
    #[inline]
    fn head(&self) -> &[u8] {
        record_head(self.buf_pos, self.buffer)
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        record_seq(self.buf_pos, self.buffer)
    }

    #[inline]
    fn full_seq(&self) -> Cow<[u8]> {
        join_lines(self.seq_lines(), self.num_seq_lines())
    }

    #[inline]
    fn full_seq_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        join_lines_given(
            self.buf_pos.seq_lines(self.buffer),
            self.buf_pos.num_seq_lines(),
            owned_fn,
        )
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        self.buf_pos.num_seq_lines()
    }

    #[inline]
    fn has_quality(&self) -> bool {
        false
    }

    #[inline]
    fn opt_qual(&self) -> Option<&[u8]> {
        None
    }

    #[inline]
    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        None
    }

    #[inline]
    fn opt_full_qual_given<'s, F>(&'s self, _: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        None
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        0
    }

    #[inline]
    fn write<W>(&self, writer: W) -> io::Result<()>
    where
        W: io::Write,
    {
        write_iter(writer, self.head(), self.seq_lines())
    }
}

impl<'a, S> Record for RefRecord<'a, S>
where
    S: SeqRecordPosition,
{
    #[inline]
    fn write_wrap<W: io::Write>(&self, mut writer: W, wrap: usize) -> io::Result<()> {
        write_wrap_iter(&mut writer, self.head(), self.seq_lines(), wrap)
    }
}

impl<'a, S> RefRecord<'a, S>
where
    S: SeqRecordPosition,
{
    #[inline]
    pub(crate) fn new(buffer: &'a [u8], buf_pos: &'a S) -> Self {
        RefRecord { buffer, buf_pos }
    }

    /// Returns an iterator over all sequence lines in the data. The exact
    /// type of the iterator depends on the generic parameter `S`.
    #[inline]
    pub fn seq_lines(&self) -> impl Iterator<Item = &'a [u8]> + DoubleEndedIterator {
        self.buf_pos.seq_lines(self.buffer)
    }

    /// Returns the sequence as owned `Vec`. **Note**: This function
    /// must be called in order to obtain a sequence that does not contain
    /// line endings (as returned by `seq()`)
    #[inline]
    pub fn owned_seq(&self) -> Vec<u8> {
        let mut seq = Vec::new();
        for segment in self.seq_lines() {
            seq.extend(segment);
        }
        seq
    }

    /// Returns a new [`OwnedRecord`](OwnedRecord) with the same
    /// data.
    #[inline]
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            head: self.head().to_vec(),
            seq: self.owned_seq(),
        }
    }

    /// Copies the data of the record into an [`OwnedRecord`](OwnedRecord)
    /// instance.
    #[inline]
    pub fn clone_into_owned(&self, rec: &mut OwnedRecord) {
        rec.head.clear();
        rec.head.extend(self.head());
        rec.seq.clear();
        self.full_seq_given(|| &mut rec.seq);
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `BaseRecord::write`
    #[inline]
    pub fn write_unchanged<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        let data =
            &self.buffer[self.buf_pos.record_start() as usize..self.buf_pos.record_end() as usize];
        writer.write_all(data)?;
        if *data.last().unwrap() != b'\n' {
            writer.write_all(&[b'\n'])?;
        }
        Ok(())
    }
}

impl_fastx_from_refrecord!(RefRecord);

/// A FASTA record that ownes its data (requiring two allocations)
/// It implements the traits [`BaseRecord`](crate::BaseRecord) and
/// [`Record`](Record).
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
}

impl BaseRecord for OwnedRecord {
    #[inline]
    fn head(&self) -> &[u8] {
        &self.head
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        &self.seq
    }

    #[inline]
    fn full_seq(&self) -> Cow<[u8]> {
        (&self.seq).into()
    }

    #[inline]
    fn full_seq_given<'s, F: FnOnce() -> &'s mut Vec<u8>>(&'s self, _: F) -> Cow<'s, [u8]> {
        (&self.seq).into()
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        1
    }

    #[inline]
    fn has_quality(&self) -> bool {
        false
    }

    #[inline]
    fn opt_qual(&self) -> Option<&[u8]> {
        None
    }

    #[inline]
    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        None
    }

    #[inline]
    fn opt_full_qual_given<'s, F>(&'s self, _: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        None
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        0
    }

    #[inline]
    fn write<W: io::Write>(&self, writer: W) -> io::Result<()> {
        write(writer, &self.head, &self.seq)
    }
}

impl Record for OwnedRecord {
    #[inline]
    fn write_wrap<W: io::Write>(&self, mut writer: W, wrap: usize) -> io::Result<()> {
        write_wrap(&mut writer, &self.head, &self.seq, wrap)
    }
}

impl_recordset!(RefRecord, SeqRecordPosition, LineStore, "fasta", "fasta");

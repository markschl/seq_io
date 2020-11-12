use std::borrow::Cow;
use std::io;
use std::iter;
use std::slice;
use std::str;

use super::{write, write_iter, Error, ErrorKind, RangeStore};
use crate::core::{
    check_lengths, join_lines, join_lines_given, record_head, record_qual, record_seq,
    QualRecordPosition,
};
use crate::{BaseRecord, ErrorPosition};
use serde::{Deserialize, Serialize};

/// FASTQ record trait implemented by both `RefRecord` and `OwnedRecord`
/// which adds more methods to [`BaseRecord`](crate::BaseRecord).
pub trait Record: BaseRecord {
    /// Return the FASTQ quality line as byte slice
    fn qual(&self) -> &[u8];

    fn full_qual(&self) -> Cow<[u8]>;

    fn full_qual_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
        Self: Sized;

    doc_record_check_lengths!(
        "use seq_io::fastq::Reader;",
        fn check_lengths(&self) -> Result<&Self, Error>;
    );
}

// TODO: necessary?

// impl<'a, R: Record> Record for &'a R {
//     fn qual(&self) -> &[u8] {
//         (**self).qual()
//     }

//     fn full_qual(&self) -> Cow<[u8]> {
//         (**self).full_qual()
//     }

//     fn full_qual_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
//     where
//         F: FnOnce() -> &'s mut Vec<u8>,
//         Self: Sized
//     {
//         (**self).full_qual_given(owned_fn)
//     }

//     fn check_lengths(&self) -> Result<&Self, Error> {
//         (**self).check_lengths()
//     }
// }

/// A FASTQ record that borrows data from a buffer
/// It implements the traits [`BaseRecord`](crate::BaseRecord) and
/// [`Record`](Record).
#[derive(Debug, Clone)]
pub struct RefRecord<'a, S = RangeStore>
where
    S: QualRecordPosition,
{
    pub(crate) buffer: &'a [u8],
    pub(crate) buf_pos: &'a S,
}

impl<'a, S> BaseRecord for RefRecord<'a, S>
where
    S: QualRecordPosition,
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
        true
    }

    #[inline]
    fn opt_qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    #[inline]
    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        Some(self.full_qual())
    }

    #[inline]
    fn opt_full_qual_given<'s, F>(&'s self, owned_fn: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        Some(self.full_qual_given(owned_fn))
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        self.buf_pos.num_qual_lines()
    }

    #[inline]
    fn write<W>(&self, writer: W) -> io::Result<()>
    where
        W: io::Write,
        Self: Sized,
    {
        write_iter(writer, self.head(), self.seq_lines(), self.qual_lines())
    }
}

impl<'a, S> Record for RefRecord<'a, S>
where
    S: QualRecordPosition,
{
    #[inline]
    fn qual(&self) -> &[u8] {
        record_qual(self.buf_pos, self.buffer)
    }

    #[inline]
    fn full_qual(&self) -> Cow<[u8]> {
        join_lines(self.qual_lines(), self.num_qual_lines())
    }

    #[inline]
    fn full_qual_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        join_lines_given(
            self.buf_pos.qual_lines(self.buffer),
            self.buf_pos.num_qual_lines(),
            owned_fn,
        )
    }

    #[inline]
    fn check_lengths(&self) -> Result<&Self, Error> {
        self._check_lengths(false)
    }
}

impl<'a, S> RefRecord<'a, S>
where
    S: QualRecordPosition,
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

    /// Returns an iterator over all quality lines in the data. The exact
    /// type of the iterator depends on the generic parameter `S`.
    #[inline]
    pub fn qual_lines(&self) -> impl Iterator<Item = &'a [u8]> + DoubleEndedIterator {
        self.buf_pos.qual_lines(self.buffer)
    }

    #[inline]
    fn _check_lengths(&self, strict: bool) -> Result<&Self, Error> {
        // TODO: remove multiilne_fastq, check by position
        check_lengths(self.buf_pos, self.buffer, strict)
            .map(|_| self)
            .map_err(|(seq, qual)| {
                let id = String::from_utf8_lossy(self.id_bytes()).into();
                let pos = ErrorPosition::new(None, None, Some(id));
                Error::new(ErrorKind::UnequalLengths { pos, seq, qual })
            })
    }

    doc_refrecord_check_lengths_strict!(
        #[inline]
        pub fn check_lengths_strict(&self) -> Result<&Self, Error> {
            self._check_lengths(true)
        }
    );

    /// Returns a new [`OwnedRecord`](OwnedRecord) with the same
    /// data.
    #[inline]
    pub fn to_owned_record(&self) -> OwnedRecord {
        OwnedRecord {
            head: self.head().to_vec(),
            seq: self.full_seq().to_vec(),
            qual: self.full_qual().to_vec(),
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
        rec.qual.clear();
        self.full_qual_given(|| &mut rec.qual);
    }

    /// Writes a record to the given `io::Write` instance
    /// by just writing the unmodified input, which is faster than `BaseRecord::write`
    #[inline]
    pub fn write_unchanged<W: io::Write>(&self, mut writer: W) -> io::Result<()> {
        let data = &self.buffer
            [self.buf_pos.record_start() as usize..self.buf_pos.record_end() as usize - 1];
        writer.write_all(data)?;
        writer.write_all(b"\n")
    }
}

impl_fastx_from_refrecord!(RefRecord);

/// A FASTQ record that ownes its data (requires allocations)
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
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
    fn full_seq_given<'s, F>(&'s self, _: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        (&self.seq).into()
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        1
    }

    #[inline]
    fn has_quality(&self) -> bool {
        true
    }

    #[inline]
    fn opt_qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    #[inline]
    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        Some(self.full_qual())
    }

    #[inline]
    fn opt_full_qual_given<'s, F>(&'s self, _: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        Some((&self.qual).into())
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        1
    }

    #[inline]
    fn write<W>(&self, writer: W) -> io::Result<()>
    where
        W: io::Write,
    {
        write(writer, self.head(), self.seq(), self.qual())
    }
}

impl Record for OwnedRecord {
    #[inline]
    fn qual(&self) -> &[u8] {
        &self.qual
    }

    #[inline]
    fn full_qual(&self) -> Cow<[u8]> {
        (&self.qual).into()
    }

    #[inline]
    fn full_qual_given<'s, F>(&'s self, _: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        (&self.qual).into()
    }

    #[inline]
    fn check_lengths(&self) -> Result<&Self, Error> {
        if self.seq.len() == self.qual.len() {
            return Ok(self);
        }
        let id = String::from_utf8_lossy(self.id_bytes()).into();
        let pos = ErrorPosition::new(None, None, Some(id));
        return Err(Error::new(ErrorKind::UnequalLengths {
            pos,
            seq: self.seq.len(),
            qual: self.qual.len(),
        }));
    }
}

impl_recordset!(RefRecord, QualRecordPosition, RangeStore, "fastq", "fastq");

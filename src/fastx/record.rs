use super::SeqFormat;
use super::{Error, ErrorKind, LineStore};
use crate::core::PositionStore;
use crate::{fasta, fastq, BaseRecord, ErrorPosition};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::io;
use std::iter;
use std::slice;
use std::str;

/// FASTX record trait implemented by both `RefRecord` and `OwnedRecord`
/// which adds more methods to [`BaseRecord`](crate::BaseRecord).
pub trait Record: BaseRecord {
    /// Write the record to the given `io::Write` instance. For FASTA,
    /// the sequence is wrapped to produce multi-line FASTA with a
    /// maximum width specified by `wrap`.
    /// Wrapping FASTQ sequence and quality lines is explicitly not supported
    /// because such files are prone to parsing problems.
    fn write_wrap<W>(&self, writer: W, wrap: usize) -> io::Result<()>
    where
        W: io::Write;

    /// Writes the record to an output given a sequence format.
    ///
    /// FASTA lines can be wrapped by specifying `wrap_fasta = Some(width)`.
    ///
    /// # Panics
    ///
    /// Panics if `SeqFormat::FASTQ` was specified, but there is no quality
    /// information available becase the input is FASTA.
    fn write_as<W>(
        &self,
        writer: W,
        format: SeqFormat,
        wrap_fasta: Option<usize>,
    ) -> io::Result<()>
    where
        W: io::Write;

    doc_record_check_lengths!(
        "use seq_io::fastx::Reader;",
        fn check_lengths(&self) -> Result<&Self, Error>;
    );
}

// impl<'a, R: Record> Record for &'a R {
//     fn write_wrap<W>(&self, writer: W, wrap: usize) -> io::Result<()>
//     where
//         W: io::Write,
//     {
//         (**self).write_wrap(writer, wrap)
//     }

//     fn write_as<W>(&self, writer: W, format: SeqFormat, wrap_fasta: Option<usize>) -> io::Result<()>
//     where
//         W: io::Write,
//     {
//         (**self).write_as(writer, format, wrap_fasta)
//     }

//     fn check_lengths(&self) -> Result<&Self, Error> {
//         (**self).check_lengths()
//     }
// }

/// A FASTX record that borrows data from a buffer
/// It implements the traits [`BaseRecord`](crate::BaseRecord) and
/// [`Record`](Record).
#[derive(Debug, Clone)]
pub struct RefRecord<'a, S = LineStore>
where
    S: PositionStore,
{
    pub(crate) buffer: &'a [u8],
    pub(crate) buf_pos: &'a S,
}

impl<'a, S> BaseRecord for RefRecord<'a, S>
where
    S: PositionStore,
{
    #[inline]
    fn head(&self) -> &[u8] {
        self.buf_pos.head(self.buffer)
    }

    #[inline]
    fn seq(&self) -> &[u8] {
        self.buf_pos.seq(self.buffer)
    }

    #[inline]
    fn full_seq(&self) -> Cow<[u8]> {
        self.buf_pos.join_seq(self.buffer)
    }

    #[inline]
    fn full_seq_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
        Self: Sized,
    {
        self.buf_pos.join_seq_given(self.buffer, owned_fn)
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        self.buf_pos.num_seq_lines()
    }

    #[inline]
    fn has_quality(&self) -> bool {
        self.buf_pos.has_qual()
    }

    #[inline]
    fn opt_qual(&self) -> Option<&[u8]> {
        if self.buf_pos.has_qual() {
            Some(self.buf_pos.qual(self.buffer))
        } else {
            None
        }
    }

    #[inline]
    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        if self.buf_pos.has_qual() {
            Some(self.buf_pos.join_qual(self.buffer))
        } else {
            None
        }
    }

    #[inline]
    fn opt_full_qual_given<'s, F>(&'s self, owned_fn: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        if self.buf_pos.has_qual() {
            Some(self.buf_pos.join_qual_given(self.buffer, owned_fn))
        } else {
            None
        }
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
        if let Some(qual_lines) = self.opt_qual_lines() {
            fastq::write_iter(writer, self.head(), self.seq_lines(), qual_lines)
        } else {
            fasta::write_iter(writer, self.head(), self.seq_lines())
        }
    }
}

impl<'a, S> Record for RefRecord<'a, S>
where
    S: PositionStore,
{
    #[inline]
    fn write_wrap<W>(&self, writer: W, wrap: usize) -> io::Result<()>
    where
        W: io::Write,
        Self: Sized,
    {
        if let Some(qual_lines) = self.opt_qual_lines() {
            fastq::write_iter(writer, self.head(), self.seq_lines(), qual_lines)
        } else {
            fasta::write_wrap_iter(writer, self.head(), self.seq_lines(), wrap)
        }
    }

    #[inline]
    fn write_as<W>(&self, writer: W, format: SeqFormat, wrap_fasta: Option<usize>) -> io::Result<()>
    where
        W: io::Write,
    {
        match format {
            SeqFormat::FASTA => match wrap_fasta {
                None => fasta::write_iter(writer, self.head(), self.seq_lines()),
                Some(w) => fasta::write_wrap_iter(writer, self.head(), self.seq_lines(), w),
            },
            SeqFormat::FASTQ => {
                let qual_lines = self
                    .opt_qual_lines()
                    .expect("no quality information present");
                fastq::write_iter(writer, self.head(), self.seq_lines(), qual_lines)
            }
        }
    }

    #[inline]
    fn check_lengths(&self) -> Result<&Self, Error> {
        self._check_lengths(false)
    }
}

impl<'a, S> RefRecord<'a, S>
where
    S: PositionStore,
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

    /// Returns an iterator over all quality lines in the data (if present).
    /// The exact type of the iterator depends on the generic parameter `S`.
    #[inline]
    pub fn opt_qual_lines(&self) -> Option<impl Iterator<Item = &'a [u8]> + DoubleEndedIterator> {
        if self.buf_pos.has_qual() {
            Some(self.buf_pos.qual_lines(self.buffer))
        } else {
            None
        }
    }

    #[inline]
    fn _check_lengths(&self, strict: bool) -> Result<&Self, Error> {
        self.buf_pos
            .check_lengths(self.buffer, strict)
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
            qual: self.opt_full_qual().map(|q| q.to_vec()),
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
        if self.buf_pos.has_qual() {
            let q = rec.qual.get_or_insert_with(|| Vec::new());
            q.clear();
            self.buf_pos.join_qual_given(self.buffer, || q);
        }
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

macro_rules! impl_fastx_from_refrecord {
    ($other:ident) => {
        impl<'a, S> From<$other<'a, S>> for crate::fastx::RefRecord<'a, S>
        where
            S: crate::PositionStore,
        {
            fn from(rec: $other<'a, S>) -> Self {
                crate::fastx::RefRecord {
                    buffer: rec.buffer,
                    buf_pos: rec.buf_pos,
                }
            }
        }

        impl<'a, S> From<crate::fastx::RefRecord<'a, S>> for $other<'a, S>
        where
            S: crate::PositionStore,
        {
            fn from(rec: crate::fastx::RefRecord<'a, S>) -> Self {
                $other {
                    buffer: rec.buffer,
                    buf_pos: rec.buf_pos,
                }
            }
        }
    };
}

/// A FASTX record that ownes its data (requires allocations)
#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OwnedRecord {
    pub head: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
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
        self.qual.is_some()
    }

    #[inline]
    fn opt_qual(&self) -> Option<&[u8]> {
        self.qual.as_ref().map(|q| q.as_slice())
    }

    #[inline]
    fn opt_full_qual_given<'s, F>(&'s self, _: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        self.qual.as_ref().map(|q| q.as_slice().into())
    }

    #[inline]
    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        self.qual.as_ref().map(|q| q.into())
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        1
    }

    #[inline]
    fn write<W: io::Write>(&self, writer: W) -> io::Result<()> {
        if let Some(q) = self.qual.as_ref() {
            fastq::write(writer, &self.head, &self.seq, q)
        } else {
            fasta::write(writer, &self.head, &self.seq)
        }
    }
}

impl Record for OwnedRecord {
    #[inline]
    fn write_wrap<W>(&self, writer: W, wrap: usize) -> io::Result<()>
    where
        W: io::Write,
    {
        if let Some(q) = self.qual.as_ref() {
            fastq::write(writer, &self.head, &self.seq, q)
        } else {
            fasta::write_wrap(writer, &self.head, &self.seq, wrap)
        }
    }

    #[inline]
    fn write_as<W>(&self, writer: W, format: SeqFormat, wrap_fasta: Option<usize>) -> io::Result<()>
    where
        W: io::Write,
    {
        match format {
            SeqFormat::FASTA => match wrap_fasta {
                None => fasta::write(writer, &self.head, &self.seq),
                Some(w) => fasta::write_wrap(writer, &self.head, &self.seq, w),
            },
            SeqFormat::FASTQ => {
                let qual = self.qual.as_ref().expect("no quality information present");
                fastq::write(writer, &self.head, &self.seq, qual)
            }
        }
    }

    #[inline]
    fn check_lengths(&self) -> Result<&Self, Error> {
        if let Some(qual) = self.qual.as_ref() {
            if self.seq.len() == qual.len() {
                return Ok(self);
            }
            let id = String::from_utf8_lossy(self.id_bytes()).into();
            let pos = ErrorPosition::new(None, None, Some(id));
            return Err(Error::new(ErrorKind::UnequalLengths {
                pos,
                seq: self.seq.len(),
                qual: qual.len(),
            }));
        }
        Ok(self)
    }
}

impl_recordset!(RefRecord, LineStore, "fastx", "fastq");

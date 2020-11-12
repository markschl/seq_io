use memchr::memchr;
use std::borrow::Cow;
use std::io;
use std::str;

pub use crate::core::{QualRecordPosition, SeqRecordPosition};

pub trait BaseRecord {
    /// Return the header line of the record as byte slice
    fn head(&self) -> &[u8];

    /// Return the record sequence as byte slice. With FASTA and multi-line
    /// FASTQ, this slice may contain line terminators if called on `RefRecord`.
    /// Use [`full_seq`](#tymethod.full_seq) or
    /// [`full_seq_given`](#tymethod.full_seq_given) instead in
    /// such cases, or iterate over lines with `fasta/q/x::RefRecord::seq_lines()`
    fn seq(&self) -> &[u8];

    /// Returns the full sequence as `Cow<[u8]>`.
    /// If the sequence consists of a single line, then the sequence will be
    /// borrowed from the underlying buffer (equivalent to calling `seq()`).
    /// If there are multiple lines, an owned copy will be created.
    fn full_seq(&self) -> Cow<[u8]>;

    /// Like `full_seq` returns the full sequence as `Cow<[u8]>`, with the
    /// difference that the vector that is used to store the contiguous
    /// sequence (in the case of multiple sequence lines) is provided in
    /// a closure, which is only executed if needed.
    /// This allows e.g. storing the necessary vectors in an arena allocator.
    fn full_seq_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>;

    /// Returns the number of sequence lines.
    /// Calling `record.seq_lines()` gives the same result as
    /// `record.seq_lines().len()`, but the latter may be slow depending on the
    /// type [`SeqRecordPosition`](crate::core::SeqRecordPosition) used since the line
    /// positions may have to be searched first.
    fn num_seq_lines(&self) -> usize;

    /// Returns a boolean specifying whether there is quality information in
    /// this record or not.
    fn has_quality(&self) -> bool;

    /// Returns the quality information or `None`.
    /// With multi-line FASTQ, this slice may contain line terminators if
    /// called on `RefRecord`.
    /// Use [`opt_full_qual`](#tymethod.opt_full_qual) or
    /// [`opt_full_qual_given`](#tymethod.opt_full_qual_given) instead in
    /// such cases, or iterate over lines with `fasta/q/x::RefRecord::[opt_]qual_lines()`
    fn opt_qual(&self) -> Option<&[u8]>;

    /// Returns the quality scores as contiguous `Cow<[u8]>`, if present. If not
    /// dealing with multi-line FASTQ, this does not involve any copying
    /// (same as if calling `opt_qual()`).
    /// If there are multiple lines, an owned copy will be created.
    fn opt_full_qual(&self) -> Option<Cow<[u8]>>;

    /// Like `opt_full_qual` returns the full quality scores as `Cow<[u8]>`,
    /// with the difference that the vector that is used to store the contiguous
    /// data (in the case of multiple quality lines) is provided in
    /// a closure, which is only executed if needed.
    fn opt_full_qual_given<'s, F>(&'s self, owned_fn: F) -> Option<Cow<'s, [u8]>>
    where
        F: FnOnce() -> &'s mut Vec<u8>;

    /// Returns the number of quality lines. Always > 0 with the current parsers.
    fn num_qual_lines(&self) -> usize;

    /// Writes the record to an output.
    ///
    /// The default for FASTA is to write one sequence line. Otherwise refer to
    /// the method [`fasta::Record::write_wrap`](crate::fasta::Record::write_wrap).
    fn write<W>(&self, writer: W) -> io::Result<()>
    where
        W: io::Write;

    /// Returns the record ID (everything before an optional space) as byte slice.
    #[inline]
    fn id_bytes(&self) -> &[u8] {
        let head = self.head();
        if let Some(pos) = memchr(b' ', head) {
            return &head[..pos];
        }
        head
    }

    /// Returns the record ID (everything before an optional space) as `&str`.
    #[inline]
    fn id(&self) -> Result<&str, str::Utf8Error> {
        str::from_utf8(self.id_bytes())
    }

    /// Returns the record description (separated from the ID by a space)
    /// as byte slice, if present.
    #[inline]
    fn desc_bytes(&self) -> Option<&[u8]> {
        let head = self.head();
        if let Some(pos) = memchr(b' ', head) {
            return Some(&head[pos + 1..]);
        }
        None
    }

    /// Returns the record description (separated from the ID by a space)
    /// as `&str` slice, if present.
    #[inline]
    fn desc(&self) -> Option<Result<&str, str::Utf8Error>> {
        self.desc_bytes().map(str::from_utf8)
    }

    /// Return both the ID and the description of the record (if present)
    /// This should be faster than calling `id()` and `desc()` separately.
    #[inline]
    fn id_desc_bytes(&self) -> (&[u8], Option<&[u8]>) {
        let head = self.head();
        if let Some(pos) = memchr(b' ', head) {
            return (&head[..pos], Some(&head[pos + 1..]));
        }
        (head, None)
    }

    /// Returns both the ID and the description of the record (if present)
    /// as `&str` slices.
    /// This should be faster than calling `id()` and `desc()` separately.
    #[inline]
    fn id_desc(&self) -> Result<(&str, Option<&str>), str::Utf8Error> {
        let (id, desc) = self.id_desc_bytes();
        Ok((str::from_utf8(id)?, desc.map(str::from_utf8).transpose()?))
    }
}

impl<'a, R> BaseRecord for &'a R
where
    R: BaseRecord,
{
    fn head(&self) -> &[u8] {
        (**self).head()
    }

    fn seq(&self) -> &[u8] {
        (**self).seq()
    }

    fn full_seq(&self) -> Cow<[u8]> {
        (**self).full_seq()
    }

    fn full_seq_given<'s, F>(&'s self, owned_fn: F) -> Cow<'s, [u8]>
    where
        F: FnOnce() -> &'s mut Vec<u8>,
    {
        (**self).full_seq_given(owned_fn)
    }

    fn num_seq_lines(&self) -> usize {
        (**self).num_seq_lines()
    }

    fn id_bytes(&self) -> &[u8] {
        (**self).id_bytes()
    }

    fn id(&self) -> Result<&str, std::str::Utf8Error> {
        (**self).id()
    }

    fn desc_bytes(&self) -> Option<&[u8]> {
        (**self).desc_bytes()
    }

    fn desc(&self) -> Option<Result<&str, std::str::Utf8Error>> {
        (**self).desc()
    }

    fn id_desc_bytes(&self) -> (&[u8], Option<&[u8]>) {
        (**self).id_desc_bytes()
    }

    fn id_desc(&self) -> Result<(&str, Option<&str>), std::str::Utf8Error> {
        (**self).id_desc()
    }

    fn has_quality(&self) -> bool {
        (**self).has_quality()
    }

    fn opt_qual(&self) -> Option<&[u8]> {
        (**self).opt_qual()
    }

    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        (**self).opt_full_qual()
    }

    fn opt_full_qual_given<'s, F: FnOnce() -> &'s mut Vec<u8>>(
        &'s self,
        owned_fn: F,
    ) -> Option<Cow<'s, [u8]>> {
        (**self).opt_full_qual_given(owned_fn)
    }

    fn num_qual_lines(&self) -> usize {
        (**self).num_qual_lines()
    }

    fn write<W: io::Write>(&self, writer: W) -> io::Result<()> {
        (**self).write(writer)
    }
}

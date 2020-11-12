// Macro for implementing a set of standard methods and traits for all Readers,
// together with the appropriate documentation
macro_rules! impl_reader {
    (
        $Reader:ident, $RefRecord:ty, $OwnedRecord:ty, $RecordSet:ty, $Error:ty,
        $RecordStoreTrait:path, $DefaultPositionStore:ty,
        $multiline_fasta:expr, $multiline_fastq:expr,
        // specific args needed for documentation
        $format:expr, ($($mod_path:expr),*), $seq1:expr, $seq2:expr,
        [$seq2_owned1:expr, $seq2_owned2:expr]
    ) => {
        _impl_reader!(
            $Reader, $RefRecord, $OwnedRecord, $RecordSet, $Error,
            $RecordStoreTrait, $DefaultPositionStore, $multiline_fasta, $multiline_fastq,
            doc_import!(($($mod_path),*), "{ReaderBuilder, Reader}"),
            doc_import!(($($mod_path),*), "Reader"),
            doc_import!(($($mod_path),*), "Reader", "RecordSet"),
            doc_import!(($($mod_path),*), "Reader", "OwnedRecord"),
            concat!("let seq_path = \"seqs.", $format, "\";"),
            concat!("let seq = b\"", $seq1, "\";"),
            concat!("let seq = b\"", $seq2, "\";"),
            concat!("vec![\n", $seq2_owned1, ",\n", $seq2_owned2, "\n]);"),
            concat!("[", stringify!($RecordStoreTrait), "](", stringify!($RecordStoreTrait), ")"),
            doc_link!(($($mod_path),*), "RefRecord"),
            doc_link!(($($mod_path),*), "Record"),
            doc_link!(($($mod_path),*), "RecordSet")
        );
    }
}

macro_rules! doc_link {
    (($_:expr, $record_path:expr), $obj:expr) => {
        concat!("[`", $obj, "`](crate::", $record_path, "::", $obj, ")")
    };
    (($record_path:expr), $obj:expr) => {
        concat!("[`", $obj, "`](crate::", $record_path, "::", $obj, ")")
    };
}

macro_rules! doc_import {
    (($module_path:expr), $name:expr, $rec_name:expr) => {
        concat!(
            "use seq_io::",
            $module_path,
            "::{",
            $name,
            ", ",
            $rec_name,
            "};"
        )
    };
    (($rdr_path:expr, $record_path:expr), $name:expr, $rec_name:expr) => {
        concat!(
            "use seq_io::",
            $rdr_path,
            "::",
            $name,
            ";\nuse seq_io::",
            $record_path,
            "::",
            $rec_name,
            ";"
        )
    };
    (($module_path:expr), $name:expr) => {
        concat!("use seq_io::", $module_path, "::", $name, ";")
    };
    (($rdr_path:expr, $_:expr), $name:expr) => {
        concat!("use seq_io::", $rdr_path, "::", $name, ";")
    };
}

macro_rules! impl_reader_builder {
    ($Reader:ident, $RecordStoreTrait:path, $DefaultPositionStore:ty,
    $import_builder:expr, $seqfile:expr, $record_store_link:expr) => {

/// Allows building a [`Reader`](Reader)
/// with various configuration options.
#[derive(Debug)]
pub struct ReaderBuilder<P: crate::policy::BufPolicy = crate::policy::StdPolicy, S: $RecordStoreTrait = $DefaultPositionStore> {
    buf_policy: P,
    _store: std::marker::PhantomData<S>,
    capacity: usize,
}

impl ReaderBuilder {
    /// Creates a new `ReaderBuilder`
    #[inline]
    pub fn new() -> Self {
        ReaderBuilder {
            buf_policy: crate::policy::StdPolicy,
            _store: std::marker::PhantomData,
            capacity: crate::core::BUFSIZE,
        }
    }
}

impl<P, S> ReaderBuilder<P, S>
where
    P: crate::policy::BufPolicy,
    S: $RecordStoreTrait,
{
    /// Creates a new [`Reader`](Reader)
    /// with the current configuration of the builder.
    #[inline]
    pub fn from_reader<'s, R>(self, reader: R) -> $Reader<R, P, S>
    where
        R: std::io::Read + 's,
        P: 's,
        S: 's,
    {
        $Reader::_new(reader, self.capacity, self.buf_policy)
    }

    /// Creates a new [`Reader`](Reader) from the given path
    /// with the current configuration of the builder.
    #[inline]
    pub fn from_path<'s, Q>(
        self,
        path: Q,
    ) -> std::io::Result<$Reader<std::fs::File, P, S>>
    where
        Q: AsRef<std::path::Path> + 's,
        P: 's,
        S: 's,
    {
        std::fs::File::open(path).map(|f| self.from_reader(f))
    }

    /// Creates a new reader with a given
    #[doc = $record_store_link]
    /// as defined in the type argument of the method. The method consumes
    /// the reader and returns a new `Reader` instance.
    ///
    /// # Example:
    ///
    /// This example makes the reader use `seq_io::fastx::LineStore` instead
    /// of the default one.
    ///
    /// ```no_run
    /// # fn main()  -> Result<(), std::io::Error> {
    /// use seq_io::prelude::*;
    #[doc = $import_builder]
    /// use seq_io::fastx::LineStore;
    ///
    #[doc = $seqfile]
    ///
    /// let mut reader = ReaderBuilder::new()
    ///     .pos_store::<LineStore>()
    ///     .from_path(seq_path)?;
    ///
    /// // or:
    ///
    /// let mut reader: Reader<_, _, LineStore> = ReaderBuilder::new()
    ///     .pos_store::<LineStore>()
    ///     .from_path(seq_path)?;
    /// // ...
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn pos_store<T: $RecordStoreTrait>(self) -> ReaderBuilder<P, T> {
        ReaderBuilder {
            buf_policy: self.buf_policy,
            _store: std::marker::PhantomData,
            capacity: self.capacity,
        }
    }

    /// Applies a [`BufPolicy`](crate::policy::BufPolicy) to the
    /// current reader. The method consumes the reader and returns a new
    /// `Reader` instance.
    #[inline]
    pub fn buf_policy<T: crate::policy::BufPolicy>(self, buf_policy: T) -> ReaderBuilder<T, S> {
        ReaderBuilder {
            buf_policy,
            _store: self._store,
            capacity: self.capacity,
        }
    }

    /// Sets a reader capacity.
    #[inline]
    pub fn capacity(mut self, cap: usize) -> Self {
        self.capacity = cap;
        self
    }
}
}
}

// Inner macro for creating Reader. Needs a few documentation lines to be
// pre-assembled.
macro_rules! _impl_reader {
    (
        $Reader:ident, $RefRecord:ty, $OwnedRecord:ty, $RecordSet:ty, $Error:ty,
        $RecordStoreTrait:path, $DefaultPositionStore:ty,
        $multiline_fasta:expr, $multiline_fastq:expr,
        // the following params are used for documentation
        $import_builder:expr, $import_rdr:expr, $import_rset:expr, $import_owned:expr,
        $seqfile:expr, $seq1:expr, $seq2:expr, $seq2_owned_vec:expr,
        $record_store_link:expr,
        $refrec_link:expr, $record_link:expr, $rset_link:expr
    ) => {

impl_records_iter!($Reader<R, P, S>, $RecordStoreTrait, $OwnedRecord, $Error);

impl_reader_builder!($Reader, $RecordStoreTrait, $DefaultPositionStore,
    $import_builder, $seqfile, $record_store_link);

impl<R> $Reader<R>
where
    R: std::io::Read,
{
    /// Creates a new reader with the default buffer size of 64 KiB
    ///
    /// # Example:
    ///
    /// ```
    /// use seq_io::prelude::*;
    #[doc = $import_rdr]
    #[doc = $seq1]
    ///
    /// let mut reader = Reader::new(&seq[..]);
    /// let record = reader.next().unwrap().unwrap();
    /// assert_eq!(record.id(), Ok("id"))
    /// ```
    #[inline]
    pub fn new(reader: R) -> Self {
        Self::with_capacity(reader, crate::core::BUFSIZE)
    }

    /// Creates a new reader with a given buffer capacity. The minimum allowed
    /// capacity is 3.
    #[inline]
    pub fn with_capacity(reader: R, capacity: usize) -> Self {
        Self::_new(reader, capacity, crate::policy::StdPolicy)
    }
}

impl $Reader<std::fs::File> {
    /// Creates a reader from a file path.
    ///
    /// # Example:
    ///
    /// ```no_run
    #[doc = $import_rdr]
    ///
    #[doc = $seqfile]
    /// let mut reader = Reader::from_path(seq_path).expect("File could not be opened.");
    ///
    /// // (... do something with the reader)
    /// ```
    #[inline]
    pub fn from_path<P: AsRef<std::path::Path>>(seq_path: P) -> std::io::Result<Reader<std::fs::File>> {
        std::fs::File::open(seq_path).map(Reader::new)
    }
}


impl<R, S> $Reader<R, crate::policy::StdPolicy, S>
where
    R: std::io::Read,
    S: $RecordStoreTrait,
{
    /// Creates a new reader with a custom
    #[doc = $record_store_link]
    /// type used for storing positional information. If more customization is
    /// required, use [`ReaderBuilder`](ReaderBuilder) instead.
    ///
    /// # Example:
    ///
    /// This example makes the reader use `seq_io::fastx::LineStore` instead
    /// of the default one.
    ///
    /// ```no_run
    /// # fn main()  -> Result<(), std::io::Error> {
    #[doc = $import_rdr]
    /// use seq_io::fastx::LineStore;
    ///
    #[doc = $seq1]
    ///
    /// let mut reader: Reader<_, _, LineStore> = Reader::with_pos_store(&seq[..]);
    /// // ...
    /// # Ok(())
    /// # }
    #[inline]
    pub fn with_pos_store(reader: R) -> Self {
        Self::_new(reader, crate::core::BUFSIZE, crate::policy::StdPolicy)
    }
}


impl<R, P> $Reader<R, P>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
{
    /// Creates a new reader with the given [`BufPolicy`](crate::policy::BufPolicy)
    /// applied. If more customization is required, use [`ReaderBuilder`](ReaderBuilder)
    /// instead.
    #[inline]
    pub fn with_buf_policy(reader: R, policy: P) -> Self {
        Self::_new(reader, crate::core::BUFSIZE, policy)
    }
}


impl<R, P, S> $Reader<R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: $RecordStoreTrait,
{
    /// Returns a reference to the underlying `BufPolicy` of the reader
    #[inline]
    pub fn buf_policy(&self) -> &P {
        self.inner.policy()
    }

    /// Creates a new Reader from an already instantiated
    /// [`BufReader`](crate::core::BufReader).
    /// This is mostly useful if doing format recognition before deciding on
    /// the exact reader type.
    #[inline]
    pub fn from_buf_reader(rdr: crate::core::BufReader<R, P>, byte_offset: usize, line_idx: u64) -> Self {
        Self::_from_buf_reader(rdr, byte_offset, line_idx)
    }

    #[inline]
    fn _next(&mut self, check_lengths: bool) -> Option<crate::fastx::Result<(&[u8], &S)>> {
        if let Some(fasta) = try_opt!(self._check_is_fasta()) {
            self.inner.next(fasta, $multiline_fasta, $multiline_fastq, false, check_lengths)
        } else {
            None
        }
    }

    /// Searches the next record and returns a
    #[doc = $refrec_link]
    /// that
    /// borrows its data from the underlying buffer of this reader.
    ///
    /// # Example:
    ///
    /// ```no_run
    /// use seq_io::prelude::*;
    #[doc = $import_rdr]
    ///
    #[doc = $seqfile]
    /// let mut reader = Reader::from_path(seq_path).unwrap();
    ///
    /// while let Some(record) = reader.next() {
    ///     let record = record.unwrap();
    ///     println!("{}", record.id().unwrap());
    /// }
    /// ```
    #[inline]
    pub fn next(&mut self) -> Option<super::Result<$RefRecord>> {
        self._next(true).map(|res| {
            res.map(|(pos, buf)| <$RefRecord>::new(pos, buf))
                .map_err(std::convert::From::from)
        })
    }

    /// Searches the next FASTQ record and returns a
    #[doc = $refrec_link]
    /// exactly like [`next()`](Reader::next), but this function does not
    /// compare the lengths of the sequence and quality scores. The check can
    /// (and should!) be done later using [`check_lengths()`](
    #[doc = $record_link]
    /// ).
    ///
    /// **Note**: not checking lengths this with multi-line FASTQ is a very bad
    /// idea, since the parser uses the sequence length as orientation when
    /// parsing the quality scores; quality lengths will never be smaller,
    /// but they can be larger than sequence lengths, and not doing a length
    /// check may lead to confusing errors later, or in the worst case
    /// "swallowing" of a record.
    #[inline]
    pub fn next_unchecked_len(&mut self) -> Option<super::Result<$RefRecord>> {
        self._next(false).map(|res| {
            res.map(|(pos, buf)| <$RefRecord>::new(pos, buf))
                .map_err(std::convert::From::from)
        })
    }

    // needs to be generic in order to allow an easy FastxReader impl
    #[inline]
    fn _read_record_set<T>(&mut self, record_set: &mut T, check_lengths: bool, n_records: Option<usize>) -> super::Result<bool>
    where
        T: crate::core::RecordSet<S>
    {
        if let Some(fasta) = self._check_is_fasta()? {
            self.inner.read_record_set(record_set, fasta, $multiline_fasta, $multiline_fastq, false, check_lengths, n_records)
            .map_err(From::from)
        } else {
            Ok(false)
        }
    }

    /// Updates a
    #[doc = $rset_link]
    /// with new data.
    /// The contents of the internal buffer are just copied over to the record
    /// set and the positions of all records are found.
    /// Old data will be erased. Returns `None` if the input reached its end.
    /// Iterating over a record set will yield
    #[doc = $refrec_link]
    /// instances.
    ///
    /// # Example:
    ///
    /// ```
    /// use seq_io::prelude::*;
    #[doc = $import_rset]
    ///
    #[doc = $seq2]
    /// // Read records the normal way (for illustration, actually
    /// // we could also use the Reader::records iterator here)
    /// let mut reader = Reader::new(&seq[..]);
    /// let mut records = vec![];
    /// while let Some(res) = reader.next() {
    ///     records.push(res.unwrap().to_owned_record());
    /// }
    ///
    /// // Use read_record_set
    /// let mut reader = Reader::new(&seq[..]);
    /// let mut record_set = RecordSet::default(); // initialize once, reuse later
    /// let mut rset_records = vec![];
    /// while reader.read_record_set(&mut record_set).unwrap() {
    ///     for record in &record_set {
    ///         rset_records.push(record.to_owned_record());
    ///     }
    /// }
    ///
    /// // The two methods should give the same result
    /// assert_eq!(records, rset_records);
    /// ```
    #[inline]
    pub fn read_record_set(&mut self, record_set: &mut $RecordSet) -> super::Result<bool> {
        self._read_record_set(record_set, true, None).map_err(From::from)
    }

    #[inline]
    pub fn read_record_set_exact(&mut self, record_set: &mut $RecordSet, n_records: usize) -> super::Result<bool> {
        self._read_record_set(record_set, true, Some(n_records)).map_err(From::from)
    }

    /// Returns a borrowed iterator over all FASTA records. The records
    /// are owned (`OwnedRecord`), this is therefore slower than using
    /// `Reader::next()`.
    ///
    /// # Example
    ///
    /// ```
    /// # fn main()  -> Result<(), std::io::Error> {
    /// use seq_io::prelude::*;
    #[doc = $import_owned]
    ///
    #[doc = $seq2]
    ///
    /// let mut reader = Reader::new(&seq[..]);
    ///
    /// let records: Result<Vec<_>, _> = reader
    ///     .records()
    ///     .collect();
    ///
    /// assert_eq!(records?,
    #[doc = $seq2_owned_vec]
    /// # Ok(())
    /// # }
    /// ```
    #[inline]
    pub fn records(&mut self) -> RecordsIter<R, P, S> {
        RecordsIter::new(self)
    }

    /// Returns an iterator over all records like `Reader::records()`,
    /// but with the difference that it owns the underlying reader.
    #[inline]
    pub fn into_records(self) -> RecordsIntoIter<R, P, S> {
        RecordsIntoIter::new(self)
    }

    /// Returns the position of the current record (found by the previous call
    /// to [`next()`](#next)). Useful with [`seek()`](#seek).
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() {
    /// use seq_io::prelude::*;
    #[doc = $import_rdr]
    /// use seq_io::Position;
    ///
    #[doc = $seq1]
    ///
    /// let mut reader = Reader::new(&seq[..]);
    ///
    /// // find first record
    /// reader.next().unwrap();
    ///
    /// assert_eq!(
    ///     &reader.position(),
    ///     Position::new().set_line(1).set_byte(1).set_record(0)
    /// );
    /// # }
    /// ```
    #[inline]
    pub fn position(&self) -> crate::Position {
        self.inner.position()
    }
}


impl<R, P, S> $Reader<R, P, S>
where
    R: std::io::Read + std::io::Seek,
    P: crate::policy::BufPolicy,
    S: $RecordStoreTrait,
{
    /// Seeks to a specified position.
    /// Keeps the underyling buffer if the seek position is found within it,
    /// otherwise it has to be discarded.
    /// If an error was returned before, seeking to that position will return
    /// the same error. The same is not always true with `None`. If there is no
    /// newline character at the end of the file, the last record will be
    /// returned instead of `None`.
    ///
    /// # Example
    ///
    /// ```
    /// # fn main() {
    /// use seq_io::Position;
    /// use seq_io::prelude::*;
    #[doc = $import_rdr]
    /// use std::io::Cursor;
    ///
    #[doc = $seq2]
    ///
    /// let mut cursor = Cursor::new(&seq[..]);
    /// let mut reader = Reader::new(cursor);
    ///
    /// // Read the first record and get its position.
    /// let record0 = reader.next().unwrap().unwrap().to_owned_record();
    /// let pos0 = reader.position();
    ///
    /// // Read the second record.
    /// reader.next().unwrap().unwrap();
    ///
    /// // Now seek to position of first record.
    /// reader.seek(&pos0);
    ///
    /// // The two records should be the same.
    /// let record1 = reader.next().unwrap().unwrap().to_owned_record();
    /// assert_eq!(record0, record1);
    /// # }
    #[inline]
    pub fn seek(&mut self, pos: &crate::Position) -> std::io::Result<()> {
        self.inner.seek(pos, $multiline_fastq)
    }
}

impl<R, P, S> crate::fastx::dynamic::FastxReader<R, P, S> for $Reader<R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: crate::core::QualRecordPosition,
{
    #[inline]
    fn next_fastx(&mut self) -> Option<crate::fastx::Result<crate::fastx::RefRecord<S>>> {
        self._next(true).map(|res| {
            res.map(|(pos, buf)| crate::fastx::RefRecord::new(pos, buf))
        })
    }

    #[inline]
    fn read_record_set_fastx(
        &mut self,
        record_set: &mut crate::fastx::RecordSet<S>,
    ) -> crate::fastx::Result<bool> {
        self._read_record_set(record_set, true, None).map_err(From::from)
    }

    #[inline]
    fn read_record_set_exact_fastx(
        &mut self,
        record_set: &mut crate::fastx::RecordSet<S>,
        n_records: usize,
    ) -> crate::fastx::Result<bool> {
        self._read_record_set(record_set, true, Some(n_records)).map_err(From::from)
    }

    #[inline]
    fn format(&self) -> Option<crate::fastx::SeqFormat> {
        self._format()
    }

    #[inline]
    fn fastx_records(&mut self) -> crate::fastx::dynamic::RecordsIter<R, P, S> {
        crate::fastx::dynamic::RecordsIter::new(self)
    }

    #[inline]
    fn into_fastx_records<'a>(self) -> crate::fastx::dynamic::RecordsIntoIter<'a, R, P, S>
    where
        R: 'a,
        P: 'a,
        S: 'a,
    {
        crate::fastx::dynamic::RecordsIntoIter::new(Box::new(self))
    }

    #[inline]
    fn position_fastx(&self) -> crate::Position {
        self.position()
    }
}

impl<R, P, S> crate::fastx::dynamic::FastxSeekReader<R, P, S> for $Reader<R, P, S>
where
    R: std::io::Read + std::io::Seek,
    P: crate::policy::BufPolicy,
    S: crate::core::QualRecordPosition,
{
    #[inline]
    fn seek_fastx(&mut self, pos: &crate::Position) -> std::io::Result<()> {
        self.seek(pos)
    }
}

}
}

use crate::PositionStore;

// used only internally
pub(crate) trait RecordSet<S>
where
    S: PositionStore,
{
    fn clear(&mut self);
    fn set_buffer(&mut self, buffer: &[u8]);
    fn set_next_pos(&mut self, pos: &S);
    fn is_empty(&self) -> bool;
    fn len(&self) -> usize;
}

macro_rules! impl_recordset {
    ($RefRecord:ident, $DefaultStore:ident, $module_path:expr, $format:expr) => {
        _impl_recordset!(
            $RefRecord,
            $DefaultStore,
            concat!("use seq_io::", $module_path, "::{Reader, RecordSet};"),
            concat!("let seq_path = \"sequences.", $format, "\";")
        );
    };
}

macro_rules! _impl_recordset {
    ($RefRecord:ident, $DefaultStore:ident, $imports:expr, $seq_path:expr) => {
        /// Set of sequence records that owns it's buffer
        /// and knows the positions of each record.
        ///
        /// # Example:
        ///
        /// ```no_run
        /// use seq_io::prelude::*;
        #[doc = $imports]
        ///
        #[doc = $seq_path]
        /// let mut reader = Reader::from_path(seq_path).unwrap();
        ///
        /// while let Some(record) = reader.next() {
        ///     let record = record.unwrap();
        ///     println!("{}", record.id().unwrap());
        /// }
        /// ```
        #[derive(Clone, Debug, Serialize, Deserialize, Default)]
        pub struct RecordSet<S = $DefaultStore>
        where
            S: crate::PositionStore,
        {
            buffer: Vec<u8>,
            positions: Vec<S>,
            npos: usize,
        }

        impl<S> RecordSet<S>
        where
            S: crate::PositionStore,
        {
            /// Returns the number of records in the record set.
            #[inline]
            pub fn len(&self) -> usize {
                self.npos
            }

            /// Returns the number of records in the record set.
            #[inline]
            pub fn is_empty(&self) -> bool {
                self.npos == 0
            }
        }

        impl<S> crate::core::RecordSet<S> for RecordSet<S>
        where
            S: crate::PositionStore,
        {
            #[inline]
            fn clear(&mut self) {
                self.npos = 0;
            }

            #[inline]
            fn set_buffer(&mut self, buffer: &[u8]) {
                self.buffer.clear();
                self.buffer.extend(buffer);
            }

            #[inline]
            fn set_next_pos(&mut self, pos: &S) {
                let i = self.npos;
                self.npos += 1;
                if let Some(p) = self.positions.get_mut(i) {
                    p.clone_from(pos);
                } else {
                    debug_assert!(self.npos == self.positions.len() + 1);
                    self.positions.push(S::default());
                    self.positions[i].clone_from(pos);
                }
            }

            #[inline]
            fn len(&self) -> usize {
                RecordSet::len(self)
            }

            #[inline]
            fn is_empty(&self) -> bool {
                RecordSet::is_empty(self)
            }
        }

        impl<'a, S> iter::IntoIterator for &'a RecordSet<S>
        where
            S: crate::PositionStore,
        {
            type Item = $RefRecord<'a, S>;
            type IntoIter = RecordSetIter<'a, S>;
            #[inline]
            fn into_iter(self) -> Self::IntoIter {
                RecordSetIter {
                    buffer: &self.buffer,
                    pos: self.positions[..self.npos].iter(),
                }
            }
        }

        /// Iterator over record sets
        pub struct RecordSetIter<'a, S = $DefaultStore>
        where
            S: crate::PositionStore,
        {
            buffer: &'a [u8],
            pos: slice::Iter<'a, S>,
        }

        impl<'a, S> Iterator for RecordSetIter<'a, S>
        where
            S: crate::PositionStore,
        {
            type Item = RefRecord<'a, S>;

            #[inline]
            fn next(&mut self) -> Option<RefRecord<'a, S>> {
                self.pos.next().map(|p| RefRecord {
                    buffer: self.buffer,
                    buf_pos: p,
                })
            }
        }
    };
}

macro_rules! impl_records_iter {
    ($Reader:ty, $OwnedRecord:ty, $Error:ty) => {
        /// Borrowed iterator of `OwnedRecord`
        pub struct RecordsIter<'a, R, P, S>
        where
            P: crate::policy::BufPolicy + 'a,
            R: std::io::Read + 'a,
            S: crate::core::PositionStore + 'a,
        {
            rdr: &'a mut $Reader,
        }

        impl<'a, R, P, S> RecordsIter<'a, R, P, S>
        where
            R: std::io::Read,
            P: crate::policy::BufPolicy,
            S: crate::core::PositionStore,
        {
            pub(crate) fn new(rdr: &'a mut $Reader) -> Self {
                RecordsIter { rdr }
            }
        }

        impl<'a, R, P, S> Iterator for RecordsIter<'a, R, P, S>
        where
            P: crate::policy::BufPolicy + 'a,
            R: std::io::Read + 'a,
            S: crate::core::PositionStore,
        {
            type Item = std::result::Result<$OwnedRecord, $Error>;
            fn next(&mut self) -> Option<Self::Item> {
                self.rdr.next().map(|rec| rec.map(|r| r.to_owned_record()))
            }
        }

        /// Iterator of `OwnedRecord` that owns the underlying reader
        pub struct RecordsIntoIter<R, P, S>
        where
            R: std::io::Read,
            P: crate::policy::BufPolicy,
            S: crate::core::PositionStore,
        {
            rdr: $Reader,
            _s: std::marker::PhantomData<S>,
        }

        impl<R, P, S> RecordsIntoIter<R, P, S>
        where
            R: std::io::Read,
            P: crate::policy::BufPolicy,
            S: crate::core::PositionStore,
        {
            pub(crate) fn new(rdr: $Reader) -> Self {
                RecordsIntoIter {
                    rdr,
                    _s: std::marker::PhantomData,
                }
            }
        }

        impl<R, P, S> Iterator for RecordsIntoIter<R, P, S>
        where
            P: crate::policy::BufPolicy,
            R: std::io::Read,
            S: crate::core::PositionStore,
        {
            type Item = std::result::Result<$OwnedRecord, $Error>;
            fn next(&mut self) -> Option<Self::Item> {
                self.rdr.next().map(|rec| rec.map(|r| r.to_owned_record()))
            }
        }
    };
}

// Macros for documenting methods common to the FASTQ and FASTX Record traits

macro_rules! doc_record_check_lengths {
    ($import_rdr:expr, $($tt:tt)*) => {
        /// Checks if the lengths of the sequence and quality lines
        /// are equal.
        /// If not, returns an error of `ErrorKind::UnequalLengths`.
        ///
        /// # RefRecord implementation
        ///
        /// The [`RefRecord`] implementation (which is also used internally
        /// when `Reader::next` is called) defaults to not checking line ends,
        /// which is faster, but may lead to wrong results in very unusual
        /// cases.
        /// More precisely, `RefRecord::check_lengths()` first compares lengths
        /// without removing the line endings.
        /// If lengths are different, the check is repeated with line endings.
        /// The downside of this strategy is that it can theoretically lead to
        /// records with different line lengths being accepted (even if it's a
        /// very unlikely scenario, see example below).
        ///
        /// If a more strict check for this corner case is required, use
        /// [`check_lengths_strict`](RefRecord::check_lengths_strict).
        ///
        /// # Example
        ///
        /// ```
        /// use seq_io::prelude::*;
        #[doc = $import_rdr]
        ///
        /// let seq = b"@id\nAGT\n+\nII\r\n";
        ///
        /// let mut reader = Reader::new(&seq[..]);
        /// let record = reader.next().unwrap().unwrap();
        ///
        /// // check not failing despite lengths being different
        /// assert!(record.check_lengths().is_ok());
        ///
        /// // more strict check fails
        /// assert!(record.check_lengths_strict().is_err());
        /// ```
        ///
        /// # Note on errors
        /// Record objects do not have any information about their
        /// position, within the file. Therefore, `Error::position()` will
        /// return `None`. `Error::record_id()` is always defined.
        /// In contrast, errors of `ErrorKind::UnequalLengths` returned by
        /// `Reader::next` will contain the correct position.
        $($tt)*
    };
}

macro_rules! doc_refrecord_check_lengths_strict {
    ($($tt:tt)*) => {
        /// Checks if the lengths of the sequence and quality lines
        /// are equal.
        /// If not, returns an error of `ErrorKind::UnequalLengths`.
        ///
        /// In contrast to  [`RefRecord::check_lengths`](RefRecord::check_lengths),
        /// this method always removes line terminators before comparing the
        /// lengths. Therefore, it will always return the correct result, even
        /// if there is a mix of `CR` and `CRLF` line ends.
        $($tt)*
    };
}

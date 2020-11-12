use std::borrow::Cow;
use std::fmt::Debug;
use std::ops::Deref;

use crate::core::trim_cr;
use crate::{LinePositionIter, LineSearchIter};

/// Indicates, which part of the sequence record is being searched. Used in
/// methods of [SeqRecordPosition](SeqRecordPosition)
#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub enum SearchPos {
    HEAD,
    SEQ,
    SEP,
    QUAL,
}

/// Trait for type constructors of line iterators.
///
/// The iterator types cannot be directly specified as associated types of
/// `SeqRecordPosition` / `QaulRecordPosition` because their lifetime would
/// pollute the trait signature.
/// Therefore, this intermediate type is used to generate the final type. For
/// more information see this article on the
/// ["family trait pattern"](http://lukaskalbertodt.github.io/2018/08/03/solving-the-generalized-streaming-iterator-problem-without-gats.html#workaround-b-hrtbs--the-family-trait-pattern).
pub trait LinesIterKind<'a> {
    type Out: Iterator<Item = &'a [u8]> + DoubleEndedIterator;
}

/// Type constructor for [`LineSearchIter`](LineSearchIter).
pub struct LineSearchKind;

impl<'a> LinesIterKind<'a> for LineSearchKind {
    type Out = LineSearchIter<'a>;
}

/// Type constructor for [`LinePositionIter`](LinePositionIter).
pub struct LinePositionIterKind;

impl<'a> LinesIterKind<'a> for LinePositionIterKind {
    type Out = LinePositionIter<'a>;
}

/// Trait for objects storing the coordinates of sequence records in the buffer.
///
/// The minimal methods required for FASTA parsing have to be implemented,
/// FASTQ methods are optional (but should be implemented if storing quality
/// data). A position store can of course choose to ignore some data by not
/// doing anything in a method call.
pub trait SeqRecordPosition: Debug + Clone + Default {
    type SeqLinesKind: for<'a> LinesIterKind<'a>;

    /// Initialize the record. Called before set_record_start(). This is the
    /// place to clear information from the previous record if necessary.
    #[inline]
    fn init(&mut self) {}

    /// Sets the start index of the record, always called first, may be
    /// called multiple times while skipping empty lines.
    fn set_record_start(&mut self, start: usize);

    /// Returns the byte index of the record start
    fn record_start(&self) -> usize;

    /// Sets the start index of the sequence.
    fn set_seq_start(&mut self, pos: usize);

    /// Adds another sequence line start index.
    fn add_seq_line_start(&mut self, pos: usize);

    /// Returns the byte index of the sequence start. If there is no sequence,
    /// this should return the start of the next line after the header, which
    /// would be the sequence start if there was one.
    fn seq_start(&self) -> usize;

    /// Returns a slice of all sequence line start indexes. If there is no
    /// sequence, this returns an empty slice.
    fn seq_starts(&self) -> &[usize];

    /// Returns the byte index of the next line after the sequence ends.
    fn seq_end(&self) -> usize;

    /// Sets the record end. Happens only if parsing was successful up to here,
    /// for FASTA it will be called after set_seq_start, [add_seq_line, ] and
    /// set_sep_pos. The index of the record end is actually the index of the
    /// actual record end (\n) + 1, which normally also is the potential start
    /// of the next record.
    /// A slice of FASTA sequence `&buffer[pos.seq_start()..pos.record_end()]`
    /// would thus include `\n` at the end, usually one would therefore use
    /// `&buffer[pos.seq_start()..pos.record_end() - 1]`. Including the newline
    /// may also fail at the end of the input if there is no newline (then, the
    /// record_end would be buffer length + 2!.
    /// - The has_line argument is usually true, it can be false if the record
    ///   end is the same as the sequence start for FASTA, meaning that there
    ///   is no sequence line. This is important to know, because the position
    ///   may assume a sequence line when set_seq_start() was called, but the
    ///   parser later finds then finds out that there actually isn't any.
    ///   For FASTQ, this is not an issue because missing lines are not
    ///   accepted.
    fn set_record_end(&mut self, end: usize, has_line: bool);

    /// Returns the byte index of the record end, which is the index of the
    /// next line after the record ends.
    fn record_end(&self) -> usize;

    // TODO: not seq. specific, anyway remove?
    /// Returns the number of lines that the whole record covers. Used for
    /// incrementing the line counter of the parser, it should thus be correct!
    fn num_lines(&self) -> usize;

    /// Returns the number of lines that the sequence has. Should be correct,
    /// used when iterating, and a number > 1 will just take a slice of the
    /// whole sequence.
    fn num_seq_lines(&self) -> usize;

    /// Returns an iterator over sequence lines
    fn seq_lines<'s>(&'s self, buffer: &'s [u8]) -> <Self::SeqLinesKind as LinesIterKind<'s>>::Out;

    /// Applies a positive or negative offset to all stored coordinates.
    ///
    /// The `search_pos` argument is only supplied if searching got stuck due to
    /// reaching the end of the buffer. It describes the part of the record,
    /// in which searching aborted.
    ///
    /// All properties that will only be defined later in the searching process
    /// may *not* be touched.
    fn apply_offset(&mut self, offset: isize, search_pos: Option<SearchPos>);

    // Returns the line offset, at which the parser is at the moment, used
    // for getting positional information for an UnexpectedEnd error.
    // Has_line has the same meaning as in set_record_end.
    fn line_offset(&self, search_pos: SearchPos, has_line: bool) -> usize;

    /// Returns whether the position is complete in the sense that
    /// `set_record_end()` was called. In some cases, this may not be easy to
    /// determine, then `true` can be returned if a call to `record_end()` will
    /// not panic.
    #[inline]
    fn maybe_complete(&self) -> bool {
        self.record_end() > self.record_start()
    }
}

/// Trait for objects storing the coordinates of sequence records in the buffer.
///
/// The minimal methods required for FASTA parsing have to be implemented,
/// FASTQ methods are optional (but should be implemented if storing quality
/// data). A position store can of course choose to ignore some data by not
/// doing anything in a method call.
pub trait QualRecordPosition: SeqRecordPosition {
    type QualLinesKind: for<'a> LinesIterKind<'a>;

    /// Sets the start index of the FASTQ separator byte. This method is *also*
    /// called for FASTA records before set_record_end() with the *same*
    /// position. The `has_line` argument has the same meaning as there.
    fn set_sep_pos(&mut self, _pos: usize, _has_line: bool);

    /// Returns the byte index of the FASTQ record separator
    fn sep_pos(&self) -> usize;

    /// Set the start of the quality information. Only ever called on FASTQ
    /// records.
    fn set_qual_start(&mut self, _pos: usize);

    /// Adds another quality line start index (only multi-line FASTQ).
    fn add_qual_line_start(&mut self, _pos: usize);

    /// Returns the byte index of the FASTQ quality line start
    fn qual_start(&self) -> usize;

    fn qual_starts(&self) -> &[usize];

    /// Should return, whether there is quality information present (e.g. test
    /// if set_qual_start() was called).
    fn has_qual(&self) -> bool;

    /// Returns the number of quality lines (1 if not multiline FASTQ, 0 if
    /// FASTA)
    fn num_qual_lines(&self) -> usize;

    /// Returns an iterator over quality lines
    fn qual_lines<'s>(
        &'s self,
        _buffer: &'s [u8],
    ) -> <Self::QualLinesKind as LinesIterKind<'s>>::Out;
}

macro_rules! make_seqpos_wrapper {
    ($name:ident) => {
        /// Internal wrapper for use in CoreReader.
        /// The Deref impl exposes the wrapped type.
        #[derive(Clone, Default, Debug)]
        pub(crate) struct $name<S: SeqRecordPosition> {
            inner: S,
        }

        impl<S> SeqRecordPosition for $name<S>
        where
            S: SeqRecordPosition,
        {
            type SeqLinesKind = S::SeqLinesKind;

            #[inline]
            fn init(&mut self) {
                self.inner.init()
            }

            #[inline]
            fn set_record_start(&mut self, pos: usize) {
                self.inner.set_record_start(pos)
            }

            #[inline]
            fn record_start(&self) -> usize {
                self.inner.record_start()
            }

            #[inline]
            fn set_seq_start(&mut self, pos: usize) {
                self.inner.set_seq_start(pos)
            }

            #[inline]
            fn add_seq_line_start(&mut self, pos: usize) {
                self.inner.add_seq_line_start(pos)
            }

            #[inline]
            fn seq_start(&self) -> usize {
                self.inner.seq_start()
            }

            #[inline]
            fn seq_starts(&self) -> &[usize] {
                self.inner.seq_starts()
            }

            #[inline]
            fn seq_end(&self) -> usize {
                self.inner.seq_end()
            }

            #[inline]
            fn set_record_end(&mut self, end: usize, has_line: bool) {
                self.inner.set_record_end(end, has_line)
            }

            #[inline]
            fn record_end(&self) -> usize {
                self.inner.record_end()
            }

            #[inline]
            fn num_lines(&self) -> usize {
                self.inner.num_lines()
            }

            #[inline]
            fn num_seq_lines(&self) -> usize {
                self.inner.num_seq_lines()
            }

            #[inline]
            fn seq_lines<'s>(
                &'s self,
                buffer: &'s [u8],
            ) -> <Self::SeqLinesKind as LinesIterKind<'s>>::Out {
                self.inner.seq_lines(buffer)
            }

            #[inline]
            fn apply_offset(&mut self, offset: isize, search_pos: Option<SearchPos>) {
                self.inner.apply_offset(offset, search_pos)
            }

            #[inline]
            fn line_offset(&self, search_pos: SearchPos, has_line: bool) -> usize {
                self.inner.line_offset(search_pos, has_line)
            }

            #[inline]
            fn maybe_complete(&self) -> bool {
                self.inner.maybe_complete()
            }
        }

        impl<S> Deref for $name<S>
        where
            S: SeqRecordPosition,
        {
            type Target = S;

            fn deref(&self) -> &S {
                &self.inner
            }
        }
    };
}

make_seqpos_wrapper!(SeqRecordPositionWrapper);

// Dummy implementation of QualRecordPosition, these methods will (and should)
// never be called, but we need this in CoreReader::find_record
impl<S> QualRecordPosition for SeqRecordPositionWrapper<S>
where
    S: SeqRecordPosition,
{
    type QualLinesKind = LinePositionIterKind;

    #[inline]
    fn set_sep_pos(&mut self, _pos: usize, _has_line: bool) {
        unreachable!()
    }

    #[inline]
    fn sep_pos(&self) -> usize {
        unreachable!()
    }

    #[inline]
    fn set_qual_start(&mut self, _pos: usize) {
        unreachable!()
    }

    #[inline]
    fn qual_start(&self) -> usize {
        unreachable!()
    }

    #[inline]
    fn add_qual_line_start(&mut self, _pos: usize) {
        unreachable!()
    }

    #[inline]
    fn qual_starts(&self) -> &[usize] {
        unreachable!()
    }

    #[inline]
    fn has_qual(&self) -> bool {
        unreachable!()
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        unreachable!()
    }

    #[inline]
    fn qual_lines<'s>(
        &'s self,
        _buffer: &'s [u8],
    ) -> <Self::QualLinesKind as LinesIterKind<'s>>::Out {
        unreachable!()
    }
}

make_seqpos_wrapper!(QualRecordPositionWrapper);

// Dummy implementation of QualRecordPosition, these methods will (and should)
// never be called, but we need this in CoreReader::find_record
impl<S> QualRecordPosition for QualRecordPositionWrapper<S>
where
    S: QualRecordPosition,
{
    type QualLinesKind = S::QualLinesKind;

    #[inline]
    fn set_sep_pos(&mut self, pos: usize, has_line: bool) {
        self.inner.set_sep_pos(pos, has_line)
    }

    #[inline]
    fn sep_pos(&self) -> usize {
        self.inner.sep_pos()
    }

    #[inline]
    fn set_qual_start(&mut self, pos: usize) {
        self.inner.set_qual_start(pos)
    }

    #[inline]
    fn qual_start(&self) -> usize {
        self.inner.qual_start()
    }

    #[inline]
    fn add_qual_line_start(&mut self, pos: usize) {
        self.inner.add_qual_line_start(pos)
    }

    #[inline]
    fn qual_starts(&self) -> &[usize] {
        self.inner.qual_starts()
    }

    #[inline]
    fn has_qual(&self) -> bool {
        self.inner.has_qual()
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        self.inner.num_qual_lines()
    }

    #[inline]
    fn qual_lines<'s>(
        &'s self,
        buffer: &'s [u8],
    ) -> <Self::QualLinesKind as LinesIterKind<'s>>::Out {
        self.inner.qual_lines(buffer)
    }
}

/// Joins lines together
#[inline]
pub(crate) fn join_lines<'a, L>(mut lines: L, num_lines: usize) -> Cow<'a, [u8]>
where
    L: Iterator<Item = &'a [u8]>,
{
    match num_lines {
        1 => lines.next().unwrap().into(),
        0 => b""[..].into(),
        _ => {
            let mut out = vec![];
            for line in lines {
                out.extend(line);
            }
            return out.into();
        }
    }
}

/// Joins lines together
#[inline]
pub(crate) fn join_lines_given<'a, L, F>(
    mut lines: L,
    num_lines: usize,
    owned_fn: F,
) -> Cow<'a, [u8]>
where
    L: Iterator<Item = &'a [u8]>,
    F: FnOnce() -> &'a mut Vec<u8>,
{
    match num_lines {
        1 => lines.next().unwrap().into(),
        0 => b""[..].into(),
        _ => {
            let output = owned_fn();
            for line in lines {
                output.extend(line);
            }
            return (&*output).into();
        }
    }
}

#[inline]
pub(crate) fn record_head<'a, P: SeqRecordPosition>(pos: &P, buffer: &'a [u8]) -> &'a [u8] {
    trim_cr(&buffer[pos.record_start() + 1..pos.seq_start() - 1])
}

#[inline]
pub(crate) fn record_seq<'a, P: SeqRecordPosition>(pos: &P, buffer: &'a [u8]) -> &'a [u8] {
    let seq = buffer
        .get(pos.seq_start() as usize..pos.seq_end() as usize - 1)
        .unwrap_or(&[]);
    trim_cr(seq)
}

#[inline]
pub(crate) fn record_qual<'a, P: QualRecordPosition>(pos: &P, buffer: &'a [u8]) -> &'a [u8] {
    let seq = buffer
        .get(pos.qual_start() as usize..pos.record_end() as usize - 1)
        .unwrap_or(&[]);
    trim_cr(seq)
}

// Compare sequence and quality score lengths
#[inline]
pub(crate) fn check_lengths<P: QualRecordPosition>(
    pos: &P,
    buffer: &[u8],
    strict: bool,
) -> Result<(), (usize, usize)> {
    if pos.has_qual() {
        let (seq_len, qual_len) = if pos.num_seq_lines() == 1 && pos.num_qual_lines() == 1 {
            if !strict && !cfg!(feature = "strict_length_check") {
                // First, we just compare line lengths, if they are equal, that's
                // enough. Improves performance.
                let seq_len = pos.seq_end() - pos.seq_start();
                let qual_len = pos.record_end() - pos.qual_start();
                if seq_len == qual_len {
                    return Ok(());
                }
            }
            // If not, check without line endings
            let seq_len = record_seq(pos, buffer).len();
            let qual_len = record_qual(pos, buffer).len();
            if seq_len == qual_len {
                return Ok(());
            }
            (seq_len, qual_len)
        } else {
            let seq_len = pos.seq_lines(buffer).fold(0, |l, line| l + line.len());
            let qual_len = pos.qual_lines(buffer).fold(0, |l, line| l + line.len());
            if seq_len == qual_len {
                return Ok(());
            }
            (seq_len, qual_len)
        };
        Err((seq_len, qual_len))
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn join_lines() {
        let mut out = vec![];
        let data = vec![
            (vec![&b"a"[..], &b"bcde"[..]], &b"abcde"[..]),
            (vec![&b""[..], &b"a"[..]], &b"a"[..]),
            (vec![&b""[..]], &b""[..]),
        ];
        for (lines, expected) in data {
            let n = lines.len();
            let joined = join_lines_given(lines.into_iter(), n, || &mut out);
            assert_eq!(&joined, &expected);
            out.clear();
        }
    }
}

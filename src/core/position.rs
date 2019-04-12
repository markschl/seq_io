use std::borrow::Cow;
use std::fmt::Debug;

use crate::core::{join_lines, join_lines_given, trim_cr};
use crate::{LinePositionIter, LineSearchIter};

/// Indicates, which part of the sequence record is being searched. Used in
/// methods of [PositionStore](PositionStore)
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
/// `PositionStore` because their lifetime would pollute the `PositionStore`
/// signature. Therefore, this intermediate type is used to generate the final
/// type. See also
/// ["family trait pattern"](http://lukaskalbertodt.github.io/2018/08/03/solving-the-generalized-streaming-iterator-problem-without-gats.html#workaround-b-hrtbs--the-family-trait-pattern).
pub trait LinesIterKind<'a> {
    type Out: Iterator<Item = &'a [u8]> + DoubleEndedIterator;
}

/// Type constructor for [`LineSearchIter`](LineSearchIter).
pub struct LinesParseKind;

impl<'a> LinesIterKind<'a> for LinesParseKind {
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
pub trait PositionStore: Debug + Clone + Default + Send + Sync {
    type SeqLinesType: for<'a> LinesIterKind<'a>;
    type QualLinesType: for<'a> LinesIterKind<'a>;

    /// Move all stored coordinates back versus the start of the buffer
    /// by `offset` bytes. Reason: searcing got stuck due to reaching the end
    /// of the buffer.
    /// `search_pos` indicates in which part of the record searching got stuck.
    /// Everything that comes after `search_pos` should not be moved because
    /// it's not yet initialized for this record.
    fn move_to_start(&mut self, search_pos: SearchPos, offset: usize);

    /// Sets the start index of the record, always called first, may be
    /// called multiple times while skipping empty lines.
    fn set_record_start(&mut self, start: usize);

    fn record_start(&self) -> usize;

    /// Sets the start index of the sequence.
    fn set_seq_start(&mut self, pos: usize);

    /// Adds another sequence line start index.
    fn add_seq_line_start(&mut self, pos: usize);

    fn seq_start(&self) -> usize;

    fn seq_starts(&self) -> &[usize];

    fn sep_pos(&self) -> usize;

    fn qual_start(&self) -> usize;

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

    fn record_end(&self) -> usize;

    /// Returns the number of lines that the whole record covers. Used for
    /// incrementing the line counter of the parser, it should thus be correct!
    fn num_lines(&self) -> usize;

    /// Returns the number of lines that the sequence has. Should be correct,
    /// used when iterating, and a number > 1 will just take a slice of the
    /// whole sequence.
    fn num_seq_lines(&self) -> usize;

    // Returns the line offset, at which the parser is at the moment, used
    // for getting positional information for an UnexpectedEnd error.
    // Has_line has the same meaning as in set_record_end.
    fn line_offset(&self, search_pos: SearchPos, has_line: bool) -> usize;

    /// Returns an iterator over sequence lines
    fn seq_lines<'s>(&'s self, buffer: &'s [u8]) -> <Self::SeqLinesType as LinesIterKind<'s>>::Out;

    /// Returns an iterator over quality lines
    fn qual_lines<'s>(
        &'s self,
        _buffer: &'s [u8],
    ) -> <Self::QualLinesType as LinesIterKind<'s>>::Out;

    /// Initialize the record. Called before set_record_start(). This is the
    /// place to clear old information if necessary.
    #[inline]
    fn init(&mut self) {}

    /// Combines init() and set_record_start()
    #[inline]
    fn init_record(&mut self, start: usize) {
        self.init();
        self.set_record_start(start);
    }

    /// Sets the start index of the FASTQ separator byte. This method is *also*
    /// called for FASTA records before set_record_end() with the *same*
    /// position. The `has_line` argument has the same meaning as there.
    #[inline]
    fn set_sep_pos(&mut self, _pos: usize, _has_line: bool) {}

    /// Set the start of the quality information. Only ever called on FASTQ
    /// records.
    #[inline]
    fn set_qual_start(&mut self, _pos: usize) {}

    /// Adds another quality line start index (only multi-line FASTQ).
    #[inline]
    fn add_qual_line_start(&mut self, _pos: usize) {}

    #[inline]
    fn qual_starts(&self) -> &[usize] {
        &[]
    }

    /// Should return, whether there is quality information present (e.g. test
    /// if set_qual_start() was called).
    #[inline]
    fn has_qual(&self) -> bool {
        false
    }

    /// Returns whether the position is complete in the sense that
    /// set_record_end() was called. In some cases, this may not be easy to
    /// determine, then true can be returned if record_end() will not panic
    #[inline]
    fn maybe_complete(&self) -> bool {
        self.record_end() > self.record_start()
    }

    /// Returns the number of quality lines (1 if not multiline FASTQ, 0 if
    /// FASTA)
    #[inline]
    fn num_qual_lines(&self) -> usize {
        0
    }

    /// Returns a slice of the head.
    #[inline]
    fn head<'b>(&self, buffer: &'b [u8]) -> &'b [u8] {
        // println!("head {:?} {:?}",buffer, self);
        trim_cr(&buffer[self.record_start() + 1..self.seq_start() - 1])
    }

    /// Returns a slice of the sequence.
    #[inline]
    fn seq<'s>(&'s self, buffer: &'s [u8]) -> &'s [u8] {
        // println!("seq {}", std::string::String::from_utf8_lossy(trim_cr(&buffer[self.seq_start()..self.sep_pos() - 1])));
        // trim_cr(&buffer[self.seq_start()..self.sep_pos() - 1])
        let seq = buffer
            .get(self.seq_start() as usize..self.sep_pos() as usize - 1)
            .unwrap_or(&[]);
        trim_cr(seq)
    }

    /// Returns a slice of the quality information.
    #[inline]
    fn qual<'s>(&'s self, buffer: &'s [u8]) -> &'s [u8] {
        // println!("qual {} {}", self.qual_start(),self.record_end());
        let qual = buffer
            .get(self.qual_start() as usize..self.record_end() as usize - 1)
            .unwrap_or(&[]);
        trim_cr(qual)
    }

    #[inline]
    fn check_lengths(&self, buffer: &[u8], strict: bool) -> Result<(), (usize, usize)> {
        check_lengths(
            self,
            buffer,
            strict,
            self.has_qual() && self.num_qual_lines() > 1,
        )
    }

    #[inline]
    fn join_seq<'a>(&'a self, buffer: &'a [u8]) -> Cow<'a, [u8]> {
        join_lines(self.seq_lines(buffer), self.num_seq_lines())
    }

    #[inline]
    fn join_qual<'a>(&'a self, buffer: &'a [u8]) -> Cow<'a, [u8]> {
        join_lines(self.qual_lines(buffer), self.num_qual_lines())
    }

    #[inline]
    fn join_seq_given<'a, F: FnOnce() -> &'a mut Vec<u8>>(
        &'a self,
        buffer: &'a [u8],
        owned_fn: F,
    ) -> Cow<'a, [u8]> {
        join_lines_given(self.seq_lines(buffer), self.num_seq_lines(), owned_fn)
    }

    #[inline]
    fn join_qual_given<'a, F: FnOnce() -> &'a mut Vec<u8>>(
        &'a self,
        buffer: &'a [u8],
        owned_fn: F,
    ) -> Cow<'a, [u8]> {
        join_lines_given(self.qual_lines(buffer), self.num_qual_lines(), owned_fn)
    }

    #[inline]
    fn from_other<S: PositionStore>(other: S) -> Self {
        // TODO: not well tested
        let mut store = Self::default();
        store.set_record_start(other.record_start());
        if other.maybe_complete() {
            let (&seq_start, more_seq) = other.seq_starts().split_first().unwrap();
            println!("seq starts {} {:?}", seq_start, more_seq);
            store.set_seq_start(seq_start);
            for &pos in more_seq {
                store.add_seq_line_start(pos);
            }
            store.set_sep_pos(other.sep_pos(), true); // TODO: always has_line?
            let (&qual_start, more_qual) = other.seq_starts().split_first().unwrap();
            store.set_qual_start(qual_start);
            for &pos in more_qual {
                store.add_qual_line_start(pos);
            }
            store.set_record_end(other.record_end(), true); // TODO: always has_line?
        }
        store
    }
}

// Compare sequence and quality score lengths
#[inline]
fn check_lengths<P: PositionStore>(
    pos: &P,
    buffer: &[u8],
    strict: bool,
    multiline_fastq: bool,
) -> Result<(), (usize, usize)> {
    if pos.has_qual() {
        let (seq_len, qual_len) = if !multiline_fastq {
            if !strict && !cfg!(feature = "strict_length_check") {
                // First, we just compare line lengths, if they are equal, that's
                // enough. Improves performance.
                let seq_len = pos.sep_pos() - pos.seq_start();
                let qual_len = pos.record_end() - pos.qual_start();
                if seq_len == qual_len {
                    return Ok(());
                }
            }
            // If not, check without line endings
            let seq_len = pos.seq(buffer).len();
            let qual_len = pos.qual(buffer).len();
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

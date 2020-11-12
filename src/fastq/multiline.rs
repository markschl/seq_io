//! FASTQ reading with multi-line FASTQ support.
//!
//! # Behaviour
//!
//! * Since the quality scores can contain `@`, this parser always compares
//!   sequence and quality score lengths, internal `@` are ignored.
//! * Sequence and quality lines are *not* optional. The following is not a
//!   valid record: `@id\n+\n`. The minimal valid (though empty) record would
//!   be `@id\n\n+\n\n`.
//!
//! Writing multi-line FASTQ is not possible. This is on purpose, since
//! multi-line FASTQ is problematic and its use discouraged by many people.
//!
//! # Example
//!
//! ```rust
//! use seq_io::prelude::*;
//! use seq_io::fastq::multiline::Reader;
//!
//! # fn main() {
//! let seq = b"@id
//! SEQU
//! ENCE
//! +
//! II
//! @EI
//! III
//! ";
//!
//! let mut reader = Reader::new(&seq[..]);
//!
//! let rec = reader.next().unwrap().unwrap();
//!
//! assert_eq!(rec.id(), Ok("id"));
//! assert_eq!(rec.seq(), b"SEQU\nENCE");
//! assert_eq!(rec.full_seq().as_ref(), b"SEQUENCE");
//! assert_eq!(rec.full_qual().as_ref(), b"II@EIIII");
//! # }
//! ```
//!
//! **Note** that even if the second line of the quality string starts with a
//! `@`, it is recognized as internal quality score because it is assumed that
//! sequence and quality lengths are the same. If they were different, it may
//! confuse the parser and lead to weird errors.

use crate::core::{LineSearchKind, QualRecordPosition, SearchPos, SeqRecordPosition};
use crate::LineSearchIter;
use serde::{Deserialize, Serialize};
use std::str;

impl_fastq_reader!(true, MultiRangeStore, ("fastq::multiline", "fastq"));

// TODO: still optimized for single-line, otherwise -> FASTX?
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct MultiRangeStore {
    pos: [usize; 5],
    n_seq_lines: usize,
    n_qual_lines: usize,
}

impl SeqRecordPosition for MultiRangeStore {
    type SeqLinesKind = LineSearchKind;

    #[inline]
    fn init(&mut self) {
        self.n_seq_lines = 0;
        self.n_qual_lines = 0;
    }

    #[inline]
    fn set_record_start(&mut self, pos: usize) {
        self.pos[0] = pos;
    }

    #[inline]
    fn record_start(&self) -> usize {
        self.pos[0]
    }

    #[inline]
    fn set_seq_start(&mut self, pos: usize) {
        self.pos[1] = pos;
        self.n_seq_lines += 1;
    }

    #[inline]
    fn add_seq_line_start(&mut self, _: usize) {
        self.n_seq_lines += 1;
    }

    #[inline]
    fn seq_start(&self) -> usize {
        self.pos[1]
    }

    #[inline]
    fn seq_starts(&self) -> &[usize] {
        &self.pos[1..2]
    }

    #[inline]
    fn seq_end(&self) -> usize {
        self.sep_pos()
    }

    #[inline]
    fn set_record_end(&mut self, end: usize, has_line: bool) {
        self.pos[4] = end;
        if !has_line && self.has_qual() {
            // one quality line too much may have been added
            self.n_qual_lines -= 1;
        }
    }

    #[inline]
    fn record_end(&self) -> usize {
        self.pos[4]
    }

    // TODO: does not work with normal FASTA (line numbers not correctly incremented)
    #[inline]
    fn num_lines(&self) -> usize {
        1 + self.n_seq_lines + 1 + self.n_qual_lines
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        self.n_seq_lines
    }

    #[inline]
    fn seq_lines<'s>(&'s self, buffer: &'s [u8]) -> LineSearchIter<'s> {
        LineSearchIter::new(
            &buffer[self.seq_start()..self.seq_end()],
            self.n_seq_lines == 1,
        )
    }

    #[inline]
    fn apply_offset(&mut self, offset: isize, search_pos: Option<SearchPos>) {
        // TODO: correct?
        let range_end = search_pos.unwrap_or(SearchPos::QUAL) as usize;
        for i in 0..=range_end {
            self.pos[i] = (self.pos[i] as isize + offset) as usize;
        }
    }

    #[inline]
    fn line_offset(&self, search_pos: SearchPos, has_line: bool) -> usize {
        self.n_seq_lines + self.n_qual_lines + (search_pos >= SearchPos::SEP) as usize
            - !has_line as usize
    }
}

impl QualRecordPosition for MultiRangeStore {
    type QualLinesKind = LineSearchKind;

    #[inline]
    fn set_sep_pos(&mut self, pos: usize, has_line: bool) {
        self.pos[2] = pos;
        if !has_line {
            self.n_seq_lines -= 1;
        }
    }

    #[inline]
    fn sep_pos(&self) -> usize {
        self.pos[2]
    }

    #[inline]
    fn set_qual_start(&mut self, pos: usize) {
        self.pos[3] = pos;
        if pos == 27 {
            panic!();
        }
        self.n_qual_lines += 1;
    }

    #[inline]
    fn add_qual_line_start(&mut self, _: usize) {
        self.n_qual_lines += 1;
    }

    #[inline]
    fn qual_start(&self) -> usize {
        self.pos[3]
    }

    #[inline]
    fn qual_starts(&self) -> &[usize] {
        &self.pos[3..4]
    }

    #[inline]
    fn has_qual(&self) -> bool {
        self.qual_start() != 0
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        self.n_qual_lines
    }

    #[inline]
    fn qual_lines<'s>(&'s self, buffer: &'s [u8]) -> LineSearchIter<'s> {
        // we need the newline at the end in order to correctly parse the lines
        // -> but records at EOF may not have a newline
        let end = if self.record_end() >= buffer.len() + 1 {
            debug_assert!(self.record_end() == buffer.len() + 1);
            buffer.len()
        } else {
            // println!("record_end {}", self.record_end());
            self.record_end()
        };
        LineSearchIter::new(&buffer[self.qual_start()..end], self.n_qual_lines == 1)
    }
}

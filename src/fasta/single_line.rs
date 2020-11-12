//! FASTA reader, which requires the sequence to be on a single line, it should
//! *not* to be wrapped to multiple lines. This is a rather unusually strict
//! definition of the format and may not be of general use. This parser was
//! mostly created for exploring performance optimizations.

use crate::core::{LinePositionIterKind, SearchPos, SeqRecordPosition};
use crate::LinePositionIter;
use serde::{Deserialize, Serialize};
use std::str;

impl_fasta_reader!(false, RangeStore, ("fasta::single_line", "fasta"));

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RangeStore {
    // All line positions of a record are stored in an array
    // This way, `move_to_start()` is easier to implement
    // than with a struct.
    pos: [usize; 3],
}

impl SeqRecordPosition for RangeStore {
    type SeqLinesKind = LinePositionIterKind;

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
    }

    #[inline]
    fn add_seq_line_start(&mut self, _: usize) {}

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
        self.record_end()
    }

    #[inline]
    fn set_record_end(&mut self, end: usize, _: bool) {
        self.pos[2] = end;
    }

    #[inline]
    fn record_end(&self) -> usize {
        self.pos[2]
    }

    // TODO: does not work with normal FASTA (line numbers not correctly incremented)
    #[inline]
    fn num_lines(&self) -> usize {
        2
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        1
    }

    #[inline]
    fn seq_lines<'s>(&'s self, buffer: &'s [u8]) -> LinePositionIter<'s> {
        // TODO: does not work with normal FASTA (lines not parsed)
        LinePositionIter::new(buffer, &self.pos[1..=2])
    }

    #[inline]
    fn apply_offset(&mut self, offset: isize, search_pos: Option<SearchPos>) {
        // TODO: correct?
        let range_end = search_pos.unwrap_or(SearchPos::SEP) as usize;
        for i in 0..=range_end {
            self.pos[i] = (self.pos[i] as isize + offset) as usize;
        }
    }

    #[inline]
    fn line_offset(&self, search_pos: SearchPos, has_line: bool) -> usize {
        search_pos as usize - !has_line as usize
    }
}

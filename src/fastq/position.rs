use crate::core::{LinePositionIterKind, PositionStore, SearchPos};
use crate::LinePositionIter;
use serde::{Deserialize, Serialize};
use std::str;

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RangeStore {
    // All line positions of a record are stored in an array
    // Allows easy implementation of move_to_start() and slicing in seq_lines()
    pos: [usize; 5],
}

impl PositionStore for RangeStore {
    type SeqLinesType = LinePositionIterKind;
    type QualLinesType = LinePositionIterKind;

    #[inline]
    fn move_to_start(&mut self, search_pos: SearchPos, offset: usize) {
        for i in 0..search_pos as usize + 1 {
            self.pos[i] -= offset;
        }
    }

    #[inline]
    fn record_start(&self) -> usize {
        self.pos[0]
    }

    #[inline]
    fn set_record_start(&mut self, start: usize) {
        self.pos[0] = start;
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
    fn set_sep_pos(&mut self, pos: usize, _: bool) {
        self.pos[2] = pos;
    }

    #[inline]
    fn sep_pos(&self) -> usize {
        self.pos[2]
    }

    #[inline]
    fn set_qual_start(&mut self, pos: usize) {
        self.pos[3] = pos;
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
    fn set_record_end(&mut self, pos: usize, _: bool) {
        self.pos[4] = pos;
    }

    #[inline]
    fn record_end(&self) -> usize {
        self.pos[4]
    }

    // TODO: does not work with multiline-FASTQ (line numbers not correctly incremented)
    #[inline]
    fn num_lines(&self) -> usize {
        4
    }

    #[inline]
    fn line_offset(&self, search_pos: SearchPos, has_line: bool) -> usize {
        search_pos as usize - !has_line as usize
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        1
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        1
    }

    #[inline]
    fn seq_lines<'s>(&'s self, buffer: &'s [u8]) -> LinePositionIter<'s> {
        // TODO: does not work with multi-line FASTQ
        //LineSearchIter::new(self.seq(buffer), false)
        LinePositionIter::new(buffer, &self.pos[1..=2])
    }

    #[inline]
    fn qual_lines<'s>(&'s self, buffer: &'s [u8]) -> LinePositionIter<'s> {
        //LineSearchIter::new(self.qual(buffer), false)
        LinePositionIter::new(buffer, &self.pos[3..=4])
    }
}

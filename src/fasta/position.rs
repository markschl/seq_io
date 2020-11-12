use std::str;

use crate::core::{LinePositionIterKind, SearchPos, SeqRecordPosition};
use crate::LinePositionIter;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct LineStore {
    /// index of '>'
    start: usize,
    /// Indicate line start, but actually it is one byte before (start - 1), which is usually
    /// the line terminator of the header (if there is one). The last index in the Vec is always
    /// the last byte of the last sequence line (including line terminator if present).
    /// Therefore, the length of this Vec should never be 0.
    seq_pos: Vec<usize>,
}

impl SeqRecordPosition for LineStore {
    type SeqLinesKind = LinePositionIterKind;

    #[inline]
    fn init(&mut self) {
        self.seq_pos.clear();
    }

    #[inline]
    fn set_record_start(&mut self, pos: usize) {
        self.start = pos;
    }

    #[inline]
    fn record_start(&self) -> usize {
        self.start
    }

    #[inline]
    fn set_seq_start(&mut self, pos: usize) {
        self.seq_pos.push(pos);
    }

    #[inline]
    fn add_seq_line_start(&mut self, pos: usize) {
        self.seq_pos.push(pos);
    }

    #[inline]
    fn seq_start(&self) -> usize {
        self.seq_pos[0]
    }

    #[inline]
    fn seq_starts(&self) -> &[usize] {
        self.seq_pos.split_last().unwrap().1
    }

    #[inline]
    fn seq_end(&self) -> usize {
        self.record_end()
    }

    #[inline]
    fn set_record_end(&mut self, end: usize, has_line: bool) {
        if has_line {
            // in order to obtain proper slices when iterating,
            // we also need the end position of the sequence
            self.seq_pos.push(end);
        }
    }

    #[inline]
    fn record_end(&self) -> usize {
        self.seq_pos.last().cloned().unwrap()
    }

    #[inline]
    fn num_lines(&self) -> usize {
        self.seq_pos.len()
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        self.seq_pos.len().checked_sub(1).unwrap()
    }

    #[inline]
    fn seq_lines<'s>(&'s self, buffer: &'s [u8]) -> LinePositionIter<'s> {
        debug_assert!(self.seq_pos.len() > 0);
        LinePositionIter::new(buffer, &self.seq_pos)
    }

    #[inline]
    fn apply_offset(&mut self, offset: isize, _: Option<SearchPos>) {
        self.start = (self.start as isize + offset) as usize;
        for s in &mut self.seq_pos {
            *s = (*s as isize + offset) as usize;
        }
    }

    #[inline]
    fn line_offset(&self, _: SearchPos, _: bool) -> usize {
        self.seq_pos.len()
    }

    #[inline]
    fn maybe_complete(&self) -> bool {
        self.seq_pos.len() > 0 && self.record_end() > self.record_start()
    }
}

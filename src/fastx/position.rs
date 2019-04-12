use std::str;

use crate::core::{LinePositionIterKind, PositionStore, SearchPos};
use crate::LinePositionIter;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct LineStore {
    start: usize,
    end: usize,
    seq_pos: Vec<usize>,
    // TODO: performance impact of storing quality positions?
    qual_pos: Vec<usize>,
}

impl PositionStore for LineStore {
    type SeqLinesType = LinePositionIterKind;
    type QualLinesType = LinePositionIterKind;

    #[inline]
    fn init(&mut self) {
        self.seq_pos.clear();
        self.qual_pos.clear();
    }

    #[inline]
    fn move_to_start(&mut self, _: SearchPos, offset: usize) {
        self.start -= offset;
        for s in &mut self.seq_pos {
            *s -= offset;
        }
        for s in &mut self.qual_pos {
            *s -= offset;
        }
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
    fn set_sep_pos(&mut self, pos: usize, has_line: bool) {
        // in order to obtain proper slices when iterating,
        // we also need the end position of the sequence
        if has_line {
            self.seq_pos.push(pos);
        }
    }

    #[inline]
    fn sep_pos(&self) -> usize {
        *self.seq_pos.last().unwrap()
    }

    #[inline]
    fn set_qual_start(&mut self, pos: usize) {
        self.qual_pos.push(pos);
    }

    #[inline]
    fn add_qual_line_start(&mut self, pos: usize) {
        self.qual_pos.push(pos);
    }

    #[inline]
    fn qual_start(&self) -> usize {
        *self.qual_pos.first().unwrap()
    }

    #[inline]
    fn qual_starts(&self) -> &[usize] {
        self.qual_pos.split_last().unwrap().1
    }

    #[inline]
    fn has_qual(&self) -> bool {
        !self.qual_pos.is_empty()
    }

    #[inline]
    fn set_record_end(&mut self, pos: usize, has_line: bool) {
        self.end = pos;
        if self.has_qual() && has_line {
            self.qual_pos.push(pos);
        }
    }

    #[inline]
    fn record_end(&self) -> usize {
        self.end
    }

    #[inline]
    fn num_lines(&self) -> usize {
        self.seq_pos.len() + self.qual_pos.len()
    }

    #[inline]
    fn line_offset(&self, _: SearchPos, has_line: bool) -> usize {
        self.seq_pos.len() + self.qual_pos.len() - !has_line as usize
    }

    #[inline]
    fn num_seq_lines(&self) -> usize {
        debug_assert!(self.seq_pos.len() > 0);
        self.seq_pos.len() - 1
    }

    #[inline]
    fn num_qual_lines(&self) -> usize {
        self.qual_pos.len().checked_sub(1).unwrap_or(0)
    }

    #[inline]
    fn seq_lines<'a>(&'a self, buffer: &'a [u8]) -> LinePositionIter<'a> {
        debug_assert!(self.seq_pos.len() > 0);
        LinePositionIter::new(buffer, &self.seq_pos)
    }

    #[inline]
    fn qual_lines<'a>(&'a self, buffer: &'a [u8]) -> LinePositionIter<'a> {
        debug_assert!(self.qual_pos.len() > 0);
        LinePositionIter::new(buffer, &self.qual_pos)
    }
}

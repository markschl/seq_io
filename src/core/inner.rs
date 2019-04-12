use super::{trim_cr, PositionStore, SearchPos};
use crate::fastx::{Error, ErrorKind, Result};
use crate::policy::{BufPolicy, StdPolicy};
use crate::{ErrorOffset, ErrorPosition};
use std::io;

/// Reader state
#[derive(Debug, Eq, PartialEq, Clone, Copy)]
enum State {
    New,
    Positioned,
    Parsing,
    Finished,
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
struct SearchPosition {
    record_pos: SearchPos,
    byte: usize,
}

impl SearchPosition {
    #[inline]
    fn new(record_pos: SearchPos, byte: usize) -> Self {
        SearchPosition { record_pos, byte }
    }
}

pub(crate) struct CoreReader<R, P, S>
where
    R: io::Read,
    P: BufPolicy,
    S: PositionStore,
{
    buf_reader: crate::core::BufReader<R, P>,
    // Position of current record within current buffer
    pos_store: S,
    // only used for multi-line FASTQ
    length_diff: isize,
    // Current search position within the record (only relevant with read_record_set)
    search_pos: Option<SearchPosition>,
    // Parsing state
    state: State,
    line_idx: u64,
    record_idx: u64,
}

impl<R, S> CoreReader<R, StdPolicy, S>
where
    R: io::Read,
    S: PositionStore,
{
    #[inline]
    pub fn with_capacity(reader: R, capacity: usize) -> Self {
        Self::_from_buf_reader(crate::core::BufReader::with_capacity(reader, capacity))
    }

    #[inline]
    fn _from_buf_reader(rdr: crate::core::BufReader<R>) -> Self {
        CoreReader {
            buf_reader: rdr,
            pos_store: S::default(),
            search_pos: None,
            length_diff: 0,
            state: State::New,
            line_idx: 0,
            record_idx: 0,
        }
        .validate_policy()
    }

    #[inline]
    pub fn from_buf_reader(
        rdr: crate::core::BufReader<R>,
        byte_offset: usize,
        line_idx: u64,
    ) -> Self {
        let mut reader = Self::_from_buf_reader(rdr);
        reader.init_pos(byte_offset, line_idx);
        reader
    }
}

impl<R, P, S> CoreReader<R, P, S>
where
    R: io::Read,
    P: BufPolicy,
    S: PositionStore,
{
    // Minimum capacity that will not cause any panics
    const MIN_CAPACITY: usize = 3;

    #[inline]
    pub fn set_policy<T: BufPolicy>(self, buf_policy: T) -> CoreReader<R, T, S> {
        CoreReader {
            buf_reader: self.buf_reader.set_policy(buf_policy),
            pos_store: self.pos_store,
            search_pos: self.search_pos,
            length_diff: self.length_diff,
            state: self.state,
            line_idx: self.line_idx,
            record_idx: self.record_idx,
        }
        .validate_policy()
    }

    #[inline]
    pub fn set_store<T: PositionStore>(self) -> CoreReader<R, P, T> {
        CoreReader {
            buf_reader: self.buf_reader,
            pos_store: T::from_other(self.pos_store),
            search_pos: self.search_pos,
            length_diff: self.length_diff,
            state: self.state,
            line_idx: self.line_idx,
            record_idx: self.record_idx,
        }
    }

    /// Set the reader position. The byte offset is relative
    /// to the buffer (not the file), the line offset is relative
    /// to the file.
    #[inline]
    pub fn init_pos(&mut self, byte: usize, line_idx: u64) {
        assert!(byte < self.buf_reader.buffer().len());
        self.pos_store.init_record(byte);
        self.line_idx = line_idx;
        self.state = State::Positioned;
    }

    #[inline]
    fn increment_pos(&mut self, reset_length: bool) {
        self.line_idx += self.pos_store.num_lines() as u64;
        self.record_idx += 1;
        // increment buffer position
        self.pos_store.init_record(self.pos_store.record_end());
        if reset_length {
            self.length_diff = 0;
        }
    }

    #[inline]
    fn validate_policy(self) -> Self {
        let cap = self.buf_reader.capacity();
        assert!(
            cap >= Self::MIN_CAPACITY,
            "Reader buffer capacity too small, should be ≥ {}",
            Self::MIN_CAPACITY
        );
        if let Some(buf_limit) = self.buf_reader.policy().limit() {
            assert!(
                cap <= buf_limit,
                "Reader buffer capacity smaller than BufPolicy::limit() ({} < {})",
                cap,
                buf_limit
            );
        }
        self
    }

    #[inline]
    pub fn buf_reader_mut(&mut self) -> &mut crate::core::BufReader<R, P> {
        &mut self.buf_reader
    }

    #[inline]
    pub fn policy(&self) -> &P {
        self.buf_reader.policy()
    }

    #[inline]
    pub fn next(
        &mut self,
        fasta: bool,
        multiline_fasta: bool,
        multiline_fastq: bool,
        guess_qual: bool,
        check_lengths: bool,
    ) -> Option<Result<(&[u8], &S)>> {
        match self._next(
            fasta,
            multiline_fasta,
            multiline_fastq,
            guess_qual,
            check_lengths,
        ) {
            Some(Ok(_)) => Some(Ok((self.buf_reader.buffer(), &self.pos_store))),
            Some(Err(e)) => {
                self.state = State::Finished;
                Some(Err(e))
            }
            None => None,
        }
    }

    #[inline]
    pub fn _next(
        &mut self,
        fasta: bool,
        multiline_fasta: bool,
        multiline_fastq: bool,
        guess_qual: bool,
        check_lengths: bool,
    ) -> Option<Result<()>> {
        //println!( "NEXT ========> f: {}, m: {}/{} check: {} {:?}", fasta, multiline_fasta, multiline_fastq, check_lengths, self.state);

        match self.state {
            State::New => {
                try_opt!(self.fill_buf());
                self.state = State::Parsing;
            }
            State::Positioned => {
                self.state = State::Parsing;
            }
            State::Finished => {
                return None;
            }
            State::Parsing => {
                self.increment_pos(guess_qual || multiline_fastq && !check_lengths);
            }
        };

        let mut search_pos = self.search_pos;
        let mut check_last_byte = false;

        loop {
            // loops until complete record is found (enlarging buffer / relocating contents if necessary)
            search_pos = try_opt!(self.find_record(
                fasta,
                multiline_fasta,
                multiline_fastq,
                guess_qual,
                search_pos,
                check_last_byte
            ));
            //println!("-> found result with {:?} check {} {:?}", search_pos, check_last_byte, self.pos_store);
            if let Some(pos) = search_pos.as_mut() {
                // complete record not found -> first check, whether EOF reached
                if !self.at_end(multiline_fasta, multiline_fastq) {
                    // not at end -> adjust buffer and try again
                    try_opt!(self.adjust_buffer(pos));
                    continue;
                } else if !check_last_byte {
                    //println!(">=>=> check last byte!");
                    check_last_byte = true;
                    // we already know that there will be no next record
                    self.state = State::Finished;
                    continue;
                } else {
                    let has_last =
                        try_opt!(self.check_end(fasta, multiline_fasta, multiline_fastq, *pos));
                    // dbg!(("found end", has_last, self.state));
                    if !has_last {
                        return None;
                    }
                }
            }

            if check_lengths && !fasta {
                try_opt!(self.check_lengths(self.buf_reader.buffer(), multiline_fastq));
            }

            //use crate::BaseRecord; println!("==> found {:?} {:?} buflen {} state {:?}", <$RefRecord>::new(self.buf_reader.buffer(), &self.pos_store).id(), &self.pos_store, self.buf_reader.buffer().len(), self.state);
            return Some(Ok(()));
        }
    }

    #[inline]
    pub fn read_record_set<T>(
        &mut self,
        record_set: &mut T,
        fasta: bool,
        multiline_fasta: bool,
        multiline_fastq: bool,
        guess_qual: bool,
        check_lengths: bool,
    ) -> Result<bool>
    where
        T: crate::core::RecordSet<S>,
    {
        self._read_record_set(
            record_set,
            fasta,
            multiline_fasta,
            multiline_fastq,
            guess_qual,
            check_lengths,
        )
        .map_err(|e| {
            self.state = State::Finished;
            e
        })
    }

    #[inline]
    pub fn _read_record_set<T>(
        &mut self,
        record_set: &mut T,
        fasta: bool,
        multiline_fasta: bool,
        multiline_fastq: bool,
        guess_qual: bool,
        check_lengths: bool,
    ) -> Result<bool>
    where
        T: crate::core::RecordSet<S>,
    {
        record_set.clear();
        //println!( "READ RECSET ========> f: {}, m: {}/{} check: {} {:?}", fasta, multiline_fasta, multiline_fastq, check_lengths, self.state);

        let adjust_buffer = match self.state {
            State::New => {
                // Parsing not yet started, or a seek() to another position
                // in the file was done, requiring the buffer to refilled.
                self.fill_buf()?;
                self.state = State::Positioned;
                false
            }
            State::Finished => {
                // End of in put reached, always return false
                return Ok(false);
            }
            State::Parsing => {
                // next() was previously called, the current record has
                // already been retrurned -> go to next record and
                // make sure it is moved to the start of the buffer if necessary.
                self.increment_pos(guess_qual || multiline_fastq && !check_lengths);
                debug_assert!(self.pos_store.record_start() > 0);
                self.state = State::Positioned;
                true
            }
            State::Positioned => {
                // The previous call to read_record_set() left the current record
                // (the last in the buffer) incomplete.
                // Or, a seek() to a position reachable from within the current
                // buffer was done.
                // In such a case, the current record needs to be moved to the
                // start of the buffer, and the remaining part needs to be
                // refilled.
                true
            }
        };

        // Restore previous search position (if any)
        let mut search_pos = self
            .search_pos
            .unwrap_or_else(|| SearchPosition::new(SearchPos::HEAD, self.pos_store.record_start()));

        // TODO: only make_room here
        //println!("new outer cap {} {:?} {:?} adj {}", self.buf_reader.capacity(), self.pos_store, search_pos, adjust_buffer);
        if adjust_buffer {
            // not at end -> adjust buffer and try again
            self.make_room(&mut search_pos)?;
            self.fill_buf()?;
            //println!("made room -> {}", self.buf_reader.capacity());
        }

        macro_rules! try_break {
            ($expr: expr, $self:ident) => {
                match $expr {
                    Ok(item) => item,
                    Err(e) => {
                        $self.state = State::Finished;
                        break Err(::std::convert::From::from(e));
                    }
                }
            };
        }

        let mut opt_search_pos = Some(search_pos);
        let mut check_last_byte = false;

        // search for records
        let retval = loop {
            opt_search_pos = try_break!(
                self.find_record(
                    fasta,
                    multiline_fasta,
                    multiline_fastq,
                    guess_qual,
                    opt_search_pos,
                    check_last_byte
                ),
                self
            );
            //println!("[RSET] -> found result with {:?} check {} {:?}", opt_search_pos, check_last_byte, self.pos_store);
            match opt_search_pos {
                None => {
                    // complete record found
                    if check_lengths && !fasta {
                        let buf = self.buf_reader.buffer();
                        try_break!(self.check_lengths(buf, multiline_fastq), self);
                    }
                    record_set.set_next_pos(&self.pos_store);
                    //println!( "[RSET] complete: {:?}, lines {} + {}", self.pos_store, self.line_idx, self.pos_store.num_lines());
                    // initiate next record
                    self.line_idx += self.pos_store.num_lines() as u64;
                    self.pos_store.init_record(self.pos_store.record_end());
                    if guess_qual || multiline_fastq && !check_lengths {
                        self.length_diff = 0;
                    }
                    // TODO: (continue)
                }
                Some(sp) => {
                    // record incomplete
                    // Check whether EOF reached
                    if self.at_end(multiline_fasta, multiline_fastq) {
                        if !check_last_byte {
                            check_last_byte = true;
                            // Last byte not yet checked, but we already know
                            // that we are finished afterwards.
                            self.state = State::Finished;
                            continue;
                        } else {
                            let has_last = try_break!(
                                self.check_end(fasta, multiline_fasta, multiline_fastq, sp),
                                self
                            );
                            // EOF reached
                            //println!("check end, has last: {}", has_last);
                            if has_last {
                                // TODO: two calls to check_length, not optimal!
                                // check in loop later? not cache friendly
                                if check_lengths && !fasta {
                                    let buf = self.buf_reader.buffer();
                                    try_break!(self.check_lengths(buf, multiline_fastq), self);
                                }
                                record_set.set_next_pos(&self.pos_store);
                            } else if record_set.is_empty() {
                                return Ok(false);
                            }
                        }
                    } else if record_set.len() == 0 {
                        // try again, enlarging the buffer / checking for end of input
                        try_break!(self.grow(), self);
                        try_break!(self.fill_buf(), self);
                        continue;
                    }
                    // TODO: where to check length??

                    // record set completed -> return
                    self.search_pos = Some(sp);
                    break Ok(true);
                }
            }
        };

        //println!("set buffer {:?}", self.buf_reader.buffer());
        record_set.set_buffer(self.buf_reader.buffer());
        self.record_idx += record_set.len() as u64;
        //println!("finished rset {:?}, lines {}",record_set.len(),self.line_idx);
        retval
    }

    #[inline]
    #[allow(unreachable_code)]
    fn find_record(
        &mut self,
        fasta: bool,
        multiline_fasta: bool,
        multiline_fastq: bool,
        guess_qual: bool,
        search_pos: Option<SearchPosition>,
        check_last_byte: bool,
    ) -> Result<Option<SearchPosition>> {
        //println!("find {:?} {:?} multi: fa {}, fq {} fasta: {}",search_pos, self.pos_store, multiline_fasta, multiline_fastq, fasta);
        let buffer = self.buf_reader.buffer();
        let mut search_buf = buffer;
        if (multiline_fasta || multiline_fastq) && !check_last_byte {
            // // make sure that it is possible to peek one byte ahead by
            // // not searching the last byte
            if buffer.is_empty() {
                return Ok(Some(SearchPosition::new(SearchPos::HEAD, 0)));
            }
            search_buf = &search_buf[..search_buf.len() - 1];
        }

        let mut lines =
            crate::core::LinesMemchr::new(search_buf, self.pos_store.record_start() as usize);

        // // only used for multi-line FASTQ
        // let mut seq_len = 0;
        // let mut qual_len = 0;
        //let search_pos = 'outer: loop {
        // Nested loops are needed for goto-style jumps, named after
        // last item that is searched in this block.
        'sep: loop {
            'seq: loop {
                'head: loop {
                    // Resume searching incomplete record
                    if let Some(pos) = search_pos {
                        lines = crate::core::LinesMemchr::new(search_buf, pos.byte);
                        match pos.record_pos {
                            SearchPos::HEAD => {}
                            SearchPos::SEQ => break 'head,
                            SearchPos::SEP => break 'seq,
                            SearchPos::QUAL => break 'sep,
                        }
                    }

                    // read header
                    loop {
                        if let Some((line, next_pos)) = lines.next() {
                            //println!("> head: {:?} at {} / {:?}",std::string::String::from_utf8_lossy(line),next_pos,buffer.get(next_pos).map(|b| *b as char));
                            if line.len() > 2 || line != b"\n" && line != b"\r\n" {
                                // check first byte of header
                                let first_byte = line[0];
                                let expected = if fasta { b'>' } else { b'@' };
                                if first_byte == expected {
                                    self.pos_store.set_seq_start(next_pos);
                                    if fasta
                                        && multiline_fasta
                                        && buffer.get(next_pos) == Some(&b'>')
                                    {
                                        // Two consequtive FASTA headers -> treat sequence as empty
                                        self.pos_store.set_sep_pos(next_pos, false);
                                        self.pos_store.set_record_end(next_pos, false);
                                        return Ok(None);
                                    }
                                    break 'head;
                                }
                                return Err(Error::new(ErrorKind::InvalidStart {
                                    pos: ErrorPosition::new(Some(self.position()), None, None),
                                    found: first_byte,
                                }));
                            }
                            // skip empty line
                            self.pos_store.set_record_start(next_pos);
                            self.line_idx += 1;
                            continue;
                        }
                        // header line incomplete
                        return Ok(Some(SearchPosition::new(SearchPos::HEAD, lines.pos())));
                    }
                } // end 'head

                // read sequence
                loop {
                    if let Some((line, next_pos)) = lines.next() {
                        //println!("seq: {:?} at {} / {:?}",std::string::String::from_utf8_lossy(line),next_pos,buffer.get(next_pos).map(|b| *b as char));
                        if fasta {
                            if !multiline_fasta || buffer.get(next_pos) == Some(&b'>') {
                                // next line has header
                                self.pos_store.set_sep_pos(next_pos, true);
                                self.pos_store.set_record_end(next_pos, true);
                                return Ok(None);
                            }
                        } else {
                            if multiline_fastq || guess_qual {
                                // !line.is_empty()
                                self.length_diff += trim_cr(&line[..line.len() - 1]).len() as isize;
                            }
                            if !multiline_fastq || buffer.get(next_pos) == Some(&b'+') {
                                self.pos_store.set_sep_pos(next_pos, true);
                                break 'seq;
                            }
                        }
                        self.pos_store.add_seq_line_start(next_pos);
                    } else {
                        // sequence line(s) incomplete
                        return Ok(Some(SearchPosition::new(SearchPos::SEQ, lines.pos())));
                    }
                }
            } // end 'seq

            // Read FASTQ separator
            // Search most likely separator line:
            // optimized for separators without ID and with LF line endings
            // If the guess fails, search for the line ending in the normal way.
            if let Some(next_pos) = lines.guess_next(b"+\n") {
                self.pos_store.set_qual_start(next_pos);
                break 'sep;
            } else if let Some((line, next_pos)) = lines.next() {
                //println!("sep: {:?} at {} / {:?}",std::string::String::from_utf8_lossy(line),next_pos,buffer.get(next_pos).map(|b| *b as char));
                // For multi-line FASTQ, it is already known that the line starts with '+',
                // for single-line we have to test
                if multiline_fastq || line.first() == Some(&b'+') {
                    self.pos_store.set_qual_start(next_pos);
                    break 'sep;
                }
                return Err(Error::new(ErrorKind::InvalidSep {
                    pos: ErrorPosition::new(
                        Some(self.position()),
                        Some(ErrorOffset {
                            line: 1 + self.pos_store.num_seq_lines() as u64,
                            byte: (self.pos_store.sep_pos() - self.pos_store.record_start()) as u64,
                        }),
                        Some(self.record_id()),
                    ),
                    found: trim_cr(&line[..line.len() - 1]).first().cloned(),
                }));
            }
            // separator line incomplete
            return Ok(Some(SearchPosition::new(SearchPos::SEP, lines.pos())));
        } // end 'sep

        // read quality scores
        loop {
            if guess_qual && !multiline_fastq {
                // This block is currently not used, it guesses the line end
                // position of the quality scores assuming that the sequence
                // and quality lines have the same length. This gives some
                // performance gain, but can be risky / lead to confusing errors
                let mut qual_end = self.pos_store.qual_start() + self.length_diff as usize;
                if qual_end + 2 < buffer.len() {
                    //println!("guess qual {} {:?}",qual_end,std::str::from_utf8(&buffer[qual_end..qual_end + 2]));
                    if buffer[qual_end] == b'\r' {
                        qual_end += 1;
                    }
                    if buffer[qual_end] == b'\n' {
                        self.pos_store.set_record_end(qual_end + 1, true);
                        return Ok(None);
                    }
                } else {
                    //println!("guess too short {} {}",lines.pos(),self.pos_store.qual_start());
                    return Ok(Some(SearchPosition::new(SearchPos::QUAL, lines.pos())));
                }
            }
            if let Some((line, next_pos)) = lines.next() {
                //println!("qual: {:?} followed by {:?} at {}",std::string::String::from_utf8_lossy(line), buffer.get(next_pos).map(|b| *b as char), next_pos);
                if !multiline_fastq {
                    self.pos_store.set_record_end(next_pos, true);
                    return Ok(None);
                }
                // multi-line FASTQ
                self.length_diff -= crate::core::trim_cr(&line[..line.len() - 1]).len() as isize;
                //}
                if buffer.get(next_pos) == Some(&b'@') {
                    // Check if quality length already >= sequence length. If not,
                    // continue reading quality.
                    if self.length_diff <= 0 {
                        self.pos_store.set_record_end(next_pos, true);
                        return Ok(None);
                    }
                }
                self.pos_store.add_qual_line_start(next_pos);
            //println!("add qual start {}", next_pos);
            } else {
                return Ok(Some(SearchPosition::new(SearchPos::QUAL, lines.pos())));
            }
        }
    }

    #[inline]
    fn at_end(&self, multiline_fasta: bool, multiline_fastq: bool) -> bool {
        // In the case of multi-line FASTA / FASTQ, the last byte of the buffer was not
        // searched. This offset makes sure that this byte doesn't get lost.
        // TODO: min_records()
        let buf = self.buf_reader.buffer();
        let offset = if multiline_fasta || multiline_fastq {
            1
        } else {
            0
        };
        //println!("at end {} {} buflen: {}, cap: {}", multiline_fasta, multiline_fastq, buf.len(), self.buf_reader.capacity());
        buf.len() + offset < self.buf_reader.capacity()
    }

    // To be called at the end of input. Checks, whether there is a complete
    // sequence record at the end. Returns true if found, otherwise false.
    #[inline]
    #[allow(unreachable_code)]
    fn check_end(
        &mut self,
        fasta: bool,
        multiline_fasta: bool,
        multiline_fastq: bool,
        pos: SearchPosition,
    ) -> Result<bool> {
        let buf = self.buf_reader.buffer();
        // It is important to know, whether there is a last line at the end
        // (without proper terminator, a single \r also counts as line)
        let has_line = pos.byte < buf.len();
        //println!("CHECK END at pos {:?}... len: {}, cap: {},  {:?}",pos,buf.len(),self.buf_reader.capacity(),self.pos_store);
        // Buffer not filled completely -> EOF reached,
        // there will be no next record.
        if pos.record_pos == SearchPos::HEAD {
            // Valid records at EOF always require the header to be terminated
            // by a newline.
            // If the 'header' line is empty, we can finish.
            let line = &buf[self.pos_store.record_start() as usize..];
            if line.is_empty() {
                return Ok(false);
            }
            // Check start byte because UnexpectedEnd errors are somewhat
            // counterintuitive if the start byte is wrong
            let expected = if fasta { b'>' } else { b'@' };
            let first_byte = line[0];
            if first_byte != expected {
                return Err(Error::new(ErrorKind::InvalidStart {
                    pos: ErrorPosition::new(Some(self.position()), None, None),
                    found: first_byte,
                }));
            }
        // Otherwise fall through to UnexpectedEnd
        } else {
            // Multi-line FASTA allows empty sequence lines
            // For FASTQ, SearchPos::QUAL is required, and line should be
            // non-empty, either a newline or some other character, even CR is
            // allowed.
            if fasta && (multiline_fasta || has_line)
                || !fasta && pos.record_pos == SearchPos::QUAL &&
                    // Make sure that quality line is non-empty or terminated.
                    // For multiline-FASTQ: at least one non-empty or terminated line required.
                    (has_line || multiline_fastq && pos.byte > self.pos_store.qual_start())
            {
                // There is a valid last record:
                // - multi-line FASTA/FASTQ: start byte of next record not found,
                //   there may be a newline at the end or not.
                // - single-line FASTA/FASTQ: no newline at end
                let end = if has_line && buf.last() == Some(&b'\n') {
                    buf.len()
                } else {
                    // If there is no newline at the end, make sure that
                    // correct slices are produced when calling seq() or qual()
                    // This end coordinate is out of bounds, but these functions will still
                    // not panic when slicing like this: &buffer[start..end-1]
                    buf.len() + 1
                };
                if fasta {
                    //println!("set fasta sep pos");
                    self.pos_store.set_sep_pos(end, has_line);
                } else if multiline_fastq {
                    // Multi-line FASTQ: consider last line for length check
                    //println!("last line multi eq {:?}", &buf[pos.byte..end - 1]);
                    self.length_diff -=
                        crate::core::trim_cr(&buf[pos.byte..end - 1]).len() as isize;
                }
                self.pos_store.set_record_end(end, has_line);
                //println!("end record {:?} has line: {}, length diff: {}", self.pos_store, has_line, self.length_diff);
                return Ok(true);
            }
        }
        // If the above checks fail, the record is regarded as truncated
        //println!("determine offset {:?} has line: {}", pos.record_pos, has_line);
        let offset = ErrorOffset {
            line: self.pos_store.line_offset(pos.record_pos, has_line) as u64,
            byte: (buf.len() - 1 - self.pos_store.record_start()) as u64,
        };
        return Err(Error::new(ErrorKind::UnexpectedEnd {
            pos: ErrorPosition::new(
                Some(self.position()),
                Some(offset),
                if pos.record_pos > SearchPos::HEAD {
                    Some(self.record_id())
                } else {
                    None
                },
            ),
        }));
    }

    // Compare sequence and quality score lengths
    #[inline]
    #[allow(unreachable_code, unused_variables)]
    fn check_lengths(&self, buffer: &[u8], multiline_fastq: bool) -> Result<()> {
        let (seq_len, qual_len) = if !multiline_fastq {
            match self
                .pos_store
                .check_lengths(self.buf_reader.buffer(), false)
            {
                Ok(()) => return Ok(()),
                Err((s, q)) => (s, q),
            }
        } else {
            let seq_len = self
                .pos_store
                .seq_lines(buffer)
                .fold(0, |l, seq| l + seq.len());
            let qual_len = (seq_len as isize - self.length_diff) as usize;
            //println!("seq diff {} {} .. {}", seq_len, qual_len, self.length_diff);
            if seq_len == qual_len {
                return Ok(());
            }
            (seq_len, qual_len)
        };
        return Err(Error::new(ErrorKind::UnequalLengths {
            pos: ErrorPosition::new(Some(self.position()), None, Some(self.record_id())),
            seq: seq_len,
            qual: qual_len,
        }));
        Ok(())
    }

    // This function is only to be called if the buffer needs to be adjusted,
    // either by making room or by growing
    #[inline]
    fn adjust_buffer(&mut self, pos: &mut SearchPosition) -> Result<()> {
        if self.pos_store.record_start() == 0 {
            // first record already incomplete -> buffer too small
            self.grow()?;
        } else {
            // not the first record -> buffer may be big enough
            self.make_room(pos)?;
        }

        self.fill_buf()?;

        Ok(())
    }

    #[inline]
    fn grow(&mut self) -> Result<()> {
        if !self.buf_reader.grow() {
            return Err(Error::new(ErrorKind::BufferLimit));
        }
        //println!("grow {:?}", self.buf_reader.capacity());
        Ok(())
    }

    #[inline]
    fn make_room(&mut self, pos: &mut SearchPosition) -> Result<()> {
        //println!("make room {:?} {:?}", self.buf_reader.buffer().len(), pos);
        // move incomplete bytes to start of buffer and retry
        let offset = self.pos_store.record_start() as usize;
        self.buf_reader.make_room(offset);
        self.pos_store.move_to_start(pos.record_pos, offset);
        pos.byte -= offset;
        Ok(())
    }

    #[inline(never)]
    fn fill_buf(&mut self) -> Result<usize> {
        self.buf_reader
            .fill_buf()
            .map_err(|e| Error::new(ErrorKind::Io(e)))
    }

    #[inline]
    pub fn position(&self) -> crate::Position {
        crate::Position {
            line: self.line_idx,
            byte: self.buf_reader.file_offset() + self.pos_store.record_start() as u64,
            record: self.record_idx,
        }
    }

    #[inline(never)]
    fn record_id(&self) -> String {
        let buffer = self.buf_reader.buffer();
        let id_bytes = self
            .pos_store
            .head(buffer)
            .split(|b| *b == b' ')
            .next()
            .unwrap();
        String::from_utf8_lossy(id_bytes).into()
    }
}

impl<R, P, S> CoreReader<R, P, S>
where
    R: io::Read + io::Seek,
    P: BufPolicy,
    S: PositionStore,
{
    /// Seeks to the specified position using the byte offset given by
    /// `Position::byte()`. If the position is reachable within the
    /// current buffer, seeking will be very fast. Otherwise,
    /// it has to be
    #[inline]
    pub fn seek(&mut self, pos: &crate::Position, multiline: bool) -> io::Result<()> {
        self.search_pos = None;
        self.line_idx = pos.line();
        self.record_idx = pos.record();
        let buf_offset = self.buf_reader.file_offset();
        let seek_to = pos.byte();
        self.length_diff = 0;
        let corr = if multiline { 1 } else { 0 };
        if seek_to >= buf_offset && seek_to + corr <= self.buf_reader.buffer().len() as u64 {
            //println!("relocate {} {}, -> {}",pos.byte(),buf_offset,pos.byte() - buf_offset);
            // position reachable within buffer -> no actual seeking necessary
            self.pos_store
                .init_record((pos.byte() - buf_offset) as usize);
            self.state = State::Positioned;
            return Ok(());
        }
        //println!("seek {}", pos.byte());
        self.pos_store.init_record(0);
        self.state = State::New;
        self.search_pos = None;
        self.buf_reader.seek_to(seek_to)
    }
}

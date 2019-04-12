use std::io::{self, Read, BufRead, BufReader};
use std::borrow::Cow;
use std::fmt::Debug;

use seq_io::{self, fastx, ErrorPosition, BaseRecord};
use seq_io::fastx::{Error, ErrorKind};


#[derive(Default, Clone, Debug)]
pub struct Record {
    head: Vec<u8>,
    seq: Vec<u8>,
    qual: Option<Vec<u8>>,
}

impl Record {
    pub fn clear(&mut self) {
        self.head.clear();
        self.seq.clear();
        if let Some(q) = self.qual.as_mut() {
            q.clear();
        }
    }
}

impl BaseRecord for Record {
    fn head(&self) -> &[u8] {
        &self.head[1..]
    }

    fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn opt_qual(&self) -> Option<&[u8]> {
        self.qual.as_ref().map(|q| q.as_slice())
    }

    fn has_quality(&self) -> bool {
        self.qual.is_some()
    }

    fn full_seq(&self) -> Cow<[u8]> {
        unimplemented!();
    }

    fn full_seq_given<'s, F: FnOnce() -> &'s mut Vec<u8>>(&'s self, _: F) -> Cow<'s, [u8]> {
        self.seq.as_slice().into()
    }

    fn opt_full_qual_given<'s, F: FnOnce() -> &'s mut Vec<u8>>(
        &'s self,
        _: F,
    ) -> Option<Cow<'s, [u8]>> {
        self.opt_full_qual()
    }

    fn opt_full_qual(&self) -> Option<Cow<[u8]>> {
        self.qual.as_ref().map(|q| q.as_slice().into())
    }

    fn num_seq_lines(&self) -> usize {
        unimplemented!();
    }

    fn num_qual_lines(&self) -> usize {
        unimplemented!();
    }

    fn write<W: io::Write>(&self, _: W) -> io::Result<()> {
        unimplemented!();
    }
}

#[derive(Debug)]
pub struct Reader<R: Read> {
    buf_reader: BufReader<R>,
    line_buf: Vec<u8>,
    record: Record,
    multiline_fasta: bool,
    multiline_fastq: bool,
    // true if at EOF without \n
    no_end: bool,
    // true if reading yielded n = 0 bytes
    finished: bool,
    // reading done, either normally or by error
    done: bool,
    previous_error: Option<fastx::Error>,
}

impl<R: Read> Reader<R> {
    pub fn new(reader: R, multiline_fasta: bool, multiline_fastq: bool) -> Self {
        Reader {
            buf_reader: BufReader::new(reader),
            line_buf: vec![],
            record: Record::default(),
            multiline_fasta,
            multiline_fastq,
            finished: false,
            no_end: false,
            done: false,
            previous_error: None
        }
    }

    #[allow(dead_code)]
    pub fn new_fasta(reader: R, multiline: bool) -> Self {
        Self::new(reader, multiline, false)
    }

    #[allow(dead_code)]
    pub fn new_fastq(reader: R, multiline: bool) -> Self {
        let mut r = Self::new(reader, false, multiline);
        r.record.qual = Some(vec![]);
        r
    }

    #[allow(dead_code)]
    pub fn new_fastx(reader: R, multiline_fasta: bool, multiline_fastq: bool) -> Self {
        let mut r = Self::new(reader, multiline_fasta, multiline_fastq);
        if let Some(start_byte) = r.go_to_start().unwrap() {
            match start_byte {
                b'>' => {},
                b'@' => r.record.qual = Some(vec![]),
                b @ _ => r.previous_error = Some(fastx::Error::new(fastx::ErrorKind::InvalidStart { 
                    pos: r.get_error_pos(),
                    found: b,
                }))
            }
        } else {
            r.done = true;
        }
        r
    }

    fn is_fasta(&self) -> bool {
        self.record.qual.is_none()
    }

    // read all empty lines until the first non-empty line is found
    fn go_to_start(&mut self) -> fastx::Result<Option<u8>> {
        loop {
            // println!("start, {:?}", self.line_buf);
            if let Some(&b) = self.line_buf.first() {
                return Ok(Some(b));
            }
            if self.finished {
                return Ok(None);
            }
            self.read_line()?;
        }
    }

    pub fn next<'a>(&'a mut self) -> Option<fastx::Result<&'a Record>> {
        // handle errors at construction with new_fastx()
        if let Some(e) = self.previous_error.take() {
            self.done = true;
            // println!("-> simple recognition err {:?}", e);
            return Some(Err(e));
        }
        
        if self.done {
            // println!("-> simple done already");
            return None;
        }

        let result = if self.is_fasta() {
            self.read_next_fasta()
        } else {
            self.read_next_fastq()
        };
        match result {
            Ok(record_found) => {
                if record_found {
                    // println!("-> simple rec {:?}", self.record);
                    return Some(Ok(&self.record));
                }
                self.done = true;
                // println!("-> simple none");
                return None;
            }
            Err(e) => {
                self.done = true;
                // println!("-> simple err {:?}", e);
                return Some(Err(e));    
            }
        }
    }

    fn read_header(&mut self, expected_start: u8) -> fastx::Result<bool> {
        // header
        let s = self.go_to_start()?;
        // println!("start {:?}", s);
        if let Some(first_byte) = s {
            // first check the first byte
            if first_byte != expected_start {
                return Err(fastx::Error::new(fastx::ErrorKind::InvalidStart { 
                    pos: self.get_error_pos(),
                    found: first_byte,
                }))
            }
            // then check if truncated
            if self.no_end {
                return Err(fastx::Error::new(fastx::ErrorKind::UnexpectedEnd { 
                    pos: self.get_error_pos(),
                }));
            }
        } else {
            // no header, empty line(s) -> finished
            return Ok(false);
        }
        self.record.head.extend_from_slice(&self.line_buf);
        self.read_line()?;
        Ok(true)
    }

    fn read_next_fasta(&mut self) -> fastx::Result<bool>  {

        self.record.clear();

        if !self.read_header(b'>')? {
            return Ok(false);
        }

        if self.multiline_fasta && (self.line_buf.first() == Some(&b'>') || self.finished) {
            // no sequence line (fine for FASTA)
            return Ok(true);
        }
        if self.finished {
            // single-line FASTA requires non-empty sequence line
            return Err(fastx::Error::new(fastx::ErrorKind::UnexpectedEnd { 
                pos: self.get_error_pos(),
            }));
        }

        // sequence
        loop {
            self.record.seq.extend_from_slice(&self.line_buf);
            // Follow oddity of seq_io to remove CR at EOF from sequence
            if self.no_end && self.record.seq.last() == Some(&b'\r') {
                self.record.seq.pop();
            }
            self.read_line()?;
            if !self.multiline_fasta || self.line_buf.first() == Some(&b'>') || self.finished {
                return Ok(true);
            }
        }
    }

    fn read_next_fastq(&mut self) -> fastx::Result<bool>  {

        self.record.clear();

        if !self.read_header(b'@')? {
            return Ok(false);
        }

        // sequence
        loop {
            if self.finished || self.no_end {
                // truncated FASTQ
                return Err(fastx::Error::new(fastx::ErrorKind::UnexpectedEnd { 
                    pos: self.get_error_pos(),
                }));
            }
            self.record.seq.extend_from_slice(&self.line_buf);
            self.read_line()?; // reads separator or next sequence line
            if !self.multiline_fastq || self.line_buf.first() == Some(&b'+') {
                break;
            }
        }

        // FASTQ separator
        // first check if truncated (separator line needs end with multi-line FASTQ)
        if self.finished || self.no_end {
            return Err(fastx::Error::new(fastx::ErrorKind::UnexpectedEnd { 
                pos: self.get_error_pos(),
            }));
        }
        // then check + byte
        let first_byte = self.line_buf.first();
        if first_byte != Some(&b'+') {
            return Err(fastx::Error::new(fastx::ErrorKind::InvalidSep { 
                pos: self.get_error_pos(),
                found: first_byte.cloned(),
            }));
        }

        self.read_line()?; // reads quality line
        if self.finished {
            // quality is allowed to have no end, but should not be empty
            return Err(fastx::Error::new(fastx::ErrorKind::UnexpectedEnd { 
                pos: self.get_error_pos(),
            }));
        }

        // FASTQ quality
        loop {
            let qual = self.record.qual.as_mut().unwrap();
            qual.extend_from_slice(&self.line_buf);
            //println!("qual {:?} finished {} no end {}", qual, self.finished, self.no_end);
            // Imitate oddity of seq_io to remove CR at EOF from sequence
            if self.no_end && qual.last() == Some(&b'\r') {
                qual.pop();
            }
            let qual_len = qual.len();
            let seq_len = self.record.seq.len();
            self.read_line()?; // next header / empty line / quality line / nothing
            if !self.multiline_fastq || self.line_buf.first() == Some(&b'@') && qual_len >= seq_len || self.finished {
                self.compare_lengths()?;
                return Ok(true)
            }
        }
    }

    fn compare_lengths(&self) -> fastx::Result<()> {
        let seq_len = self.record.seq.len();
        let qual_len = self.record.qual.as_ref().unwrap().len();
        if seq_len == qual_len {
            return Ok(());
        }
        let seq_len = self.record.seq.len();
        let qual_len = self.record.qual.as_ref().unwrap().len();
        if seq_len == qual_len {
            return Ok(());
        }
        return Err(fastx::Error::new(fastx::ErrorKind::UnequalLengths {
            pos: self.get_error_pos(),
            seq: seq_len,
            qual: qual_len,
        }))
    }

    // TODO: document CR behaviour
    fn read_line(&mut self) -> fastx::Result<()> {
        // println!("read line {:?}", self.line_buf);
        self.line_buf.clear();
        let n = self.buf_reader.read_until(b'\n', &mut self.line_buf)?;
        if n == 0 {
            self.finished = true;
            return Ok(());
        }
        if self.line_buf.last() == Some(&b'\n') {
            self.line_buf.pop();
            if self.line_buf.last() == Some(&b'\r') {
                self.line_buf.pop();
            }
        } else {
            self.no_end = true;
        }
        //println!("=> {:?}  (finished {} is fasta: {})", std::string::String::from_utf8_lossy(&self.line_buf), self.finished, self.is_fasta());
        Ok(())
    }

    // TODO: not yet implemented
    fn get_error_pos(&self) -> ErrorPosition {
        ErrorPosition::default()
    }
}


pub fn compare_simple<X, R, P, S>(mut reader: X, mut simple: Reader<R>, match_format: bool)
where
    R: std::io::Read,
    X: seq_io::fastx::dynamic::FastxReader<R, P, S>,
    P: seq_io::policy::BufPolicy,
    S: seq_io::PositionStore,
{
    let mut buf1 = vec![];
    let mut buf2 = vec![];
    let mut match_format = match_format;
    while let Some(result) = reader.next_fastx() {
        // println!("seq_io: {:?}", result);
        let simple_next = simple.next();
        // println!("simple: {:?}", simple_next);
        if match_format && simple_next.is_none() {
            continue;
        }
        let simple_result = simple_next.expect("Simple reader has no next record");
        if match_format && (simple_result.is_err() || result.is_err()) {
            continue;
        }
        match result {
            Ok(rec) => {
                let rec2 = simple_result.unwrap_or_else(|e| panic!("Simple reader returned error: {:?}", e));
                if match_format {
                    if rec.has_quality() == rec2.has_quality() {
                        //println!("format: {} {}", rec.has_quality(), rec2.has_quality());
                        match_format = false;
                    } else {
                        return;
                    }
                }
                compare_records(rec, rec2, &mut buf1, &mut buf2);
            }
            Err(e) => {
                let simple_err = simple_result.err().expect("Simple reader should return error");
                compare_errors(e, simple_err);
            }
        }
    }
    if !match_format {
        assert!(simple.next().is_none());
    }
}


pub fn compare_readers<X, Y, R, P1, P2, S1, S2>(mut reader1: X, mut reader2: Y)
where
    R: std::io::Read,
    X: seq_io::fastx::dynamic::FastxReader<R, P1, S1>,
    Y: seq_io::fastx::dynamic::FastxReader<R, P2, S2>,
    P1: seq_io::policy::BufPolicy,
    P2: seq_io::policy::BufPolicy,
    S1: seq_io::PositionStore,
    S2: seq_io::PositionStore,
{
    let mut buf1 = vec![];
    let mut buf2 = vec![];
    while let Some(result1) = reader1.next_fastx() {
        let result2 = reader2.next_fastx().expect("Reader 2 has no next record");
        //println!("rdr 1: {:?}", result1);
        //println!("rdr 2: {:?}", result2);
        match result1 {
            Ok(rec1) => {
                let rec2 = result2.unwrap_or_else(|e| panic!("Reader 2 returned error: {:?}", e));
                compare_records(rec1, rec2, &mut buf1, &mut buf2);
            }
            Err(err1) => {
                let err2 = result2.err().expect("Reader 2 should return error");
                compare_errors(err1, err2);
            }
        }
    }
    assert!(reader2.next_fastx().is_none(), "Reader 2 should have no next record");
}


pub fn compare_recset<X, Y, R, P, S>(mut reader1: X, mut reader2: Y)
where
    R: std::io::Read,
    X: seq_io::fastx::dynamic::FastxReader<R, P, S>,
    Y: seq_io::fastx::dynamic::FastxReader<R, P, S>,
    P: seq_io::policy::BufPolicy,
    S: seq_io::PositionStore,
{
    let mut buf1 = vec![];
    let mut buf2 = vec![];
    let mut recset = fastx::dynamic::RecordSet::default();
    loop {
        let result1 = reader1.read_record_set_fastx(&mut recset);
        // compare complete records
        //println!("recset {:?}", recset);
        for rec1 in &recset {
            //println!("... -> rec {:?}", rec1);
            let result2 = reader2.next_fastx().expect("Reader 2 has no next record");
            //println!("... -> direcot result {:?}", result2);
            let rec2 = result2.unwrap_or_else(|e| panic!("Reader 2 returned error: {:?}", e));
            compare_records(rec1, rec2, &mut buf1, &mut buf2);
        }
        match result1 {
            Ok(parsing) => {
                if !parsing {
                    break;
                }
            }
            Err(err1) => {
                // next record should be errorneous
                let result2 = reader2.next_fastx().expect("Reader 2 has no next record");
                let err2 = result2.err().expect("Reader 2 should return error");
                compare_errors(err1, err2);
            }
        }
    }
    assert!(reader2.next_fastx().is_none(), "Reader 2 should have no next record");
}



fn compare_records<R1, R2>(rec1: R1, rec2: R2, mut buf1: &mut Vec<u8>, mut buf2: &mut Vec<u8>)
where
    R1: BaseRecord + Debug,
    R2: BaseRecord + Debug,
{
    //println!("rec 1: {:?}", rec1);
    //println!("rec 2: {:?}", rec2);
    assert_eq!(rec1.head(), rec2.head());
    buf1.clear();
    buf2.clear();
    let seq1 = rec1.full_seq_given(|| &mut buf1);
    let seq2 = rec1.full_seq_given(|| &mut buf2);
    assert_eq!(seq1, seq2);
    buf1.clear();
    buf2.clear();
    let q1 = rec1.opt_full_qual_given(|| &mut buf1);
    let q2 = rec1.opt_full_qual_given(|| &mut buf2);
    assert_eq!(q1, q2);
}

fn compare_errors(e1: Error, e2: Error) {
    match (e1.kind(), e2.kind()) {
        (ErrorKind::Io(_), ErrorKind::Io(_)) | 
            (ErrorKind::BufferLimit, ErrorKind::BufferLimit) => {},
        (ErrorKind::InvalidStart { pos: _, found: f1 },
            ErrorKind::InvalidStart { pos: _, found: f2 }) => assert_eq!(f1, f2),
        (ErrorKind::InvalidSep { pos: _, found: f1 },
            ErrorKind::InvalidSep { pos: _, found: f2 }) => assert_eq!(f1, f2),
        (ErrorKind::UnexpectedEnd { pos: _},
            ErrorKind::UnexpectedEnd { pos: _}) => {},
        (ErrorKind::UnequalLengths { pos: _, seq: s1, qual: q1 },
            ErrorKind::UnequalLengths { pos: _, seq: s2, qual: q2 }) => {
                assert_eq!(s1, s2);
                assert_eq!(q1, q2);
            },
        _ => {}
    }
}

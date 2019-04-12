use seq_io::Position;
use std::str::Utf8Error;

pub struct ExpectedRecord {
    pub id: Result<&'static str, Utf8Error>,
    pub desc: Option<Result<&'static str, Utf8Error>>,
    pub head: &'static [u8],
    pub seq: &'static [u8],
    pub qual: Option<&'static [u8]>,
    pub seq_lines: Vec<&'static [u8]>,
    pub qual_lines: Option<Vec<&'static [u8]>>,
    pub pos: Position,
}

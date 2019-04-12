use crate::seq::ExpectedRecord;
use seq_io::Position;

pub const FASTA: &[u8] = b"
>id desc\r
ACCGTAGGCT

CCGTAGGCTG
CGTAGGCTGA\r
\r
GTAGGCTGAA\r
CCCC

>id2
ATTGTTGTTT\r
ATTGTTGTTT

ATTGTTGTTT\r
GGGG
>id3  

>id4
>id5\r
";

lazy_static! {
    pub static ref FASTA_EXPECTED: [ExpectedRecord; 5] = [
        ExpectedRecord {
            id: Ok("id"),
            desc: Some(Ok("desc")),
            head: &b"id desc"[..],
            seq: b"ACCGTAGGCTCCGTAGGCTGCGTAGGCTGAGTAGGCTGAACCCC",
            qual: None,
            seq_lines: vec![
                &b"ACCGTAGGCT"[..],
                &b""[..],
                &b"CCGTAGGCTG"[..],
                &b"CGTAGGCTGA"[..],
                &b""[..],
                &b"GTAGGCTGAA"[..],
                &b"CCCC"[..],
                &b""[..]
            ],
            qual_lines: None,
            pos: Position::new()
                .set_line(1)
                .set_byte(1)
                .set_record(0)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id2"),
            desc: None,
            head: &b"id2"[..],
            seq: b"ATTGTTGTTTATTGTTGTTTATTGTTGTTTGGGG",
            qual: None,
            seq_lines: vec![
                &b"ATTGTTGTTT"[..],
                &b"ATTGTTGTTT"[..],
                &b""[..],
                &b"ATTGTTGTTT"[..],
                &b"GGGG"[..]
            ],
            qual_lines: None,
            pos: Position::new()
                .set_line(10)
                .set_byte(66)
                .set_record(1)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id3"),
            desc: Some(Ok(" ")),
            head: &b"id3  "[..],
            seq: b"",
            qual: None,
            seq_lines: vec![&b""[..]],
            qual_lines: None,
            pos: Position::new()
                .set_line(16)
                .set_byte(112)
                .set_record(2)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id4"),
            desc: None,
            head: &b"id4"[..],
            seq: b"",
            qual: None,
            seq_lines: vec![],
            qual_lines: None,
            pos: Position::new()
                .set_line(18)
                .set_byte(120)
                .set_record(3)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id5"),
            desc: None,
            head: &b"id5"[..],
            seq: b"",
            qual: None,
            seq_lines: vec![],
            qual_lines: None,
            pos: Position::new()
                .set_line(19)
                .set_byte(125)
                .set_record(4)
                .clone()
        },
    ];
}

pub const FASTA_SINGLE: &[u8] = b"

>id desc
ACCGTAGGCTCCGTAGGCTG\r

>id2\r
ATTGTTGTTTATTGTTGTTT
>id3  
A
>id4

";

lazy_static! {
    pub static ref FASTA_SINGLE_EXPECTED: [ExpectedRecord; 4] = [
        ExpectedRecord {
            id: Ok("id"),
            desc: Some(Ok("desc")),
            head: &b"id desc"[..],
            seq: b"ACCGTAGGCTCCGTAGGCTG",
            qual: None,
            seq_lines: vec![&b"ACCGTAGGCTCCGTAGGCTG"[..]],
            qual_lines: None,
            pos: Position::new()
                .set_line(2)
                .set_byte(2)
                .set_record(0)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id2"),
            desc: None,
            head: &b"id2"[..],
            seq: b"ATTGTTGTTTATTGTTGTTT",
            qual: None,
            seq_lines: vec![&b"ATTGTTGTTTATTGTTGTTT"[..]],
            qual_lines: None,
            pos: Position::new()
                .set_line(5)
                .set_byte(34)
                .set_record(1)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id3"),
            desc: Some(Ok(" ")),
            head: &b"id3  "[..],
            seq: b"A",
            qual: None,
            seq_lines: vec![&b"A"[..]],
            qual_lines: None,
            pos: Position::new()
                .set_line(7)
                .set_byte(61)
                .set_record(2)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id4"),
            desc: None,
            head: &b"id4"[..],
            seq: b"",
            qual: None,
            seq_lines: vec![&b""[..]],
            qual_lines: None,
            pos: Position::new()
                .set_line(9)
                .set_byte(70)
                .set_record(3)
                .clone()
        },
    ];
}

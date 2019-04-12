use crate::seq::ExpectedRecord;
use seq_io::Position;

pub const FASTQ: &[u8] = b"
@id desc
ATGC
+ id
IIII\r
@id2
CGAT\r
+
IHII\r


@id3 \r

+ id3

";

lazy_static! {
    pub static ref FASTQ_EXPECTED: [ExpectedRecord; 3] = [
        ExpectedRecord {
            id: Ok("id"),
            desc: Some(Ok("desc")),
            head: &b"id desc"[..],
            seq: b"ATGC",
            qual: Some(b"IIII"),
            seq_lines: vec![&b"ATGC"[..]],
            qual_lines: Some(vec![&b"IIII"[..]]),
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
            seq: b"CGAT",
            qual: Some(b"IHII"),
            seq_lines: vec![&b"CGAT"[..]],
            qual_lines: Some(vec![&b"IHII"[..]]),
            pos: Position::new()
                .set_line(5)
                .set_byte(26)
                .set_record(1)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id3"),
            desc: Some(Ok("")),
            head: &b"id3 "[..],
            seq: b"",
            qual: Some(b""),
            seq_lines: vec![&b""[..]],
            qual_lines: Some(vec![&b""[..]]),
            pos: Position::new()
                .set_line(11)
                .set_byte(47)
                .set_record(2)
                .clone()
        },
    ];
}

pub const FASTQ_MULTILINE: &[u8] = b"
@id desc
AT

GC
+ id
II\r
I
I\r
@id2
CGA\r
T\r
+
IHII\r


@id3 \r


+ id3


@id4
ATG
+

@@

@
@id5

+

";

lazy_static! {
    pub static ref FASTQ_MULTI_EXPECTED: [ExpectedRecord; 5] = [
        ExpectedRecord {
            id: Ok("id"),
            desc: Some(Ok("desc")),
            head: &b"id desc"[..],
            seq: b"ATGC",
            qual: Some(b"IIII"),
            seq_lines: vec![&b"AT"[..], &b""[..], &b"GC"[..]],
            qual_lines: Some(vec![&b"II"[..], &b"I"[..], &b"I"[..]]),
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
            seq: b"CGAT",
            qual: Some(b"IHII"),
            seq_lines: vec![&b"CGA"[..], &b"T"[..]],
            qual_lines: Some(vec![&b"IHII"[..], &b""[..], &b""[..]],),
            pos: Position::new()
                .set_line(9)
                .set_byte(31)
                .set_record(1)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id3"),
            desc: Some(Ok("")),
            head: &b"id3 "[..],
            seq: b"",
            qual: Some(b""),
            seq_lines: vec![&b""[..], &b""[..]],
            qual_lines: Some(vec![&b""[..], &b""[..]]),
            pos: Position::new()
                .set_line(16)
                .set_byte(54)
                .set_record(2)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id4"),
            desc: None,
            head: &b"id4"[..],
            seq: b"ATG",
            qual: Some(b"@@@"),
            seq_lines: vec![&b"ATG"[..]],
            qual_lines: Some(vec![&b""[..], &b"@@"[..], &b""[..], &b"@"[..]]),
            pos: Position::new()
                .set_line(22)
                .set_byte(71)
                .set_record(3)
                .clone()
        },
        ExpectedRecord {
            id: Ok("id5"),
            desc: None,
            head: &b"id5"[..],
            seq: b"",
            qual: Some(b""),
            seq_lines: vec![&b""[..]],
            qual_lines: Some(vec![&b""[..]]),
            pos: Position::new()
                .set_line(29)
                .set_byte(89)
                .set_record(4)
                .clone()
        },
    ];
}

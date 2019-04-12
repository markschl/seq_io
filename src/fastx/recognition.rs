use super::Result;
use crate::core::BufReader;
use crate::policy::BufPolicy;
use crate::{ErrorPosition, Position};
use serde::{Deserialize, Serialize};
use std::io::Read;

#[derive(Debug, Clone, Copy, Eq, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum SeqFormat {
    FASTA,
    FASTQ,
}

/// Returns the sequence format guess based on the first character:
/// `>` for FASTA, `@` for FASTA. Empty input results in `None`, an invalid
/// character results in an error of kind `ErrorKind::InvalidStart`.
/// The tuple at the second position contains the byte offset of the first
/// character in the buffer and the line offset in the file.
/// Returns None if the file is empty.
#[inline]
pub fn recognize_format<R, P>(
    reader: &mut BufReader<R, P>,
) -> Result<Option<(SeqFormat, (usize, u64))>>
where
    R: Read,
    P: BufPolicy,
{
    assert!(reader.capacity() >= 2);
    let mut line_offset = 0;
    let mut byte_offset = 0;
    let ret = 'outer: loop {
        if reader.fill_buf()? == 0 {
            reader.make_room(byte_offset);
            // println!("make room {}", byte_offset);
            byte_offset = 0;
            if reader.fill_buf()? == 0 {
                let buf = reader.buffer();
                byte_offset += buf.len();
                debug_assert!(buf.len() <= 1);
                let last = buf.get(0).cloned();
                if last == Some(b'\n') {
                    return Ok(None);
                } else {
                    break 'outer last;
                }
            }
        }
        let buf = reader.buffer();
        let mut windows = buf.windows(2);
        while let Some(bytes) = windows.next() {
            // println!(".. {:?}, byte: {}, line: {}", bytes, byte_offset, line_offset);
            byte_offset += 1;
            if bytes == b"\r\n" {
                // println!("skip");
                windows.next();
                byte_offset += 1;
                line_offset += 1;
            } else if bytes[0] == b'\n' {
                line_offset += 1;
            } else {
                break 'outer Some(bytes[0]);
            }
        }
    };
    ret.map(|start_byte| {
        let format = match start_byte {
            b'>' => SeqFormat::FASTA,
            b'@' => SeqFormat::FASTQ,
            found @ _ => {
                let mut pos = Position::new();
                pos.set_byte(reader.file_offset() + byte_offset as u64 - 1)
                    .set_line(line_offset);
                let pos = ErrorPosition::new(Some(pos), None, None);
                return Err(crate::fastx::Error::new(
                    crate::fastx::ErrorKind::InvalidStart { pos, found },
                ));
            }
        };
        Ok((format, (byte_offset - 1, line_offset)))
    })
    .transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::BufReader;
    use crate::fastx::ErrorKind;
    use crate::{ErrorPosition, Position};
    // TODO: take care about policy
    #[test]
    fn recognize_empty() {
        let mut rdr = BufReader::with_capacity(&b"\n\r\n\n\r\n"[..], 3);
        let fmt = recognize_format(&mut rdr).unwrap();
        assert_eq!(fmt, None);
    }

    #[test]
    fn recognize_fasta1() {
        let mut rdr = BufReader::with_capacity(&b"\n\r\n\n\r\n>"[..], 3);
        let fmt = recognize_format(&mut rdr).unwrap();
        assert_eq!(fmt, Some((SeqFormat::FASTA, (0, 4))));
    }

    #[test]
    fn recognize_fasta2() {
        let mut rdr = BufReader::with_capacity(&b"\n\n\r\n>\n\n"[..], 3);
        let fmt = recognize_format(&mut rdr).unwrap();
        assert_eq!(fmt, Some((SeqFormat::FASTA, (0, 3))));
    }

    #[test]
    fn recognize_fasta3() {
        let mut rdr = BufReader::with_capacity(&b"\n\r\n\n>id"[..], 3);
        let fmt = recognize_format(&mut rdr).unwrap();
        assert_eq!(fmt, Some((SeqFormat::FASTA, (1, 3))));
    }

    #[test]
    fn recognize_fastq1() {
        let mut rdr = BufReader::with_capacity(&b"\n\r\n\n@\n"[..], 3);
        let fmt = recognize_format(&mut rdr).unwrap();
        assert_eq!(fmt, Some((SeqFormat::FASTQ, (1, 3))));
    }

    #[test]
    fn recognize_fastq2() {
        let mut rdr = BufReader::with_capacity(&b"\n\r\n\n@"[..], 3);
        let fmt = recognize_format(&mut rdr).unwrap();
        assert_eq!(fmt, Some((SeqFormat::FASTQ, (0, 3))));
    }

    #[test]
    fn recognize_invalid1() {
        let mut rdr = BufReader::with_capacity(&b"\n\r>id"[..], 3);
        let res = recognize_format(&mut rdr);
        assert!(res.is_err());
        let err = res.err().unwrap();
        match err.kind() {
            ErrorKind::InvalidStart { pos, found: b'\r' } => {
                let exp_pos = ErrorPosition::new(
                    Some(Position::new().set_byte(1).set_line(1).clone()),
                    None,
                    None,
                );
                assert_eq!(pos, &exp_pos);
            }
            e @ _ => panic!("Wrong error: {:?}", e),
        }
    }
}

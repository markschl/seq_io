use crate::HeadWriter;
use std::borrow::Borrow;
use std::io;

/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTA format.
///
/// The header line may be supplied as either a `&[u8]` / `&str` header line,
/// or as separate ID and description parts in the form `(&[u8], Option<&[u8]>)`
/// or `(&str, Option<&str>)`.
#[inline]
pub fn write<W, H>(writer: W, head: H, seq: &[u8]) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
{
    write_iter(writer, head, Some(seq))
}

/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTA format. In contrast to [`write`](write()), this
/// function accepts a sequence iterator.
#[inline]
pub fn write_iter<W, H, S, L>(mut writer: W, head: H, seq: S) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
    S: IntoIterator<Item = L>,
    L: Borrow<[u8]>,
{
    write_head(&mut writer, head)?;
    write_seq_iter(writer, seq)
}

/// Writes data to the FASTA format. Wraps the sequence to produce multi-line FASTA
/// with a maximum width specified by the `wrap` parameter.
#[inline]
pub fn write_wrap<W, H>(mut writer: W, head: H, seq: &[u8], wrap: usize) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
{
    write_head(&mut writer, head)?;
    write_wrap_seq(writer, seq, wrap)
}

/// Writes data to the FASTA format. Wraps the sequence to produce multi-line FASTA
/// with a maximum width specified by the `wrap` parameter. Accepts a sequence
/// iterator.
#[inline]
pub fn write_wrap_iter<'a, W, H, S>(mut writer: W, head: H, seq: S, wrap: usize) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
    S: IntoIterator<Item = &'a [u8]>,
{
    write_head(&mut writer, head)?;
    write_wrap_seq_iter(writer, seq, wrap)
}

/// Writes only the sequence line.
#[inline]
pub fn write_head<W, H>(writer: W, head: H) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
{
    head.write_head(writer, b'>')
}

/// Writes only the sequence line.
#[inline]
pub fn write_seq<W>(writer: W, seq: &[u8]) -> io::Result<()>
where
    W: io::Write,
{
    write_seq_iter(writer, Some(seq))
}

/// Writes the sequence line, and wraps the output to a maximum width specified by `wrap`.
#[inline]
pub fn write_wrap_seq<W>(mut writer: W, seq: &[u8], wrap: usize) -> io::Result<()>
where
    W: io::Write,
{
    assert!(wrap > 0);
    for chunk in seq.chunks(wrap) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// Writes the sequence line from an iterator of lines.
#[inline]
pub fn write_seq_iter<'a, W, S, L>(mut writer: W, seq: S) -> io::Result<()>
where
    W: io::Write,
    S: IntoIterator<Item = L>,
    L: Borrow<[u8]>,
{
    for line in seq {
        writer.write_all(line.borrow())?;
    }
    writer.write_all(b"\n")
}

/// Writes the sequence line from an iterator (such as `SeqLines`) and wraps the output
/// to a maximum width specified by `wrap`.
#[inline]
pub fn write_wrap_seq_iter<'a, W, S, L>(mut writer: W, seq: S, wrap: usize) -> io::Result<()>
where
    W: io::Write,
    S: IntoIterator<Item = L>,
    L: Borrow<[u8]>,
{
    assert!(wrap > 0);
    let mut n_line = 0;
    for line in seq {
        let mut chunk = line.borrow();
        loop {
            let remaining = wrap - n_line;
            if chunk.len() <= remaining {
                writer.write_all(chunk)?;
                n_line += chunk.len();
                break;
            }
            // chunk longer than line -> break
            let (line, rest) = chunk.split_at(remaining);
            chunk = rest;
            //  println!("write {:?}", line);
            writer.write_all(line)?;
            writer.write_all(b"\n")?;
            n_line = 0;
        }
    }
    writer.write_all(b"\n")
}

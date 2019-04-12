use crate::HeadWriter;
use std::borrow::Borrow;
use std::io;

/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTQ format.
///
/// The header line may be supplied as either a `&[u8]` / `&str` header line,
/// or as separate ID and description parts in the form `(&[u8], Option<&[u8]>)`
/// or `(&str, Option<&str>)`.
#[inline]
pub fn write<W, H>(writer: W, head: H, seq: &[u8], qual: &[u8]) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
{
    write_iter(writer, head, Some(seq), Some(qual))
}

/// Helper function for writing data (not necessarily stored in a `Record` instance)
/// to the FASTQ format. In contrast to [`write`](write()), this
/// function allows specifying sequence and quality iterators.
#[inline]
pub fn write_iter<W, H, S, Ls, Q, Lq>(mut writer: W, head: H, seq: S, qual: Q) -> io::Result<()>
where
    W: io::Write,
    H: HeadWriter,
    S: IntoIterator<Item = Ls>,
    Ls: Borrow<[u8]>,
    Q: IntoIterator<Item = Lq>,
    Lq: Borrow<[u8]>,
{
    head.write_head(&mut writer, b'@')?;
    for line in seq {
        writer.write_all(line.borrow())?;
    }
    writer.write_all(b"\n+\n")?;
    for line in qual {
        writer.write_all(line.borrow())?;
    }
    writer.write_all(b"\n")
}

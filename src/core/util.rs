use std::borrow::Cow;
use std::mem::replace;

macro_rules! try_opt {
    ($expr: expr) => {
        match $expr {
            Ok(item) => item,
            Err(e) => return Some(Err(::std::convert::From::from(e))),
        }
    };
}

/// Remove line ends ('\r\n' or '\n' or '\r') from a byte slice
#[inline]
pub(crate) fn trim_end(line: &[u8]) -> &[u8] {
    if let Some((&byte, remaining)) = line.split_last() {
        match byte {
            b'\n' => return trim_cr(remaining),
            b'\r' => return remaining,
            _ => {}
        }
    }
    line
}

/// Remove a final '\r' from a byte slice
#[inline]
pub(crate) fn trim_cr(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', remaining)) = line.split_last() {
        remaining
    } else {
        line
    }
}

/// Joins lines together
#[inline]
pub(crate) fn join_lines<'a, L>(mut lines: L, num_lines: usize) -> Cow<'a, [u8]>
where
    L: Iterator<Item = &'a [u8]>,
{
    match num_lines {
        1 => lines.next().unwrap().into(),
        0 => b""[..].into(),
        _ => {
            let mut out = vec![];
            for line in lines {
                out.extend(line);
            }
            // println!("joined {} lines: {:?}", num_lines, out);
            return out.into();
        }
    }
}

/// Joins lines together
#[inline]
pub(crate) fn join_lines_given<'a, L, F>(
    mut lines: L,
    num_lines: usize,
    owned_fn: F,
) -> Cow<'a, [u8]>
where
    L: Iterator<Item = &'a [u8]>,
    F: FnOnce() -> &'a mut Vec<u8>,
{
    match num_lines {
        1 => lines.next().unwrap().into(),
        0 => b""[..].into(),
        _ => {
            let output = owned_fn();
            for line in lines {
                output.extend(line);
            }
            return (&*output).into();
        }
    }
}

/// Iterator over lines in a byte string tailored for the
/// needs of this crate.
/// It can start searching in a text at any `byte_start`.
///
/// The iterator yields:
/// - the line itself
/// - the byte offset of the next line start (NOT the current line end)
///
/// It yields ONLY lines terminated by '\n', if the last byte in the text
/// is not '\n', the last line will not be returned
#[derive(Debug)]
pub(crate) struct LinesMemchr<'a> {
    text: &'a [u8],
    position: usize,
}

impl<'a> LinesMemchr<'a> {
    #[inline]
    pub fn new(text: &[u8], byte_start: usize) -> LinesMemchr {
        LinesMemchr {
            text: &text[byte_start..],
            position: byte_start,
        }
    }

    #[inline]
    pub fn pos(&self) -> usize {
        self.position
    }

    #[inline]
    pub fn guess_next(&mut self, guess: &'a [u8]) -> Option<usize> {
        if guess.len() <= self.text.len() {
            let (cmp, rest) = self.text.split_at(guess.len());
            if cmp == guess {
                self.position += guess.len();
                self.text = rest;
                return Some(self.position);
            }
        }
        None
    }
}

impl<'a> Iterator for LinesMemchr<'a> {
    type Item = (&'a [u8], usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        memchr::memchr(b'\n', self.text).map(|pos| {
            let pos = pos + 1;
            self.position += pos;
            let (line, text) = self.text.split_at(pos);
            self.text = text;
            (line, self.position)
        })
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.text.len()))
    }
}

/// Iterator over text lines separated by '\n' or '\r\n';
/// also removes the line ends
/// Returns empty lines or incomplete lines (without separator at end)
/// Note: Orphan \r at the end are trimmed as well!
#[derive(Debug)]
pub(crate) struct SimpleLines<'a> {
    text: &'a [u8],
    is_new: bool,
}

impl<'a> SimpleLines<'a> {
    #[inline]
    pub fn new(text: &[u8]) -> SimpleLines {
        // println!("new iter {:?}", text);
        SimpleLines { text, is_new: true }
    }
}

impl<'a> Iterator for SimpleLines<'a> {
    type Item = &'a [u8];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(pos) = memchr::memchr(b'\n', self.text) {
            let pos = pos + 1;
            let (line, text) = self.text.split_at(pos);
            self.text = text;
            // println!(">line {:?}", line);
            // println!(".line {:?}", trim_cr(&line[..line.len() - 1]));
            return Some(trim_cr(&line[..line.len() - 1]));
        }
        if !self.text.is_empty() {
            // println!(">rest {:?}", self.text);
            // println!(".rest {:?}", trim_cr(&self.text));
            let rest = trim_cr(replace(&mut self.text, &[][..]));
            return Some(rest);
        }
        None
    }

    #[inline]
    fn size_hint(&self) -> (usize, Option<usize>) {
        (0, Some(self.text.len()))
    }
}

impl<'a> DoubleEndedIterator for SimpleLines<'a> {
    #[inline]
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.is_new {
            // remove trailine newline if present
            if let Some((last, before)) = self.text.split_last() {
                if last == &b'\n' {
                    self.text = before;
                    self.is_new = false;
                }
            }
        }
        if let Some(pos) = memchr::memrchr(b'\n', self.text) {
            let (text, line) = self.text.split_at(pos);
            self.text = text;
            return Some(trim_cr(&line[1..]));
        }
        if !self.text.is_empty() {
            return Some(replace(&mut self.text, &[][..]));
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_lines() {
        let text = b"a\nbbb\nc\n\nd\n";
        let lines: Vec<_> = SimpleLines::new(text).collect();
        let expected = vec![&b"a"[..], b"bbb", b"c", b"", b"d"];
        assert_eq!(lines, expected);
        let lines: Vec<_> = SimpleLines::new(text).rev().collect();
        let exp_rev: Vec<&[u8]> = expected.into_iter().rev().collect();
        assert_eq!(lines, exp_rev);
    }

    #[test]
    fn join_lines() {
        let mut out = vec![];
        let data = vec![
            (vec![&b"a"[..], &b"bcde"[..]], &b"abcde"[..]),
            (vec![&b""[..], &b"a"[..]], &b"a"[..]),
            (vec![&b""[..]], &b""[..]),
        ];
        for (lines, expected) in data {
            let n = lines.len();
            let joined = join_lines_given(lines.into_iter(), n, || &mut out);
            assert_eq!(&joined, &expected);
            out.clear();
        }
    }
}

macro_rules! impl_fastx_reader {
    ($multiline_fastq:expr, $DefaultPositionStore:ty, ($($mod_path:expr),*)) => {

impl_reader!(
    Reader, crate::fastx::RefRecord<S>, crate::fastx::OwnedRecord,
    crate::fastx::RecordSet<S>, crate::fastx::Error,
    crate::core::QualRecordPosition, $DefaultPositionStore, true, $multiline_fastq,
    "fastx", ($($mod_path),*), "\n@id\nACGT\n+\nIIII\n", "\n@id1\nACGT\n+\nIIII\n@id2\nTGCA\n+\nIIII\n",
    ["OwnedRecord {head: b\"id1\".to_vec(), seq: b\"ACGT\".to_vec(), qual: Some(b\"IIII\".to_vec())}",
     "OwnedRecord {head: b\"id2\".to_vec(), seq: b\"TGCA\".to_vec(), qual: Some(b\"IIII\".to_vec())}"]
);

/// FASTX parser
pub struct Reader<R, P = crate::policy::StdPolicy, S = $DefaultPositionStore>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: crate::core::QualRecordPosition,
{
    inner: crate::core::CoreReader<R, P, crate::core::QualRecordPositionWrapper<S>, S>,
    fasta: Option<Option<bool>>,
}


impl<R, P, S> Reader<R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: crate::core::QualRecordPosition,
{
    #[inline]
    fn _new(reader: R, capacity: usize, policy: P) -> Self {
        Reader {
            inner: crate::core::CoreReader::new(reader, capacity, policy),
            fasta: None,
        }
    }

    #[inline]
    fn _from_buf_reader(rdr: crate::core::BufReader<R, P>, byte_offset: usize, line_idx: u64) -> Self {
        Reader {
            inner: crate::core::CoreReader::from_buf_reader(rdr, byte_offset, line_idx),
            fasta: None,
        }
    }

    #[inline]
    fn _check_is_fasta(&mut self) -> super::Result<Option<bool>> {
        if let Some(fasta) = self.fasta {
            Ok(fasta)
        } else {
            let fmt = match crate::fastx::recognize_format(self.inner.buf_reader_mut()) {
                Ok(f) => f,
                Err(e) => {
                    self.fasta = Some(None);
                    return Err(e);
                }
            };
            let is_fasta = fmt.map(|(fmt, (byte, line))| {
                self.inner.init_pos(byte, line);
                match fmt {
                    crate::fastx::SeqFormat::FASTA => true,
                    crate::fastx::SeqFormat:: FASTQ => false,
                }
            });
            self.fasta = Some(is_fasta);
            Ok(is_fasta)
        }
    }

    #[inline]
    fn _format(&self) -> Option<crate::fastx::SeqFormat> {
        self.fasta
            .and_then(|f| f)
            .map(|f| if f {
                crate::fastx::SeqFormat::FASTA
            } else {
                crate::fastx::SeqFormat::FASTQ
            })
    }
}

}
}

impl_fastx_reader!(false, super::LineStore, ("fastx"));

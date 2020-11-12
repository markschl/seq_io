macro_rules! impl_fasta_reader {
    ($multiline:expr, $DefaultPositionStore:ty, ($($mod_path:expr),*)) => {

impl_reader!(
    Reader, crate::fasta::RefRecord<S>, crate::fasta::OwnedRecord,
    crate::fasta::RecordSet<S>, crate::fasta::Error,
    crate::core::SeqRecordPosition, $DefaultPositionStore, $multiline, false,
    "fasta", ($($mod_path),*), "\n>id\nSEQUENCE\n", ">id1\nACGT\n>id2\nTGCA\n",
    ["OwnedRecord {head: b\"id1\".to_vec(), seq: b\"ACGT\".to_vec()}",
     "OwnedRecord {head: b\"id2\".to_vec(), seq: b\"TGCA\".to_vec()}"]
);

/// FASTA parser
pub struct Reader<R, P = crate::policy::StdPolicy, S = $DefaultPositionStore>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: crate::core::SeqRecordPosition,
{
    inner: crate::core::CoreReader<R, P, crate::core::SeqRecordPositionWrapper<S>, S>,
}


impl<R, P, S> Reader<R, P, S>
where
    R: std::io::Read,
    P: crate::policy::BufPolicy,
    S: crate::core::SeqRecordPosition,
{
    #[inline]
    fn _new(reader: R, capacity: usize, policy: P) -> Self {
        Reader {
            inner: crate::core::CoreReader::new(reader, capacity, policy),
        }
    }

    #[inline]
    fn _from_buf_reader(rdr: crate::core::BufReader<R, P>, byte_offset: usize, line_idx: u64) -> Self {
        Reader {
            inner: crate::core::CoreReader::from_buf_reader(rdr, byte_offset, line_idx),
        }
    }

    #[inline]
    fn _check_is_fasta(&mut self) -> super::Result<Option<bool>> {
        Ok(Some(true))
    }

    #[inline]
    fn _format(&self) -> Option<crate::fastx::SeqFormat> {
        Some(crate::fastx::SeqFormat::FASTA)
    }
}

}
}

// apply the macro
impl_fasta_reader!(true, super::LineStore, ("fasta"));

//! FASTX reading with multi-line FASTQ support.
//!
//! [`fastx::multiline_qual::Reader`](crate::fastx::multiline_qual::Reader)
//! has the same behaviour as [`Reader`](crate::fastx::Reader)
//! in [`fastx`](crate::fastx),
//! but it also supports multi-line FASTQ (see [`fastq::multiline`](crate::fastq::multiline)).

impl_fastx_reader!(true, super::LineStore, ("fastx::multiline_qual", "fastx"));

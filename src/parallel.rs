//! Functions for parallel processing of record sets and records.
//!
//! Sequences are read and processed in batches (`RecordSet`) because sending
//! data across channels has a performance impact. The process works as follows:
//!
//! * Sequence parsing is done in a background thread
//! * Record sets are sent to worker threads, where expensive operations take
//!   place (e.g. sequence analysis).
//! * The results are sent to the main thread along with the record sets.
//! * The record sets are recycled by sending them back to the background
//!   reader.
//!
//! # Per-record processsing
//!
//! The easiest to use are the functions, which operate directly on sequence
//! records without having to deal with record sets:
//!
//! * [`read_process_fasta_records`](read_process_fasta_records)
//! * [`read_process_fastq_records`](read_process_fastq_records)
//! * [`read_process_fastx_records`](read_process_fastx_records)
//!
//! They are specific for the given sequence format, but it is possible to
//! generate functions for other types using the
//! [`parallel_record_impl`](parallel_record_impl) macro.
//!
//! ## Example
//!
//! This example filters sequences by the occurrence of a pattern:
//!
//! ```no_run
//! use seq_io::prelude::*;
//! use seq_io::fastq::{Reader,Record};
//! use seq_io::parallel::read_process_fastq_records;
//! use std::fs::File;
//! use std::io::BufWriter;
//!
//! let reader = Reader::from_path("seqs.fastq").unwrap();
//! let mut writer = BufWriter::new(File::create("filtered.fastq").unwrap());
//!
//! read_process_fastq_records(reader, 4, 2,
//!     |record, found| { // runs in worker
//!         *found = record.seq().windows(3).position(|s| s == b"AAA").is_some();
//!     },
//!     |record, found| { // runs in main thread
//!         if *found {
//!             record.write(&mut writer).unwrap();
//!         }
//!         // Some(value) will stop the reader, and the value will be returned.
//!         // In the case of never stopping, we need to give the compiler a hint about the
//!         // type parameter, thus the special 'turbofish' notation is needed.
//!         None::<()>
//! }).unwrap();
//! ```
//!
//! # Record set processing
//!
//! It is still possible to directly work with record sets using the following
//! generic functions:
//!
//! * [`read_process_recordsets`](read_process_recordsets)
//! * [`read_process_recordsets_init`](read_process_recordsets_init)
//!
//! ## Example
//!
//! This example searches for the first occurrence of a sequence pattern and
//! then stops the parser.
//!
//! ```no_run
//! use seq_io::prelude::*;
//! use seq_io::fastq;
//! use seq_io::parallel::read_process_recordsets;
//!
//! let reader = fastq::Reader::from_path("seqs.fastq").unwrap();
//!
//! read_process_recordsets(reader, 4, 2,
//!     |record_set, position| {
//!         // This function does the heavy work.
//!         // The code is not necessarily very efficient, just for demonstration.
//!         for (i, record) in record_set.into_iter().enumerate() {
//!             if let Some(pos) = record.seq().windows(3).position(|s| s == b"AAA") {
//!             *position = Some((i, pos));
//!             }
//!         }
//!         *position = None;
//!     }, |mut record_sets| {
//!         // This function runs in the main thread. It provides a streaming iterator over
//!         // record sets and the corresponding return values from the worker function
//!         // (not necessarily in the same order as in the file)
//!         while let Some(result) = record_sets.next() {
//!             let (record_set, position) = result?;
//!             if let Some(&(i, pos)) = position.as_ref() {
//!                 let record = record_set.into_iter().nth(i).unwrap();
//!                 println!("Found AAA in record {} at position {}", record.id().unwrap(), pos);
//!                 return Ok(());
//!             }
//!         }
//!         // Here, we need to give the compiler a type hint about the returned
//!         // result, since it is not smart enough to infer it.
//!         // In real-world programs, this may be less of an issue because the
//!         // returned result type is often known.
//!         Ok::<_, fastq::Error>(())
//!     }
//! ).expect("FASTQ reading error");
//! ```

use crate::core::{QualRecordPosition, SeqRecordPosition};
use crate::{fasta, fastq, fastx};
use std::sync::mpsc;

/// A simple trait required to be implemented for readers fed into the
/// functions in this module.
pub trait RecordSetReader {
    type RecordSet: Send;
    type Err: Send;
    fn fill_data(&mut self, record: &mut Self::RecordSet) -> Result<bool, Self::Err>;
}

/// This function reads record sets and processes them in parallel threads.
///
/// * It takes a [`RecordSetReader`](RecordSetReader), which reads data into
///   record sets in a background thread.
/// * These are then sent to `n_workers` worker threads, where the heavy work
///   is done in the `work` closure.
/// * Once ready, the record sets and work results are sent to the main thread
///   and provided to the `func` closure. The won't necessarily arrive in the
///   same order as they were read.
pub fn read_process_recordsets<R, W, F, O, Out>(
    reader: R,
    n_workers: u32,
    queue_len: usize,
    work: W,
    func: F,
) -> Out
where
    R: RecordSetReader + Send,
    R::RecordSet: Default + Send,
    O: Default + Send,
    W: Send + Sync + Fn(&mut R::RecordSet, &mut O),
    F: FnOnce(ParallelDataSets<R::RecordSet, R::Err, O>) -> Out,
{
    read_process_recordsets_init(|| Ok::<_, ()>(reader), n_workers, queue_len, work, func).unwrap()
}

/// Like [`read_process_recordsets`](read_process_recordsets), but additionally
/// allows initiating the reader in the background thread using a closure
/// (`reader_init`).
/// This is useful for readers, which don't implement `Send`.
/// The `reader_init` closure has to return a result. Errors are returned from
/// the main function witout being mixed with reading errors. This may lead to
/// nested `Result` being returned if the `func` closure returns `Result`.
pub fn read_process_recordsets_init<R, Ri, Ei, W, F, O, Out>(
    reader_init: Ri,
    n_workers: u32,
    queue_len: usize,
    work: W,
    func: F,
) -> Result<Out, Ei>
where
    R: RecordSetReader,
    Ri: Send + FnOnce() -> Result<R, Ei>,
    R::RecordSet: Default + Send,
    O: Default + Send,
    W: Send + Sync + Fn(&mut R::RecordSet, &mut O),
    F: FnOnce(ParallelDataSets<R::RecordSet, R::Err, O>) -> Out,
    Ei: Send,
{
    let (done_send, done_recv) = mpsc::sync_channel(queue_len);
    let (empty_send, empty_recv) = mpsc::sync_channel(queue_len);

    crossbeam::scope(|scope| {
        let handle = scope.spawn::<_, Result<(), Ei>>(move |_| {
            let mut reader = reader_init()?;

            let mut pool = scoped_threadpool::Pool::new(n_workers);

            pool.scoped(|pool_scope| {
                let work = &work;

                loop {
                    // recycle an old RecordSet sent back after use
                    let (mut data, mut out) = if let Ok(r) = empty_recv.recv() {
                        r
                    } else {
                        // ParallelDataSets dropped -> stop
                        return;
                    };

                    let done_send = done_send.clone();

                    match reader.fill_data(&mut data) {
                        Ok(has_data) => {
                            if !has_data {
                                break;
                            }
                            // expensive work carried out by func()
                            pool_scope.execute(move || {
                                work(&mut data, &mut out);
                                done_send.send(Some(Ok((data, out)))).ok();
                            });
                        }
                        Err(e) => {
                            done_send.send(Some(Err(e))).ok();
                            break;
                        }
                    }
                }

                pool_scope.join_all();

                done_send.send(None).ok();
            });
            Ok(())
        });

        for _ in 0..queue_len {
            if empty_send
                .send((R::RecordSet::default(), O::default()))
                .is_err()
            {
                break;
            }
        }

        let dsets: ParallelDataSets<R::RecordSet, R::Err, O> = ParallelDataSets {
            empty_send,
            done_recv,
            current_recordset: (R::RecordSet::default(), O::default()),
        };

        let out = func(dsets);

        handle.join().unwrap()?;
        Ok(out)
    })
    .unwrap()
}

pub struct ParallelDataSets<D, E, O = ()>
where
    D: Send,
    E: Send,
    O: Send,
{
    empty_send: mpsc::SyncSender<(D, O)>,
    done_recv: mpsc::Receiver<Option<Result<(D, O), E>>>,
    current_recordset: (D, O),
}

impl<D, E, O> ParallelDataSets<D, E, O>
where
    D: Send,
    E: Send,
    O: Send,
{
    /// Returns a tuple of the next processed record set and processing results,
    /// if present.
    pub fn next(&mut self) -> Option<Result<(&mut D, &mut O), E>> {
        self.done_recv.recv().unwrap().map(move |result| {
            match result {
                Ok(d) => {
                    let prev_rset = std::mem::replace(&mut self.current_recordset, d);
                    self.empty_send.send(prev_rset).ok(); // error: channel closed is not a problem, happens after calling stop()
                    Ok((&mut self.current_recordset.0, &mut self.current_recordset.1))
                }
                Err(e) => Err(e),
            }
        })
    }
}

/// Allows generating functions equivalent to the `read_process_xy_records`
/// functions in this crate for your own types. This is rather a workaround
/// because the generic approach ([`read_process_records_init`](read_process_records_init))
/// does currently not work.
///
/// * `$format`: String specifying the name of the sequence format
///   (for genrated documentation)
/// * `$name`: name of the generated function
/// * `$name_init`: name of another generated function, which takes a closure
///   initializing the readers in the background thread.
/// * <X: Trait, ...>: Optional set of trait bounds to be added to the
///   functions. If none are to be added, specify `<>`.
/// * `$RecordSet`: record set type (see [`RecordSetReader::RecordSet`](RecordSetReader::RecordSet)).
///   In addition to the trait requirements, `&$RecordSet` needs to implement
///   `IntoIterator<Item=$Record>`.
/// * `$Record`: record type returned by the record set iterator.
/// * `$Error`: reading error type ([`RecordSetReader::Err`](RecordSetReader::Err))
#[macro_export]
macro_rules! parallel_record_impl {
    ($format:expr, $name:ident, $name_init:ident,
        ( $($bounds:tt)* ),
        $RecordSet:ty, $Record:ty, $Error:ty) => {
        _parallel_record_impl!(
            $format,
            $name,
            $name_init,
            ($($bounds)*),
            $RecordSet,
            $Record,
            $Error,
            concat!("[`", stringify!($name), "`](", stringify!($name), ")")
        );
    };
}

macro_rules! _parallel_record_impl {
    ($format:expr, $name:ident, $name_init:ident,
        ( $($bounds:tt)* ),
        $RecordSet:ty, $Record:ty, $Error:ty,
        $name_link:expr) => {

        /// This function wraps [`read_process_recordsets`](read_process_recordsets),
        /// hiding the complexity related to record sets and allowing it to
        /// directly work on
        #[doc = $format]
        /// sequence records.
        ///
        /// Apart from this, the process is similar:
        ///
        /// * The records are read (as part of record sets) in a background
        ///   thread.
        /// * Then they are sent to `n_workers` worker threads. Work is done
        ///   in the `work` closure supplied to this function.
        /// * Once ready, records an results are sent to the main thread,
        ///   where they are supplied to the `func` closure. The order of the
        ///   records may be different.
        pub fn $name<R, $($bounds)*, W, F, O, Out>(
            reader: R,
            n_workers: u32,
            queue_len: usize,
            work: W,
            func: F,
        ) -> Result<Option<Out>, R::Err>
        where
            R: RecordSetReader<RecordSet = $RecordSet, Err = $Error> + Send,
            O: Default + Send,
            W: Send + Sync + Fn($Record, &mut O),
            F: FnMut($Record, &mut O) -> Option<Out>,
        {
            let out: Result<_, $Error> = $name_init(
                || Ok(reader), n_workers, queue_len, work, func
            );
            out
        }

        /// Like
        #[doc = $name_link]
        ///, but instead of a [`RecordSetReader`](RecordSetReader), it takes a
        /// closure (`reader_init`) returning an `RecordSetReader` instance.
        /// This allows using readers that don't  implement `Send`.
        /// `reader_init` should return a result. The error type needs to
        /// implement `From<RecordSetReader::Err>`
        ///
        pub fn $name_init<R, Ri, $($bounds)*, W, F, O, Out, E>(
            reader_init: Ri,
            n_workers: u32,
            queue_len: usize,
            work: W,
            mut func: F,
        ) -> Result<Option<Out>, E>
        where
            R: RecordSetReader<RecordSet = $RecordSet, Err = $Error>,
            Ri: Send + FnOnce() -> Result<R, E>,
            O: Default + Send,
            W: Send + Sync + Fn($Record, &mut O),
            F: FnMut($Record, &mut O) -> Option<Out>,
            E: Send + From<R::Err>,
        {
            read_process_recordsets_init(
                reader_init,
                n_workers,
                queue_len,
                |rset: &mut $RecordSet, output: &mut Vec<O>| {
                    let mut record_iter = rset.into_iter();
                    for (out, record) in output.iter_mut().zip(&mut record_iter) {
                        work(record, out);
                    }
                    for record in record_iter {
                        output.push(O::default());
                        work(record, output.last_mut().unwrap());
                    }
                },
                |mut records| {
                    while let Some(result) = records.next() {
                        let (rset, out) = result?;
                        for (record, o) in rset.into_iter().zip(out.iter_mut()) {
                            if let Some(out) = func(record, o) {
                                return Ok(Some(out));
                            }
                        }
                    }
                    Ok(None)
                },
            ).and_then(From::from)
        }
    };
}

parallel_record_impl!(
    "FASTA",
    read_process_fasta_records,
    read_process_fasta_records_init,
    (S: SeqRecordPosition + Send + Sync),
    fasta::RecordSet<S>,
    fasta::RefRecord<S>,
    fasta::Error
);

parallel_record_impl!(
    "FASTQ",
    read_process_fastq_records,
    read_process_fastq_records_init,
    (S: QualRecordPosition + Send + Sync),
    fastq::RecordSet<S>,
    fastq::RefRecord<S>,
    fastq::Error
);

parallel_record_impl!(
    "FASTX",
    read_process_fastx_records,
    read_process_fastx_records_init,
    (S: QualRecordPosition + Send + Sync),
    fastx::RecordSet<S>,
    fastx::RefRecord<S>,
    fastx::Error
);

/// Using this function currently does not work due to a
/// [compiler bug](https://github.com/rust-lang/rust/issues/62529).
///
/// [`read_process_fasta_records`](read_process_fasta_records),
/// [`read_process_fastq_records`](read_process_fastq_records) and
/// [`read_process_fastx_records`](read_process_fastx_records)
///  provide the same functionality for now
/// (implemented using [`parallel_record_impl`](parallel_record_impl) macro).
pub fn read_process_records_init<R, Ri, W, F, O, Out, E>(
    reader_init: Ri,
    n_workers: u32,
    queue_len: usize,
    work: W,
    mut func: F,
) -> Result<Option<Out>, E>
where
    R: RecordSetReader,
    Ri: Send + FnOnce() -> Result<R, E>,
    R::RecordSet: Default + Send,
    for<'a> &'a R::RecordSet: IntoIterator + Send,
    O: Default + Send,
    W: Send + Sync + Fn(<&R::RecordSet as IntoIterator>::Item, &mut O),
    F: FnMut(<&R::RecordSet as IntoIterator>::Item, &mut O) -> Option<Out>,
    E: From<<R as RecordSetReader>::Err> + Send,
{
    read_process_recordsets_init(
        reader_init,
        n_workers,
        queue_len,
        |rset, out: &mut Vec<O>| {
            let mut record_iter = rset.into_iter();
            for mut d in (&mut record_iter).zip(out.iter_mut()) {
                work(d.0, &mut d.1);
            }
            for record in record_iter {
                out.push(O::default());
                work(record, out.last_mut().unwrap());
            }
        },
        |mut records| {
            while let Some(result) = records.next() {
                let (rset, out) = result?;
                for (record, o) in rset.into_iter().zip(out.iter_mut()) {
                    if let Some(out) = func(record, o) {
                        return Ok(Some(out));
                    }
                }
            }
            Ok(None)
        },
    )
    .and_then(From::from)
}

macro_rules! impl_parallel_reader {
    ($($l:lifetime)?; $SeqReader:ty, $RecordPositionTrait:path, $RecordSet:ty, $Error:ty, $read_fn:ident) => {
        impl<$($l,)? R, P, S> RecordSetReader for $SeqReader
        where
            R: std::io::Read,
            P: crate::policy::BufPolicy + Send,
            S: $RecordPositionTrait + Send + Sync
        {
            type RecordSet = $RecordSet;
            type Err = $Error;
            fn fill_data(&mut self, rset: &mut $RecordSet) -> Result<bool, $Error> {
                self.$read_fn(rset)
            }
        }
    }
}

impl_parallel_reader!(; fasta::Reader<R, P, S>, SeqRecordPosition, fasta::RecordSet<S>, fasta::Error, read_record_set);
impl_parallel_reader!(; fasta::single_line::Reader<R, P, S>, SeqRecordPosition, fasta::RecordSet<S>, fasta::Error, read_record_set);
impl_parallel_reader!(; fastq::Reader<R, P, S>, QualRecordPosition, fastq::RecordSet<S>, fastq::Error, read_record_set);
impl_parallel_reader!(; fastq::multiline::Reader<R, P, S>, QualRecordPosition, fastq::RecordSet<S>, fastq::Error, read_record_set);
impl_parallel_reader!(; fastx::Reader<R, P, S>, QualRecordPosition, fastx::RecordSet<S>, fastx::Error, read_record_set);
impl_parallel_reader!(; fastx::multiline_qual::Reader<R, P, S>, QualRecordPosition, fastx::RecordSet<S>, fastx::Error, read_record_set);
impl_parallel_reader!('a ; &'a mut (dyn fastx::dynamic::FastxReader<R, P, S> + Send), QualRecordPosition, fastx::RecordSet<S>, fastx::Error, read_record_set_fastx);
impl_parallel_reader!('a ; Box<dyn fastx::dynamic::FastxReader<R, P, S> + Send + 'a>, QualRecordPosition, fastx::RecordSet<S>, fastx::Error, read_record_set_fastx);


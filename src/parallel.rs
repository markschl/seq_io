//! Experiments with parallel processing
//!
//! The provided functions focus on the possibility of
//! returning results while the parser proceeds. Sequences are processesd in
//! batches (`RecordSet`) because sending across channels has a performance
//! impact. FASTA/FASTQ records can be accessed in both the 'worker' function and
//! (after processing) a function running in the main thread.
//!
//! # Search first occurrence of a sequence pattern
//!
//! ```no_run
//! use seq_io::fastq::{Reader,Record};
//! use seq_io::parallel::read_parallel;
//!
//! let reader = Reader::from_path("seqs.fastq").unwrap();
//!
//! read_parallel(reader, 4, 2, |record_set| {
//!     // this function does the heavy work
//!     for (i, record) in record_set.into_iter().enumerate() {
//!         // this is not very efficient code, just for demonstration
//!         if let Some(pos) = record.seq().windows(3).position(|s| s == b"AAA") {
//!             return Some((i, pos));
//!         }
//!     }
//!     None
//! }, |record_sets| {
//!     // This function runs in the main thread. It provides a streaming iterator over
//!     // record sets and the corresponding return values from the worker function
//!     // (not necessarily in the same order as in the file)
//!     while let Some(result) = record_sets.next() {
//!         let (record_set, found) = result.unwrap();
//!         if let Some((i, pos)) = found {
//!             let record = record_set.into_iter().nth(i).unwrap();
//!             println!("Found AAA in record {} at position {}", record.id().unwrap(), pos);
//!              // this will also stop the worker threads, although with some delay
//!             return;
//!         }
//!     }
//! });
//! ```
//!
//! # Per-record processsing
//! The `parallel_fasta` / `parallel_fastq` functions are designed to efficiently pass
//! results for **each record** to the main thread without having to care about record sets.
//! This example filters sequences by the occurrence of a pattern:
//!
//! ```no_run
//! use seq_io::fastq::{Reader,Record};
//! use seq_io::parallel::parallel_fastq;
//! use std::fs::File;
//! use std::io::BufWriter;
//!
//! let reader = Reader::from_path("seqs.fastq").unwrap();
//! let mut writer = BufWriter::new(File::create("filtered.fastq").unwrap());
//!
//! parallel_fastq(reader, 4, 2,
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


use std::io;

extern crate crossbeam;
extern crate scoped_threadpool;

use std::sync::mpsc;
use std::marker::PhantomData;



pub trait Reader: Send {
    type DataSet: Default + Send;
    type Err: Send;
    fn fill_data(&mut self, record: &mut Self::DataSet) -> Option<Result<(), Self::Err>>;
}


pub fn read_parallel<P, O, W, F, Out>(mut reader: P, n_threads: u32, queue_len: usize, work: W, func: F) -> Out
    where P: Reader,
          O: Default + Send,
          W: Send + Sync,
          W: Fn(&mut P::DataSet) -> O,
          F: FnOnce(&mut ParallelRecordsets<P::DataSet, P::Err, O>) -> Out,
{

    let (done_send, done_recv) = mpsc::sync_channel(queue_len);
    let (empty_send, empty_recv): (mpsc::SyncSender<Option<P::DataSet>>, _) = mpsc::sync_channel(queue_len);

    crossbeam::scope(|scope| {

        scope.spawn(move || {

            let mut pool = scoped_threadpool::Pool::new(n_threads);

            pool.scoped(|pool_scope| {

                let work = &work;

                loop {
                    // recycle an old DataSet sent back after use by the streaming iterator
                    let mut data =
                        if let Ok(Some(r)) = empty_recv.recv() {
                            r
                        } else {
                            // 'ParallelRecordsets::stop()' called
                            return;
                        };

                    // each time, we need a new reference
                    let done_send = done_send.clone();

                    if let Some(res) = reader.fill_data(&mut data) {
                        match res {
                            Ok(_) => {
                                // expensive work carried out by func()
                                pool_scope.execute(move || {

                                    let out = work(&mut data);

                                    done_send.send(Some(Ok((data, out)))).ok();
                                });
                            }
                            Err(e) => {
                                done_send.send(Some(Err(e))).ok();
                                break;
                            }
                        }
                    } else {
                        break;
                    }
                }

                // make sure that the 'done' signal is only sent after everything else is done
                pool_scope.join_all();

                done_send.send(None).ok();
            });
        });

        for _ in 0..queue_len  {
            empty_send.send(Some(P::DataSet::default())).unwrap();
        }

        let mut records = ParallelRecordsets {
            empty_send: empty_send,
            done_recv: done_recv,
            current_recordset: P::DataSet::default(),
        };

        let out = func(&mut records);

        records.stop();

        out
    })
}



pub struct ParallelRecordsets<R, E, O>
    where R: Default + Send,
          E: Send,
          O: Send,
{
    empty_send: mpsc::SyncSender<Option<R>>,
    done_recv: mpsc::Receiver<Option<Result<(R, O), E>>>,
    current_recordset: R,
}

impl<R, E, O> ParallelRecordsets<R, E, O>
    where R: Default + Send,
          E: Send,
          O: Send,
 {

    pub fn next<'b>(&'b mut self) -> Option<Result<(&'b R, O), E>> {
        if let Some(result) = self.done_recv.recv().unwrap() {
            match result {
                Ok((r, o)) => {
                     let prev_rset = ::std::mem::replace(&mut self.current_recordset, r);
                     self.empty_send.send(Some(prev_rset)).ok(); // error: channel closed is not a problem, happens after calling stop()
                    return Some(Ok((&self.current_recordset, o)));
                }
                Err(e) => return Some(Err(e))
            }
        } else {
            // 'done' signal received
            return None;
        }
    }

    // has to be called before object goes out of scope
    // currently, the signal is delayed (has to wait until it is popped from the
    // 'empty' queue)
    // TODO: should an AtomicBool be used?
    fn stop(self) {
        self.empty_send.send(None).ok();
    }
}


#[macro_export]
macro_rules! parallel_record_impl {
    ($name:ident, $io_r:tt, $rdr:ty, $record:ty, $err:ty) => {

        pub fn $name<$io_r, O, W, F, Out>(parser: $rdr, n_threads: u32, queue_len: usize, work: W, mut func: F) -> Result<Option<Out>, $err>
            where $io_r: io::Read + Send,
                  O: Default + Send,
                  W: Send + Sync,
                  W: Fn($record, &mut O),
                  F: FnMut($record, &O) -> Option<Out>,
        {
            let reader = $crate::parallel::ReusableReader::new(parser);

            $crate::parallel::read_parallel(reader, n_threads, queue_len, |d| {
                let mut iter = d.0.into_iter();
                let out: &mut Vec<O> = &mut d.1;
                for mut x in out.iter_mut().zip(&mut iter) {
                    work(x.1, &mut x.0);
                }
                for i in iter {
                    out.push(O::default());
                    work(i, out.last_mut().unwrap())
                }

            }, |records| {
                while let Some(result) = records.next() {
                    let (r, _) = result?;
                    for x in r.0.into_iter().zip(&r.1) {
                        if let Some(out) = func(x.0, x.1) {
                            return Ok(Some(out));
                        }
                    }
                }
                Ok(None)
            })
        }
     };
}

parallel_record_impl!(parallel_fasta, R, fasta::Reader<R>, fasta::RefRecord, fasta::ParseError);

parallel_record_impl!(parallel_fastq, R, fastq::Reader<R>, fastq::RefRecord, fastq::ParseError);


/// Wrapper for `parallel::Reader` instances allowing
/// the output to be reused in order to save allocations.
/// Used by `parallel_fasta`/`parallel_fastq`
///
/// ```no_run
/// use seq_io::fastq::{Reader,Record,RecordSet};
/// use seq_io::parallel::{read_parallel,ReusableReader};
///
/// let inner = Reader::from_path("seqs.fastq").unwrap();
/// let reader = ReusableReader::new(inner);
///
/// read_parallel(reader, 4, 2, |&mut (ref record_set, ref mut out): &mut (RecordSet, Vec<bool>)| {
///     out.clear();
///     for record in record_set {
///         let found = record.seq().windows(3).position(|s| s == b"AAA").is_some();
///         out.push(found);
///     }
/// }, |record_sets| {
///     while let Some(result) = record_sets.next() {
///         let &(ref record_set, ref out) = result.unwrap().0;
///         for (record, found) in record_set.into_iter().zip(out) {
///             // ...
///         }
///     }
/// });
/// ```
pub struct ReusableReader<P, O>(P, PhantomData<O>);


impl<P, O> ReusableReader<P, O> {
    pub fn new(p: P) -> ReusableReader<P, O> {
        ReusableReader(p, PhantomData)
    }
}

impl<P, O> Reader for ReusableReader<P, O>
    where P: Reader,
          O: Default + Send {

    type DataSet = (P::DataSet, O);
    type Err = P::Err;
    fn fill_data(&mut self, data: &mut Self::DataSet) -> Option<Result<(), P::Err>> {
        self.0.fill_data(&mut data.0)
    }
}

/// Using this function currently does not work due to a
/// [compiler bug](https://github.com/rust-lang/rust/issues/42950).
///
/// `parallel_fasta`/`parallel_fastq` provide the same functionality for now
/// (implemented using `parallel_record_impl` macro)
pub fn parallel_records<P, O, W, F, Out>(parser: P, n_threads: u32, queue_len: usize, work: W, mut func: F) -> Result<Option<Out>, P::Err>
    where P: Reader,
          for<'a> &'a P::DataSet: IntoIterator,
          O: Default + Send,
          W: Send + Sync,
          W: Fn(<&P::DataSet as IntoIterator>::Item, &mut O),
          F: FnMut(<&P::DataSet as IntoIterator>::Item, &O) -> Option<Out>,
{
    let reader = ReusableReader(parser, PhantomData);

    read_parallel(reader, n_threads, queue_len, |d| {
        let mut iter = d.0.into_iter();
        let out: &mut Vec<O> = &mut d.1;
        for mut x in out.iter_mut().zip(&mut iter) {
            work(x.1, &mut x.0);
        }
        for i in iter {
            out.push(O::default());
            work(i, out.last_mut().unwrap())
        }

    }, |records| {
        while let Some(result) = records.next() {
            let (r, _) = result?;
            for x in r.0.into_iter().zip(&r.1) {
                if let Some(out) = func(x.0, x.1) {
                    return Ok(Some(out));
                }
            }
        }
        Ok(None)
    })
}


// trait impls

use super::fasta;

impl<R, S> Reader for fasta::Reader<R, S>
    where R: io::Read + Send,
          S: super::strategy::BufGrowStrategy + Send
{
    type DataSet = fasta::RecordSet;
    type Err = fasta::ParseError;
    fn fill_data(&mut self, rset: &mut fasta::RecordSet) -> Option<Result<(), fasta::ParseError>> {
        self.read_record_set(rset)
    }
}


use super::fastq;

impl<R, S> Reader for fastq::Reader<R, S>
    where R: io::Read + Send,
          S: super::strategy::BufGrowStrategy + Send
{
    type DataSet = fastq::RecordSet;
    type Err = fastq::ParseError;
    fn fill_data(&mut self, rset: &mut fastq::RecordSet) -> Option<Result<(), fastq::ParseError>> {
        self.read_record_set(rset)
    }
}

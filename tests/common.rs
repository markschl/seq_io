#[macro_export]
macro_rules! validate_record {
    ($record:expr, $exp:expr) => {
        assert_eq!($record.id(), $exp.id, "id mismatch");
        assert_eq!($record.desc(), $exp.desc, "description mismatch");
        assert_eq!($record.head(), $exp.head, "head mismatch");
        assert_eq!($record.full_seq(), $exp.seq, "full sequence mismatch");
        let qual = $record.opt_full_qual().map(|q| q.to_vec());
        assert_eq!(
            qual.as_ref().map(|q| q.as_slice()),
            $exp.qual,
            "full quality mismatch"
        );
    };
}

#[macro_export]
macro_rules! validate_position {
    ($position:expr, $expected:expr) => {
        assert_eq!($position.line(), $expected.line(), "line mismatch");
        assert_eq!($position.byte(), $expected.byte(), "byte offset mismatch");
        assert_eq!(
            $position.record(),
            $expected.record(),
            "record index mismatch"
        );
    };
}

#[macro_export]
macro_rules! impl_common_tests {
    ($input:expr, $expected:expr, $RecordSet:ident, $test_reader:ident, $validate_ref:path,
    $next_fn:ident, $read_record_set_fn:ident, $seek_fn:ident, $parallel_reader_func:path) => {
        #[test]
        fn reader() {
            // try different initial capacities to test
            // buffer growing feature
            $test_reader!($input, reader, {
                let mut exp_iter = $expected.iter();
                while let Some(exp) = exp_iter.next() {
                    let record = reader.$next_fn().unwrap().unwrap();
                    // println!("{:?}",record.id());
                    $validate_ref!(record, exp);
                    validate_record!(record.to_owned_record(), exp);
                    validate_position!(reader.position(), exp.pos);
                }
                assert!(reader.$next_fn().is_none());
            });
        }

        #[test]
        fn empty() {
            let input = &b""[..];
            $test_reader!(input, reader, {
                assert!(reader.$next_fn().is_none());
            });
        }

        #[test]
        fn record_set() {
            $test_reader!($input, reader, {
                let mut rsets = vec![];
                loop {
                    let mut rset = $RecordSet::default();
                    if !reader.$read_record_set_fn(&mut rset).unwrap() {
                        break;
                    }
                    rsets.push(rset);
                }
                let mut rset_iter = rsets.iter().flat_map(|r| r.into_iter());

                for exp in $expected.into_iter() {
                    let rec = rset_iter.next().unwrap();
                    $validate_ref!(rec, exp);
                    validate_record!(rec.to_owned_record(), exp);
                }
            });
        }

        #[test]
        fn parallel() {
            $test_reader!($input, par_reader, {
                let mut i = 0;
                $parallel_reader_func(
                    par_reader,
                    1,
                    2,
                    |rec, out| {
                        // runs in worker
                        *out = rec.id().unwrap().to_owned();
                    },
                    |rec, out| {
                        // runs in main thread
                        let exp = &$expected[i];
                        $validate_ref!(rec, exp);
                        validate_record!(rec.to_owned_record(), exp);
                        assert_eq!(rec.id(), Ok(out.as_str()));
                        i += 1;
                        None::<()>
                    },
                )
                .unwrap();
            });
        }

        #[test]
        fn seek() {
            use Action::*;
            #[derive(Debug)]
            enum Action {
                Seek(usize),
                Read,
                ReadRecordset,
            }

            let actions = &[
                Read,
                ReadRecordset,
                Read,
                Seek(0),
                ReadRecordset,
                Read,
                Seek(2),
                Read,
                ReadRecordset,
                Read,
                Seek(1),
                Read,
                Read,
                Read,
                Read,
                ReadRecordset,
                Seek(0),
                Read,
                Seek(2),
                Read,
                Seek(1),
                Read,
                Read,
                Read,
            ];
            $test_reader!(std::io::Cursor::new($input), reader, {
                let mut i = 0;
                let mut record_set = $RecordSet::default();
                for action in actions {
                    // println!("{:?}", action);
                    match action {
                        Read => {
                            let has_next = i < $expected.len();
                            if let Some(res) = reader.$next_fn() {
                                assert!(has_next, "There should be no record at index {}", i);
                                let record = res.unwrap();
                                let exp = &$expected[i];
                                $validate_ref!(record, exp);
                                validate_record!(record.to_owned_record(), exp);
                                validate_position!(reader.position(), exp.pos);
                                i += 1;
                            } else {
                                assert!(!has_next, "There should be a record at index {}", i);
                            }
                        }
                        ReadRecordset => {
                            if reader.read_record_set(&mut record_set).unwrap() {
                                let n = record_set.len();
                                // println!("rcord set i {} n {}", i, n);
                                assert!(i <= $expected.len());
                                for (record, exp) in record_set.into_iter().zip(&$expected[i..]) {
                                    validate_record!(record.to_owned_record(), exp);
                                    //validate_position!(reader.position(), &exp.pos);
                                    // let pos = reader.position();
                                    // assert_eq!(pos.byte(), exp.pos.byte());
                                    // assert_eq!(pos.line(), exp.pos.line());
                                    // assert_eq!(pos.record(), exp.pos.record());
                                }
                                i += n;
                            } else {
                                // at the end
                                assert!(i == $expected.len());
                            }
                        }
                        Seek(to) => {
                            i = *to;
                            reader.$seek_fn(&$expected[i].pos).unwrap();
                        }
                    }
                }
            });
        }
    };
}

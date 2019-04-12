#![no_main]
#[macro_use] extern crate libfuzzer_sys;

include!("_fasta.rs");

fuzz_target!(|data: &[u8]| {
    evaluate(data);
});

mod fasta {
    include!("../../fuzz_targets/_fasta.rs");
}

mod fastq {
    include!("../../fuzz_targets/_fastq.rs");
}

use std::env::args;
use std::fs::File;
use std::io::Read;

fn main() {
    let mut data = vec![];
    let filename = args().skip(1).next().unwrap().as_str().to_string();
    File::open(&filename)
        .unwrap()
        .read_to_end(&mut data)
        .expect("could not open file");
    let data = data.as_slice();
    println!(
        "data: {:?}\n{:?}'",
        data,
        String::from_utf8(data.to_owned())
    );

    fasta::evaluate(data);
    fastq::evaluate(data);
}

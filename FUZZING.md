Fuzzing was very helpful for finding bugs, especially when implementing everything
from scratch for `v0.4`.

[cargo-fuzz](https://rust-fuzz.github.io/book/cargo-fuzz/tutorial.html) is used,
which in turn uses [libFuzzer](http://llvm.org/docs/LibFuzzer.html).

In order to catch as many problematic edge cases, a simple parser was built,
which can read FASTA, FASTQ and FASTX and shows the same behaviour as the
`seq_io` parsers. The code 
[is found here](fuzz/fuzz_targets/simple_reader.rs).

The random input is then parsed by the *seq\_io* readers *and* the "simple"
implementation, and the resulting records / errors are compared.

Features still not tested:

* `fastx::dynamic` readers
* `position()`

# Setup

Install cargo-fuzz:

```sh
cargo install cargo-fuzz
```

# Running

Two fuzzing targets were created, one for FASTA and FASTQ. The FASTX readers
are tested with both targets. These commands may be run in separate terminals:

```sh
cargo fuzz run --release fasta -- -only_ascii=1
cargo fuzz run --release fastq -- -only_ascii=1
```

# Debugging

If a problem was found, it often needed to be minified, here an example:

```sh
cargo fuzz tmin --release fasta fuzz/artifacts/fasta/crash-0a7cc920d077cd5454a397fbe8fd5833509c5086
```

Then, the minified input was analysed using a separate test crate. There may
be an easier way that I'm unaware of, but this works:

```sh
cd fuzz/fuzz_debug

cargo run ../artifacts/fasta/minimized-from-a7d26fe5cfaea49fb041c2f4fa8ca2e811c362ca
```

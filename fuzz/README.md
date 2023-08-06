[Cargo-fuzz](https://rust-fuzz.github.io/book/cargo-fuzz/tutorial.html) is used
to find more bugs.

# Setup

Install cargo-fuzz:

```sh
cargo install cargo-fuzz
```

# Running

The following example runs the fuzzing targets in four separate processes.

```sh
# simple parsers, checking only for panics
cargo fuzz run --release -j4 fasta -- -only_ascii=1
cargo fuzz run --release -j4 fastq -- -only_ascii=1
# comparison of seq_io::fasta::Reader with a simple FASTA reader
cargo fuzz run --release -j4 cmp_fasta_rdr -- -only_ascii=1
```

# Debugging

If a problem is found, it can be minified for debugging:

```sh
cargo fuzz tmin --release fasta fuzz/artifacts/fasta/crash-0a7cc920d077cd5454a397fbe8fd5833509c5086
```

Then, the minified input was analysed using a separate test crate. There may
be an easier way that I'm unaware of, but this works:

```sh
cd fuzz/fuzz_debug

cargo run ../artifacts/fasta/minimized-from-a7d26fe5cfaea49fb041c2f4fa8ca2e811c362ca
```

# Change Log

## v0.4.0 [...]

This release includes a rewrite of the whole library. The API didn't change
fundamentally, mostly new features were added:

* FASTX parsing
* Multi-line FASTQ parsing

Still, there are a few API changes (list may be incomplete):

* A set of new traits was introduced. The easiest is to include them by adding
  `use seq_io::prelude::*` on top.
* The behaviour of readers changed slightly in some cases. Refer to the `fasta`
  and `fastq` module docs for details on the new behaviour.
* The writing functions were renamed and their number reduced by combining 
  functionality.
* Error types now have an associated `ErrorKind`, which is similar in all
  formats.
* The `BufPolicy` trait was changed, and the standard policy now limits buffer
  growth to 1 GiB.
* `Reader::position()` now returns an owned copy of `Position`, not a reference.

**[Still not done:]**

* Review the API
* Evaluate, whether adding a `min_records()` method to `BufPolicy` will help.
  In cases where only a few records are in the buffer, frequent relocations are
  necessary, which could hurt performance -> buffer should be enlarged.
* Make it easier to work with paired-end files, especially in parallel
  processing.
* Evaluate even more code with fuzzing
* Resolve TODOs, write more tests, remove println! in comments.

## v0.3.1 [2020-10-18]

* Updated dependencies (fixing #4) and silenced warnings

#!/bin/bash

set -e

cargo bench "$@"

Rscript scripts/bench_analysis.R target/criterion

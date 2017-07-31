#!/bin/bash

# sudo cpupower frequency-set -u 2.2G

mkdir -p bench_results
scripts/bench.sh > bench_results/bench.txt

Rscript scripts/bench_analysis.R bench_results/bench.txt


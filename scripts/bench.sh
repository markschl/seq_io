#!/bin/bash

echo 'format method seqlen other owned parallel reader ns ns_dev mb_per_s'
cargo bench --benches | grep -e "^test " | grep -v "test result:" \
  | sed -E 's|test ([^ ]+)[^0-9]+([0-9,]+) ns/iter \(\+/\- ([0-9,]+)\) = ([0-9]+) MB/s|\1 \2 \3 \4|g' \
  | tr -d ',' | tr '_' ' '


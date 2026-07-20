#!/usr/bin/env bash
set -e

# Run FUSOR-hibff build
LD_LIBRARY_PATH=/home/okcut/FUSOR-hibff/jemalloc_install/lib /home/okcut/FUSOR-hibff/build/main/fusor build \
    --input-file /home/okcut/hibff-index/fusor_input.tsv \
    --input-sequence-dir /home/okcut/hibff-index/flat_genomes \
    --output-filename /home/okcut/FUSOR-hibff/fusor_hibff.hixf \
    --threads 4 \
    --kmer-size 22 \
    --bff-arity 4 \
    --largest-max-stash 100000 \
    --regular-max-stash 0 \
    --crypto

echo "Done building. Running verify_fusor..."
LD_LIBRARY_PATH=/home/okcut/FUSOR-hibff/jemalloc_install/lib /home/okcut/FUSOR-hibff/build/main/verify_fusor /home/okcut/FUSOR-hibff/fusor_hibff.hixf /home/okcut/hibff-index/binning.out > /home/okcut/FUSOR-hibff/verify_out.txt

grep FAIL /home/okcut/FUSOR-hibff/verify_out.txt || echo "ALL SUCCESS"

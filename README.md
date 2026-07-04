# FUSOR-hibff

FUSOR-hibff is an optimized pipeline and indexer for building Hybrid Hierarchical Interleaved Binary Fuse Filter (HHIBFF) indexes over large sets of FASTA/FASTQ genomic references. It wraps FUSOR.

## Requirements

- CMake (>= 3.16)
- GCC (>= 10)
- Linux / WSL environment

## Build

```bash
mkdir build && cd build
cmake ../src
make -j$(nproc)
```

The compiled binaries will be located in `build/main/`.

## Input Format

The `build_hibff` wrapper takes a File-Of-Files (FOF) input containing the absolute or relative paths to your target genomes:

```text
path/to/genome1.fa
path/to/genome2.fa.gz
```

Empty lines and lines starting with `#` are ignored. Supported sequence inputs include FASTA, FASTQ, and gzip-compressed FASTA/FASTQ.

## Building the Index

Use the `build_hibff` wrapper to construct an optimized HIBFF index.

```bash
build/main/build_hibff \
  --fof genomes.fof \
  --k 22 \
  --threads 4 \
  --interleaved \
  --optimize-memory \
  --output database.hixf
```

Important build options:

- `--fof`: File containing paths to reference genomes.
- `-k`, `--k`: Canonical k-mer length.
- `--interleaved`: Forces an interleaved HIBFF layout mapping technical bins directly to reference genomes. Leaving it out builds a hierarchical interleaved index.
- `--optimize-memory`: Forces the layout calculation algorithm to pick the layout that minimizes the index. The default behavior, without this flag, is to minimize the expected query time of the index. 
- `--fast-layout`: Optionally disables exhaustive layout calculation for a faster build.

## Searching the Index

You can query sequences against the generated `.hixf` index using the core `fusor search` binary.

```bash
build/main/fusor search \
  --index-file database.hixf \
  --query-file reads.fna \
  --output-file results.tsv \
  --threads 4
```

## Some notes
determine_best_number_of_technical_bins in src/main/taxor_build.cpp is responsible for calculating the optimal layout for the HHIBFF.
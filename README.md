# FUSOR-hibff

FUSOR-hibff is an optimized pipeline and indexer for building Hierarchical Interleaved Binary Fuse Filter (HIBFF) indexes over large sets of FASTA/FASTQ genomic references. It builds upon FUSOR to enable low false-positive rate k-mer matching with significantly reduced memory footprints and improved scaling for large datasets.

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
- `--interleaved`: Forces an interleaved HIBFF layout mapping technical bins directly to reference genomes.
- `--optimize-memory`: Employs an advanced layout calculation algorithm prioritizing the absolute minimum RAM and disk footprint over theoretical query speed.
- `--fast-layout`: Optionally disables exhaustive layout calculation for faster graph compilation.

## Searching the Index

You can query sequences against the generated `.hixf` index using the core `fusor search` binary.

```bash
build/main/fusor search \
  --index-file database.hixf \
  --query-file reads.fna \
  --output-file results.tsv \
  --threads 4
```

## Testing and Verification

You can verify the accuracy of the generated index using the provided `verify_fusor` utility. It compares index query results against an expected layout or binning output.

```bash
build/main/verify_fusor database.hixf expected_binning.out > verify_out.txt
grep FAIL verify_out.txt || echo "ALL SUCCESS"
```

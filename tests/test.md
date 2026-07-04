# Testing Interleaved Binary Fuse Filters

This directory contains rigorous testing logic for the 3-way and 4-way interleaved binary fuse filters.

## Compilation

To compile `test_interleaved_bff.cpp`, you need to link the necessary include directories for SeqAn3 and SDSL. Run the following command from the root `FUSOR` project directory:

```bash
g++ -std=c++20 -O3 -march=native \
    -I. \
    -Isrc \
    -Isrc/external_libs \
    -Ibuild/seqan/seqan3/seqan3-src/submodules/sdsl-lite/include \
    -Ibuild/seqan/seqan3/seqan3-src/include \
    -o test_interleaved_bff \
    tests/test_interleaved_bff.cpp
```

## Running the Tests

Once compiled, you can run the test suite. By default, the suite runs exactly one test iteration, testing both the 3-way and 4-way filters. Each test verifies false positive rates (~2^-8) and ensures zero false negatives while checking stash capacity logic.

```bash
./test_interleaved_bff
```

### Multiple Test Iterations

You can specify the `-n` parameter to run multiple continuous test iterations. For example, to run 100 consecutive tests verifying 100 different IBFF constructions (with different unique items to ensure the stash is rigorously validated 100 times):

```bash
./test_interleaved_bff -n 100
```

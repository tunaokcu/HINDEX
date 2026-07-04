# External Library Modifications

This directory contains modified versions of external libraries used by Taxor.

## SeqAn3 Modifications

**File:** `seqan3/search/dream_index/interleaved_xor_filter.hpp`

### Changes Made:

1. **Added `get_seed()` method** (line ~900)
   - Purpose: Allow debugging and verification that seeds are changing during construction
   - Returns the current seed value used for hashing

2. **Increased `max_stash` from 1 to 10** (line ~124)
   - Purpose: Allow more elements to be stashed during binary fuse filter construction
   - Helps with construction success for large bins (e.g., 447K+ elements)

### Reason for Local Copy:

These modifications are necessary for Taxor's hierarchical XOR filter construction to work reliably with large viral genome datasets. The changes were made to fix persistent construction failures when building indexes with bins containing hundreds of thousands of elements.

### Original Source:

SeqAn3 library downloaded via CMake from: https://github.com/seqan/seqan3

### Maintenance:

When updating SeqAn3 version, these modifications must be reapplied to the new version.

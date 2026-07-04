// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <utility>
#include <ankerl/unordered_dense.h>
#include <chrono>
#include <utility>

#include <cereal/archives/binary.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>

namespace taxor
{

//! Type alias for bin location: (ixf_idx, bin_idx)
using bin_location = std::pair<uint32_t, uint32_t>;

//! Type alias for stash map: kmer -> vector of bin locations
using stash_map_type = std::unordered_map<uint64_t, std::vector<bin_location>>;

//! Type alias for runtime inverted stash map: ixf_idx -> vector of (kmer, bin_idx)
using inverted_stash_type = ankerl::unordered_dense::map<int64_t, std::vector<std::pair<uint64_t, uint32_t>>>;

/**
 * @brief Invert the stash map from kmer -> bin locations to ixf_idx -> (kmer, bin_idx)
 */
inline inverted_stash_type invert_stash_map(stash_map_type const& stash_map)
{
    inverted_stash_type inverted;
    for (auto const& [kmer, locations] : stash_map)
    {
        for (auto const& loc : locations)
        {
            // loc.first = ixf_idx, loc.second = bin_idx
            inverted[loc.first].emplace_back(kmer, loc.second);
        }
    }
    return inverted;
}

/**
 * @brief Store the stash map (kmer -> bin locations) to a binary file.
 * @param path Path to the output file (e.g., "index.stashmap")
 * @param stash_map The stash map to serialize
 */
inline void store_stash_map(std::filesystem::path const& path,
                            stash_map_type const& stash_map)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(stash_map);
}

/**
 * @brief Load the stash map (kmer -> bin locations) from a binary file.
 * @param path Path to the input file (e.g., "index.stashmap")
 * @param stash_map The stash map to populate
 * @param io_time [out] Time spent on I/O in seconds (added to existing value)
 */
inline void load_stash_map(std::filesystem::path const& path,
                           stash_map_type& stash_map,
                           double& io_time)
{
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    auto start = std::chrono::high_resolution_clock::now();
    iarchive(stash_map);
    auto end = std::chrono::high_resolution_clock::now();

    io_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
}

/**
 * @brief Load the stash map without timing.
 */
inline void load_stash_map(std::filesystem::path const& path,
                           stash_map_type& stash_map)
{
    double unused_time = 0.0;
    load_stash_map(path, stash_map, unused_time);
}

} // namespace taxor


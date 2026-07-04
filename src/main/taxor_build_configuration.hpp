#pragma once

#include <filesystem>

namespace taxor::build
{

struct configuration
{
    std::string input_file_name{}; // provided by user
    std::vector<std::string> input_files{};
    std::string input_sequence_folder{};
    std::vector<std::string> input_folders{};
    std::string output_file_name{}; // provided by user
    std::string layout_file{}; // optional pre-existing layout file
    int threads{1u};
    int kmer_size{20u};
    int window_size{20u};
    int syncmer_size{10u};
    int scaling{1u};
    int max_stash{1u}; // maximum stash size for XOR filter construction
    bool output_verbose_statistics{false};
    bool debug{false};
    bool use_syncmer{false};
    bool use_xor{false};
    bool fast_layout{false};
    bool interleaved_only{false};
    bool optimize_memory{false};
    bool build_taxor{false};
    bool build_ganon{false};
    std::string level{"species"}; // default taxonomy level
    uint8_t bff_arity{3}; // Binary Fuse Filter arity: 3 or 4
};

} // namespace taxor::build

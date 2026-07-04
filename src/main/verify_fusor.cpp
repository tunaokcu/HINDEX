#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>

#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "index.hpp"
#include "store_stash_map.hpp"
#include <build/adjust_seed.hpp>
#include <build/dna4_traits.hpp>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <fusor.hixf> <binning.out>\n";
        return 1;
    }

    std::string index_path = argv[1];
    std::string binning_path = argv[2];

    std::cout << "Loading index from " << index_path << "..." << std::endl;
    taxor::taxor_index<taxor::hixf_t> index;
    {
        std::ifstream is(index_path, std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(index);
    }

    // Try to load stash map if it exists
    std::string stashmap_path = index_path.substr(0, index_path.find_last_of('.')) + ".stashmap";
    taxor::stash_map_type stash_map{};
    if (std::filesystem::exists(stashmap_path)) {
        std::cout << "Loading stash map from " << stashmap_path << "..." << std::endl;
        std::ifstream is(stashmap_path, std::ios::binary);
        cereal::BinaryInputArchive archive(is);
        archive(stash_map);
    }


    uint8_t kmer_size = index.kmer_size();
    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{kmer_size}},
                                                      seqan3::window_size{kmer_size},
                                                      seqan3::seed{hixf::adjust_seed(kmer_size)});
    using traits_type = hixf::dna4_traits;
    auto agent = index.ixf().membership_agent();
    
    std::ifstream binning_file(binning_path);
    std::string line;
    while (std::getline(binning_file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string genomes_str;
        if (!(iss >> genomes_str)) continue;
        
        std::stringstream ss(genomes_str);
        std::string single_genome_path;
        while(std::getline(ss, single_genome_path, ';')) {
            std::filesystem::path fname {single_genome_path};
            if (std::filesystem::is_symlink(fname)) {
                fname = std::filesystem::read_symlink(fname);
            }
            std::string resolved_path = fname.string();
            
            // Find expected user_bin_id for this genome by scanning index.ixf().user_bins
            int64_t expected_user_bin = -1;
            for (size_t i = 0; i < index.ixf().user_bins.num_user_bins(); ++i) {
                if (index.ixf().user_bins.filename_of_user_bin(i) == resolved_path) {
                    expected_user_bin = i;
                    break;
                }
            }
            if (expected_user_bin == -1) {
                std::cout << "Could not find mapping for " << resolved_path << " in FUSOR user bins!" << std::endl;
                continue;
            }
            // std::cout << "Genome " << resolved_path << " -> FUSOR user_bin_id " << expected_user_bin << std::endl;
            
            seqan3::sequence_file_input<traits_type, seqan3::fields<seqan3::field::seq>> fin{single_genome_path};
            size_t fn_count = 0;
            size_t kmer_count = 0;
            for (auto & record : fin) {
                auto hash_view = record.sequence() | hash_adaptor;
                for (uint64_t hash : hash_view) {
                    kmer_count++;
                    std::vector<uint64_t> query{hash};
                    auto const & result = agent.bulk_contains(query, 1);
                    bool found = false;
                    for (auto const & pair : result) {
                        if (pair.first == expected_user_bin) { found = true; break; }
                    }
                    if (!found) {
                        // Check stash map
                        auto stash_it = stash_map.find(hash);
                        if (stash_it != stash_map.end()) {
                            for (uint64_t u_bin : stash_it->second) {
                                if (u_bin == expected_user_bin) {
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                    if (!found) { fn_count++; }
                }
            }
            if (fn_count > 0) {
                std::cout << "FAIL: Genome " << single_genome_path << " has " << fn_count << "/" << kmer_count << " false negatives!\n";
            } else {
                std::cout << "SUCCESS: " << single_genome_path << " all kmers found!\n";
            }
        }
    }
    return 0;
}

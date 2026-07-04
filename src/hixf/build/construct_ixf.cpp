
#include <lemon/list_graph.h> /// Must be first include.

#include <seqan3/search/dream_index/interleaved_3way_binary_fuse_filter.hpp>
#include <seqan3/search/dream_index/interleaved_4way_binary_fuse_filter.hpp>

#include "bin_size_in_bits.hpp"
#include "construct_ixf.hpp"
#include "insert_into_ixf.hpp"
#include "temp_hash_file.hpp"

namespace hixf
{

/*
Does not even have max_stash????? TODO homogenize
*/
hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf(std::vector<ankerl::unordered_dense::set<size_t>> &node_hashes, bool use_xor, uint8_t bff_arity, uint32_t max_stash)
{
    std::vector<std::vector<size_t>> tmp{};

    for (const auto& hash_bin : node_hashes)
    {
        std::vector<size_t> c{};
        c.reserve(hash_bin.size());
        std::ranges::copy(hash_bin, std::back_inserter(c));
        tmp.emplace_back(std::move(c));
    }
    
    size_t max_bin_size = 0;
    for (const auto& bin : tmp)
        if (bin.size() > max_bin_size) max_bin_size = bin.size();

    hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t ixf;
    size_t floor_threshold = (bff_arity == 4) ? BFF4_FLOOR : BFF3_FLOOR;
    try {
        if (max_bin_size < floor_threshold || use_xor)
        {
            ixf = seqan3::interleaved_xor_filter<>{tmp, max_stash};
        }
        else
        {
            if (bff_arity == 4)
                ixf = seqan3::interleaved_4way_binary_fuse_filter<>{tmp, max_stash};
            else
                ixf = seqan3::interleaved_3way_binary_fuse_filter<>{tmp, max_stash};
        }
    } catch (std::exception const& e) {
        std::cerr << "Exception in construct_ixf(node_hashes): " << e.what() << std::endl;
        throw;
    }

    return std::move(ixf);
}


hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               bool is_second,
                                               size_t const & current_node_ixf_pos,
                                               uint32_t max_stash,
                                               bool use_xor,
                                               uint8_t bff_arity)
{
    auto &current_node_data = data.node_map[current_node];

    // Current guess for max_bin_size
    size_t max_bin_size = current_node_data.max_bin_hashes;
    bool new_max_bin_size = false;
    bool using_IBFF = false;

    auto create_filter = [&]() -> hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t {
        size_t floor_threshold = (bff_arity == 4) ? BFF4_FLOOR : BFF3_FLOOR;
        if (max_bin_size < floor_threshold || use_xor)
        {
            using_IBFF = false;
            return seqan3::interleaved_xor_filter<>{current_node_data.number_of_technical_bins, max_bin_size, max_stash};
        }
        else
        {
            using_IBFF = true;
            if (bff_arity == 4) {
                auto filter = seqan3::interleaved_4way_binary_fuse_filter<>{current_node_data.number_of_technical_bins, max_bin_size};
                filter.set_max_stash(max_stash);
                return filter;
            } else {
                auto filter = seqan3::interleaved_3way_binary_fuse_filter<>{current_node_data.number_of_technical_bins, max_bin_size};
                filter.set_max_stash(max_stash);
                return filter;
            }
        }
    };

    // create empty IXF/IBFF based on number of technical bins and max number of hashes per bin
    hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t ixf = create_filter();
    // first iterate over all child IXFs 
    
    bool success{false};

    std::map<size_t, int64_t> bins{};
    for (lemon::ListDigraph::OutArcIt arc_it(data.ixf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
    {
            auto child = data.ixf_graph.target(arc_it);
            auto& child_node_data = data.node_map[child];
            int64_t child_ixf_pos = ixf_positions[child_node_data.parent_bin_index];
            bins.insert(std::make_pair(child_node_data.parent_bin_index, child_ixf_pos));
    }

    ankerl::unordered_dense::set<std::string> tmp_files{};
    ankerl::unordered_dense::set<size_t> hashset{};

    bool has_failed = false;

    while (!success)
    {
        new_max_bin_size = false;
        success = true;
        size_t failed_bin_id = 0;
        for (std::map<size_t, int64_t>::iterator it = bins.begin(); it != bins.end(); ++it)
        {

            std::vector<size_t> hashes{};
            // read in hashes of child IXF and add all hashes to corresponding bin of current node IXF
            read_from_temp_hash_file((*it).second, hashes, tmp_files);

            // Update current max_bin_size guess
            if (hashes.size() > max_bin_size)
            {
                max_bin_size = hashes.size(); 
                new_max_bin_size = true;
            }

            success = std::visit([&](auto& f) { return f.add_bin_elements((*it).first, hashes); }, ixf);
            if(!success) {
                failed_bin_id = (*it).first;
                break;
            }
            
            if (is_second)
            {
                for (size_t hash : hashes)
                    hashset.insert(hash);
            }

        }
        // reset seed if adding bin to IXF was not successful
        if (!success)
        {   
            std::string filter_name = using_IBFF ? "IBFF" : "IXF";
            if (new_max_bin_size){
                ixf = create_filter();
            } else {
                std::visit([](auto& f) { f.clear(); }, ixf);
            }
            std::visit([](auto& f) { f.set_seed(); }, ixf);
            hashset.clear();
            has_failed = true;
            std::cerr << " [Failed on (" << filter_name << ", child bin): (" << current_node_ixf_pos << ", " << failed_bin_id << ")";
            if (new_max_bin_size)
                std::cerr << ", NEW BIN SIZE DISCOVERED!";
            std::cerr << ")] ... " << std::flush;
            continue;
        }

        // iterate over new hashes
        // add hashes of bins for newly computed hashes on that level
        for (uint16_t bin_idx = 0; bin_idx <= current_node_data.number_of_technical_bins; ++bin_idx)
        {
        
            if (bins.contains(bin_idx))
                continue;
            std::vector<size_t> c{};
            read_from_temp_hash_file(current_node_ixf_pos, bin_idx, c, tmp_files);

            // Update current max_bin_size guess
            if (c.size() > max_bin_size) {
                max_bin_size = c.size(); 
                new_max_bin_size = true;
            }
            
            if (c.size() == 0)
                continue;
            //std::ranges::copy(hash_bin, std::back_inserter(c));
            if ( c.size() > current_node_data.max_bin_hashes)
                std::cerr << "False max number of bin hashes: " << c.size() << "\t" << current_node_data.max_bin_hashes << std::endl;

            //std::cerr << bin_idx << "\t-\t" << c.size() << "\t" << current_node_data.max_bin_hashes << std::endl << std::flush;
            success = std::visit([&](auto& f) { return f.add_bin_elements(bin_idx, c); }, ixf);
            if(!success){
                failed_bin_id = bin_idx;
                break;
            }

            if (is_second)
            {
                for (size_t hash : c)
                    hashset.insert(hash);
            }

            //bin_idx++;
        }
        
        // reset seed if adding bin to IXF was not successful
        if (!success)
        {   
            std::string filter_name = using_IBFF ? "IBFF" : "IXF";
            if (new_max_bin_size){
                ixf = create_filter();
            } else {
                std::visit([](auto& f) { f.clear(); }, ixf);
            }
            std::visit([](auto& f) { f.set_seed(); }, ixf);
            hashset.clear();
            has_failed = true;
            std::cerr << " [Failed on (" << filter_name << ", technical bin): (" << current_node_ixf_pos << ", " << failed_bin_id << ")";
            if (new_max_bin_size)
                std::cerr << ", NEW BIN SIZE DISCOVERED!";
            std::cerr << ")] ... " << std::flush;
            continue;
        }
        
    }

    if (has_failed) {
        std::cerr << "Success!" << std::endl;
    }

    if (is_second)
    {
        current_node_data.number_of_hashes = hashset.size();

        create_temp_hash_file(current_node_ixf_pos, hashset);
        hashset.clear();
    }
    
    for (auto &file : tmp_files)
    {
        if (std::filesystem::exists(file))
            std::filesystem::remove(file);
    }

    return std::move(ixf);
}
}
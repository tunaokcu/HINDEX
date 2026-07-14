
#include <lemon/list_graph.h> /// Must be first include.

#include <seqan3/search/dream_index/interleaved_3way_binary_fuse_filter.hpp>
#include <seqan3/search/dream_index/interleaved_4way_binary_fuse_filter.hpp>

#include <atomic>
#include <omp.h>

#include "bin_size_in_bits.hpp"
#include "construct_ixf.hpp"
#include "insert_into_ixf.hpp"
#include "temp_hash_file.hpp"

namespace hixf
{

// Global pool to reuse large vectors and avoid millions of page faults
static std::vector<std::vector<size_t>> g_vector_pool;
static std::mutex g_vector_pool_mutex;

static std::vector<size_t> get_reusable_vector() {
    std::lock_guard<std::mutex> lock(g_vector_pool_mutex);
    if (!g_vector_pool.empty()) {
        auto v = std::move(g_vector_pool.back());
        g_vector_pool.pop_back();
        return v;
    }
    return {};
}

static void return_reusable_vector(std::vector<size_t>&& v) {
    v.clear();
    std::lock_guard<std::mutex> lock(g_vector_pool_mutex);
    g_vector_pool.push_back(std::move(v));
}

// Helper: number of 8-bit fingerprint bins that share a single uint64_t word
// For uint8_t fingerprints: 64 / 8 = 8 bins per word
static constexpr size_t BINS_PER_WORD = 8;

/*
Overload 1: construct from in-memory node_hashes
Used for lower-level IXFs (not root/second).
Uses sized constructor + parallel word-aligned batch add_bin_elements.
*/
hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf(std::vector<ankerl::unordered_dense::set<size_t>> &node_hashes, bool use_xor, uint8_t bff_arity, uint32_t max_stash, uint8_t threads, bool use_crypto_hash)
{
    size_t num_bins = node_hashes.size();

    // Convert hash sets to vectors and find max bin size
    std::vector<std::vector<size_t>> tmp{};
    tmp.reserve(num_bins);

    size_t max_bin_size = 0;
    for (const auto& hash_bin : node_hashes)
    {
        std::vector<size_t> c = get_reusable_vector();
        c.reserve(hash_bin.size());
        std::ranges::copy(hash_bin, std::back_inserter(c));
        if (c.size() > max_bin_size) max_bin_size = c.size();
        tmp.emplace_back(std::move(c));
    }

    if (max_bin_size == 0 || num_bins == 0)
    {
        // Edge case: return empty filter via batch constructor
        hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t ixf;
        size_t floor_threshold = (bff_arity == 4) ? BFF4_FLOOR : BFF3_FLOOR;
        try {
            if (max_bin_size < floor_threshold || use_xor)
                ixf = seqan3::interleaved_xor_filter<>{tmp, max_stash};
            else if (bff_arity == 4) {
                auto filter = seqan3::interleaved_4way_binary_fuse_filter<>{tmp, max_stash};
                filter.use_crypto_hash = use_crypto_hash;
                ixf = std::move(filter);
            } else {
                ixf = seqan3::interleaved_3way_binary_fuse_filter<>{tmp, max_stash};
            }
        } catch (std::exception const& e) {
            std::cerr << "Exception in construct_ixf(node_hashes): " << e.what() << std::endl;
            throw;
        }
        return std::move(ixf);
    }

    // Determine filter type and create sized filter
    bool using_IBFF = false;
    size_t floor_threshold = (bff_arity == 4) ? BFF4_FLOOR : BFF3_FLOOR;

    auto create_filter = [&]() -> hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t {
        if (max_bin_size < floor_threshold || use_xor)
        {
            using_IBFF = false;
            return seqan3::interleaved_xor_filter<>{num_bins, max_bin_size, max_stash};
        }
        else
        {
            using_IBFF = true;
            if (bff_arity == 4) {
                auto filter = seqan3::interleaved_4way_binary_fuse_filter<>{num_bins, max_bin_size};
                filter.set_max_stash(max_stash);
                filter.use_crypto_hash = use_crypto_hash;
                return filter;
            } else {
                auto filter = seqan3::interleaved_3way_binary_fuse_filter<>{num_bins, max_bin_size};
                filter.set_max_stash(max_stash);
                return filter;
            }
        }
    };

    hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t ixf = create_filter();

    // Compute word-aligned batches
    size_t num_batches = (num_bins + BINS_PER_WORD - 1) / BINS_PER_WORD;
    size_t actual_threads = std::min(static_cast<size_t>(threads), num_batches);

    bool success = false;
    bool has_failed = false;

    while (!success)
    {
        success = true;
        std::atomic<bool> any_failed{false};

        #pragma omp parallel for num_threads(actual_threads) schedule(dynamic)
        for (size_t batch = 0; batch < num_batches; ++batch)
        {
            if (any_failed.load(std::memory_order_relaxed)) continue;

            size_t bin_start = batch * BINS_PER_WORD;
            size_t bin_end = std::min(bin_start + BINS_PER_WORD, num_bins);

            for (size_t bin_idx = bin_start; bin_idx < bin_end; ++bin_idx)
            {
                if (any_failed.load(std::memory_order_relaxed)) break;
                if (tmp[bin_idx].empty()) continue;

                bool ok = std::visit([&](auto& f) { return f.add_bin_elements(bin_idx, tmp[bin_idx]); }, ixf);
                if (!ok)
                {
                    any_failed.store(true, std::memory_order_relaxed);
                    break;
                }
            }
        }

        if (any_failed)
        {
            success = false;
            has_failed = true;
            std::string filter_name = using_IBFF ? "IBFF" : "IXF";
            std::cerr << " [Failed on " << filter_name << " batch construction, retrying...] " << std::flush;
            std::visit([](auto& f) { f.clear(); }, ixf);
            std::visit([](auto& f) { f.set_seed(); }, ixf);
        }
    }

    if (has_failed) {
        std::cerr << "Success!" << std::endl;
    }

    for (auto& vec : tmp) {
        return_reusable_vector(std::move(vec));
    }

    return std::move(ixf);
}


/*
Overload 2: construct from temp hash files (root/second-level IXFs)
Restores low memory profile by streaming bins from disk in word-aligned batches
and constructing the filter incrementally.
*/
hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               bool is_second,
                                               size_t const & current_node_ixf_pos,
                                               uint32_t max_stash,
                                               bool use_xor,
                                               uint8_t bff_arity,
                                               uint8_t threads,
                                               bool use_crypto_hash)
{
    auto &current_node_data = data.node_map[current_node];

    // Current guess for max_bin_size
    size_t max_bin_size = current_node_data.max_bin_hashes;
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
                filter.use_crypto_hash = use_crypto_hash;
                return filter;
            } else {
                auto filter = seqan3::interleaved_3way_binary_fuse_filter<>{current_node_data.number_of_technical_bins, max_bin_size};
                filter.set_max_stash(max_stash);
                return filter;
            }
        }
    };

    hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t ixf = create_filter();

    // Collect child IXF bin positions
    std::map<size_t, int64_t> child_bins{};
    for (lemon::ListDigraph::OutArcIt arc_it(data.ixf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
    {
            auto child = data.ixf_graph.target(arc_it);
            auto& child_node_data = data.node_map[child];
            int64_t child_ixf_pos = ixf_positions[child_node_data.parent_bin_index];
            child_bins.insert(std::make_pair(child_node_data.parent_bin_index, child_ixf_pos));
    }

    ankerl::unordered_dense::set<std::string> tmp_files{};
    ankerl::unordered_dense::set<size_t> hashset{};

    bool has_failed = false;
    bool new_max_bin_size = false;
    bool success = false;

    struct bin_entry {
        size_t bin_idx;
        std::vector<size_t> hashes;
    };

    // Calculate number of batches total
    size_t num_tech_bins = current_node_data.number_of_technical_bins;
    size_t num_batches = (num_tech_bins + BINS_PER_WORD - 1) / BINS_PER_WORD;
    size_t batches_per_chunk = std::max(static_cast<size_t>(1), static_cast<size_t>(threads));

    while (!success)
    {
        new_max_bin_size = false;
        success = true;
        size_t failed_bin_id = 0;

        for (size_t chunk_start = 0; chunk_start < num_batches; chunk_start += batches_per_chunk)
        {
            size_t chunk_end = std::min(chunk_start + batches_per_chunk, num_batches);
            
            // loaded_chunk contains multiple word-aligned batches
            std::vector<std::vector<bin_entry>> loaded_chunk(chunk_end - chunk_start);

            // Load all batches for this chunk
            for (size_t batch = chunk_start; batch < chunk_end; ++batch)
            {
                size_t batch_idx_in_chunk = batch - chunk_start;
                size_t bin_start = batch * BINS_PER_WORD;
                size_t bin_end = std::min(bin_start + BINS_PER_WORD, num_tech_bins);

                for (size_t bin_idx = bin_start; bin_idx < bin_end; ++bin_idx)
                {
                    std::vector<size_t> hashes = get_reusable_vector();
                    if (child_bins.contains(bin_idx))
                    {
                        read_from_temp_hash_file(child_bins[bin_idx], hashes, tmp_files);
                    }
                    else
                    {
                        read_from_temp_hash_file(current_node_ixf_pos, bin_idx, hashes, tmp_files);
                        if (hashes.empty()) {
                            return_reusable_vector(std::move(hashes));
                            continue;
                        }
                    }

                    if (hashes.size() > max_bin_size)
                    {
                        max_bin_size = hashes.size();
                        new_max_bin_size = true;
                    }

                    loaded_chunk[batch_idx_in_chunk].push_back({bin_idx, std::move(hashes)});
                }
            }

            // Parallel insert for this chunk of batches
            // Each thread takes one full word-aligned batch to avoid word tearing
            std::atomic<bool> any_failed{false};
            size_t local_failed_bin_id = 0;

            size_t actual_threads = std::min(static_cast<size_t>(threads), loaded_chunk.size());

            #pragma omp parallel for num_threads(actual_threads) schedule(dynamic)
            for (size_t b = 0; b < loaded_chunk.size(); ++b)
            {
                if (any_failed.load(std::memory_order_relaxed)) continue;

                // Process the 8 bins within this batch sequentially to prevent data races in sdsl::int_vector
                for (auto& entry : loaded_chunk[b])
                {
                    if (any_failed.load(std::memory_order_relaxed)) break;

                    bool ok = std::visit([&](auto& f) { return f.add_bin_elements(entry.bin_idx, entry.hashes); }, ixf);
                    
                    if (!ok)
                    {
                        any_failed.store(true, std::memory_order_relaxed);
                        local_failed_bin_id = entry.bin_idx;
                        break;
                    }
                }
            }

            if (any_failed)
            {
                success = false;
                failed_bin_id = local_failed_bin_id;
                for (auto& batch_bins : loaded_chunk) {
                    for (auto& entry : batch_bins) {
                        return_reusable_vector(std::move(entry.hashes));
                    }
                }
                break; // Break chunk loop, trigger filter retry
            }

            // Successfully inserted, collect hashes for is_second if needed
            if (is_second)
            {
                for (auto& batch_bins : loaded_chunk)
                {
                    for (auto& entry : batch_bins)
                    {
                        for (size_t hash : entry.hashes)
                            hashset.insert(hash);
                    }
                }
            }

            // Return memory to pool
            for (auto& batch_bins : loaded_chunk) {
                for (auto& entry : batch_bins) {
                    return_reusable_vector(std::move(entry.hashes));
                }
            }
        }

        // Retry handler
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
            std::cerr << " [Failed on (" << filter_name << ", child/technical bin): (" << current_node_ixf_pos << ", " << failed_bin_id << ")";
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

hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf_two_pass(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               bool is_second,
                                               size_t const & current_node_ixf_pos,
                                               uint32_t max_stash,
                                               bool use_xor,
                                               uint8_t bff_arity,
                                               uint8_t threads,
                                               bool use_crypto_hash)
{
    auto &current_node_data = data.node_map[current_node];

    // Current guess for max_bin_size
    size_t max_bin_size = current_node_data.max_bin_hashes;
    bool using_IBFF = false;

    // Collect child IXF bin positions
    std::map<size_t, int64_t> child_bins{};
    for (lemon::ListDigraph::OutArcIt arc_it(data.ixf_graph, current_node); arc_it != lemon::INVALID; ++arc_it)
    {
            auto child = data.ixf_graph.target(arc_it);
            auto& child_node_data = data.node_map[child];
            int64_t child_ixf_pos = ixf_positions[child_node_data.parent_bin_index];
            child_bins.insert(std::make_pair(child_node_data.parent_bin_index, child_ixf_pos));
    }

    size_t num_tech_bins = current_node_data.number_of_technical_bins;

    // First pass: definitively find max_bin_size
    #pragma omp parallel for num_threads(threads) schedule(dynamic) reduction(max:max_bin_size)
    for (size_t bin_idx = 0; bin_idx < num_tech_bins; ++bin_idx)
    {
        ankerl::unordered_dense::set<std::string> local_tmp_files{};
        std::vector<size_t> hashes = get_reusable_vector();
        if (child_bins.contains(bin_idx))
        {
            read_from_temp_hash_file(child_bins.at(bin_idx), hashes, local_tmp_files);
        }
        else
        {
            read_from_temp_hash_file(current_node_ixf_pos, bin_idx, hashes, local_tmp_files);
        }
        
        if (hashes.size() > max_bin_size)
        {
            max_bin_size = hashes.size();
        }
        
        return_reusable_vector(std::move(hashes));
    }

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
                filter.use_crypto_hash = use_crypto_hash;
                return filter;
            } else {
                auto filter = seqan3::interleaved_3way_binary_fuse_filter<>{current_node_data.number_of_technical_bins, max_bin_size};
                filter.set_max_stash(max_stash);
                return filter;
            }
        }
    };

    hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t ixf = create_filter();

    ankerl::unordered_dense::set<std::string> tmp_files{};
    ankerl::unordered_dense::set<size_t> hashset{};

    bool has_failed = false;
    bool success = false;

    struct bin_entry {
        size_t bin_idx;
        std::vector<size_t> hashes;
    };

    // Calculate number of batches total
    size_t num_batches = (num_tech_bins + BINS_PER_WORD - 1) / BINS_PER_WORD;
    size_t batches_per_chunk = std::max(static_cast<size_t>(1), static_cast<size_t>(threads));

    while (!success)
    {
        success = true;
        size_t failed_bin_id = 0;

        for (size_t chunk_start = 0; chunk_start < num_batches; chunk_start += batches_per_chunk)
        {
            size_t chunk_end = std::min(chunk_start + batches_per_chunk, num_batches);
            
            // loaded_chunk contains multiple word-aligned batches
            std::vector<std::vector<bin_entry>> loaded_chunk(chunk_end - chunk_start);

            // Load all batches for this chunk
            for (size_t batch = chunk_start; batch < chunk_end; ++batch)
            {
                size_t batch_idx_in_chunk = batch - chunk_start;
                size_t bin_start = batch * BINS_PER_WORD;
                size_t bin_end = std::min(bin_start + BINS_PER_WORD, num_tech_bins);

                for (size_t bin_idx = bin_start; bin_idx < bin_end; ++bin_idx)
                {
                    std::vector<size_t> hashes = get_reusable_vector();
                    if (child_bins.contains(bin_idx))
                    {
                        read_from_temp_hash_file(child_bins.at(bin_idx), hashes, tmp_files);
                    }
                    else
                    {
                        read_from_temp_hash_file(current_node_ixf_pos, bin_idx, hashes, tmp_files);
                        if (hashes.empty()) {
                            return_reusable_vector(std::move(hashes));
                            continue;
                        }
                    }

                    loaded_chunk[batch_idx_in_chunk].push_back({bin_idx, std::move(hashes)});
                }
            }

            // Parallel insert for this chunk of batches
            // Each thread takes one full word-aligned batch to avoid word tearing
            std::atomic<bool> any_failed{false};
            size_t local_failed_bin_id = 0;

            size_t actual_threads = std::min(static_cast<size_t>(threads), loaded_chunk.size());

            #pragma omp parallel for num_threads(actual_threads) schedule(dynamic)
            for (size_t b = 0; b < loaded_chunk.size(); ++b)
            {
                if (any_failed.load(std::memory_order_relaxed)) continue;

                // Process the 8 bins within this batch sequentially to prevent data races in sdsl::int_vector
                for (auto& entry : loaded_chunk[b])
                {
                    if (any_failed.load(std::memory_order_relaxed)) break;

                    bool ok = std::visit([&](auto& f) { return f.add_bin_elements(entry.bin_idx, entry.hashes); }, ixf);
                    
                    if (!ok)
                    {
                        any_failed.store(true, std::memory_order_relaxed);
                        local_failed_bin_id = entry.bin_idx;
                        break;
                    }
                }
            }

            if (any_failed)
            {
                success = false;
                failed_bin_id = local_failed_bin_id;
                for (auto& batch_bins : loaded_chunk) {
                    for (auto& entry : batch_bins) {
                        return_reusable_vector(std::move(entry.hashes));
                    }
                }
                break; // Break chunk loop, trigger filter retry
            }

            // Successfully inserted, collect hashes for is_second if needed
            if (is_second)
            {
                for (auto& batch_bins : loaded_chunk)
                {
                    for (auto& entry : batch_bins)
                    {
                        for (size_t hash : entry.hashes)
                            hashset.insert(hash);
                    }
                }
            }

            // Return memory to pool
            for (auto& batch_bins : loaded_chunk) {
                for (auto& entry : batch_bins) {
                    return_reusable_vector(std::move(entry.hashes));
                }
            }
        }

        // Retry handler
        if (!success)
        {   
            std::string filter_name = using_IBFF ? "IBFF" : "IXF";
            std::visit([](auto& f) { f.clear(); }, ixf);
            std::visit([](auto& f) { f.set_seed(); }, ixf);
            hashset.clear();
            has_failed = true;
            std::cerr << " [Failed on (" << filter_name << ", child/technical bin): (" << current_node_ixf_pos << ", " << failed_bin_id << ")";
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

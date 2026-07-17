import sys

with open('/home/okcut/FUSOR-hibff/src/hixf/build/construct_ixf.cpp', 'r') as f:
    content = f.read()

# We want to insert our new function right before the last closing brace.
# Find the last '}\n'
last_brace_idx = content.rfind('}\n')

new_func = """
hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf_two_pass(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               bool is_second,
                                               size_t const & current_node_ixf_pos,
                                               uint32_t largest_max_stash,
                                               uint32_t regular_max_stash,
                                               bool use_xor,
                                               uint8_t bff_arity,
                                               uint8_t threads)
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
            return seqan3::interleaved_xor_filter<>{current_node_data.number_of_technical_bins, max_bin_size, largest_max_stash, regular_max_stash};
        }
        else
        {
            using_IBFF = true;
            if (bff_arity == 4) {
                auto filter = seqan3::interleaved_4way_binary_fuse_filter<>{current_node_data.number_of_technical_bins, max_bin_size};
                filter.set_max_stash(largest_max_stash, regular_max_stash);
                return filter;
            } else {
                auto filter = seqan3::interleaved_3way_binary_fuse_filter<>{current_node_data.number_of_technical_bins, max_bin_size};
                filter.set_max_stash(largest_max_stash, regular_max_stash);
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

                    bool ok = std::visit([&](auto& f) { 
                        uint32_t current_max_stash = (entry.hashes.size() == max_bin_size) ? largest_max_stash : regular_max_stash;
                        return f.add_bin_elements(entry.bin_idx, entry.hashes, current_max_stash); 
                    }, ixf);
                    
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

"""

new_content = content[:last_brace_idx] + new_func + content[last_brace_idx:]

with open('/home/okcut/FUSOR-hibff/src/hixf/build/construct_ixf.cpp', 'w') as f:
    f.write(new_content)


// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/JensUweUlrich/Taxor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


#include <lemon/list_graph.h> /// Must be first include.

#include "create_ixfs_from_chopper_pack.hpp"
#include "hierarchical_build.hpp"
#include "read_chopper_pack_file.hpp"
#include <seqan3/search/dream_index/interleaved_xor_filter.hpp>

#include <seqan3/search/views/minimiser_hash.hpp>
#include <build/adjust_seed.hpp>
#include <build/dna4_traits.hpp>
#include <syncmer.hpp>

namespace hixf
{

void create_ixfs_from_chopper_pack(build_data& data, build_arguments const & arguments)
{   

    read_chopper_pack_file(data, arguments.bin_file);

    lemon::ListDigraph::Node root = data.ixf_graph.nodeFromId(0); // root node = high level IXF node
    ankerl::unordered_dense::set<size_t> root_hashes{};

    // TODO t_max is likely dead code, remove it
    size_t const t_max{data.node_map[root].number_of_technical_bins};

    hierarchical_build(root_hashes, root, data, arguments, true, false, false);
}

} 

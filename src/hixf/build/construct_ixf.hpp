// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

//#include "robin_hood.h"
#include <ankerl/unordered_dense.h>
#include "build_arguments.hpp"
#include "build_data.hpp"

namespace hixf
{

constexpr size_t BFF3_FLOOR = 17000; //16,802
constexpr size_t BFF4_FLOOR = 6800;  //6,779


hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf(build_data & data, 
                                               lemon::ListDigraph::Node const & current_node,
                                               std::vector<int64_t> & ixf_positions,
                                               bool is_second,
                                               size_t const & current_node_ixf_position,
                                               uint32_t max_stash = 1,
                                               bool use_xor = false,
                                               uint8_t bff_arity = 3);

hixf::hierarchical_interleaved_xor_filter<uint8_t>::ixf_t construct_ixf(std::vector<ankerl::unordered_dense::set<size_t>> &node_hashes, bool use_xor = false, uint8_t bff_arity = 3, uint32_t max_stash = 1);


} // namespace raptor::hibf

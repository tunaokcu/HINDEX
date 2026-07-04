#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>
#include <seqan3/search/dream_index/interleaved_3way_binary_fuse_filter.hpp>
#include <seqan3/search/dream_index/interleaved_4way_binary_fuse_filter.hpp>

std::vector<size_t> generate_unique_keys(size_t count, std::mt19937_64& gen) {
    std::vector<size_t> keys;
    keys.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        keys.push_back(gen());
    }
    return keys;
}

template <typename FilterType>
void test_filter(int arity, std::mt19937_64& gen) {
    std::cout << "Testing " << arity << "-way interleaved binary fuse filter..." << std::endl;
    
    // We use 64 bins (1 batch word) and enough elements to force victim selection
    const size_t num_bins = 64;
    const size_t elements_per_bin = 15000;
    
    uint32_t max_stash = 200;
    
    // We will retry with different seeds until the stash is actually used
    size_t total_stash = 0;
    std::unique_ptr<FilterType> filter_ptr;
    std::vector<std::vector<size_t>> bins(num_bins);
    
    while (total_stash == 0) {
        for (size_t i = 0; i < num_bins; ++i) {
            bins[i] = generate_unique_keys(elements_per_bin, gen);
        }
        
        filter_ptr = std::make_unique<FilterType>(bins, max_stash);
        total_stash = filter_ptr->stash_usage();
        
        if (total_stash == 0) {
            std::cout << "Stash was not used, retrying with new random keys..." << std::endl;
        }
    }
    
    FilterType& filter = *filter_ptr;
    std::cout << "Total stash usage: " << total_stash << std::endl;
    std::cout << "SUCCESS: Stash is used. Testing correctness..." << std::endl;
    
    // 1. Verify NO false negatives
    std::cout << "Verifying false negatives (expecting 0)..." << std::endl;
    auto agent = filter.membership_agent();
    size_t false_negatives = 0;
    for (size_t b = 0; b < num_bins; ++b) {
        for (size_t key : bins[b]) {
            auto const & result = agent.bulk_contains(key);
            if (!result[b]) {
                // If it's not in the main filter array, check the stash
                if (!filter.stash_contains(b, key)) {
                    false_negatives++;
                }
            }
        }
    }
    std::cout << "False negatives: " << false_negatives << std::endl;
    if (false_negatives != 0) {
        std::cerr << "test failed! Filter should have NO false negatives! Count: " << false_negatives << std::endl;
        std::exit(1);
    }
    
    // 2. Verify false positive rate
    std::cout << "Verifying false positive rate..." << std::endl;
    size_t test_keys = 100000;
    size_t false_positives = 0;
    for (size_t i = 0; i < test_keys; ++i) {
        size_t key = gen();
        auto const & result = agent.bulk_contains(key);
        for (size_t b = 0; b < num_bins; ++b) {
            if (result[b]) {
                false_positives++;
            }
        }
    }
    
    double expected_fpr = std::pow(2, -8); // uint8_t fingerprint -> ~1/256
    double actual_fpr = static_cast<double>(false_positives) / (test_keys * num_bins);
    std::cout << "Expected FPR: " << expected_fpr << " (" << expected_fpr * 100 << "%)" << std::endl;
    std::cout << "Actual FPR:   " << actual_fpr << " (" << actual_fpr * 100 << "%)" << std::endl;
    
    // 20% tolerance on FPR
    double margin = 0.2;
    if (actual_fpr < expected_fpr * (1.0 - margin) || actual_fpr > expected_fpr * (1.0 + margin)) {
        std::cerr << "test failed! Actual FPR (" << actual_fpr << ") deviates too much from expected FPR (" << expected_fpr << ")." << std::endl;
        std::exit(1);
    }
    
    std::cout << arity << "-way test passed!\n" << std::endl;
}

int main(int argc, char* argv[]) {
    int num_tests = 1;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-n" && i + 1 < argc) {
            num_tests = std::stoi(argv[++i]);
        }
    }

    std::mt19937_64 gen(1337);
    for (int i = 0; i < num_tests; ++i) {
        std::cout << "\n=== Test Iteration " << (i + 1) << " / " << num_tests << " ===" << std::endl;
        test_filter<seqan3::interleaved_3way_binary_fuse_filter<uint8_t>>(3, gen);
        test_filter<seqan3::interleaved_4way_binary_fuse_filter<uint8_t>>(4, gen);
    }
    
    std::cout << "All " << num_tests << " test iterations passed successfully!" << std::endl;
    return 0;
}

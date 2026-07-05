
#include "temp_hash_file.hpp"
#include <filesystem>
#include <fstream>

namespace hixf
{

// Helper function to get a unique temporary directory path
// Returns hixf_tmp, hixf_tmp1, hixf_tmp2, etc. depending on availability
static std::filesystem::path get_tmp_directory()
{
    static std::filesystem::path cached_tmp_dir = []() {
        auto base_dir = std::filesystem::current_path();
        
        auto tmp_dir = base_dir / "hixf_tmp";
        if (!std::filesystem::exists(tmp_dir))
        {
            std::filesystem::create_directory(tmp_dir);
            return tmp_dir;
        }
        
        for (int i = 1; i < 1000; ++i)
        {
            tmp_dir = base_dir / ("hixf_tmp" + std::to_string(i));
            if (!std::filesystem::exists(tmp_dir))
            {
                std::filesystem::create_directory(tmp_dir);
                return tmp_dir;
            }
        }
        
        return base_dir / "hixf_tmp";
    }();
    
    return cached_tmp_dir;
}

void create_temp_hash_file(size_t const ixf_pos, ankerl::unordered_dense::set<size_t> &node_hashes)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_pos) + ".tmp";
    auto tmp_dir = get_tmp_directory();
    auto tmp_file = tmp_dir / ixf_tmp_name;
    std::ofstream tmp_stream(tmp_file, std::ios::out | std::ios::trunc | std::ios::binary);
    
    if (!node_hashes.empty())
    {
        std::vector<size_t> buffer(node_hashes.begin(), node_hashes.end());
        tmp_stream.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(size_t));
    }
    
    tmp_stream.close();
}

void create_temp_hash_file(size_t const ixf_pos, size_t const bin_index, ankerl::unordered_dense::set<size_t> &node_hashes)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_pos) + "_" + std::to_string(bin_index) + ".tmp";
    auto tmp_dir = get_tmp_directory();
    auto tmp_file = tmp_dir / ixf_tmp_name;
    std::ofstream tmp_stream(tmp_file, std::ios::out | std::ios::trunc | std::ios::binary);
    
    if (!node_hashes.empty())
    {
        std::vector<size_t> buffer(node_hashes.begin(), node_hashes.end());
        tmp_stream.write(reinterpret_cast<const char*>(buffer.data()), buffer.size() * sizeof(size_t));
    }

    tmp_stream.close();
}


void read_from_temp_hash_file(int64_t & ixf_position,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_position) + ".tmp";
    auto tmp_dir = get_tmp_directory();
    auto tmp_file = tmp_dir / ixf_tmp_name;
    if (!std::filesystem::exists(tmp_file))
    {
        std::cerr << ixf_tmp_name << "does not exist!" << std::endl << std::flush;
        return;
    }
    
    auto file_size = std::filesystem::file_size(tmp_file);
    if (file_size > 0)
    {
        size_t num_hashes = file_size / sizeof(size_t);
        size_t current_size = node_hashes.size();
        node_hashes.resize(current_size + num_hashes);
        
        std::ifstream tmp_stream(tmp_file, std::ios::in | std::ios::binary);
        tmp_stream.read(reinterpret_cast<char*>(node_hashes.data() + current_size), file_size);
        tmp_stream.close();
    }
    
    tmp_files.insert(tmp_file.string());
}

void read_from_temp_hash_file(size_t const ixf_position,
                              uint16_t const bin_index,
                              std::vector<size_t> &node_hashes,
                              ankerl::unordered_dense::set<std::string>& tmp_files)
{
    std::string ixf_tmp_name = "interleavedXOR_" + std::to_string(ixf_position) + "_" + std::to_string(bin_index) + ".tmp";
    auto tmp_dir = get_tmp_directory();
    auto tmp_file = tmp_dir / ixf_tmp_name;
    if (!std::filesystem::exists(tmp_file))
    {
        return;
    }
    
    auto file_size = std::filesystem::file_size(tmp_file);
    if (file_size > 0)
    {
        size_t num_hashes = file_size / sizeof(size_t);
        size_t current_size = node_hashes.size();
        node_hashes.resize(current_size + num_hashes);
        
        std::ifstream tmp_stream(tmp_file, std::ios::in | std::ios::binary);
        tmp_stream.read(reinterpret_cast<char*>(node_hashes.data() + current_size), file_size);
        tmp_stream.close();
    }
    
    tmp_files.insert(tmp_file.string());
}

}
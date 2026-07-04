#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

int main(int argc, char** argv) {
    std::string fof_file = "";
    int kmer_size = 20;
    int threads = 1;
    bool interleaved = false;
    bool fast_layout = true;
    std::string output_file = "index.hixf";
    bool optimize_memory = false;

    // Parse arguments matching the original hibff-index build_hibff interface
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-fof" || arg == "--fof") {
            if (i + 1 < argc) fof_file = argv[++i];
        } else if (arg == "-k" || arg == "--kmer") {
            if (i + 1 < argc) kmer_size = std::stoi(argv[++i]);
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) threads = std::stoi(argv[++i]);
        } else if (arg == "-o" || arg == "--out") {
            if (i + 1 < argc) output_file = argv[++i];
        } else if (arg == "--interleaved") {
            interleaved = true;
        } else if (arg == "--optimize-memory") {
            optimize_memory = true;
        } else if (arg == "--fast-layout") {
            if (i + 1 < argc) {
                std::string val = argv[++i];
                if (val == "false" || val == "0") fast_layout = false;
                else fast_layout = true;
            }
        } else if (arg == "-m" || arg == "--memory") {
            if (i + 1 < argc) ++i; // Safely ignored, handled natively by FUSOR disk streaming
        }
    }

    if (fof_file.empty()) {
        std::cerr << "Error: -fof is required\n";
        std::cerr << "Usage: " << argv[0] << " -fof <files.txt> -o <index.hixf> -k <kmer_size> -t <threads> [--fast-layout true|false] [--interleaved]\n";
        return 1;
    }

    // 1. Generate TSV file for FUSOR (Format: SpeciesName \t DummyRank \t FilePath)
    std::string tsv_file = output_file + ".fusor_input.tsv";
    std::ifstream in_fof(fof_file);
    std::ofstream out_tsv(tsv_file);
    
    if (!in_fof || !out_tsv) {
        std::cerr << "Error opening files for TSV generation\n";
        return 1;
    }

    std::string line;
    int id = 1;
    while (std::getline(in_fof, line)) {
        if (line.empty()) continue;
        
        // If FOF has "TaxID /path/to/file", extract just the path, or use the taxID if we want
        std::string tax_id = "dummy_" + std::to_string(id);
        std::string filepath = line;
        
        size_t space_pos = line.find_first_of(" \t");
        if (space_pos != std::string::npos) {
            tax_id = line.substr(0, space_pos);
            size_t path_start = line.find_first_not_of(" \t", space_pos);
            if (path_start != std::string::npos) {
                filepath = line.substr(path_start);
            }
        }
        
        out_tsv << tax_id << "\t" << id << "\t" << filepath << "\n";
        id++;
    }
    in_fof.close();
    out_tsv.close();

    // 2. Call FUSOR build
    // Resolving absolute path to fusor to avoid issues if run from elsewhere
    std::string fusor_bin = std::string(argv[0]);
    size_t last_slash = fusor_bin.find_last_of('/');
    if (last_slash != std::string::npos) {
        fusor_bin = fusor_bin.substr(0, last_slash) + "/fusor";
    } else {
        fusor_bin = "./fusor";
    }

    std::string cmd = fusor_bin + " build "
                      "--input-file " + tsv_file + " " +
                      "--output-filename " + output_file + " " +
                      "--threads " + std::to_string(threads) + " " +
                      "--kmer-size " + std::to_string(kmer_size) + " " +
                      "--window-size " + std::to_string(kmer_size) + " " +
                      "--bff-arity 4";
                      
    if (fast_layout) {
        cmd += " --fast-layout";
    }
    
    if (interleaved) {
        cmd += " --interleaved";
    }

    if (optimize_memory) {
        cmd += " --optimize-memory";
    }

    std::cout << "Wrapper generated TSV and is now running FUSOR:\n" << cmd << std::endl;
    int ret = std::system(cmd.c_str());
    
    if (ret != 0) {
        std::cerr << "FUSOR build failed with exit code " << ret << "!\n";
        return ret;
    }
    
    std::cout << "Index successfully built at " << output_file << std::endl;
    return 0;
}

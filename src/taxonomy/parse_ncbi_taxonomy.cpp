
#include "parse_ncbi_taxonomy.hpp"
#include "taxutil.hpp"
#include <stdexcept>

namespace taxor::taxonomy
{

    /**
     * @brief Parses a RefSeq taxonomy TSV file.
     * 
     * The file is expected to be a tab-separated values (TSV) file where each line 
     * represents a species/genome with the following columns:
     *   - Column 0: Accession ID (e.g., "NZ_CP012345.1")
     *   - Column 1: Taxonomy ID (e.g., "562")
     *   - Column 2: File path or stem (e.g., "/path/to/NZ_CP012345.1.fna.gz") [Optional, defaults to Accession ID]
     *   - Column 3: Organism Name (e.g., "Escherichia coli") [Optional]
     *   - Column 4: Taxonomy Name Lineage (e.g., "cellular organisms;Bacteria;...") [Optional]
     *   - Column 5: Taxonomy ID Lineage (e.g., "131567;2;1224;...") [Optional]
     * 
     * Example input line:
     * NZ_CP012345.1	562	ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz	Escherichia coli str. K-12 substr. MG1655	cellular organisms;Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli	131567;2;1224;1236;91347;543;561;562
     */
    std::vector<Species> parse_refseq_taxonomy_file(std::string const filepath)
    {
        std::vector<std::vector<std::string> > tax_file_lines{};
	    read_tsv(filepath, tax_file_lines);
	
		std::vector<Species> org_list{};
		org_list.reserve(tax_file_lines.size());
		for (auto & line : tax_file_lines)
		{
			Species sp{};
			sp.accession_id = std::move(line[0]);
			sp.taxid        = std::move(line[1]);

			sp.organism_name = "";
			sp.taxnames_string = "";
			sp.taxid_string = "";

			if (line.size() > 3)
				sp.organism_name = std::move(line[3]);
			if (line.size() > 4)
				sp.taxnames_string = std::move(line[4]);
			if (line.size() > 5)
				sp.taxid_string = std::move(line[5]);

			std::string stem_source = line.size() > 2 ? std::move(line[2]) : sp.accession_id;
			
			// Remove any trailing slashes to prevent extracting an empty string
			while (!stem_source.empty() && (stem_source.back() == '/' || stem_source.back() == '\\'))
				stem_source.pop_back();
				
            // KEEP ABSOLUTE PATH
			sp.file_stem = std::move(stem_source);

			// strip common extensions
			if (sp.file_stem.ends_with(".gz"))    sp.file_stem.resize(sp.file_stem.size() - 3);
			if (sp.file_stem.ends_with(".fna"))   sp.file_stem.resize(sp.file_stem.size() - 4);
			if (sp.file_stem.ends_with(".fasta")) sp.file_stem.resize(sp.file_stem.size() - 6);
			if (sp.file_stem.ends_with(".fa"))    sp.file_stem.resize(sp.file_stem.size() - 3);

			if (sp.file_stem.empty() || sp.file_stem == " ")
				throw std::runtime_error{"No file name found for " + sp.accession_id + " !!!"};

			org_list.emplace_back(std::move(sp));
		}
		return org_list;
    }
    
} // namespace taxor::taxonomy



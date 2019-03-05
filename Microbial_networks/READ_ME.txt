## Software ##
R (>= 3.5)

## Dependencies ##
To run this analysis, users need to install the following R packages:
- MetaboSignal (Bioconductor)
- KEGGREST (Bioconductor)
- taxize (CRAN)
- RCurl (CRAN)

## Instructions ##
Three R scripts are required to run this analysis:
* "Compilation_gut_microbial_species.R" within the subfolder ("Microbial_networks/Gut_microbial_species")
This script allows compiling a list of gut microbial species from published literature or from the HMP Catalog

* "01_Build_tryptophan_reaction_network_and_get_indole_genes_with_FASTA" within the subfolder ("Microbial_networks/Indole_metabolism/1_Build_tryptophan_network")
In this script, the MetaboSignal package is used to build a indole reaction network. The reactions of these network are then linked to the corresponding protein FASTA sequences

* "02_Analyse_BLAST_results" within the subfolder ("Microbial_networks/Indole_metabolism/2_Analyse_BLAST_results")
This script is used to analyse the BLAST outputs, and predict novel bacterial species involved in indole metabolism

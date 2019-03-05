################################################################################
############################ PRELIMINARIES #####################################
################################################################################

## Set working directory ##
Microbial_networks_path <- "M:/PhD_work/Microbial_networks" ## Change as required (based on the location of the folder "Microbial_networks")
setwd(Microbial_networks_path)

## Load packages ##
library(taxize)

## Import list of gut microbiota species (HMP + culturomics)
gut_microbiota <- read.csv("./Gut_microbial_species/all_gut_microbial_species/gut_microbiota.csv") ## Species from HMP and culturomics

## Import KEGG input ##
kegg_input <- read.csv("./Indole_metabolism/1_Build_tryptophan_network/FASTA_files/KEGG_FASTA.csv")
# Match kegg species included in gut microbiota dataset (by taxon, by short name, or by full name)
kegg_input$gut1 <- sapply(kegg_input$taxon, function(x) x %in% gut_microbiota$taxon)
kegg_input$gut2 <- sapply(kegg_input$organism_name_short, function(x) x %in% gut_microbiota$name)
kegg_input$gut3 <- sapply(kegg_input$organism_name, function(x) x %in% gut_microbiota$name)
kegg_input$gut <- kegg_input$gut1 + kegg_input$gut2 + kegg_input$gut3
kegg_input$gut <- ifelse(kegg_input$gut == 0, FALSE, TRUE) 

## Define variable names
names_blast_out <- c("sseqid", "qseqid", "pident", "length", "mismatch", "gapopen", 
                     "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                     "staxid", "ssciname","scomname","sblastname", "sskingdom", 
                     "staxids", "sscinames", "scomnames", "sblastnames","sskingdoms")
names_blast_filtered <- c("staxid", "ssciname","scomname","sblastname", "sskingdom")

## Get list of blast files
blast_files <- list.files("./Indole_metabolism/2_Analyse_blast_results/BLAST_output")


################################################################################
################################ ANALYSIS ######################################
################################################################################

## Loop to compare kegg input with blast output for each orthology gene ## 
res_file <- NULL
blast_filtered <- vector(mode = "list", length = length(blast_files))
names(blast_filtered) <- gsub(".csv", "", blast_files)

for (file in blast_files) {
    print(file)
    gene_name <- gsub(".csv", "", file)
    gene_subset <- subset(kegg_input, orthology_ID == gene_name)
    kegg_taxon <- as.character(gene_subset$taxon)
    kegg_name <- as.character(gene_subset$organism_name)
    kegg_name_short <- as.character(gene_subset$organism_name_short)
    
    ## Import BLAST output ##
    blast_out <- na.omit(read.csv(paste("./Indole_metabolism/2_Analyse_blast_results/BLAST_output/", file, sep = ""), header = FALSE))
    colnames(blast_out) <- names_blast_out
    
    ## Check whether all KEGG queries have been matched ##
    queries <- unique(sort(as.character(blast_out$qseqid)))
    ans_query <- sapply(queries, function(x) 
        max(as.numeric(as.character(subset(blast_out, qseqid == x)$pident))))
    
    message(paste("Number of queries in BLAST file:", length(queries)))
    message(paste("Number of queries in KEGG file:", length(unique(as.character(gene_subset$FASTA)))))
    
    if(!all(ans_query == 100)) { # This could be possible, but should be very close to 100% (>99.5). It could also be due to lines incorrectly read because of "," separator
        message("Some queries where not self-matched")
    }
    
    ## Filter BLAST output by identify >= 95 %
    blast_out <- blast_out[which(as.numeric(as.character(blast_out$pident)) > 95), ]
    # Notice that some lines have additional "," in the column sscinames and therefore will be cut; producing NAs.
    # However this is not a problem here, as we look only at columns which are before sscinames (names_blast_filtered)
    blast_out_filtered <- unique(subset(blast_out)[, c(names_blast_filtered)]) ## Unique organisms (after removing sequence ids, etc)
    blast_out_filtered$ssciname <- gsub("[[]", "", blast_out_filtered$ssciname)
    blast_out_filtered$ssciname <- gsub("[]]", "", blast_out_filtered$ssciname)
    blast_out_filtered$name_short <- unlist(lapply(strsplit(as.character(blast_out_filtered$ssciname), " "), 
                                                   function(x) paste(x[1], tolower(x[2]))))
    blast_out_filtered <- blast_out_filtered[!grepl("NA", blast_out_filtered$name_short), ]  # Remove multispecies (ssciname only has one word, therefore the second element of name_short is NA)
    
    ## Check if source organism of sequences is included in gut (match by organism name, organism name short, or by taxon id)
    blast_out_filtered$gut1 <- sapply(as.character(blast_out_filtered$staxid), 
                                      function (x) x %in% as.character(gut_microbiota$taxon))
    blast_out_filtered$gut2  <- sapply(as.character(blast_out_filtered$ssciname), 
                                       function (x) x %in% as.character(gut_microbiota$name))
    blast_out_filtered$gut3  <- sapply(as.character(blast_out_filtered$name_short), 
                                       function (x) x %in% as.character(gut_microbiota$name))
    blast_out_filtered$gut <- blast_out_filtered$gut1 + blast_out_filtered$gut2 + blast_out_filtered$gut3
    blast_out_filtered$gut <- ifelse(blast_out_filtered$gut == 0, FALSE, TRUE)
    
    ## Check if source organism of sequences is included in KEGG (match by organism name or by taxon id)
    blast_out_filtered$kegg1 <- sapply(as.character(blast_out_filtered$staxid), 
                                       function (x) x %in% kegg_taxon)
    blast_out_filtered$kegg2  <- sapply(as.character(blast_out_filtered$ssciname), 
                                        function (x) x %in% kegg_name)
    blast_out_filtered$kegg3  <- sapply(as.character(blast_out_filtered$name_short), 
                                       function (x) x %in% kegg_name_short)
    blast_out_filtered$kegg <- blast_out_filtered$kegg1 + blast_out_filtered$kegg2
    blast_out_filtered$kegg <- ifelse(blast_out_filtered$kegg == 0, FALSE, TRUE)
    
    blast_filtered[[gene_name]] <- cbind(blast_out_filtered, kegg_gene = gene_name)

    ## Compare species from kegg input and blast output at organism level (org) and species level (sp)
    all_kegg_org <- length(unique(kegg_taxon)) ## number of taxon_ids reported in KEGG for gene_name
    # Noticed that one taxon_id can be associated to several protein sequences (i.e. several genes)
    # Also, several taxon_ids can be associated to the same protein sequence. Therefore, all_kegg
    # will not necessarily match the number of queries in the KEGG input file (or BLAST output file)
    all_kegg_sp <- length(unique(kegg_name_short))
    
    all_kegg_org_gut <- length(intersect(kegg_taxon, as.character(subset(kegg_input, gut)$taxon))) ## number of taxon_ids reported in KEGG & gut
    all_kegg_sp_gut <- length(intersect(kegg_name_short, as.character(subset(kegg_input, gut)$organism_name_short)))
    
    in_blast_org <- sum(!blast_out_filtered$kegg) ## Total number of new hits
    in_blast_sp <- length(unique(subset(blast_out_filtered, !kegg3)$name_short))
    
    in_blast_org_gut <- sum(!blast_out_filtered$kegg & blast_out_filtered$gut) ## Total number of new hits included in gut
    in_blast_sp_gut <- length(unique(subset(blast_out_filtered, !kegg3 & gut)$name_short))
    
    res_file <- rbind(res_file, c(gene_name, all_kegg_org, all_kegg_sp, 
                                  in_blast_org, in_blast_sp,
                                  all_kegg_org_gut, all_kegg_sp_gut,
                                  in_blast_org_gut, in_blast_sp_gut))
}

colnames(res_file) <- c("kegg_gene", "N_organisms_kegg_input", "N_species_kegg_input", 
                        "N_organisms_new_blast", "N_species_new_blast",
                        "N_organisms_kegg_input_gut", "N_species_kegg_input_gut", 
                        "N_organisms_new_blast_gut", "N_species_new_blast_gut")

write.csv(res_file, "./Indole_metabolism/2_Analyse_blast_results/Post_BLAST_files/BLAST_hits.csv")

################################################################################
################################################################################

## Tryptophanase ## K01667
trp_blast <- unique(as.character(subset(blast_filtered[["K01667"]], gut)$name_short)) # Species level (i.e not sub-species)
trp_kegg <- unique(as.character(subset(kegg_input, orthology_ID == "K01667" & gut)$organism_name_short))
all_trp <- unique(c(trp_blast, trp_kegg))
phylum <- tax_name(all_trp, ask = FALSE, db = "ncbi", get = c("phylum"))$phylum

trp_final <- data.frame(organism = all_trp, phylum = phylum)
trp_final <- trp_final[order(trp_final$phylum), ]
trp_final$new_blast <- sapply(trp_final$organism, function(x) x %in% setdiff(trp_blast, trp_kegg))

write.csv(trp_final, "./Indole_metabolism/2_Analyse_blast_results/Post_BLAST_files/trypthophanase_network.csv")

################################################################################
################################################################################

## Production of I3A ## 
I3A_genes <- c("K00466", "K01593", "K01426", "K04103", "K00274", "K00128")
blast_filteredM <- as.data.frame(do.call(rbind, blast_filtered))
rownames(blast_filteredM) <- NULL
I3A_data <- subset(blast_filteredM, kegg_gene %in% I3A_genes)[, c("name_short", "gut", "kegg3", "kegg_gene")]
I3A_data <- unique(I3A_data)

I3A_hm <- matrix(0, nrow = length(unique(I3A_data$name_short)), ncol = length(I3A_genes))
colnames(I3A_hm) <- I3A_genes
rownames(I3A_hm) <- unique(as.character(I3A_data$name_short))

for(gene in colnames(I3A_hm)) {
    print(gene)
    my_subset <- subset(I3A_data, kegg_gene == gene)
    I3A_hm[as.character(my_subset$name_short), gene] <- 1
}
I3A_hm <- as.data.frame(I3A_hm)

## Split K00128 into K00128a and K00128b as this gene is involved into two different paths producing I3A
colnames(I3A_hm) <- gsub("K00128", "K00128a", colnames(I3A_hm))
I3A_hm$K00128b <- I3A_hm$K00128a

## Add sum of reactions and gut columns
I3A_hm$sum <- rowSums(I3A_hm)
I3A_hm$gut <- ifelse(rownames(I3A_hm) %in% as.character(subset(I3A_data, gut)$name_short), TRUE, FALSE)

## Score paths
path_1 <- c("K00466", "K01426") # score = 1/2
path_2 <- c("K04103", "K00128a") # score = 1/2
path_3 <- c("K01593", "K00274", "K00128b") # score = 1/3

ind_0 <- which(I3A_hm == 0, arr.ind = T) ## Keep track of 0 values in I3A_hm

I3A_hm[, path_1] <- 1/length(path_1)
I3A_hm[, path_2] <- 1/length(path_2)
I3A_hm[, path_3] <- 1/length(path_3)

I3A_hm[ind_0] <- 0 ## Restore the 0 values from the original data frame
I3A_hm$sum_score <- rowSums(I3A_hm[, grep("K", names(I3A_hm))])

## Select gut microbial species
I3A_hm_gut <- subset(I3A_hm, gut == 1)

write.csv(I3A_hm_gut, "./Indole_metabolism/2_Analyse_blast_results/Post_BLAST_files/I3A_heatmap.csv")



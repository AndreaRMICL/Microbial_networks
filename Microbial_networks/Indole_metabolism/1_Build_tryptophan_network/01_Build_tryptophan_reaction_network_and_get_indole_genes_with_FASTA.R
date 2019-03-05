################################################################################
############################ PRELIMINARIES #####################################
################################################################################

## Set working directory ##
Microbial_networks_path <- "M:/PhD_work/Microbial_networks" ## Change as required (based on the location of the folder "Microbial_networks")
setwd(Microbial_networks_path)

## Load packages ##
library(MetaboSignal)
library(taxize)
library(RCurl)
library(KEGGREST)

## Functions to be used in this script ## 
# kegg_fasta #
kegg_fasta <- function(gene) {
    print(gene)
    file <- paste("http://rest.kegg.jp/get/", gene, "/aaseq", sep = "")
    lines = try(readLines(file))
    if(grepl("Error", lines[1]) == FALSE) {
        fasta = paste(lines[-1], collapse = "")
    } else {
        fasta = "not_found"
    }
    return(fasta)
}

# convertTable #
convertTable <- function(res) {
    if (nchar(res) == 0) {
        print("no result")
        result <- NULL
    } else {
        rows <- strsplit(res, "\n")
        rows.len <- length(rows[[1]])
        result <- matrix(unlist(lapply(rows, strsplit, "\t")), nrow = rows.len,
                         byrow = TRUE)
    }
    return(result)
}

################################################################################
################################ ANALYSIS ######################################
################################################################################

## Build tryptophan reaction network ## 
tryp <- MS_reactionNetwork("rn00380")

## Export tryptophan reaction network as .txt file to be imported into Cytoscape ##
setwd("./Indole_metabolism/1_Build_tryptophan_network/Tryptophan_network")
network <- MS_exportCytoscape(tryp, organism_code = "hsa", file_name = "tryptophan_network") # warnings are expected
setwd(Microbial_networks_path) # reset working directory 

# Select indoles (end-products) of interest: 5-Methoxyindoleacetate (cpd:C05660), 
# Indoxyl (cpd:C05658), 3-Methylindolepyruvate (cpd:C05644), 
# Indolelactate (cpd:C02043), Indole-3-acetate (cpd:C00954), 3-Indoleglycolaldehyde (cpd:C03230)
# The entry of indolelactate seems to have been temporarily removed from KEGG (last checked on 03/03/19)
indole_cpd <- c("cpd:C05660", "cpd:C05658", "cpd:C05644", "cpd:C02043", "cpd:C00954", "cpd:C03230")
distances_indole <- MS_distances(tryp, mode = "out", target_metabolites = indole_cpd) # distances from all reactions to all compounds
distances_indole[distances_indole == "Inf"] <- 1000
indole_reactions <- rownames(distances_indole)[apply(distances_indole, 1, min) < 1000] # 36 reactions 
# Notice that: R01971 seems to have been temporarily removed from KEGG (last checked on 03/03/19)

## Link indole reactions to orthologies (ko) ##
file_reaction <- "http://rest.kegg.jp/link/reaction/ko"
response_reaction <- getURL(file_reaction)
reaction_table <- convertTable(response_reaction) ## all ko_rn KEGG pairs

reaction_ko <- NULL
for(reaction in indole_reactions) {
    print(reaction)
    ind <- which(reaction_table[, 2] == reaction) 
    if (length(ind) > 0) {
        ans <- cbind(rep(reaction, length(ind)), reaction_table[ind, 1])
        reaction_ko <- rbind(reaction_ko, ans)
    }
}
reaction_ko <- gsub("ko:", "", reaction_ko)
setdiff(indole_reactions, reaction_ko) ## list of reactions involved in indole production but not linked to an orthology

## Link orthologies to specific genes ## 
genesM_indoles <- NULL
for (i in 1:length(reaction_ko[, 2])) {
    print(reaction_ko[i, 2])
    genes <- keggGet(reaction_ko[i, 2])[[1]]$GENES
    genesM <- do.call(rbind, strsplit(genes, ": "))
    genesM[, 1] <- tolower(genesM[, 1])
    ans <- cbind(reaction_ko[i, 1], reaction_ko[i, 2], genesM)
    for (z in 1:nrow(ans)) {
        row <- ans[z, ]
        all_entrez <- unlist(strsplit(row[4], " "))
        if(length(all_entrez) > 1) {
            #print(row)
        }
        ans2 <- cbind(row[1], row[2], row[3], all_entrez)
        genesM_indoles <- rbind(genesM_indoles, ans2)
    }
}
colnames(genesM_indoles) <- c("reaction_ID", "orthology_ID", "organism_code", "gene_ID")

## Get list of prokaryotic organisms from KEGG
org_prok <- as.matrix(read.csv("./Indole_metabolism/1_Build_tryptophan_network/KEGG_bacteria/bacteria_KEGG.csv", 
                               header = TRUE, row.names = 1)) ## List of prokaryotic organisms 
short_name <- unlist(lapply(strsplit(org_prok[, 1], " "), function(x) paste(x[1], tolower(x[2]))))
#org_prok <- org_prok[!grepl("sp.", short_name), ] # remove the ones where name is just one word + .sp

## Link prokaryotic organisms from KEGG to taxon (takes a while)
org_aa <- rep(NA, nrow(org_prok))
names(org_aa) <- rownames(org_prok)
org_tax <- org_aa
for (org in row.names(org_prok)) {
    print(org)
    lines <- readLines(paste("https://www.genome.jp/kegg-bin/show_organism?org=", org, sep = ""))
    info_aa <-  strsplit(strsplit(strsplit(strsplit(lines[62], "Assembly")[[1]][2], 
                                           "<br/>")[[1]][1], "[\">]")[[1]][4], "<")[[1]][1]
    info_tax <- strsplit(strsplit(lines[60], "Info&id=")[[1]][2], "[\">]")[[1]][1]
    org_aa[org] <- info_aa
    org_tax[org] <- info_tax
}

org_prok <- cbind(org_prok, Taxon = org_tax, Assembly = org_aa)
write.csv(org_prok, "./Indole_metabolism/1_Build_tryptophan_network/KEGG_bacteria/KEGG_prokaryotic_organisms.csv")

## Filter by prokaryotic organisms ## 
genes_prokaryote <- genesM_indoles[genesM_indoles[, 3] %in% rownames(org_prok), ]
table(genes_prokaryote[, "orthology_ID"]) ## see prokaryotic orthologies involved in indole production
table(genes_prokaryote[, "reaction_ID"]) ## see prokaryotic reactions involved in indole production

## Clean gene IDs (remove "()", eg. EcWSU1_00526(betB))
gene_clean <- lapply(strsplit(genes_prokaryote[, "gene_ID"], "[(]"), function(x) x[1])
genes_prokaryote[, "gene_ID"] <- paste(genes_prokaryote[, "organism_code"], 
                                       unlist(gene_clean), sep = ":")

## Get organism name and phylum
org_name <- org_prok[genes_prokaryote[, "organism_code"], 1]
org_name <- gsub("subsp. ", "", org_name)
short_org_name <- unlist(lapply(strsplit(org_name, " "), 
                                function(x) paste(x[1], tolower(x[2]))))
phylum <- tax_name(unique(short_org_name), ask = FALSE, db = "ncbi", get = c("phylum"))
rownames(phylum) <- phylum$query

## Get FASTA (protein) for each gene ##
fasta_seq <- sapply(genes_prokaryote[, "gene_ID"], kegg_fasta)
genes_fasta <- cbind(genes_prokaryote[, 1:3], 
                     organism_name = org_prok[genes_prokaryote[, "organism_code"], "Organism.name"],
                     organism_name_short = short_org_name, 
                     assembly_accesion = org_prok[genes_prokaryote[, "organism_code"], "Assembly"], 
                     taxon = org_prok[genes_prokaryote[, "organism_code"], "Taxon"], 
                     organism_phylum = phylum[short_org_name, "phylum"], 
                     gene_ID = genes_prokaryote[, 4], 
                     FASTA = fasta_seq)

write.csv(genes_fasta, "./Indole_metabolism/1_Build_tryptophan_network/FASTA_files/KEGG_FASTA.csv")

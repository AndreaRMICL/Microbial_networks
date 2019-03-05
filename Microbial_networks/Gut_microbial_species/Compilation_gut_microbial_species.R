library(taxize)
setwd("M:/PhD_work/Microbial_networks") # Change as required, based on the location of the folder "Microbial_networks"

## HMP ##
hmp <- read.csv("./Gut_microbial_Species/HMP/HMRGD.csv")[, 1:2]
hmp_gut <- subset(hmp, Body.Site == "Gastrointestinal_tract")
hmp_gut_sp <- as.character(hmp_gut[, 1])
hmp_gut_sp <- hmp_gut_sp[nchar(hmp_gut_sp) > 1]
hmp_gut_sp <- gsub("[[]", "", hmp_gut_sp)
hmp_gut_sp <- gsub("[]]", "", hmp_gut_sp)
hmp_gut_sp[grepl("Bifidobacterium longum subsp. infantis ATCC 15697", hmp_gut_sp)] <- 
    "Bifidobacterium longum subsp. infantis ATCC 15697"
hmp_gut_sp[grepl("Lactobacillus plantarum subsp. plantarum ATCC 14917", hmp_gut_sp)] <- 
    "Lactobacillus plantarum subsp. plantarum ATCC 14917"
hmp_gut_sp <- unique(hmp_gut_sp[!grepl("Gastrointestinal_tract", hmp_gut_sp)])
tax_hmp <- get_ids(hmp_gut_sp, db = "ncbi")
pl_hmp <- tax_name(hmp_gut_sp, ask = FALSE, db = "ncbi", get = c("phylum"))
hmp_gut_tax <- cbind(name = hmp_gut_sp, taxon = as.character(tax_hmp$ncbi), 
                     phylum = pl_hmp[, 3], source = "HMP_2018")

## Browne ##
br <- as.matrix(read.csv("./Gut_microbial_Species/Browne/all_bacteria_Browne.csv"))
br_l <- strsplit(br[, 2], "[[]")
br_length <- sapply(br_l, function(x) length(x))
table(br_length)
br_1 <- br_l[br_length == 2] # take last one
br_2 <- br_l[br_length == 3] # take last one
br_3 <- br_l[br_length == 4] # take last one

br_sp <- unique(unlist(lapply(strsplit(br[, 2], "[[]"), function(x) x[length(x)])))
br_sp <- gsub("]", "",br_sp)
br_sp <- trimws(br_sp)
br_sp <- br_sp[!grepl("uncultured", br_sp)]
br_clean <- unique(unlist(lapply(strsplit(br_sp, " strain"), function(x) x[1])))
tax_br <- get_ids(br_clean, db = "ncbi")
pl_br<- tax_name(br_clean, ask = FALSE, db = "ncbi", get = c("phylum"))
br_tax <- cbind(name = br_clean, taxon = as.character(tax_br$ncbi), phylum = pl_br[, 3], source = "Browne_2016")

## Lagier ## 
lag <- as.matrix(read.csv("./Gut_microbial_Species/Lagier/all_bacteria_Lagier.csv"))
lag <- gsub("Zimmermanella faecalis", "Zimmermannella faecalis", lag)
tax_lag <- get_ids(lag[, 1], db = "ncbi")
lag_tax <- cbind(name = lag[, 1], taxon = as.character(tax_lag$ncbi), phylum = lag[, "Phylum"], source = "Lagier_2016")

## Rajilic-Stojanovic ## (Review)
rs <- as.matrix(read.csv("./Gut_microbial_Species/Rajilic-Stojanovic/gut_microbial_species_24861948.csv"))
tax_rs <- get_ids(rs[, "Species"], db = "ncbi")
rs_tax <- cbind(name = rs[, "Species"], taxon = as.character(tax_rs$ncbi), phylum = rs[, "Phylum"], source = "Rajilic-Stojanovic_2014")

## Merge datasets ##
gut_microbiota <- rbind(hmp_gut_tax, br_tax, lag_tax, rs_tax)
write.csv(gut_microbiota, "./Gut_microbial_Species/all_gut_microbial_species/gut_microbiota.csv", 
          row.names = FALSE)

setdiff(br_tax[, "taxon"], lag_tax[, "taxon"]) # most of Browne are in Lagier
length(setdiff(hmp_gut_tax[, "taxon"], lag_tax[, "taxon"])) # most HMP are not in Lagier
sum(is.na(lag_tax[, "taxon"])) ## 95 missing values
length(unique(lag_tax[, "taxon"]))


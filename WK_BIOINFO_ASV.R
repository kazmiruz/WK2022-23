library(tidyverse)

##############################################
##### Create an ASV table for William Kennedy samples

#
### Load the general ASV tables for 16S and 18S

## The tables are samples (rows) vs 16S/18S sequences (columns)
# These are raw: for many different samples and not annotated yet 
seqtab.16S = readRDS("WK_DNA/515F-806RB/seqtab_final.rds")
rownames(seqtab.16S) = gsub(pattern = "716S",replacement = "7-16S",x = rownames(seqtab.16S)) # fix 716S to be 7-16S in rownames

seqtab.18S = readRDS("WK_DNA/Euk1391F-EukBr/seqtab_final.rds")

## Explore using subset b/c the files are large
# subset <- seqtab.18S[, 1:10]
# rownames(seqtab.18S)

## Split rownames into columns, make a plate.geno with all upper cases
sepnames.16S <- rownames(seqtab.16S) %>%
  as_tibble() %>%
  separate(value, into = c("project", "primer", "plate", "well"), sep = "\\|", fill = "right") %>%
  separate(plate, into = c("p1", "p2", "p3", "p4", "p5", "p6"), sep = "-", fill = "right") %>%
  mutate(plate.geno = toupper(paste(p2, well, sep = "_")))

sepnames.18S <- rownames(seqtab.18S) %>%
  as_tibble() %>%
  separate(value, into = c("project", "primer", "plate", "well"), sep = "\\|", fill = "right") %>%
  separate(plate, into = c("p1", "p2", "p3", "p4", "p5", "p6"), sep = "-", fill = "right") %>%
  mutate(plate.geno = toupper(paste(p2, well, sep = "_")))

#
### Find William Kennedy sample names (pattern WK22 and WK23) and combine them

## Identify William Kennedy samples and combine them
sam.wk.16S = grep(pattern = "WK22|WK23",x = sepnames.16S$plate.geno, value=T)
sam.wk.18S = grep(pattern = "WK22|WK23",x = sepnames.18S$plate.geno, value=T)
sam.wk.both = intersect(sam.wk.16S, sam.wk.18S)

## Locate where the William Kennedy sample are in the sepnames table, and hence the main ASV table - aka order ID
sam.16S.pick = match(sam.wk.both,sepnames.16S$plate.geno)                 # identify the location in the sepnames table
sam.18S.pick = match(sam.wk.both,sepnames.18S$plate.geno)

## Create ASV table for 16S/18S William Kennedy samples
proptab.small <- cbind(prop.table(seqtab.16S[sam.16S.pick,],margin=1),      # make a proportional table by rows (margin=1) for 16S
                      prop.table(seqtab.18S[sam.18S.pick,],margin=1)) %>%   # do the same for 18S and combine the columns
  replace(is.na(.), 0) %>%                                                  # set NA to 0
  .[rowSums(. > 0) > 2,                                                     # select samples with more than 2 species
    colSums(. > 0) > 2] %>%                                                 # select species present in more than 2 samples
  .[, colMeans(.) > 0.0001]                                                 # remove rare species (<= 0.01% sample relab)


##########################
##### Taxonomic annotation

#
### Load taxa

# load cleaned up taxa
# see WK_BIOINFO_TAXA to see how the data was processed
tax.16S <- readRDS("OUTPUT/tax.16S.rds")
tax.18S <- readRDS("OUTPUT/tax.18S.rds")

# explore
subset <- tax.16S[1:1000,]
colnames(tax.16S)

#
### Rename sequences to taxa

## Identify which sequences belong to 16S and 18S
sequence_names <- colnames(proptab.small)

## Create a lookup table with original sequences and their taxonomy names
tax_lookup <- tibble(sequence = sequence_names) %>%
  mutate(
    taxon_raw = case_when(
      sequence %in% rownames(tax.16S) ~ paste(tax.16S[sequence, "Phylum"], 
                                              tax.16S[sequence, "Order"], 
                                              tax.16S[sequence, "Genus"], sep = "; "),
      sequence %in% rownames(tax.18S) ~ paste(tax.18S[sequence, "Division"], 
                                              tax.18S[sequence, "Order"], 
                                              tax.18S[sequence, "Genus"], sep = "; "),
      TRUE ~ "Unknown")) %>%                                                   # Handle cases where a sequence is not found in tax tables
  rowwise() %>%  # Process each row independently
  mutate(                                                                      # Process to remove redundant taxonomy names
    taxon_name = {
      tax_levels <- unlist(strsplit(taxon_raw, ";\\s*"))                       # Split by "; " and clean spaces
      unique_taxon <- unique(tax_levels)                                       # Remove redundant levels
      paste(unique_taxon, collapse = "; ")}) %>%                               # Recombine into final taxonomy string
  ungroup() %>%                                                                # Exit rowwise mode to ensure proper numbering
  mutate(taxon_name = paste0(taxon_name, "; ASV_", row_number()))              # Correct ASV numbering


## Create a taxonomy-annotated table using the lookup table 
proptab.small.taxan <- proptab.small %>%
  `colnames<-`(tax_lookup$taxon_name)


####################################
##### Sample annotation and matching

#
### Shorten rownames to plate.geno

# Convert row names into a tibble and extract plate.geno
sepnames.small <- rownames(proptab.small.taxan) %>%
  tibble(value = .) %>%
  separate(value, into = c("project", "primer", "plate", "well"), sep = "\\|", fill = "right") %>%
  separate(plate, into = c("p1", "p2", "p3", "p4", "p5", "p6"), sep = "-", fill = "right") %>%
  mutate(plate.geno = toupper(paste(p2, well, sep = "_"))) %>%
  select(plate.geno)  # Keep only new row names

# Assign sepnames.small as new row names (the order is preserved)
proptab.small.taxan.shortrow <- proptab.small.taxan %>%
  as.data.frame() %>%
  `rownames<-`(sepnames.small$plate.geno)

#
### Load metadata specific to the William Kennedy samples analysed and reorder to match the ASV table (72 of both 16S and 18S minus 1 sample)
meta.wk.small <- read.csv("OUTPUT/meta.wk.complete.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
  filter(plate.geno %in% sepnames.small$plate.geno) %>%                                # Keep only matching rows
  arrange(match(plate.geno, sepnames.small$plate.geno))                                # Reorder to match sepnames.small


#
### Build sample names for the William Kennedy samples

## Extract the sample name part that is informed by the metadata
wk_names <- meta.wk.small %>%
  mutate(Coast = paste0(round(`Distance to coast [km]`), " km"),             # Distance to coast, km
         Depth = paste0(round(`Depth [m]`), " m")) %>%                       # Add " m" to Depth
  transmute(plate.geno,                                                      # Keep plate.geno for matching
            new_colname = paste(Cruise, Station, Coast, Depth, sep = "; "))  # Construct the new column name

#
### Assign "wk_names; sepnames.small" as new row names (the order is still preserved)
# Thus, both meta and genohub info is part of the sample names
proptab.small.taxan.saman <- proptab.small.taxan %>%
  as.data.frame() %>%
  `rownames<-`(sepnames.small %>%
                 left_join(wk_names, by = "plate.geno") %>%
                 mutate(final_rowname = paste(new_colname, plate.geno, sep = "; ")) %>%
                 pull(final_rowname))  # Extract combined row names and assign directly


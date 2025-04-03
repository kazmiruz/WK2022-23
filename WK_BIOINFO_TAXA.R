
##
## 16S
# tax.16S is sequences vs Domain/Phylum/Class/Order/Family/Genus
# boot.16S is bootstrap confidence scores for the tax.16S table
# bootstrapping is a resampling technique used to estimate confidence intervals / variability by repeatedly sampling with replacement
tax.16S = data.frame(readRDS("WK_DNA/515F-806RB/taxout.rds")$`ref_dada2_silva_v138-2_train_set_uniq.fasta.gz`)
boot.16S = data.frame(readRDS("WK_DNA/515F-806RB/bootout.rds")$`ref_dada2_silva_v138-2_train_set_uniq.fasta.gz`)

# Explore using subset b/c the files are large
# subset <- tax.16S[1:100,]
# colnames(tax.16S)
# subset <- boot.16S[1:10,]
# colnames(tax.16S)

# Ensure that both tables contain matching entries and have the same row order
boot.16S = boot.16S[rownames(tax.16S),]

# Fix taxonomic names
for(i in 1:nrow(tax.16S)) {                     # Iterate over each ASV (row)
  last_known = NA                               # A fall back in case all classifications are unreliable
  last_known_level = NA
  
  for(j in 1:ncol(tax.16S)) {                   # Iterate through Domain to Genus
    if(is.na(tax.16S[i, j]) |
       tax.16S[i, j] == "Incertae Sedis" |
       boot.16S[i, j] < 80) {                   # Adjust depending on how high of a confidence you want
      
      tax.16S[i, j] = NA}                       # Assign NAs if unknown species or low confidence
    
    else {
      last_known = tax.16S[i, j]                # If non-NA, store the classification
      # Gets replaced if the next column also provides a non-NA - hence, a more precise taxa
      last_known_level = colnames(tax.16S)[j]}  # Store its taxonomic level (column name)
  }
  
  if(all(is.na(tax.16S[i, ]))) {                # If the entire row is NA (even Domain is unknown), assign "Unknown" to all cells
    tax.16S[i, ] <- "Unknown"}
  else {
    for(j in 1:ncol(tax.16S)) {                 # Otherwise, fill all NA values in the row with "Best_knows_lvl Best_known_name"
      if(is.na(tax.16S[i, j])) {
        tax.16S[i, j] = paste0(last_known_level, " ", last_known)}}
  }
}

saveRDS(tax.16S, "OUTPUT/tax.16S.rds")

###################

##
## 18S - do the same as for 16S, notice how there are more taxonomic classifications (columns)
tax.18S = data.frame(readRDS("WK_DNA/Euk1391F-EukBr/taxout.rds")$ref_dada2_pr2_version_5.0.0_SSU.fasta.gz)
boot.18S = data.frame(readRDS("WK_DNA/Euk1391F-EukBr/bootout.rds")$ref_dada2_pr2_version_5.0.0_SSU.fasta.gz)

# Explore using subset b/c the files are large
# subset <- tax.18S[1:1000,]
# colnames(tax.18S)
# subset <- boot.18S[1:10,]
# colnames(tax.18S)

# Ensure that both tables contain matching entries and have the same row order
boot.18S = boot.18S[rownames(tax.18S),]

# Fix taxonomic names
for(i in 1:nrow(tax.18S)) {                     # Iterate over each ASV (row)
  last_known = NA                               # A fall back in case all classifications are unreliable
  last_known_level = NA
  
  for(j in 1:ncol(tax.18S)) {                   # Iterate through Domain to Genus
    if(is.na(tax.18S[i, j]) |
       tax.18S[i, j] == "Incertae Sedis" |
       boot.18S[i, j] < 80) {                   # Adjust depending on how high of a confidence you want
      
      tax.18S[i, j] = NA}                       # Assign NAs if unknown species or low confidence
    
    else {
      last_known = tax.18S[i, j]                # If non-NA, store the classification
      # Gets replaced if the next column also provides a non-NA - hence, a more precise taxa
      last_known_level = colnames(tax.18S)[j]}  # Store its taxonomic level (column name)
  }
  
  if(all(is.na(tax.18S[i, ]))) {                # If the entire row is NA (even Domain is unknown), assign "Unknown" to all cells
    tax.18S[i, ] <- "Unknown"}
  else {
    for(j in 1:ncol(tax.18S)) {                 # Otherwise, fill all NA values in the row with "Best_knows_lvl Best_known_name"
      if(is.na(tax.18S[i, j])) {
        tax.18S[i, j] = paste0(last_known_level, " ", last_known)}}
  }
}

saveRDS(tax.18S, "OUTPUT/tax.18S.rds")


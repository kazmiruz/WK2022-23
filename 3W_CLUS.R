# 3------------------------------------------------------------------------ [can be run as is w/out running 1W_META and 2W_TAX first]
# Purpose:
#   Cluster samples based on taxonomic composition (RelAb of non-rare ASVs), and identify bioindicators of clusters and cosmopolitans.
#   Specifically, assess cluster stability (silhouette / bootstrap Jaccard), build a Ward‑D2 dendrogram (rotated by latitude),
#   perform indicator species analysis (IndVal.g), select bioindicators and cosmopolitans, produce a ComplexHeatmap and a pure dendrogram,
#   run TITAN for community‑level change points along salinity and latitude, and check for flowthrough contamination.
#
# Required external files:
#   EXT/meta.wk.rds                         (cleaned up metadata; can also produce by running 1W_META and 2W_TAX)
#   EXT/proptab.wk.rds                      (cleaned up ASV table; can also produce by running 1W_META and 2W_TAX)
#   EXT/indicator_res.rds                   (can reproduce using commented code)
#   EXT/titan_offshore_surface_salinity.rds (can reproduce using commented code)
#   EXT/titan_offshore_surface_latitude.rds (can reproduce using commented code)
#
# Outputs:
#   OUT/3W_CLUS/wk_simple_hm.pdf                           (exploration: diagnostic heatmap to check the transformation and clustering)
#   OUT/3W_CLUS/wk_complex_hm.pdf / .png                   (Figure S1: annotated ComplexHeatmap)
#   OUT/3W_CLUS/dendro_pure.pdf / .png                     (cluster dendrogram shown alongside the maps in Figure 2)
#   OUT/3W_CLUS/titan_offshore_surface_latitude.pdf / .csv (TITAN for latitude: change points in .csv reported in Section 3.3; .pdf for exploration)
#   OUT/3W_CLUS/titan_offshore_surface_salinity.pdf / .csv (TITAN for salinity: change points in .csv reported in Section 3.3; .pdf for exploration)
# -------------------------------------------------------------------------

library(tidyverse)
library(gplots)                          # diagnostic heatmap.2
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)                  # annotated heatmap and pure dendrogram
library(patchwork)                       # combine TITAN panels
library(viridis)                         # heatmap colours
library(cluster)                         # silhouette widths
library(dendextend)                      # rotate and colour dendrograms
library(vegan)                           # Bray-Curtis dissimilarity
library(fpc)                             # bootstrap cluster stability
library(indicspecies)                    # indicator taxa analysis
library(TITAN2)                          # threshold indicator taxa analysis
options("lifecycle_verbosity"="warning") # always get warnings if there are (not just once every 8 hours)
# -------------------------------------------------------------------------

### Things to keep in mind before before diving into the analysis here and in the following scripts ###

#   Since microbes are the main drivers of biogeochemical cycles, community composition data can directly explain
#   and predict environmental shifts better than biochemical parameters alone. This mirrors medical diagnostics:
#   doctors would rather test the patient to identify the parasite than assign treatment based on symptoms alone.
#   Moreover, microbes act as sensitive environmental sensors - much like clinical biomarkers for early cancer detection.
#   Finally, microbes can be used as tracers of water masses, complementing traditional parameters like T and S.

# - EXPLORE the data with simple figures and tables before deciding which parameters contribute to the story.
# - Prioritize NOVEL or UNEXPECTED patterns over well-established differences.
# - Focus on findings with IMPLICATIONS for researchers, modelers, and northern communities.
# - Use each figure to COMMUNICATE a clear message and support the overall manuscript narrative.
# -------------------------------------------------------------------------

# upload cleaned up metadata and ASV table for WK samples
meta.wk <- readRDS("EXT/meta.wk.rds")
proptab.wk <- readRDS("EXT/proptab.wk.rds")

#
### What does your ASV table look like in terms of data distribution? What transformation should you apply?
# Apply one of these transformation for optimal clustering/grouping, adjust based on visual inspection of histograms and heatmaps.
# Fourth root (^.25) and square root (^.5) transformations are more common and easier to defend.

hist(as.matrix(proptab.wk))        # looks pretty bad
hist(as.matrix(proptab.wk)^.5)     # looks better
hist(as.matrix(proptab.wk)^.25)    # looks good: a lot of 0s but normal distribution to the left
hist(as.matrix(proptab.wk)^.125)   # looks good - use if need more pronounced clusters

proptab.trans <- as.matrix(proptab.wk)^.25  # transformed ASV table for grouping: row/col order is the same


############################################
### Select the number of sample clusters ###
# Use silhouette analysis for your Bray-Curtis dissimilarity; don't use WCSS as it is for Euclidean distance.
# Silhouette width measures cluster separation, whereas bootstrap Jaccard measures cluster stability.
# Among solutions with all cluster-specific Jaccard values >0.6, select the k with the highest mean silhouette width.

eval_k <- function(x, ks = 5:9, B = 1000) {
  set.seed(161193)
  d <- vegdist(x, method = "bray")
  tree <- hclust(d, method = "ward.D2")
  
  purrr::map_dfr(ks, \(k) {
    groups <- cutree(tree, k = k)
    sil <- cluster::silhouette(groups, d)[, 3]
    jacc <- fpc::clusterboot(d, B = B, bootmethod = "boot", clustermethod = fpc::disthclustCBI,
                             method = "ward.D2", k = k, count = FALSE)$bootmean
    
    tibble::tibble(k = k, sil_mean = mean(sil), sil_min = min(tapply(sil, groups, mean)),
                   jacc_mean = mean(jacc), jacc_min = min(jacc))
  })
}

res <- eval_k(proptab.trans)

ggplot(res, aes(k)) +
  geom_line(aes(y = sil_mean, colour = "Mean silhouette")) +
  geom_point(aes(y = sil_mean, colour = "Mean silhouette")) +
  geom_line(aes(y = jacc_min, colour = "Minimum bootstrap Jaccard")) +
  geom_point(aes(y = jacc_min, colour = "Minimum bootstrap Jaccard")) +
  labs(x = "Number of clusters", y = "Score", colour = NULL)

res %>% filter(jacc_min > 0.6) %>% slice_max(sil_mean, n = 1, with_ties = FALSE) %>% pull(k)  # 9


#########################
### Build dendrograms ###
#########################

k <- 9
cl_labels <- c("oHB", "BI22", "cHB", "dBI", "dJB", "EST", "MR", "nJB", "sJB")

# Sample dendrogram
dist_row <- vegdist(proptab.trans, method = "bray")
tree_row <- hclust(dist_row, method = "ward.D2")
dend_row <- as.dendrogram(tree_row)

# ASV dendrogram
dist_col <- vegdist(t(proptab.trans), method = "bray")
tree_col <- hclust(dist_col, method = "ward.D2")
dend_col <- as.dendrogram(tree_col)


##############################################################
### Rotate sample branches from highest to lowest latitude ###
##############################################################

# # Default
# dend_row_order <- dend_row

# # reorder() rotation
# reorder() recursively orders sister branches at each dendrogram split.
# With agglo.FUN = mean, it compares the mean latitude of the samples beneath each child branch.
# The negative sign places branches with higher mean latitude earlier in the leaf order.
# This produces a north-to-south tendency based on actual latitude values, but does not guarantee a monotonic leaf order.
# Only entire sister subtrees are flipped; the dendrogram topology, branch heights and cluster memberships are unchanged.
# dend_row_order <- reorder(dend_row, wts=(-meta.wk$`Latitude [deg]`), agglo.FUN = mean)

# rotate() rotation
# First, list the samples in the desired north-to-south order.
# rotate() converts that sequence into rank weights and uses them to orient sister branches, so that
# the leaf order resembles the requested sequence as closely as the existing topology permits. 
# Like reorder(), it does not change branch heights, clustering topology, or cluster membership.
# While reorder() places sister branches according to their average latitude,
# rotate() directly attempts to reproduce the requested order, which more closely matches my stated visual objective.
ord_row <- rownames(meta.wk)[order(meta.wk$`Latitude [deg]`, decreasing = TRUE)] # sample names ordered from high to low latitude
dend_row_order <- dendextend::rotate(dend_row, ord_row)                # make the dendrogram appear in a SIMILAR order by “rotating” branches

# Visual top-to-bottom sample order after rotation
row_ordered <- labels(dend_row_order)


########################################
### Label and colour sample clusters ###
########################################

# cutree() cluster numbers are arbitrary; relabel them 1:k by their first top-to-bottom appearance in the rotated dendrogram
clus <- cutree(tree_row, k = k)[row_ordered] %>%
  (\(x) factor(x, levels = unique(x), labels = as.character(seq_len(k))))()

cluster_palette <- colorRampPalette(c("darkviolet", "magenta2", "darkgoldenrod", "red2", "coral",
                                      "green2", "cyan3", "blue", "dodgerblue2"))(k)
names(cluster_palette) <- levels(clus)

dend_row_order <- dend_row_order %>% color_branches(k = k, col = cluster_palette) %>% set("branches_lwd", 3)
meta.wk$Cluster <- clus[rownames(meta.wk)]


#####################################################
### Rotate ASV branches by weighted mean latitude ###
#####################################################

# ASVs can be ordered based on the mean latitude among the samples where they occur.
# taxon_latitude_pa <- apply(proptab.trans > 0, 2, \(x) mean(meta.wk$`Latitude [deg]`[x]))

# However, it is better to order ASVs based on the mean sample latitude weighted by fourth-root-transformed ASV abundance.
# Samples in which the ASV is more abundant contribute more strongly to its display-order latitude.
# apply(proptab.trans, 2, …) - take each column (taxon's RelAb in different samples) and transfer to the function as x
# weighted.mean(y, w=x) computes sum(y*w)/sum(w) aka multiply Lat of each sample by the RelAb of the taxon in that sample ->
# sum it -> divide  by the sum of RelAb of the taxon
taxon_latitude <- apply(proptab.trans, 2, \(x) weighted.mean(meta.wk$`Latitude [deg]`, w = x))

# Now, let's proceed with rotating branches to follow the desired order.

# Default
# dend_col_order <- dend_col

# reorder() rotation: based on mean weighted latitude within sister subtrees
# dend_col_order <- reorder(dend_col, wts = -taxon_latitude, agglo.FUN = mean)

# rotate() rotation: toward the complete high-to-low weighted-latitude order
ord_col <- names(sort(taxon_latitude, decreasing = TRUE))
dend_col_order <- dendextend::rotate(dend_col, ord_col)
col_ordered <- colnames(proptab.trans)[order.dendrogram(dend_col_order)]

# Optional: cut the ASV dendrogram into labelled communities for column splitting and decoration.
# tax_k <- 10
# com_labels <- LETTERS[seq_len(tax_k)]
# com <- cutree(tree_col, k = tax_k)[col_ordered] %>%
#   (\(x) factor(x, levels = unique(x), labels = com_labels))()
# community <- data.frame(ASV = names(com), Community = as.character(com))
# col_slice_labels <- paste0("Network ", levels(com))
# dend_col_order <- dend_col_order %>% set("branches_lwd", 3)


##########################
### Diagnostic heatmap ###
##########################
# Exploration: check the transformation and clustering before the final heatmap
pdf("OUT/3W_CLUS/wk_simple_hm.pdf", width = 35, height = 20)
heatmap.2(proptab.trans, Rowv = dend_row_order, Colv = dend_col_order, trace = "none", key = FALSE,
          margins = c(20, 30), col = viridis(256),
          RowSideColors = unname(cluster_palette[as.character(clus[rownames(proptab.trans)])]))
dev.off()


###############################
### Indicator taxa analysis ###
###############################

# Use meta.wk$Cluster because its sample names align with proptab.wk
# Run on untransformed relative abundances, allow single clusters and cluster pairs, and save the result for later sessions.
set.seed(161193)
# indicator_res <- multipatt(proptab.wk, meta.wk$Cluster, func = "IndVal.g",
#                            duleg = FALSE, # If TRUE, Cluster combos aren't considered, only the original Clusters (1 to k, Dufrêne & Legendre, 1997)
#                            max.order = 2, control = how(nperm = 9999))
# write_rds(indicator_res, "EXT/indicator_res.rds")       # 9 clusters w/ pairs
indicator_res <- read_rds("EXT/indicator_res.rds")

# # Explore the indicspecies output to understand the concept
# indicator_res$cluster   # sample clusters
# indicator_res$comb      # candidate single-cluster and cluster-pair patterns
#
# # IndVal = A*B, where "A" indicates how strongly the ASV is assoc-d with the target Cluster
# # A = "specificity" = [Mean RelAb in the target Cluster samples]/[SUM of Mean RelAbs in ALL samples]
# # B = "fidelity" = [# of samples in the target Cluster where the ASV exists]/[overall # of samples in the target Cluster]
# hist(indicator_res$A)   # "A" component of IndVal (= A * B) of each ASV for each cluster or cluster combination: "specificity" of ASV
# hist(indicator_res$B)
#
# # $sign gives some important values: s.1, ... - group preferred by each ASV; index - cluster comb-n; stat - sq.root(IndVal); p.value - sign-ce of stat
# indicator_res$sign
# plot(colMeans(proptab.wk), indicator_res$sign$stat, log = 'xy') # for each ASV, plot Mean RelAb (cols = ASVs) vs stat (same order of ASVs)
#
# # The permutation test shuffles the group labels 9,999 times across samples while holding ASV abundances fixed. It tests if the microbe-group link is real.
# # A p-value <0.01 means that there is <1% chance that the link happened by accident. 
# hist(indicator_res$sign$p.value, breaks = 100)                                 # most are <0.1
#
# The statistics we are interested in are A, B, and p.value
# summary(indicator_res, alpha = 0.01, At = 0.6, Bt = 0.6, indvalcomp = TRUE)


###########################################
### Indicator selection + cosmopolitans ###
###########################################

rn <- rownames(indicator_res$A); cn <- colnames(indicator_res$A); idx <- indicator_res$sign$index
valid <- !is.na(idx) & !is.na(indicator_res$sign$p.value) & indicator_res$sign$p.value < 0.01 &
  indicator_res$A[cbind(seq_along(rn), idx)] >= 0.6 & indicator_res$B[cbind(seq_along(rn), idx)] >= 0.6

specific_asvs <- rn[valid]
chosen_group <- setNames(rep(NA_character_, length(rn)), rn)
chosen_group[valid] <- cn[idx[valid]]

# For pair-associated indicators, assign the ASV for display and interpretation to the member cluster with the highest mean relative abundance.
single_cols <- cn[!str_detect(cn, "[-+]")]
sample_cluster <- as.character(meta.wk$Cluster[match(rownames(proptab.wk), rownames(meta.wk))])
mean_rel_ab <- sapply(single_cols, \(g) colMeans(proptab.wk[sample_cluster == g, , drop = FALSE], na.rm = TRUE))

chosen_cluster <- map2_chr(chosen_group, rn, \(group, asv) {
  if (is.na(group) || !str_detect(group, "[-+]")) return(group)
  pair <- intersect(str_split(group, "[-+]")[[1]], single_cols)
  if (!length(pair)) return(NA_character_)
  pair[which.max(mean_rel_ab[asv, pair])]
}) %>% set_names(rn)

thr_euk_indic <- 0.015; thr_prok_indic <- 0.015
proptab.indic <- proptab.wk[, specific_asvs, drop = FALSE]
proptab.not.indic <- proptab.wk[, setdiff(colnames(proptab.wk), specific_asvs), drop = FALSE]

# Cosmopolitans are non-indicator ASVs present in at least 60% of samples in every non-riverine cluster.
cosmo_prev_cutoff <- 0.6
presence <- as.matrix(proptab.not.indic > 0)
keep <- !is.na(sample_cluster); presence <- presence[keep, , drop = FALSE]
presence_cluster <- factor(sample_cluster[keep], levels = single_cols)
present_by_clus <- rowsum(presence * 1L, group = presence_cluster)
prev_by_clus <- sweep(present_by_clus, 1, as.numeric(table(presence_cluster)[rownames(present_by_clus)]), "/")
nonriver_clus <- setdiff(single_cols, as.character(which(cl_labels == "MR")))
cosmo_asvs <- colnames(prev_by_clus)[colSums(prev_by_clus[nonriver_clus, , drop = FALSE] >= cosmo_prev_cutoff) == length(nonriver_clus)]
proptab.cosmo <- proptab.not.indic[, cosmo_asvs, drop = FALSE]

keep_asvs <- union(specific_asvs, cosmo_asvs)
proptab.indic.cosmo <- proptab.wk[, keep_asvs, drop = FALSE]
proptab.not.indic.cosmo <- proptab.wk[, setdiff(colnames(proptab.wk), keep_asvs), drop = FALSE]


################################
### Indicator summary tables ###
################################

clean_tax <- function(x) x %>% str_remove(";\\s*ASV_\\d+$") %>% str_remove_all("\\s*\\[\\d+\\]") %>% str_squish()

build_tax_tbl <- function(asv) {
  tx <- clean_tax(asv); parts <- str_split_fixed(tx, ";\\s*", 4)
  tibble(ASV = asv, Taxonomy = tx, Domain = na_if(parts[, 1], ""), Phylum = na_if(parts[, 2], ""),
         Order = na_if(parts[, 3], ""), Genus = na_if(parts[, 4], ""))
}

indic_asv_tbl <- build_tax_tbl(specific_asvs) %>% mutate(Clusters = unname(chosen_cluster[ASV]))
indic_genus_tbl <- indic_asv_tbl %>% group_by(Domain, Phylum, Order, Genus) %>%
  summarise(Clusters = str_c(sort(unique(na.omit(Clusters))), collapse = ", "), Num_ASVs = n(), .groups = "drop")


###############################################
### ASV-level information for later figures ###
###############################################

# Overall max-abundance cluster for each ASV; use only where that interpretation is appropriate.
asv_max_cluster_vec <- setNames(single_cols[max.col(mean_rel_ab[rn, , drop = FALSE], ties.method = "first")], rn)

tax_for_sel <- tibble(ASV = colnames(proptab.indic.cosmo)) %>%
  mutate(Domain = map_chr(str_split(clean_tax(ASV), ";\\s*"), 1),
         Domain2 = if_else(Domain %in% c("Bacteria", "Archaea"), "Prokaryote", "Eukaryote"))

samp_clus <- meta.wk %>% rownames_to_column("Sample") %>% transmute(Sample, Cluster = as.character(Cluster))

long0 <- proptab.indic.cosmo %>% as_tibble(rownames = "Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
  left_join(samp_clus, by = "Sample") %>% left_join(tax_for_sel %>% select(ASV, Domain2), by = "ASV") %>%
  left_join(indic_asv_tbl %>% transmute(ASV, Aff = as.character(Clusters)), by = "ASV") %>%
  mutate(use = is.na(Aff) | Cluster == Aff)

mean_aff <- long0 %>% group_by(ASV, Domain2) %>% summarise(meanRelAb_aff = mean(RelAb[use], na.rm = TRUE), .groups = "drop")
asv_in_dotplot <- mean_aff %>% filter(meanRelAb_aff >= if_else(Domain2 == "Eukaryote", thr_euk_indic, thr_prok_indic)) %>% pull(ASV)



######################
### ComplexHeatmap ###
######################

draw_wk_complex_hm <- function() {
  ht_opt$DIMNAME_PADDING <- unit(0.2, "mm")
  ht_opt$ROW_ANNO_PADDING <- unit(2.5, "mm")
  Heatmap(proptab.trans,
          name = "wk_complex_hm",
          col = viridis(256),
          cluster_rows = dend_row_order,
          row_split = k,
          cluster_columns = dend_col_order,
          # column_split = tax_k,
          right_annotation = rowAnnotation(
            # bar1 = anno_points(log1p(meta.wk$`Distance to coast [km]`), size = unit(3, "mm"), width = unit(2.6, "cm"),
            #                    axis_param = list(at = log1p(c(0, 25, 250)), labels = c(0, 25, 250), labels_rot = 0, gp = gpar(fontsize = 10))),
            # bar2 = anno_points(log10(meta.wk$`Depth [m]`), size = unit(3, "mm"), width = unit(2.6, "cm"),
            #                    axis_param = list(at = log10(c(1, 10, 75)), labels = c(1, 10, 75), labels_rot = 0, gp = gpar(fontsize = 10))),
            bar1 = anno_points(meta.wk$`Temperature [deg C]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar2 = anno_points(meta.wk$`Salinity`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar3 = anno_points(meta.wk$`Chlorophyll-a fluorescence [mg/m^3]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar4 = anno_points(meta.wk$`CDOM fluorescence [mg/m^3]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar5 = anno_points(meta.wk$`Beam attenuation [1/m]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            gap = unit(3, "mm"), show_annotation_name = FALSE),
          
          row_labels = paste0("  ", rownames(proptab.trans)),  # move rownames a bit to the right
          
          top_annotation = {
            asv_ids <- colnames(proptab.trans)     # ASVs in heatmap (columns)
            mc <- chosen_cluster[asv_ids]          # Member of the tested association selected for each bioindicator
            is_indic <- asv_ids %in% specific_asvs # bioindicators by name
            is_cosmo <- asv_ids %in% cosmo_asvs    # cosmopolitans by name
            pal_annot <- c(cluster_palette, cosmo = "black")      # Palette for annotation bars: clusters + cosmopolitans (black)
            
            ## --- TOP BAR: all indicators + all cosmopolitans in the heatmap ---
            top_vec <- rep(NA_character_, length(asv_ids))
            names(top_vec) <- asv_ids
            top_vec[is_indic] <- mc[is_indic]      # indicators: color by selected associated cluster
            top_vec[is_cosmo] <- "cosmo"           # cosmopolitans: black
            
            ## --- BOTTOM BAR: only ASVs that appear in the INDIC+COSMO dotplot ---
            bottom_vec <- rep(NA_character_, length(asv_ids))
            names(bottom_vec) <- asv_ids
            dot_common <- intersect(asv_in_dotplot, asv_ids)
            cosmo_dot  <- intersect(cosmo_asvs, dot_common)         # cosmopolitans in dotplot
            indic_dot  <- setdiff(dot_common, cosmo_dot)            # indicators in dotplot
            bottom_vec[indic_dot] <- mc[match(indic_dot, asv_ids)]  # indicators in dotplot: selected associated cluster
            bottom_vec[cosmo_dot] <- "cosmo"                        # cosmopolitans in dotplot: black
            
            # Factors with explicit levels: clusters 1..k + "cosmo"
            top    <- factor(top_vec,    levels = c(as.character(1:k), "cosmo"))
            bottom <- factor(bottom_vec, levels = c(as.character(1:k), "cosmo"))
            
            columnAnnotation("Bioindicators & cosmopolitans" = top, "Most abundant of the above" = bottom,
                             col = list("Bioindicators & cosmopolitans" = pal_annot, "Most abundant of the above" = pal_annot),
                             show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
                             show_legend = FALSE)
          },
          
          heatmap_height = unit(40, "cm"), heatmap_width = unit(68, "cm"),
          row_dend_width = unit(4, "cm"), column_dend_height = unit(2, "cm"),
          row_title = "Samples", row_title_gp = gpar(fontsize = 34, fontface = "bold"),
          column_title = "ASVs", column_title_gp = gpar(fontsize = 34, fontface = "bold"),
          row_names_gp = gpar(fontsize = 13, fontface = "bold"),
          row_names_max_width = unit(7.6, "cm"),
          column_names_gp = gpar(fontsize = 1, fontface = "bold"), column_names_rot = 80,
          
          show_heatmap_legend = FALSE,
          heatmap_legend_param = list(
            title = NULL,
            # title = "Fourth-root transformed relative abundances\n",
            title_gp = gpar(fontsize = 30, fontface = "bold"),
            legend_direction = "horizontal", legend_width = unit(8, "cm"),
            color_bar = "continuous")) %>%
    ComplexHeatmap::draw(., heatmap_legend_side = "top", padding = unit(c(10, 10, 10, 10), "mm"))  # top, right, bottom, left)
  
  # Decorate slices (colors first, then labels on top)
  for (i in 1:k) {
    decorate_row_names("wk_complex_hm", {
      grid.rect(gp = gpar(col = cluster_palette[i], lwd = 3.2, fill = NA),
                width = unit(7.6, "cm"), x = unit(1, "npc") + unit(-0.4, "cm"), just = "center")
    }, slice = i)
  }
  for (i in 1:k) {
    decorate_row_names("wk_complex_hm", {
      grid.text(cl_labels[i], x = unit(2, "npc") + unit(-4, "cm"), just = "left",
                gp = gpar(fontsize = 27, fontface = "bold"))
    }, slice = i)
  }
  # for (i in seq_along(com_labels)) {
  #   decorate_column_names("wk_complex_hm", {
  #     grid.text(com_labels[i], y = unit(1, "npc") - unit(2.5, "cm"), just = "bottom",
  #               gp = gpar(fontsize = 25, fontface = "bold.italic"))
  #   }, slice = i)
  # }
  
  annotation_labels <- list(
    list(name = "Temp", unit = expression(bold("(" * degree * "C)"))),
    list(name = "Sal", unit = NULL),
    list(name = "fChl-a", unit = expression(bold("(mg m"^{-3} * ")"))),
    list(name = "fCDOM", unit = expression(bold("(mg m"^{-3} * ")"))),
    list(name = "BA", unit = expression(bold("(m"^{-1} * ")")))
  )
  
  for (i in seq_along(annotation_labels)) {
    decorate_annotation(paste0("bar", i), {
      if (is.null(annotation_labels[[i]]$unit)) {
        grid.text(annotation_labels[[i]]$name,
                  y = unit(1, "npc") + unit(5.5, "mm"), just = "bottom",
                  gp = gpar(fontsize = 18, fontface = "bold"))
      } else {
        grid.text(annotation_labels[[i]]$name,
                  y = unit(1, "npc") + unit(11.0, "mm"), just = "bottom",
                  gp = gpar(fontsize = 18, fontface = "bold"))
        
        grid.text(annotation_labels[[i]]$unit,
                  y = unit(1, "npc") + unit(1.5, "mm"), just = "bottom",
                  gp = gpar(fontsize = 18, fontface = "bold"))
      }
    })
  }
}

cairo_pdf(file = "OUT/3W_CLUS/wk_complex_hm.pdf", width = 31, height = 17)
draw_wk_complex_hm()
dev.off()

png("OUT/3W_CLUS/wk_complex_hm.png", width = 31, height = 17, units = "in", res = 1000)
draw_wk_complex_hm()                       # calls the function which includes all decorations
dev.off()


######################################
### Pure Dendrogram + Sample Names ###
######################################

# Create a completely empty heatmap structure
empty_mat <- matrix(0, nrow = nrow(proptab.trans), ncol = 0)
rownames(empty_mat) <- rownames(proptab.trans)

# Create the heatmap w/ a name
ht <- Heatmap(empty_mat,
              name = "dendro_plot",
              cluster_rows = dend_row_order,
              row_split = k,
              row_labels = rownames(proptab.trans),
              row_names_side = "right",
              row_dend_side = "left",
              show_row_names = TRUE,
              show_column_names = FALSE,
              show_heatmap_legend = FALSE,
              row_title = NULL,
              column_title = NULL,
              heatmap_height = unit(55, "cm"),
              row_dend_width = unit(10, "cm"),
              row_names_gp = gpar(fontsize = 20, fontface = "bold"))

cairo_pdf(file = "OUT/3W_CLUS/dendro_pure.pdf", width = 12, height = 24)
draw(ht)
dev.off()

png("OUT/3W_CLUS/dendro_pure.png", width = 12, height = 24, units = "in", res = 300)
draw(ht)
dev.off()


#############
### TITAN ###
#############

# Offshore surface samples and ASVs occurring in at least three of them.
off_surf <- meta.wk %>% rownames_to_column("Sample") %>% filter(`Distance to coast [km]` > 10, between(`Depth [m]`, 1, 5))
txa_off_surf <- proptab.wk[off_surf$Sample, , drop = FALSE] %>% .[, colSums(. > 0) >= 3, drop = FALSE]
subset_titan <- function(ti, keep) { ti$sppmax <- ti$sppmax[keep, , drop = FALSE]; ti$taxa <- ti$taxa[, keep, drop = FALSE]; ti }

# Environmental gradients analyzed with the same workflow.
titan_specs <- c(salinity = "Salinity", latitude = "Latitude [deg]")
titan_results <- list()

for (nm in names(titan_specs)) {
  variable <- titan_specs[[nm]]
  
  # Recalculate and save only when needed.
  # set.seed(161193)
  # ti <- titan(off_surf[[variable]], txa_off_surf)
  # saveRDS(ti, paste0("EXT/titan_offshore_surface_", nm, ".rds"))
  
  ti <- readRDS(paste0("EXT/titan_offshore_surface_", nm, ".rds"))
  titan_results[[nm]] <- ti
  
  print(ti$sumz.cp[c("sumz-", "sumz+", "fsumz-", "fsumz+"), ])
  print(range(ti$sumz.cp[c("sumz-", "sumz+"), "cp"]))
  
  keep_euk <- rownames(ti$sppmax) %in% intersect(asv_in_dotplot, specific_asvs) &
    grepl("^Eukaryota", rownames(ti$sppmax)) & ti$sppmax[, "zscore"] >= 5
  keep_prok <- rownames(ti$sppmax) %in% intersect(asv_in_dotplot, specific_asvs) &
    grepl("^(Bacteria|Archaea)", rownames(ti$sppmax)) & ti$sppmax[, "zscore"] >= 5
  
  p_dens <- plot_sumz_density(ti)
  p_euk <- plot_taxa_ridges(subset_titan(ti, keep_euk), n_ytaxa = sum(keep_euk), printspp = FALSE, xlabel = variable)
  p_prok <- plot_taxa_ridges(subset_titan(ti, keep_prok), n_ytaxa = sum(keep_prok), printspp = FALSE, xlabel = variable)
  
  # PDF (exploration): top panel shows the bootstrap distributions of community-level negative and positive change points;
  # Lower panels show taxon-specific change-point distributions for eukaryotic and prokaryotic bioindicator ASVs
  # retained in asv_in_dotplot and having a TITAN z-score >=5 (blue = decreasing; red = increasing along the gradient).
  pdf(paste0("OUT/3W_CLUS/titan_offshore_surface_", nm, ".pdf"), width = 12, height = 24)
  print(p_dens / p_euk / p_prok)
  dev.off()
  
  # CSV (change points reported in Section 3.3): community change-point estimates and bootstrap intervals calculated using all TITAN taxa (sumz-/sumz+)
  # and using only taxa with "pure" and "reliable" bootstrap responses (fsumz-/fsumz+).
  # "pure": the taxon consistently responds in the same direction across bootstraps;
  # "reliable": its response is repeatedly statistically supported across bootstraps.
  write.csv(ti$sumz.cp, paste0("OUT/3W_CLUS/titan_offshore_surface_", nm, ".csv"))
}


###########################################
### Check for flowthrough contamination ###
###########################################

# Test whether any offshore surface ASV is specifically associated with flow-through rather than rosette collection.
# Separate from TITAN: offshore surface samples (>10 km from the coast, 1-5 m depth) and all ASVs present in any of them.
ft_samples <- meta.wk %>% rownames_to_column("Sample") %>% filter(`Distance to coast [km]` > 10, between(`Depth [m]`, 1, 5))
txa_ft <- proptab.wk[ft_samples$Sample, , drop = FALSE] %>% .[, colSums(.) > 0, drop = FALSE]
group <- ifelse(ft_samples$Sampling == "FT", "FT", "nonFT")
table(group)   # FT = 25, nonFT = 14

set.seed(161193)
ind <- multipatt(txa_ft, group, func = "IndVal.g", control = how(nperm = 9999))

ft_idx <- match("FT", colnames(ind$A))
ft_tested <- !is.na(ind$sign$index) & ind$sign$index == ft_idx & ind$sign$p.value < 0.01
if (any(ft_tested)) plot(ind$A[ft_tested, ft_idx], ind$B[ft_tested, ft_idx], xlab = "A", ylab = "B")
keep_ft <- ft_tested & ind$A[, ft_idx] >= 0.6 & ind$B[, ft_idx] >= 0.6
ft_asvs <- rownames(ind$A)[keep_ft]
if (length(ft_asvs)) {
  cbind(FT = colMeans(txa_ft[group == "FT", ft_asvs, drop = FALSE]),
        nonFT = colMeans(txa_ft[group == "nonFT", ft_asvs, drop = FALSE]))} else {"No FT indicators"}

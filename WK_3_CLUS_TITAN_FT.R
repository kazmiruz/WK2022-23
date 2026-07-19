#######################################################################################################
### Things to keep in mind before before diving into the analysis here and in the following scripts ###

#   Since microbes are the main drivers of biochem-l cycles, knowing the state of microbial communities provides
#   more and better leverages to explain env-l changes, as opposed to just looking at "symptoms" (e.g.,chl a, CO2, O2, etc.),
#   just how doctors would rather know the parasite causing the disease than blindly treat the symptoms.
#   Moreover, microbes can serve as versatile sensors reacting to env-l changes, similar to diagnostic biomarkers for detecting cancer.
#   For example, specific taxa of SAR11, can be used 

# - EXPLORE the data: use simple figures/tables to see what data is available for samples and if certain parameters contribute to the story.
# - Finding out something is nice, but focus on the NOVEL and/or UNEXPECTED findings/patterns
#   (e.g., surface vs subsurface MICB being different is already well known, but Synechococcus being abundant in Hudson Bay is novel)?
# - From novel findings, what are some w/ IMPLICATIONS for target audiences (e.g., modelers, other researchers, local population)
# - How can you SHOW these novel important findings with figures? What MESSAGE are you trying to deliver with each figures?
# - Finally, what is the overall STORY you are trying to tell in the manuscript?


# WK_3---------------------------------------------------------------------
# Purpose:
#   Assess cluster stability (silhouette / bootstrap Jaccard), build a Ward‑D2 dendrogram (rotated by latitude),
#   perform indicator species analysis (IndVal.g), select bioindicators and cosmopolitans, produce a ComplexHeatmap and a pure dendrogram,
#   run TITAN for community‑level change points along salinity and latitude, and check for flow‑through contamination.
#
# Prerequisites:
#   The objects `proptab.small.taxan.saman`, `meta.wk.small`, and `ctd` must be in the workspace (created by WK_2_TAX.R and WK_1_META.R).
#
# Required external files - you can also reproduce them all yourself using the commented out code:
#   OUTPUT/indicator_res.rds
#   OUTPUT/titan_offshore_surface_salinity.rds
#   OUTPUT/titan_offshore_surface_lat.rd
#
# Output:
#   OUTPUT/wk_simple_hm.pdf                      (diagnostic heatmap)
#   OUTPUT/wk_complex_hm.pdf / .png              (annotated ComplexHeatmap)
#   OUTPUT/dendro_pure.pdf / .png                (dendrogram with cluster labels)
#   OUTPUT/titan_offshore_surface_salinity.pdf   (TITAN salinity ridges)
#   OUTPUT/titan_offshore_surface_salinity.csv   (change‑point summary)
#   OUTPUT/titan_offshore_surface_lat.pdf        (TITAN latitude ridges)
#   OUTPUT/titan_offshore_surface_lat.csv        (change‑point summary)
# -------------------------------------------------------------------------


library(tidyverse)
library(gplots)         # for simple heatmaps (heatmap.2)
library(ComplexHeatmap) # for well-annotated heatmaps
library(patchwork)
library(viridis)        # colouring, viridis()
library(viridisLite)    # for turbo coloring (b/c better than jet/rainbow for colour distinguishing and colour blindness), turbo()
library(scales)         # for alpha
library(cluster)        # for silhoutte plot
library(dendextend)     # for cutree and rotate
library(vegan)          # vegdist
library(fpc)            # cluster 
library(indicspecies)   # for indicator species analysis
library(TITAN2)

set.seed(161193)
options("lifecycle_verbosity"="warning") # always get warnings if there are (not just once every 8 hours)

#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
# Apply one of these transformation for optimal clustering/grouping, adjust based on visual inspection of historgrams and heatmaps
# fourth root (^.25) and square root (^.5) transformations are more acceptable / easier to explain to reviewers

hist(as.matrix(proptab.small.taxan.saman))        # looks pretty bad
hist(as.matrix(proptab.small.taxan.saman)^.5)     # looks better
hist(as.matrix(proptab.small.taxan.saman)^.25)    # looks good: a lot of 0s but normal distribution to the left
hist(as.matrix(proptab.small.taxan.saman)^.125)   # looks good - use if need more pronounced clusters

proptab.trans <- as.matrix(proptab.small.taxan.saman)^.25  # your transformed ASV table for grouping: row/col order is the same


##########################################################
### What is the optimal number of sample/taxon groups? ###
# use silhouette analysis for your Bray-Curtis dissimilarity; don't use WCSS as it is for Euclidean distance

eval_k <- function(x, ks = 5:9, B = 1000) {
  d <- vegdist(x, "bray"); tr <- hclust(d, "ward.D2")
  purrr::map_dfr(ks, \(k) {
    sil <- cluster::silhouette(cutree(tr, k), d)[, 3]
    jacc <- fpc::clusterboot(d, B = B, bootmethod = "boot",
                             clustermethod = fpc::disthclustCBI,
                             method = "ward.D2", k = k)$bootmean
    tibble::tibble(k = k,
                   sil_mean = mean(sil), sil_min = min(tapply(sil, cutree(tr, k), mean)),
                   jacc_mean = mean(jacc), jacc_min = min(jacc))
  })
}

res <- eval_k(proptab.trans)

ggplot(res, aes(k)) +
  geom_line(aes(y = sil_mean, color = "Mean silhouette")) +
  geom_point(aes(y = sil_mean, color = "Mean silhouette")) +
  geom_line(aes(y = jacc_min, color = "Minimum bootstrap Jaccard")) +
  geom_point(aes(y = jacc_min, color = "Minimum bootstrap Jaccard"))

res %>% filter(jacc_min > 0.6) %>% slice_max(sil_mean, n = 1, with_ties = FALSE) %>% pull(k) # 9


# sil_widths <- \(x, mode = c("rows", "columns"), k_range = 3:20) {
#   mode <- match.arg(mode)
#   if (mode == "columns") x <- t(x)
#   d <- vegdist(x, method = "bray")
#   tree <- hclust(d, method = "ward.D2")
#   sapply(k_range, \(k) mean(silhouette(cutree(tree, k), d)[, 3]))
# }
# ggplot(data.frame(k = 3:20, width = sil_widths(proptab.trans, mode = "rows", k_range = 3:20)),
#        aes(k, width)) + geom_line() + geom_point() + labs(title = "samples") |
# ggplot(data.frame(k = 3:20, width = sil_widths(proptab.trans, mode = "columns", k_range = 3:20)),
#        aes(k, width)) + geom_line() + geom_point() + labs(title = "ASVs")

# library(NbClust) # can also do silhoutte with NbClust: same result
# ggplot(data.frame(k = 3:20, width = NbClust(proptab.trans, diss = vegdist(proptab.trans, method = "bray"),             # for samples
#                                             distance = NULL, min.nc = 3, max.nc = 20, method = "ward.D2", index = "silhouette")$All.index),
#        aes(k, width)) + geom_line() + geom_point() + labs (title = "samples") |
# ggplot(data.frame(k = 3:20, width = NbClust(t(proptab.trans), diss = vegdist(t(proptab.trans), method = "bray"),       # for ASVs
#                                             distance = NULL, min.nc = 3, max.nc = 20, method = "ward.D2", index = "silhouette")$All.index),
#        aes(k, width)) + geom_line() + geom_point() + labs (title = "ASVs")


#
## Test how your selected number of groups (sample clusters and microbial communities) looks on a heatmap
pdf("OUTPUT/wk_simple_hm.pdf", width = 35, height = 20)

hm <- heatmap.2(proptab.trans, trace = "none", col = viridis,
                distfun = \(x) vegdist(x, "bray"),      # \(x) is function(x)
                hclustfun = \(x) hclust(x, "ward.D2"))
row_gr <- cutree(hm$rowDendrogram, k = 7)
col_gr <- cutree(hm$colDendrogram, k = 5)

heatmap.2(proptab.trans, trace = "none", key = FALSE, margins = c(20, 30), col = viridis,    # replot with group's colors
          distfun = \(x) vegdist(x, "bray"),
          hclustfun = \(x) hclust(x, "ward.D2"),
          RowSideColors = turbo(7)[row_gr],      # the visual order of groups does NOT correspond to its numbering (e.g., 312 from top)
          ColSideColors = turbo(5)[col_gr])

dev.off()


#####################
### Cluster build ###

#
## Basic build of groups (clusters/communities)

# How many groups and do we call them?
k <- 9         # number of sample clusters (rows)
# tax_k <- 10     # number of taxon clusters (columns), not used, but kept
cl_labels <- c("oHB", "BI22", "cHB", "dBI", "dJB", "EST", "MR", "nJB", "sJB")
# com_labels <- LETTERS[1:tax_k]

# Basic row clustering
dist_row <- vegdist(proptab.trans, method = "bray")
tree_row <- hclust(dist_row, method = "ward.D2")
dend_row <-  as.dendrogram(tree_row)                   # row/sample tree

# Basic column clustering
dist_col <- vegdist(t(proptab.trans), method = "bray")
tree_col <- hclust(dist_col, method = "ward.D2")
dend_col <-  as.dendrogram(tree_col)                   # taxa/column tree

#
## Reorder rows based on salinity from highest to lowest

# Default ordering
# dend_row_order <- dend_row

# Local branch rotation
# dend_row_order <- reorder(dend_row, wts=-(meta.wk.small$`Salinity [PSU]`), agglo.FUN = mean)
# ?reorder/dendrogram
# dendrogram leaves are reordered to be in an order as consistent as possible with the wts (note: meta.wk.small is in the same order as proptab.trans)
# at each split, reorder() computes a sum wt for each child branch by applying agglo.FUN to the leaves beneath it,
# then flips (or not) the node to order by incr wt
# here, it works by: what's the median salinity at each branch?
# only entire subtrees are flipped and the underlying tree topology is not broken

# Global branch rotation (still calculates at each node)
ord_row <- rownames(meta.wk.small)[order(meta.wk.small$`Latitude [deg]`, decreasing = TRUE)] # sample names ordered from high to low latitude
dend_row_order <- dendextend::rotate(dend_row, ord_row)               # force the dendrogram appear in a SIMILAR order by “rotating” branches
# here, ord_row - from high to low S - is "the order which the user would like to see the tree rotated by"
# and rotate "attempts to rotate the tree to fit that order"
# hence, rotate() flips any left-right splits so that, as closely as possible, the leaves read out in the order asked
# only entire subtrees are flipped and the underlying tree topology is not broken
# ord_row
# labels(dend_row_order)[ order.dendrogram(dend_row_order) ]

# Just like reorder(), rotate() only ever flips two sister branches at a time and never rewires the tree.
# But! While rotate() keeps flipping ALL branches to achieve the OVERALL decreasing in S order (w/out breaking the tree ofc),
# reorder() only compares the median S of the 2 child subtrees,
# hence it makes changes more LOCALLY and thus preserves the tree’s Ward‐clustering structure, but doesn't order by decreasing S as well

#
## Cluster labeling and coloring

# Get visual top-to-bottom order of samples from dendrogram
row_ordered <- rownames(proptab.trans)[order.dendrogram(dend_row_order)]

# Cut dendrogram tree into k clusters using that order
clus <- cutree(tree_row, k = k)[row_ordered] %>%
  (\(x) factor(x, levels = unique(x), labels = as.character(1:k)))()        # assign 1 to k labels

# Set color palette after final cluster labels (aligned with top-to-bottom order)
# cluster_palette <- c("#F233F2", "#0000F2", "#00F2F2", "#00F200", "#A6734D", "#D9D900", "#F20000")
cluster_palette <- colorRampPalette(c("darkviolet","magenta2",
                                      "darkgoldenrod",
                                      # "yellow3",
                                      "red2", "coral",
                                      # "olivedrab",
                                      "green2", "cyan3",
                                      "blue", "dodgerblue2"))(k)

names(cluster_palette) <- levels(clus)  # e.g., "1", "2", ..., "k"

# Apply coloring/length to the Cluster dendrogram
dend_row_order <- color_branches(dend_row_order, k = k, col = cluster_palette) %>% set("branches_lwd", 3)

# Adding the cluster info to the meta (make sure to match the samples) and creating Cluster labels
meta.wk.small$Cluster <- clus[rownames(meta.wk.small)]
row_slice_labels <- paste0("Cluster ", levels(factor(clus)))

#
## Find a good way to find which latitude each taxon prefers the most

# Presence/Absence
# latwted = colMeans(sweep((proptab.trans>0), STATS = meta.wk.small$`Latitude [deg]`, MARGIN = 1, FUN = "*"))
# (proptab.trans>0) turns it into presence/absence (1/0), then with sweep() each row/sample is multiplied by its latitude (0->0; 1->S),
# then get mean Lat per column/taxon
# all.equal(colMeans(sweep((proptab.trans>0), STATS = meta.wk.small$`Latitude [deg]`, MARGIN = 1, FUN = "*")),
#           colMeans((proptab.trans>0) * meta.wk.small$`Latitude [deg]`))    # even simpler, don't need sweep, just multiply your matrix by Lat

# use weighted RelAB rather than presence/absence, b/c even if the taxon doesn't like the given Lat, it may still be technically "present", hence = 1
# e.g., the taxon is 0.001 at 20, 0.01 at 25, and 0.1 at 30 -> pres/abs will tell you that its optimal Lat is 25, but it is 30, given no other factors
# how big the RelAb at the favorite Lat doesn't really matter - only its value to the rest of the samples

# Weighted RelAb
# apply(proptab.trans, 2, …) - take each column (taxon's RelAb in different samples) and transfer to the function as x
# weighted.mean(y, w=x) computes sum(y*w)/sum(w) aka multiply Lat of each sample by the RelAb of the taxon in that sample ->
# sum it -> divide  by the sum of RelAb of the taxon
latwted <- apply(proptab.trans, 2, function(x) {weighted.mean(meta.wk.small$`Latitude [deg]`, w = x)})  # each column

#
## Reorder cols based on salinity from highest to lowest

# Default ordering
# dend_col_order <- dend_col

# Local branch rotation
# dend_col_order <- reorder(dend_col, wts=-(latwted), agglo.FUN = median)

# Global branch rotation
ord_col <- names(sort(latwted, decreasing = TRUE))
dend_col_order <- dendextend::rotate(dend_col, ord_col)

#
## Extract column order (and assign community labels if you want them)

# Get visual left-to-right order of taxa from dendrogram
col_ordered <- colnames(proptab.trans)[order.dendrogram(dend_col_order)]

# # Cut dendrogram tree into tax_k clusters using that order
# com <- cutree(tree_col, k = tax_k)[col_ordered] %>%
#   (\(x) factor(x, levels = unique(x), labels = com_labels))()           # com_labels are letters in the alphabetical order

# # Apply length to the Community dendrogram
# dend_col_order <- dend_col_order %>% set("branches_lwd", 3)

# # Adding the community info to a data frame and creating Community labels
# community <- data.frame(ASV = names(com), Community = as.character(com))
# col_slice_labels <- paste0("Network ", levels(com))

##################################
### Indicator species analysis ###
# use meta.wk.small$Cluster b/c while clus DOES have the correct sample-cluster assignment, the sample order ISN'T the same as in proptab.small.taxan.saman
#
# indicator_res <- multipatt(proptab.small.taxan.saman, # use the original ASV table b/c real abundance matters
#                            meta.wk.small$Cluster,
#                            func = "IndVal.g",
#                            duleg = FALSE,  # If TRUE, Cluster combos aren't considered, only the orig Clusters (1 to k, Dufrêne & Legendre, 1997)
#                            # restcomb <- c(1:10, 11, 28, 41, 50, 55), # 10 clusters, 1+2, 3+4, 5+6, 7+8, 9+10
#                            # restcomb <- c(1:55), # 10 cluster + all pairs (45); max.order = 2 should also work
#                            max.order = 2,
#                            control = how(nperm = 999))
# write_rds(indicator_res, "OUTPUT/indicator_res.rds")       # 9 clusters w/ pairs; takes 5 min to run
indicator_res <- read_rds("OUTPUT/indicator_res.rds")

# # Explore the indicspecies output
# indicator_res$cluster # clusters of 72 samples
# indicator_res$comb    # cluster combinations each sample belongs to - I consider only singles here
# indicator_res$A       # "A" component of indicator values (IndVal) for each cluster (or cluster combination) for each ASV: "specificity" of ASV;
# # indicates how strongly the ASV is assoc-d with the target Cluster; Mean RelAb in the target Cluster samples / SUM of Mean RelAbs in ALL samples
# indicator_res$B       # "B" component of IndVal: "fidelity" of each ASV; for the target Cluster: the # of samples where the ASV is found / overall # of samples
# indicator_res$sign    # s.1, s.2, etc. - which clus is pref by the ASV; index - comb-n; stat - sq.root of IndVal; p.value - sign-ce of stat based on 999 perm-ns
# indicator_res$sign$stat      # stat = sq.root(IndVal), where IndVal = A*B

# hist(indicator_res$sign$stat, breaks=100) # look at the distribution of stat; most are 0.5-0.8
# plot(proptab.small.taxan.saman %>% colMeans(), indicator_res$sign$stat, log = 'xy') # plot Mean RelAb (for each col aka ASV) vs stat for each ASV; same order
# mean(indicator_res$A) # 0.12
# mean(indicator_res$B) # 0.25
# hist(indicator_res$A)
# hist(indicator_res$B)

# 999 permutations: shuffle the meta.wk.small$Cluster labels across samples while keeping the ASV abundances fixed -> recompute that ASV’s IndVal.g for the same Cl.
# This simulates what values of IndVal.g you’d expect if there is no real link between the ASV's distribution and and the Cluster.
# Goal: fewer than 1% of my random‐label runs ever produced an IndVal.g as large as the real one.
# hist(indicator_res$sign$p.value, breaks=100)   # most are <0.1

# Select only Indicators for each specific Cluster (no Indicators Cluster combinations; don't need [, 1:k] if max.order = 1)
# A: the Mean RelAb of the ASV in the Cluster samples divided by the SUM of the Mean RelAbs in each Cluster (e.g., if 8 C, then 8 Mean RelAbs summed up)
# B: at least 50% of the samples in the Cluster contain the ASV
# p.value: the relationship is significant

###########################################
### Indicator selection + cosmopolitans ###
rn <- rownames(indicator_res$A); cn <- colnames(indicator_res$A); ng <- ncol(indicator_res$A)
idx <- indicator_res$sign$index                                             # column of the single best pattern actually permutation-tested per ASV
valid <- !is.na(idx) & indicator_res$sign$p.value < .01 &                   # significance + A/B must hold for THAT tested pattern...
  indicator_res$A[cbind(seq_along(rn), idx)] >= .6 &                        # ...not some other A-maximizing combination whose p-value was never tested
  indicator_res$B[cbind(seq_along(rn), idx)] >= .6
pick <- ifelse(valid, idx, NA_integer_)                                     # keep the tested pattern (no re-pick by max.col over A)
specific_asvs <- rn[!is.na(pick)]
chosen_group <- setNames(rep(NA_character_, length(rn)), rn); chosen_group[!is.na(pick)] <- cn[pick[!is.na(pick)]]

single_cols <- cn[!str_detect(cn, "[-+]")]
chosen_cluster <- map2_chr(chosen_group, rn, ~{
  if (is.na(.x) || !str_detect(.x, "[-+]")) return(.x)
  clus <- intersect(str_split(.x, "[-+]")[[1]], single_cols); if (!length(clus)) return(NA_character_)
  clus[which.max(indicator_res$A[match(.y, rn), clus])]
}) %>% set_names(rn)

thr_euk_indic <- 0.015; thr_prok_indic <- 0.015

proptab.indic <- proptab.small.taxan.saman[, specific_asvs, drop = FALSE]
proptab.not.indic <- proptab.small.taxan.saman %>% select(-all_of(specific_asvs))

cosmo_clus_cutoff <- 0.85; cosmo_prev_cutoff <- 0.6
presence <- (proptab.not.indic > 0) %>% as.matrix()
clus <- meta.wk.small$Cluster[match(rownames(presence), rownames(meta.wk.small))]
keep <- !is.na(clus); presence <- presence[keep, , drop = FALSE]; clus <- factor(clus[keep])
n_per <- as.numeric(table(clus)[levels(clus)])
prev_by_clus <- sweep(rowsum(presence * 1L, group = clus), 1, n_per, "/")
cosmo_asvs <- names(colMeans(prev_by_clus >= cosmo_prev_cutoff))[colMeans(prev_by_clus >= cosmo_prev_cutoff) >= cosmo_clus_cutoff]
proptab.cosmo <- proptab.not.indic[, cosmo_asvs, drop = FALSE]

keep_asvs <- union(specific_asvs, cosmo_asvs)
proptab.indic.cosmo <- proptab.small.taxan.saman[, keep_asvs, drop = FALSE]
proptab.not.indic.cosmo <- proptab.small.taxan.saman[, setdiff(colnames(proptab.small.taxan.saman), keep_asvs), drop = FALSE]

## =========================
## Indicator summary table, if you want to look up which ASVs or genera characterize which clusters
clean_tax <- function(x) x %>% str_remove(";\\s*ASV_\\d+$") %>% str_remove_all("\\s*\\[\\d+\\]") %>% str_squish()
build_tax_tbl <- function(asv) { tx <- clean_tax(asv); parts <- str_split_fixed(tx, ";\\s*", 4)
tibble(ASV = asv, Taxonomy = tx, Domain = na_if(parts[,1],""), Phylum = na_if(parts[,2],""), Order = na_if(parts[,3],""), Genus = na_if(parts[,4],"")) }

indic_asv_tbl <- bind_rows(build_tax_tbl(specific_asvs) %>% mutate(Clusters = chosen_cluster[ASV]))

indic_genus_tbl <- indic_asv_tbl %>% group_by(Domain, Phylum, Order, Genus) %>%
  summarise(Clusters = str_c(sort(unique(na.omit(Clusters))), collapse = ", "), Num_ASVs = n(), .groups="drop")

##############################################
### ASV-level info for heatmap annotations ###
# careful: it assigns cluster to all based on highest RelAb, but we should use it only for bioindicators!
asv_max_cluster_vec <- setNames(single_cols[max.col(indicator_res$A[, single_cols, drop = FALSE], ties.method = "first")], rn)

tax_for_sel <- tibble(ASV = colnames(proptab.indic.cosmo)) %>%
  mutate(Domain = map_chr(str_split(clean_tax(ASV), ";\\s*"), 1),
         Domain2 = if_else(Domain %in% c("Bacteria","Archaea"), "Prokaryote", "Eukaryote"))

samp_clus <- meta.wk.small %>% rownames_to_column("Sample") %>% select(Sample, Cluster)

long0 <- proptab.indic.cosmo %>% as_tibble(rownames = "Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
  left_join(samp_clus, by = "Sample") %>% left_join(tax_for_sel %>% select(ASV, Domain2), by = "ASV") %>%
  left_join(indic_asv_tbl %>% transmute(ASV, Aff = as.character(Clusters)), by = "ASV") %>%
  mutate(use = is.na(Aff) | Cluster == Aff)

mean_aff <- long0 %>% group_by(ASV, Domain2) %>% summarise(meanRelAb_aff = mean(RelAb[use], na.rm = TRUE), .groups = "drop")
asv_in_dotplot <- mean_aff %>% filter(meanRelAb_aff >= if_else(Domain2 == "Eukaryote", thr_euk_indic, thr_prok_indic)) %>% pull(ASV)


######################
### ComplexHeatmap ###
pdf(file = "OUTPUT/wk_complex_hm.pdf", width = 31, height = 17)
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
            # bar1 = anno_points(log1p(meta.wk.small$`Distance to coast [km]`), size = unit(3, "mm"), width = unit(2.6, "cm"),
            #                    axis_param = list(at = log1p(c(0, 25, 250)), labels = c(0, 25, 250), labels_rot = 0, gp = gpar(fontsize = 10))),
            # bar2 = anno_points(log10(meta.wk.small$`Depth [m]`), size = unit(3, "mm"), width = unit(2.6, "cm"),
            #                    axis_param = list(at = log10(c(1, 10, 75)), labels = c(1, 10, 75), labels_rot = 0, gp = gpar(fontsize = 10))),
            bar1 = anno_points(meta.wk.small$`Temperature [deg C]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar2 = anno_points(meta.wk.small$`Salinity [PSU]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar3 = anno_points(meta.wk.small$`Chlorophyll-a fluorescence [mg/m^3]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar4 = anno_points(meta.wk.small$`CDOM fluorescence [mg/m^3]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            bar5 = anno_points(meta.wk.small$`Beam attenuation [1/m]`, size = unit(3, "mm"), width = unit(2.6, "cm"),
                               axis_param = list(labels_rot = 0, gp = gpar(fontsize = 10))),
            gap = unit(3, "mm"), show_annotation_name = FALSE),
          
          row_labels = paste0("  ", rownames(proptab.trans)),  # move rownames a bit to the right

          top_annotation = {
            asv_ids <- colnames(proptab.trans)     # ASVs in heatmap (columns)
            mc <- asv_max_cluster_vec[asv_ids]     # Cluster where each ASV is most abundant (character "1".."11")
            is_indic <- asv_ids %in% specific_asvs # Indicators and cosmopolitans by name
            is_cosmo <- asv_ids %in% cosmo_asvs
            pal_annot <- c(cluster_palette, cosmo = "black")      # Palette for annotation bars: clusters + cosmopolitans (black)
            
            ## --- TOP BAR: all indicators + all cosmopolitans in the heatmap ---
            top_vec <- rep(NA_character_, length(asv_ids))
            names(top_vec) <- asv_ids
            top_vec[is_indic] <- mc[is_indic]      # indicators: color by max-abundance cluster
            top_vec[is_cosmo] <- "cosmo"           # cosmopolitans: black
            
            ## --- BOTTOM BAR: only ASVs that appear in the INDIC+COSMO dotplot ---
            bottom_vec <- rep(NA_character_, length(asv_ids))
            names(bottom_vec) <- asv_ids
            dot_common <- intersect(asv_in_dotplot, asv_ids)
            cosmo_dot  <- intersect(cosmo_asvs, dot_common)         # cosmopolitans in dotplot
            indic_dot  <- setdiff(dot_common, cosmo_dot)            # indicators in dotplot
            bottom_vec[indic_dot] <- mc[match(indic_dot, asv_ids)]  # indicators in dotplot: cluster of max abundance
            bottom_vec[cosmo_dot] <- "cosmo"                        # cosmopolitans in dotplot: black
            
            # Factors with explicit levels: clusters 1..k + "cosmo"
            top    <- factor(top_vec,    levels = c(as.character(1:k), "cosmo"))
            bottom <- factor(bottom_vec, levels = c(as.character(1:k), "cosmo"))
            
            columnAnnotation("Bioindicators & cosmopolitans" = top, "Most abundant of the above" = bottom,
                             col = list("Bioindicators & cosmopolitans" = pal_annot, "Most abundant of the above" = pal_annot),
              show_annotation_name = TRUE, annotation_name_side = "left", annotation_name_gp   = gpar(fontsize = 16, fontface = "bold"), show_legend = FALSE)
          },

          heatmap_height = unit(40, "cm"), heatmap_width = unit(68, "cm"),
          row_dend_width = unit(4, "cm"), column_dend_height = unit(2, "cm"),
          row_title = "Samples", row_title_gp = gpar(fontsize = 34, fontface = "bold"),
          column_title = "ASVs", column_title_gp = gpar(fontsize = 34, fontface = "bold"),
          row_names_gp = gpar(fontsize = 13, fontface = "bold"),
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
    grid.rect(gp = gpar(col = cluster_palette[i], lwd = 3.2, fill = NULL),
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
  
draw_wk_complex_hm()

dev.off()

png("OUTPUT/wk_complex_hm.png", width = 31, height = 17, units = "in", res = 300, type = "cairo")
draw_wk_complex_hm() # Calls the function which includes all decorations
dev.off()


######################################
### Pure Dendrogram + Sample Names ###
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

# Decorations are painted into the CURRENT device (not stored in ht), so replay them for each device
draw_dendro_pure <- function() {
  draw(ht)
  # Add colored borders to sample names
  for (i in 1:k) {
    decorate_row_names("dendro_plot", {
      grid.rect(gp = gpar(col = cluster_palette[i], lwd = 3.2, fill = NULL),
                width = unit(11.4, "cm"), x = unit(1, "npc") + unit(-0.5, "cm"), just = "center")
    }, slice = i)
  }
  # Add cluster labels
  for (i in 1:k) {
    decorate_row_names("dendro_plot", {
      grid.text(cl_labels[i], x = unit(1, "npc") + unit(1, "cm"), just = "left", gp = gpar(fontsize = 35, fontface = "bold"))
    }, slice = i)
  }
}

pdf(file = "OUTPUT/dendro_pure.pdf", width = 12, height = 24)
draw_dendro_pure()
dev.off()

png("OUTPUT/dendro_pure.png", width = 12, height = 24, units = "in", res = 300)
draw_dendro_pure()
dev.off()


#############
### TITAN ###

## select offshore surface
off_surf <- meta.wk.small %>% rownames_to_column("Sample") %>% filter(`Distance to coast [km]` > 10, between(`Depth [m]`, 1, 5))
txa_off_surf <- proptab.small.taxan.saman[off_surf$Sample, , drop = FALSE] %>% .[, colSums(. > 0) >= 3, drop = FALSE]

## perform analysis vs salinity
# ti_off_surf <- titan(off_surf$`Salinity [PSU]`, txa_off_surf) # 9 min
# saveRDS(ti_off_surf, "OUTPUT/titan_offshore_surface_salinity.rds")
ti_off_surf <- readRDS("OUTPUT/titan_offshore_surface_salinity.rds")
# look at cp = change point; consider sumz; if sumz is the same as fsumz means threshold is robust to filtering
ti_off_surf$sumz.cp[c("sumz-", "sumz+", "fsumz-", "fsumz+"), ]
range(ti_off_surf$sumz.cp[c("sumz-", "sumz+"), "cp"])
p_dens <- plot_sumz_density(ti_off_surf) # shows the density of bootstrap change points; the peak is the community‑level threshold
# blue - relab decreases with salinity; red - relab increases with salinity
# zenv.cp for each taxa - change pint where the repsonse is strongest
# zscore - magnitude is the strongest
# most taxa are blue aka adapted to low salinity
# plot_taxa_ridges(ti_off_surf)
keep_euk  <- rownames(ti_off_surf$sppmax) %in% asv_in_dotplot &
  grepl("^Eukaryota", rownames(ti_off_surf$sppmax)) & ti_off_surf$sppmax[, "zscore"] >= 5
keep_prok <- rownames(ti_off_surf$sppmax) %in% asv_in_dotplot &
  grepl("^(Bacteria|Archaea)", rownames(ti_off_surf$sppmax)) & ti_off_surf$sppmax[, "zscore"] >= 5
p_euk  <- plot_taxa_ridges({ti <- ti_off_surf; ti$sppmax <- ti$sppmax[keep_euk, , drop = FALSE];  ti$taxa <- ti$taxa[, keep_euk, drop = FALSE];  ti},
                           n_ytaxa = sum(keep_euk),  printspp = FALSE, xlabel = "Salinity")
p_prok <- plot_taxa_ridges({ti <- ti_off_surf; ti$sppmax <- ti$sppmax[keep_prok, , drop = FALSE]; ti$taxa <- ti$taxa[, keep_prok, drop = FALSE]; ti},
                           n_ytaxa = sum(keep_prok), printspp = FALSE, xlabel = "Salinity")
pdf("OUTPUT/titan_offshore_surface_salinity.pdf", width = 12, height = 24); print(p_dens + p_euk + p_prok + plot_layout(ncol = 1)); dev.off()
write.csv(ti_off_surf$sumz.cp, "OUTPUT/titan_offshore_surface_salinity.csv")

# ## perform analysis vs temperature - not necessary b/c here, T and S are indirectly related but, unlike S, T overlaps between Hudson Bay and James Bay
# # ti_off_surf_t <- titan(off_surf$`Temperature [deg C]`, txa_off_surf) # 9 min
# # saveRDS(ti_off_surf_t, "OUTPUT/titan_offshore_surface_temp.rds")
# ti_off_surf_t <- readRDS("OUTPUT/titan_offshore_surface_temp.rds")
# ti_off_surf_t$sumz.cp[c("sumz-", "sumz+", "fsumz-", "fsumz+"), ]
# range(ti_off_surf_t$sumz.cp[c("sumz-", "sumz+"), "cp"])
# p_dens <- plot_sumz_density(ti_off_surf_t) # shows the density of bootstrap change points; the peak is the community‑level threshold
# # plot_taxa_ridges(ti_off_surf_t)
# keep_euk  <- rownames(ti_off_surf_t$sppmax) %in% asv_in_dotplot &
#   grepl("^Eukaryota", rownames(ti_off_surf_t$sppmax)) & ti_off_surf_t$sppmax[, "zscore"] >= 5
# keep_prok <- rownames(ti_off_surf_t$sppmax) %in% asv_in_dotplot &
#   grepl("^(Bacteria|Archaea)", rownames(ti_off_surf_t$sppmax)) & ti_off_surf_t$sppmax[, "zscore"] >= 5
# p_euk  <- plot_taxa_ridges({ti <- ti_off_surf_t; ti$sppmax <- ti$sppmax[keep_euk, , drop = FALSE];  ti$taxa <- ti$taxa[, keep_euk, drop = FALSE];  ti},
#                            n_ytaxa = sum(keep_euk),  printspp = FALSE, xlabel = "Temperature")
# p_prok <- plot_taxa_ridges({ti <- ti_off_surf_t; ti$sppmax <- ti$sppmax[keep_prok, , drop = FALSE]; ti$taxa <- ti$taxa[, keep_prok, drop = FALSE]; ti},
#                            n_ytaxa = sum(keep_prok), printspp = FALSE, xlabel = "Temperature")
# pdf("OUTPUT/titan_offshore_surface_temp.pdf", width = 12, height = 24); print(p_dens + p_euk + p_prok + plot_layout(ncol = 1)); dev.off()
# write.csv(ti_off_surf_t$sumz.cp, "OUTPUT/titan_offshore_surface_temp.csv")

## perform analysis vs latitude
# ti_off_surf_l <- titan(off_surf$`Latitude [deg]`, txa_off_surf) # 9 min
# saveRDS(ti_off_surf_l, "OUTPUT/titan_offshore_surface_lat.rds")
ti_off_surf_l <- readRDS("OUTPUT/titan_offshore_surface_lat.rds")
ti_off_surf_l$sumz.cp[c("sumz-", "sumz+", "fsumz-", "fsumz+"), ]
range(ti_off_surf_l$sumz.cp[c("sumz-", "sumz+"), "cp"])
p_dens <- plot_sumz_density(ti_off_surf_l) # shows the density of bootstrap change points; the peak is the community‑level threshold
# plot_taxa_ridges(ti_off_surf_l)
keep_euk  <- rownames(ti_off_surf_l$sppmax) %in% asv_in_dotplot &
  grepl("^Eukaryota", rownames(ti_off_surf_l$sppmax)) & ti_off_surf_l$sppmax[, "zscore"] >= 5
keep_prok <- rownames(ti_off_surf_l$sppmax) %in% asv_in_dotplot &
  grepl("^(Bacteria|Archaea)", rownames(ti_off_surf_l$sppmax)) & ti_off_surf_l$sppmax[, "zscore"] >= 5
p_euk  <- plot_taxa_ridges({ti <- ti_off_surf_l; ti$sppmax <- ti$sppmax[keep_euk, , drop = FALSE];  ti$taxa <- ti$taxa[, keep_euk, drop = FALSE];  ti},
                           n_ytaxa = sum(keep_euk),  printspp = FALSE, xlabel = "Latitude")
p_prok <- plot_taxa_ridges({ti <- ti_off_surf_l; ti$sppmax <- ti$sppmax[keep_prok, , drop = FALSE]; ti$taxa <- ti$taxa[, keep_prok, drop = FALSE]; ti},
                           n_ytaxa = sum(keep_prok), printspp = FALSE, xlabel = "Latitude")
pdf("OUTPUT/titan_offshore_surface_lat.pdf", width = 12, height = 24); print(p_dens + p_euk + p_prok + plot_layout(ncol = 1)); dev.off()
write.csv(ti_off_surf_l$sumz.cp, "OUTPUT/titan_offshore_surface_lat.csv")


####################################
### Check if FT was contaminated ###
# Indicator ASVs for FT in offshore surface samples: NONE
group <- ifelse(off_surf$Code == "FT", "FT", "nonFT")
ind <- multipatt(txa_off_surf, group, func = "IndVal.g", control = how(nperm = 999))
plot(ind$A[ind$sign$p.value < 0.01, "FT"], ind$B[ind$sign$p.value < 0.01, "FT"])
ft_idx <- match("FT", colnames(ind$A))                                                          # only accept ASVs whose permutation-tested pattern IS the FT column...
keep <- !is.na(ind$sign$index) & ind$sign$index == ft_idx & ind$sign$p.value < 0.01 &            # ...not one whose p-value belongs to a nonFT/combined pattern
  ind$A[cbind(seq_len(nrow(ind$A)), ind$sign$index)] >= 0.6 & ind$B[cbind(seq_len(nrow(ind$A)), ind$sign$index)] >= 0.6
ft_asvs <- rownames(ind$A)[keep]
tryCatch(cbind(FT = colMeans(txa_off_surf[group == "FT", ft_asvs]), nonFT = colMeans(txa_off_surf[group == "nonFT", ft_asvs])), error = function(e) "No FT indicators")


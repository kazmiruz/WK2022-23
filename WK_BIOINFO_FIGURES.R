# for simple heatmaps
library(gplots)

library(viridis)
library(viridisLite)  # for turbo coloring

library(vegan)        # for RDA, standardization, etc.

library(dendextend)   # for cutting tree

library("ComplexHeatmap") # for complexheatmap

library(ggrepel)      # for labeling with geom_text_repel
library(ggforce)      # for geom_mark_ellipse

library(gridExtra)    # for grid.draw and grid.table
library(GGally)       # for ggpairs

library(ggtree)       # for dendrogram drawing
library(ape)          # for as.phylo

# build William Kennedy ASV by running WK_BIOINFO_ASV
# source("WK_BIOINFO_ASV.R")

#
### filter data

################################################################################

#
### What does your ASV table look like in terms of data distribution? What normalization should you apply?
# Apply one of these transformation for optimal clustering, choose based on visual inspection of the elbow/silhouette and heatmaps
hist(as.matrix(proptab.small.taxan.saman))        # looks pretty bad
hist(as.matrix(proptab.small.taxan.saman)^.5)     # looks better
hist(as.matrix(proptab.small.taxan.saman)^.25)    # looks good: a lot of 0s but normal distribution to the left
hist(as.matrix(proptab.small.taxan.saman)^.2)     # looks good
hist(as.matrix(proptab.small.taxan.saman)^.125)   # looks good


#############################################################
### Statistically discover the optimal number of clusters ###

kr <- 3:20

# --------------------------------------------
# Step 1: Compute Bray-Curtis distance matrix
# --------------------------------------------
bray_dist <- vegdist(as.matrix(proptab.small.taxan.saman)^0.2, method = "bray")

# Convert to matrix for calculations
bray_matrix <- as.matrix(bray_dist)

# --------------------------------------------
# Step 2: Hierarchical clustering with ward.D2
# --------------------------------------------
hclust_result <- hclust(bray_dist, method = "ward.D2")

# --------------------------------------------
# Step 3: Elbow Plot (WCSS)
# --------------------------------------------
# Calculate Within-Cluster Sum of Squares (WCSS)
wcss <- sapply(kr, function(k) {
  clusters <- cutree(hclust_result, k)
  sum(sapply(1:k, function(i) {
    cluster_indices <- which(clusters == i)
    if (length(cluster_indices) > 0) {
      sum(bray_matrix[cluster_indices, cluster_indices]^2) / (2 * length(cluster_indices))
    } else 0
  }))
})

# Plot elbow curve
elbow_plot <- ggplot(data.frame(k = kr, wcss = wcss), aes(k, wcss)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "steelblue", size = 3) +
  labs(title = "Elbow Plot for Optimal Clusters",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares (WCSS)") +
  theme_minimal()

# --------------------------------------------
# Step 4: Silhouette Analysis
# --------------------------------------------
# Calculate average silhouette width
library(cluster)
sil_width <- sapply(kr, function(k) {
  clusters <- cutree(hclust_result, k)
  sil <- silhouette(clusters, bray_dist)
  mean(sil[, 3])
})

# Plot silhouette widths
silhouette_plot <- ggplot(data.frame(k = kr, sil = sil_width), aes(k, sil)) +
  geom_line(color = "firebrick", size = 1.2) +
  geom_point(color = "firebrick", size = 3) +
  labs(title = "Silhouette Analysis",
       x = "Number of Clusters (k)",
       y = "Average Silhouette Width") +
  theme_minimal()

# --------------------------------------------
# Combine Plots and Print out optimal number of clusters
# --------------------------------------------
library(patchwork)
combined_plot <- (elbow_plot | silhouette_plot)
print(combined_plot)


######################
### Simple Heatmap ### 
# build a simple heatmap to explore all the data and visually assess how many clusters you want

## open pdf
pdf(file="OUTPUT/wk_simple_hm.pdf", width=35, height=20)

## simple heatmap
hm = heatmap.2(as.matrix(proptab.small.taxan.saman)^.2, trace="none", col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"))

## see how many clusters you may have for samples (rows)
clusnum=7
oldclus = cutree(hm$rowDendrogram, k=clusnum)
oldorder = unname(rle(oldclus[as.hclust(hm$rowDendrogram)$order])$values)
neworder = (1:clusnum)
names(neworder) = oldorder
clus = unname(neworder[as.character(oldclus)])
names(clus) = names(oldclus)

## redraw the map with clusters colored
hm = heatmap.2(as.matrix(proptab.small.taxan.saman)^.2, trace="none", key=F, margins=c(20,30), col=viridis,
               distfun=function(x) vegdist(x, method="bray"),
               hclustfun=function(x) hclust(x, method="ward.D2"),
               RowSideColors = rainbow(max(clus))[clus])

## close the pdf
dev.off()


#################################
### Cluster # + color palette ###
k = 7   # number of clusters you decided on (here, 5 is optimal according to elbow/silhouette)
        # but one of the clusters is split into two b/c too big and has a distinguishing microbial block

# cluster_palette <- setNames(turbo(k), as.character(1:k))   # how you will color your clusters from now on
# # turbo: an improved Jet rainbow color map for reducing false detail, banding and color blindness ambiguity

# cluster_palette <- c("#C832FF", "#002AFF", "#00FFFF", "#00FF00", "#FFD700", "#FF1900") %>%           # custom palette for 6 clusters
#   setNames(as.character(1:k))

cluster_palette <- c("#F233F2", "#0000F2", "#00F2F2", "#00F200", "#A6734D", "#D9D900", "#F20000") %>%  # custom palette for 7 clusters
  setNames(as.character(1:k))


######################
### ComplexHeatmap ###

colnames(meta.wk.small) # see which meta parameters you want to investigate

# Compute Bray-Curtis distance and hierarchical clustering
bray_dist <- function(x) vegdist(x, method = "bray")

tdist <- vegdist(as.matrix(proptab.small.taxan.saman)^.2, method = "bray")
ttree <- hclust(tdist, method = "ward.D2")

# Make a dendrogram variable and color branches with Turbo
row_dend <- as.dendrogram(ttree) %>%
  dendextend::color_branches(., k = k, col = cluster_palette)

#
## open pdf
pdf(file="OUTPUT/wk_complex_hm.pdf", width=33, height=20)

wk_complex_hm <- function() {

ht_opt$ROW_ANNO_PADDING = unit(2, "mm")                                          # Padding between the heatmap and annotations

# Create heatmap with dendrogram
wk_complex_hm <- Heatmap(as.matrix(proptab.small.taxan.saman)^.2,
                           name = "wk_complex_hm",
                           col = viridis(256),
                           
                           # Row clustering with pre-computed dendrogram
                           cluster_rows = row_dend %>%
                             set("branches_lwd", 3),                             # set branch line width for better visibility
                           row_split = k,
                           
                           # Column clustering
                           clustering_distance_columns = bray_dist,
                           clustering_method_columns = "ward.D2",
                           column_dend_reorder = TRUE,
                           
                           # Annotations
                           right_annotation = rowAnnotation(
                             bar1 = anno_points(log10(meta.wk.small$"Depth [m]"),
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(
                                                  at = log10(c(1, 10, 100)),
                                                  labels = c(1, 10, 100),
                                                  gp = gpar(fontsize = 10))),
                             bar2 = anno_points(meta.wk.small$"Temperature [deg C]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar3 = anno_points(meta.wk.small$"Salinity [PSU]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar4 = anno_points(meta.wk.small$"CDOM [mg/m^3]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar5 = anno_points(meta.wk.small$"Fluorescence [mg/m^3]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"),
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             bar6 = anno_points(meta.wk.small$"Oxygen [umol/kg]",
                                                size = unit(3, "mm"),
                                                width = unit(3, "cm"), 
                                                axis_param = list(gp = gpar(fontsize = 10))),
                             gap = unit(2, "mm"),
                             show_annotation_name = FALSE
                           ),
                           
                           # Rest of your parameters
                           heatmap_height = unit(40, "cm"),
                           heatmap_width = unit(65, "cm"),
                           
                           row_dend_width = unit(6, "cm"),
                           column_dend_height = unit(2, "cm"),
                           
                           row_title = "Samples",
                           row_title_gp = gpar(fontsize = 25, fontface = "bold"),
                           
                           column_title = "Taxa",
                           column_title_gp = gpar(fontsize = 25, fontface = "bold"),
                           
                           row_names_gp = gpar(fontsize = 11, fontface = "bold"),
                           
                           column_names_gp = gpar(fontsize = 1.2, fontface = "bold"), 
                           # column_names_rot = 45,
                           
                           heatmap_legend_param = list(
                             title = "Fifth-root transformed relative abundances of microbial taxa during the 2022 and 2023 James Bay expeditions\n",
                             title_gp = gpar(fontsize = 30, fontface = "bold"),
                             legend_direction = "horizontal",
                             legend_width = unit(2, "cm"),
                             color_bar = "continuous")
)

#
## Draw the complex heatmap
draw(wk_complex_hm, heatmap_legend_side = "top")

#
## Add decorations

# Fix annotation bars
annotation_labels <- c( "Depth",                                              # Define parameter labels in order (bar1, bar2, etc.)
                        "Temp", "Sal",
                        "CDOM", "Fluo", "O2")
for (i in 1:6) {
  decorate_annotation(paste0("bar", i), {
    grid.text(annotation_labels[i],
              y = unit(1, "npc") + unit(2, "mm"), just = "bottom",               # Position label above the annotation
              gp = gpar(fontsize = 17, fontface = "bold"))
  })
}

# Color the cluster row labels (samples) from top/first slice (later named as "Cluster 1") to the bottom/last (later named"Cluster k")
for (i in 1:k) {                                    # loop over cluster slices
  decorate_row_names("wk_complex_hm", {
    grid.rect(gp = gpar(fill = alpha(cluster_palette[i], 0.5)),    # transparency set at 0.5
              width = unit(13.3, "cm"),
              x = unit(1, "npc") + unit(0.6, "cm"),
              just = "center")
  }, slice = i)
}
}

wk_complex_hm()

## close the pdf
dev.off()


########################################
### Cluster order + add to metadata ####
# reorder clusters in accordance with the heatmap order
tclus_reordered <- cutree(ttree, k = k) %>%
  setNames(rownames(proptab.small.taxan.saman)) %>%  # Assign sample names
  (\(tclus) tclus %>%
     factor(levels = unique(tclus[order.dendrogram(row_dend)]), labels = 1:k) %>%
     as.numeric() %>%
     setNames(names(tclus)))()

meta.wk.small$Cluster <- tclus_reordered[match(meta.wk.small$plate.geno, sub(".*; ", "", names(tclus_reordered)))]


#
## quick exploration
# plot(y = meta.wk.small$`Latitude [deg]`, x = (meta.wk.small$`Longitude [deg]`*-1), col = cluster_palette[meta.wk.small$Cluster])
# plot(y = meta.wk.small$`Salinity [PSU]`, x = meta.wk.small$`Temperature [deg C]`, col = cluster_palette[meta.wk.small$Cluster],
#      pch = as.numeric(as.factor(meta.wk.small$Year)))
# plot(x = meta.wk.small$`Longitude [deg]`, y = meta.wk.small$`Latitude [deg]`, col = cluster_palette[meta.wk.small$Cluster],
#      pch = as.numeric(as.factor(meta.wk.small$Year)))
# plot(meta.wk.small$`Fluorescence [mg/m^3]`, meta.wk.small$Cluster)


###########################
### Circular Dendrogram ###

# Create cluster data frame for annotation
cluster_data <- data.frame(label = names(tclus_reordered),
                           cluster = tclus_reordered,
                           stringsAsFactors = FALSE)
phylo <- as.phylo(ttree)

# Open PDF with optimized dimensions
pdf(file = "OUTPUT/wk_circ_dend.pdf", width = 8, height = 8)

# Create base circular tree with compact layout
p <- ggtree(phylo, layout = 'circular', ladderize = TRUE, size = 0.8) %<+% cluster_data +
  geom_tree(aes(color = factor(cluster))) +                           # color branches by cluster
  geom_tiplab(aes(color = factor(cluster)),                           # label spacing
              size = 2,                        # label size
              offset = 0.5,                      # pushing labels outward
              align = TRUE,                      # prevent overlap
              show.legend = TRUE) +
  scale_color_manual(name = "Cluster", values = cluster_palette,
                     guide = guide_legend(override.aes = list(size = 5))) +     # your predefined colors for clusters
  expand_limits(x = 5)                           # expands - smaller number means greater exppansion

print(p)
dev.off()


#####################################################################################################################
### db-RDA (distance-based Redundancy Analysis; in this case, we will base it on the Bray Curtis Distance matrix) ###
# to identify how environmental variables influence microbial relative abundance
# the vector of the arrow points in the direction of its greatest influence
# PCA is an unconstrained form of ordination - searches for any variables best explaining the species composition
# RDA is constrained - searches for the best explanatory (environmental) variables

# Info that helped me to put together a db-RDA pipeline below:
# https://youtu.be/Wqa7_18jf5A
# https://youtu.be/_bpQ16uUqs4
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# https://wiki.qcbs.ca/r_workshop10
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/rda-and-dbrda/
# https://sites.ualberta.ca/~ahamann/teaching/renr690/seminars/dbRDA.pdf
# https://vegandevs.github.io/vegan/reference/goodness.cca.html

#
## Select and transform microbial and environmental 

# Load microbial ASV data (response variables)
microbial_data <- as.matrix(proptab.small.taxan.saman)^0.2

# Check how your ASV data looks like
# If the distribution is not normal - go back and modify.
# Here we have already done prop.table (by sample) and ^0.25 transformations on the original ASV matrix
hist(microbial_data)    # looks optimal: a lot of 0s but normal distribution to the left

# Load metadata (explanatory variables)
metadata <- meta.wk.small

# Select numeric environmental variables correctly
env_vars <- metadata %>%
  select ("Depth [m]", 
          "Salinity [PSU]", "Temperature [deg C]",
          # "Oxygen [umol/kg]",                                     # Most oxygen is NA
          # "CDOM [mg/m^3]", "Fluorescence [mg/m^3]",               # CDOM strongly correlates with salinity and Fluo is not significant
          # "Beam Attenuation [1/m]",                               # Useless if we already have Fluo, and too few values
          "Distance to coast [km]") %>%
  rename_with(~ str_replace_all(., "\\s*\\[.*?\\]", "")) %>%      # Remove brackets
  rename_with(~ gsub(" ", "_", .))                                # Replace spaces with underscores

# Identify the variables that strongly correlate with each other according to Pearson cor. coef. (R >= ~0.9) - remove one of them
# If you leave two strongly correlating environmental variables, it does not add much to the total variation explained b/c 
# the corr. var-le appearing in the model formula first already explains the variance the other corr. var-le also tries to explain
# It leads to inflated variance inflation factors (VIF) for both corr. var-les: it is hard to know which var-le drives the variation
ggpairs(env_vars, 
        lower = list(continuous = wrap("smooth", method = "lm", se = TRUE, color = "blue", size = 0.8)), 
        diag = list(continuous = wrap("barDiag", bins = 20, color = "cyan")),  # Histogram on the diagonal
        upper = list(continuous = wrap("cor", size = 3, color = "black")))     # Correlation in the upper panel

# Leave only the data where all parameters are non-NA (complete cases), unless you can find a way to harmlessly interpolate
# Check if the data was evenly removed
complete_cases <- complete.cases(env_vars)                           # 70/72
microbial_data <- microbial_data[complete_cases, , drop = FALSE]
env_vars <- env_vars[complete_cases, , drop = FALSE]

# The environmental/explanatory variables are all in different units and need to be standardized prior to computing distance measures to perform most ordination analyses
env_stand <- decostand(env_vars, method="standardize")  # subtract the mean and divide by the standard deviation
# this centers each variable around 0 (mean = 0) and scales it so that its standard deviation is 1
apply(env_stand, 2, mean) # check if the data are now centered (means~0)
apply(env_stand, 2, sd)   # check if the data are now scaled (standard deviations=1)

#
## Build db-RDA model and remove/add environmental variables above + rerun, as needed based on the results of the checks below

# Before applying db-RDA, you can check which dissimilarity index works best for you
# Usually, we choose Bray anyway because it works best with microbial data, calculates distances better
rankindex(env_stand, microbial_data, indices = c("euc", "man", "gow", "bra", "kul"),
          stepacross= FALSE, method = "pearson")

# Perform db-RDA (distance-based redundancy analysis) with all variables
rda_model <- dbrda(microbial_data ~ ., data = env_stand, distance = "bray", dfun = vegdist)

# Check if any variables significantly correlate (VIF>10~20) in your constrained ordination (e.g., RDA)
# If yes, remove the ones that makes sense to remove above and rerun
vif.cca(rda_model) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%  # Convert row names to a proper column
  dplyr::rename(VIF = ".") %>%               # Rename the unnamed column to "VIF"
  mutate(VIF = round(VIF, 0)) %>%     # Round to 0 decimal places
  arrange(desc(VIF)) %>%              # Sort in descending order for easy review)
  print()

# Perform ANOVA checks for the model, environmental variables an axes.
# anova.cca performs an ANOVA like permutation test for cca, rda, or dbRDA to assess the significance of constraints
# Also, extract and look at environmental variable scores, % variance explained, R²

# Print significance of the model (has to be p < 0.001)
anova.cca(rda_model, permutations = 1000)

# Extract scores for environmental variables (arrows)
# Scaling = 1 if interested in clustering of samples (distance between points)
# Scaling = 2 if interested in correlations between env. variables (arrows) and how they influence community structure
# Look at the variables that impact each of your dbRDA axis
# interpret how environmental variables influence community structure
sc_bp <- scores(rda_model, display = "bp", choices = c(1:ncol(env_stand)), scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  print()

# Print significance of your environmental variables
# Leave only significant environmental variables (high F value and p < 0.01** or p < 0.001***), remove the rest above and rerun
# High F value: the values within each group do not differ significantly from one another, while those across different groups do
anova_per_var <- anova.cca(rda_model, permutations = 1000, by = "term") %>%
  rownames_to_column("Variable") %>%                        # Convert rownames to a column
  select(Variable, F, `Pr(>F)`) %>%                         # Keep relevant columns
  dplyr::rename(`P-value` = `Pr(>F)`) %>%                   # Rename Pr(>F) to p-value
  filter(Variable != "Residual") %>%                        # Remove "Residual" row
  mutate(
    F = round(F, 2),  
    `P-value` = sprintf("%.6f", as.numeric(`P-value`)),     # Ensure 6 decimal places
    `Significance code` = case_when(
      as.numeric(`P-value`) <= 0.001 ~ "***",
      as.numeric(`P-value`) <= 0.01  ~ "**",
      as.numeric(`P-value`) <= 0.05  ~ "*",
      as.numeric(`P-value`) <= 0.1   ~ ".",
      TRUE                           ~ "")) %>%             # Assign stars based on p-value significance
  arrange(desc(F)) %>%                                      # Sort by descending F-value
  as.data.frame() %>%                                       # Convert to data frame
  print()                                                   # Explicitly print the table

# Print marginal effects of terms: how well each variable adds to explaining the total variance
# Yet another check as to whether or not you should remove some variables due to their strong correlation
# While "terms" assesses the sequential effects of var-les (order matters), "margin" evaluates each var-le independently
anova.cca(rda_model, permutations = 1000, by = "margin") %>%
  arrange(desc(F))                                          # Sort by descending F-value

# Print significance of the first two dbRDA axis (has to be p < 0.001)
anova.cca(rda_model, permutations = 1000, by = "axis")

# Print variance explained by axes
perc <- round(100 * (summary(rda_model)$cont$importance[2, 1:2]), 2) %>%
  print()

# Print proportion of variance explained by environmental variables
variance_partitioning <- (summary(rda_model)$constr.chi / summary(rda_model)$tot.chi) %>%
  {paste("variance explained:", round(. * 100, 2), "%")} %>%
  print()

# Print final adjusted R²
rda_adj_r2 <- RsquareAdj(rda_model)$adj.r.squared %>%
  {paste("adjusted R²:", round(. * 100, 1), "%")} %>%
  print()

#
## Prepare data for plotting

# Add ALL metadata parameters to plot by size
# metadata_filtered <- meta.wk.small %>%
#   mutate(Cluster = as.factor(tclus_reordered[match(ordered_samples, names(tclus_reordered))]))    # assign corrected clusters

# Extract RDA scores for plotting and add metadata to it
sc_si <- scores(rda_model, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  mutate(plate.geno = sapply(strsplit(Sample, "; "), tail, 1)) %>%
  left_join(metadata, by = "plate.geno") %>%
  mutate(Cluster = as.factor(Cluster))

# Define cluster colors
cluster_levels <- levels(as.factor(metadata$Cluster))  # get unique cluster numbers
cluster_colors <- cluster_palette[cluster_levels]      # assign cluster colors

# Compute dynamic arrow scaling factor
arrow_scaling_factor <- max(abs(c(sc_bp$dbRDA1, sc_bp$dbRDA2)), na.rm = TRUE) * 2

# Define proper axis limits
xlim_vals <- range(c(sc_si$dbRDA1, sc_bp$dbRDA1 * arrow_scaling_factor), na.rm = TRUE) * 1.3
ylim_vals <- range(c(sc_si$dbRDA2, sc_bp$dbRDA2 * arrow_scaling_factor), na.rm = TRUE) * 1.3

## Open PDF
pdf(file="OUTPUT/wk_db-RDA.pdf", width = 35, height = 20)

# Loop through all numeric variables in `sc_si` to generate additional db-RDA plots
numeric_vars <- meta.wk.small %>% select(where(is.numeric), -Cluster) %>% colnames()  # size by different variables

for (param in numeric_vars) {
  
  param_min <- min(sc_si[[param]], na.rm = TRUE)
  param_max <- max(sc_si[[param]], na.rm = TRUE)
  
  param_plot <- ggplot() +
    
    geom_point(
      data = sc_si,
      aes(x = dbRDA1, y = dbRDA2,
          color = Cluster,
          size = !!sym(param)),
      alpha = 0.5, stroke = 3  # Strong outline for better visibility
    ) +

    geom_segment(data = sc_bp, 
                 aes(x = 0, y = 0, xend = dbRDA1 * arrow_scaling_factor, yend = dbRDA2 * arrow_scaling_factor), 
                 arrow = arrow(length = unit(0.5, "cm")), size = 1.5, color = "black") +
    
    geom_text_repel(data = sc_bp, 
                    aes(x = dbRDA1 * arrow_scaling_factor, 
                        y = dbRDA2 * arrow_scaling_factor, 
                        label = Variable),
                    color = "black", size = 10, fontface = 'bold',
                    max.overlaps = Inf, box.padding = 1.5, point.padding = 0.7,
                    segment.color = "grey50") + 
    
    scale_color_manual(name = "Cluster", values = cluster_colors) +
    scale_size(range = scales::rescale(range(sc_si[[param]], na.rm = TRUE), to = c(2, 20))) +
    
    # Create custom legend for filtration type
    guides(
      color = guide_legend(
        # title = "Cluster",
        override.aes = list(size = 15)),
      size = guide_legend(order = 2)) +
    
    labs(title = paste0(
      "db-RDA plot for the 2022 and 2023 James Bay expeditions: ", param,
      "\n(", variance_partitioning, "; ", rda_adj_r2, ")"),
      x = paste0("dbRDA1 (", perc[1], "%)"),
      y = paste0("dbRDA2 (", perc[2], "%)")) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.text = element_text(size = 16),
      axis.title.x = element_text(size = 35, face = "bold"),
      axis.title.y = element_text(size = 35, face = "bold"),
      plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 35, face = "bold"),
      legend.key.size = unit(2, "cm"),    # Reduce legend key size
      legend.spacing.y = unit(0.1, "cm")  # Reduce vertical spacing between legend items
    )
  
  print(param_plot)
}

# Print a table with ANOVA results for the significance of environmental variables
grid.newpage()
grid.draw(tableGrob(anova_per_var, rows = NULL,
                    theme = ttheme_minimal(
                      base_size = 60,
                      padding = unit(c(10, 10), "mm"),
                      core = list(bg_params = list(fill = c("white", "grey95"), col = "black")), 
                      colhead = list(bg_params = list(fill = "grey80", col = "black")))))

# Close PDF
dev.off()


############################################
### Alpha-diversity (Shannon) by cluster ###
# Shannon index: richness and evenness
# The Shannon index formula inherently normalizes counts into proportions internally
# Hence, you can use your relative abundances instead of raw ASV values - it makes no difference

meta.wk.small$`Shannon diversity` <- vegan::diversity(proptab.small.taxan.saman, index = "shannon")

# Test for differences in richness as a function of cluster
# The Kruskal-Wallis test can be thought of as the non-parametric equivalent to ANOVA
kruskal <- kruskal.test(meta.wk.small$`Shannon diversity` ~ meta.wk.small$Cluster) %>% print()

shannon_fig <- ggplot(meta.wk.small, aes(x = factor(Cluster), y = `Shannon diversity`, fill = factor(Cluster))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  
  scale_fill_manual(values = cluster_palette) +    # use predefined cluster_palette colors
  labs(title = paste0("Shannon diversity by Cluster (Kruskal-Wallis test p-value = ", formatC(kruskal$p.value, format = "e", digits = 2), ")"),
       x = "Cluster",
       y = "Shannon Diversity") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 35, face = "bold"),
    axis.title.y = element_text(size = 35, face = "bold"),
    plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
    legend.position = "none"
  )


###################################
### Environmental cluster stats ###
# look at the parameters found to be significant as a result of the db-RDA analysis in defining microbial community clusters
# can include some additional parameters, for example if they were excluded due to unavailability for a lot of samples
# make a summary table for clusters

# Define the columns to summarize, also give NA to them if >50% of the values are NAs
cols_to_check <- c("Depth [m]", "Distance to coast [km]",
                   "Temperature [deg C]", "Salinity [PSU]",
                   "CDOM [mg/m^3]", "Fluorescence [mg/m^3]",
                   "Shannon diversity")

cluster_summary <- proptab.small.taxan.saman %>%
  as.data.frame() %>%
  cbind(meta.wk.small) %>%
  mutate(Cluster = as.character(Cluster)) %>%   # Convert Cluster to character
  group_by(Cluster) %>%
  summarize("Code" = names(which.max(table(Code))),
            across(all_of(cols_to_check), 
                   ~ if (mean(is.na(.x)) > 0.5) NA_character_                                    # 50% or more of data should be present
                   else paste0(round(mean(.x, na.rm = TRUE), 1), 
                               " ± ", 
                               round(sd(.x, na.rm = TRUE), 1)),
                   .names = "{.col}")) %>%
  bind_rows(
    proptab.small.taxan.saman %>%
      as.data.frame() %>%
      cbind(meta.wk.small) %>%
      summarize(
        "Cluster" = "Kruskal-Wallis p",   # Will match type now
        "Code" = NA,
        across(all_of(cols_to_check), 
               ~ if (mean(is.na(.x)) > 0.5) NA_character_                                        # 50% or more of data should be present
               else formatC(kruskal.test(.x ~ meta.wk.small$Cluster)$p.value, format = "e", digits = 2),
               .names = "{.col}")
      )
  ) %>%
  setNames(str_replace_all(names(.), "(.*?) \\[(.*?)\\]", "\\1\n(\\2)"))

print(cluster_summary)


########################
### Cluster analysis ###

# Compute Spearman correlation matrix between ASVs
# Spearman is useful because it uses ranks rather than absolute values
adjacency_matrix <- cor(proptab.small.taxan.saman, method = "spearman")
adjacency_matrix[abs(adjacency_matrix) < 0.6] <- 0            # Replace weak correlations with zero (set threshold)
diag(adjacency_matrix) <- 0                                   # Remove self-connections

# Create a network graph object
library(igraph)
network <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE)
E(network)$width <- E(network)$weight * 0.5                   # Set edge widths

library(Rtsne)
network_layout <- as.matrix(Rtsne(t(proptab.small.taxan.saman)^0.2, perplexity=30, dims=2, check_duplicates=FALSE, 
                                  pca=TRUE, max_iter=2500, verbose=FALSE)$Y)

# Aggregate ASV abundances by cluster
asv_clus <- aggregate(proptab.small.taxan.saman, by=list(meta.wk.small$Cluster), FUN = mean) # mean ASV abundance by cluster
rownames(asv_clus) = asv_clus$Group.1                                       # row names set to cluster labels
asv_clus$Group.1 = NULL

# For each ASV, determine in which cluster it is the most dominant (aka the highest mean abundance)
# This helps determine which cluster each ASV is most associated with - the nodes (ASVs) are colored accordingly later
asv_dominant_clus <- rownames(asv_clus)[apply(asv_clus, 2, which.max)]

# Define colors for cluster
vertex_colors <- cluster_palette[asv_dominant_clus]                      # assign cluster colors
vertex_colors[is.na(vertex_colors)] <- "gray"                            # default color for unassigned ASVs

# Out of the ASVs most associated with each cluster, extract those with the highest relative abundant (e.g., >1%)
# Notice: this is different from the most abundant ASVs per cluster!
# Some ASVs may be very abundant in the cluster but NOT associated with it as they are more abundant in other clusters
# Such ASVs may have bigger nods than Cluster-spec. ASVs, but the color is dif. - it is of the cluster they are more associated w/
top_asv_clus <- asv_clus %>%                                             # mean ASV abundance by cluster
  as.data.frame() %>%
  rownames_to_column("Cluster") %>%                                      # cluster as column
  pivot_longer(-Cluster, names_to = "ASV", values_to = "Abundance") %>%  # 2 columns: ASV and its mean Abundance per cluster
  left_join(
    tibble(ASV = colnames(asv_clus),                                     # a tibble mapping ASV names to their dominant cluster
           Cluster = as.character(asv_dominant_clus),                    # assigns dominant cluster for each ASV
           ASV_Number = gsub(".*; ASV_", "ASV_", colnames(asv_clus))),   # extract ASV numbers from names
    by = "ASV"                                                           
  ) %>%
  filter(Cluster.x == Cluster.y,                                         # keep only ASVs that are assigned to the correct cluster
         !grepl("Unknown|Domain", ASV),                                  # remove ASVs with "Unknown" or "Domain" in the name
         Abundance >= 0.02) %>%                                          # retain only ASVs with at least 2% relative abundance
  select(Cluster = Cluster.x, ASV, Abundance) %>%                        # remove intermediate columns
  mutate(Abundance = signif(Abundance, 2)) %>%                           # round to 2 significant figures
  arrange(Cluster, desc(Abundance)) %>%                                  # sort clusters & abundances in order
  print()

#
## Open PDF
pdf(file="OUTPUT/wk_cluster_analysis.pdf", width = 35, height = 20)

#
## Regenerate the complex heatmap
wk_complex_hm()

#
## Generate the network map

# Loop through clusters 1 to k
# Circle/nod size - relative abundance of the ASV
# All ASVs are included in the layout, regardless of their primary cluster
# If sizable nods of other colors appear in the ASV network graph for a specific cluster, it is due to cross-cluster relationships
# If the cluster nodes are closely connected, it is highly self-contained, but if dif. colors are mixed - cross-cluster association

for (i in 1:k) {
  node_size <- sqrt(asv_clus[as.character(i),])

  # Extract only part of the taxonomic name
  # asv_labels <- gsub(".*; (ASV_\\d+)", "\\1", colnames(asv_clus))          # Extracts ASV_#"
  asv_labels <- gsub(".*; ([^;]+); ASV_\\d+$", "\\1", colnames(asv_clus))    # Extract Genus

  plot(network,
       layout = network_layout,
       vertex.size = node_size * 30,
       vertex.color = vertex_colors,
       # vertex.frame.color = NA,
       edge.width = E(network)$width,
       edge.color = "gray50",
       vertex.label = asv_labels,      # show only part of the name
       vertex.label.cex = 0.2)         # make the name smaller

  title(main = paste("Square-root transformed relative abundances of microbial taxa in Cluster", i, "(ASVs are colored by Cluster they are most associated with)"), cex.main = 3)
}

#
## Add alpha-diversity by Cluster (Shannon index)
print(shannon_fig)

#
## Add the table with the most abundant ASVs specific to each cluster
grid.arrange(textGrob("Most abundant identified Cluster-specific ASVs\n(>2% mean relative abundance)", 
                      gp = gpar(fontsize = 30, fontface = "bold")),
             grid.table(top_asv_clus, 
                        rows = NULL,                                                        # remove row numbers
                        theme = ttheme_default(
                          core = list(fg_params = list(cex = 2)),                           # cell text size
                          colhead = list(fg_params = list(cex = 2, fontface = "bold")))),   # column label size
             nrow = 2, heights = c(0.2, 1))


# Add the table with environmental parameter values specific to each cluster
grid.arrange(textGrob("Cluster-specific values of defining environmental parameters",
                      gp = gpar(fontsize = 35, fontface = "bold"),
                      vjust = 15),                                                          # move down the title
             grid.table(cluster_summary,
                        rows = NULL,                                                        # remove row numbers
                        theme = ttheme_default(
                          core = list(fg_params = list(cex = 2.5)),                           # cell text size
                          colhead = list(fg_params = list(cex = 2.5, fontface = "bold")))),   # column label size
             nrow = 2, heights = c(0.2, 1))

## close PDF
dev.off()


#########################################
#
## Facet plot

# --- Prepare long-format data ---
cluster_plot_data <- proptab.small.taxan.saman %>%
  as.data.frame() %>%
  cbind(meta.wk.small) %>%
  select(Cluster, all_of(cols_to_check)) %>%
  pivot_longer(cols = -Cluster, names_to = "Parameter", values_to = "Value") %>%
  mutate(
    Cluster = as.factor(Cluster),
    Parameter = factor(Parameter, levels = cols_to_check)  # Enforce row order
  )

# --- Calculate Kruskal-Wallis p-values ---
param_labels <- cluster_plot_data %>%
  group_by(Parameter) %>%
  summarize(p_value = formatC(kruskal.test(Value ~ Cluster)$p.value, format = "e", digits = 2)) %>%
  mutate(facet_label = paste0(Parameter, "\n(Kruskal p = ", p_value, ")"))

# --- Merge labels ---
cluster_plot_data <- cluster_plot_data %>%
  left_join(param_labels, by = "Parameter") %>%
  mutate(facet_label = factor(facet_label, levels = param_labels$facet_label))  # Explicit control!

# --- Plot safely inside a function ---
plot_cluster_params <- function() {
  ggplot(cluster_plot_data, aes(x = Cluster, y = Value)) +
    geom_boxplot(fill = NA, outlier.shape = NA, color = "black") +
    geom_jitter(aes(color = Cluster), width = 0.2, alpha = 0.8, size = 2) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = cluster_palette) +
    labs(x = "Cluster", y = NULL, title = "Parameter Distribution by Cluster") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      strip.text = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
      panel.spacing = unit(1, "cm"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

# --- Save PDF ---
pdf(file = "OUTPUT/wk_cluster_param.pdf", width = 18, height = 27)
print(plot_cluster_params())   # Print ggplot inside PDF safely
dev.off()


#########################
### Indicator species ###

#
### Determine indicator species for your clusters ###

# Extract sample names (plate.geno part)
asv_sample_names <- rownames(proptab.small.taxan.saman)
asv_sample_names_clean <- trimws(sub(".*;\\s*", "", asv_sample_names))  # Extract AFTER last semicolon & space

# Rename ASV rownames to match plate.geno
proptab_for_indicator <- proptab.small.taxan.saman
rownames(proptab_for_indicator) <- asv_sample_names_clean

# Filter & align metadata
meta_clean <- meta.wk.small %>%
  filter(plate.geno %in% rownames(proptab_for_indicator)) %>%
  select(plate.geno, Cluster) %>%
  column_to_rownames("plate.geno")

# Filter ASV table
asv_rel_abund <- proptab_for_indicator[rownames(meta_clean), ] %>%
  sweep(1, rowSums(.), "/")  # Convert to relative abundances

# Confirm match
print(table(meta_clean$Cluster))

# Run indicator species analysis
library(indicspecies)
indicator_res <- multipatt(asv_rel_abund, meta_clean$Cluster, control = how(nperm = 999))

# Extract, clean, and filter:
significant_asvs <- indicator_res$sign %>%
  rownames_to_column("ASV") %>%
  filter(stat >= 0.9) %>%                            # Only high indicator value (adjustable)
  filter(!is.na(p.value) & p.value <= 0.001) %>%     # Significant only
  mutate(Assigned_Cluster = apply(select(., starts_with("s.")), 1, 
                                  function(x) paste(which(x == 1), collapse = ","))) %>% 
  select(ASV, Assigned_Cluster, stat, p.value)

# Sort by significance
significant_asvs <- significant_asvs %>%
  arrange(p.value, desc(stat))

# Output to CSV
write.csv(significant_asvs, "OUTPUT/important_indicator_species_boot80.csv", row.names = FALSE)


#
### Aggregate abundant indicator genera and orders for ODV ###

# Function to extract taxonomic levels with filtering rules
extract_taxa <- function(taxonomy, level) {
  parts <- strsplit(taxonomy, ";\\s*")[[1]]
  if (length(parts) >= level && !grepl("Chloroplast|Mitochondria", taxonomy)) {
    phylum <- parts[1]
    order <- parts[2]
    genus <- if (length(parts) >= 3) parts[3] else NA
    if (level == 2) {
      return(paste(phylum, order, sep = "; "))
    } else if (level == 3 && phylum != "Cyanobacteriota") {
      return(paste(phylum, order, genus, sep = "; "))
    }
  }
  return(NA)
}

# Extract significant Genera and Orders
significant_genera <- significant_asvs$ASV %>%
  sapply(extract_taxa, level = 3) %>%
  na.omit() %>%
  unique() %>%
  .[!str_detect(., "ASV_")]
significant_orders <- significant_asvs$ASV %>%
  sapply(extract_taxa, level = 2) %>%
  na.omit() %>%
  unique() %>%
  .[!str_detect(., "ASV_")]

# Extract taxonomic levels for all ASVs in the abundance table
asv_genus <- sapply(colnames(proptab.small.taxan.saman), extract_taxa, level = 3)
asv_order <- sapply(colnames(proptab.small.taxan.saman), extract_taxa, level = 2)

process_abundance <- function(taxa_list, taxa_level) {
  sapply(taxa_list, function(taxa) {
    col_indices <- which(if (taxa_level == "Genus") asv_genus == taxa else asv_order == taxa)
    if (length(col_indices) > 0) {
      rowSums(proptab.small.taxan.saman[, col_indices, drop = FALSE])
    } else {
      rep(0, nrow(proptab.small.taxan.saman))
    }
  }) %>%
    as.data.frame() %>%
    { .[, colSums(. >= 0.1) > 0, drop = FALSE] } %>%                      # Filter for taxa with ≥10% abundance in at least one sample
    rownames_to_column("rowname") %>%
    mutate(plate.geno = sub(".*;\\s*", "", rowname)) %>%
    select(-rowname)
}

# Process Genus and Order in one call each
genus_abundance <- process_abundance(significant_genera, "Genus")
order_abundance <- process_abundance(significant_orders, "Order")

# Merge to metadata
meta.wk.indic <- meta.wk.small %>%
  left_join(genus_abundance, by = "plate.geno") %>%
  left_join(order_abundance, by = "plate.geno")



########################################
### ODV for the James Bay expedition ###

write.csv(meta.wk.indic, "OUTPUT/ODV.csv", row.names = FALSE)



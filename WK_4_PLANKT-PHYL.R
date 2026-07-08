# WK_4---------------------------------------------------------------------
# Purpose:
#   Attach pre‑computed Shannon diversity (can compute them here yourself) to the metadata, build a unified taxonomy lookup,
#   compute average relative abundances per cluster for each major group (usually at the phylum level),
#   and produce a combined figure (Shannon boxplots + domain‑wise stacked bar charts).
#   Statistical tests on Shannon diversity are also performed (Spearman correlations, Kruskal–Wallis, Dunn post‑hoc).
#
# Prerequisites (objects must be in the workspace):
#   proptab.small.taxan.saman   – ASV table with relative abundances (WK_2_ASV.R)
#   meta.wk.small               – metadata with cluster assignments (WK_1_META.R and WK_3_CLUS_TITAN_FT.R)
#   asv16, asv18                – ASV taxonomy tables (WK_2_ASV.R)
#   cl_labels, cluster_palette  – cluster names and colours (WK_3_CLUS_TITAN_FT.R)
#
# Required external files - you can also reproduce them all yourself using the commented out code:
#   OUTPUT/shannon_diversity.rds  – pre‑computed Shannon indices
#
# Outputs:
#   OUTPUT/figure_community_diversity.pdf / .png
#   Console: Shannon medians, Spearman for Shannon vs Env-l param-s, Kruskal + Dunn for Shannon vs Cluster
# -------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(dunn.test)
library(broom)     # for tidy to summarize stats
library(patchwork)

################################
### SHANNON INDEX by cluster ###
# [rerun the commented out lines if you want to compute it yourself]
# Shannon index (α-diversity) : richness and evenness - positively correlated to both
#
# # 1. Decide which samples have very low number of reads and remove them
# hist(rowSums(tab.16S), breaks = 100)
# hist(rowSums(tab.18S), breaks = 100)
# samples_to_remove <- c("2022: Stn 18 (1 m)", "2022: Stn 2 (2 m)", "2022: Stn 4 (2 m)", "2022: Stn 6 (2 m)",  "2022: Stn 19 (52 m)")
# 
# # 2. Remove the samples with low reads count and ensure the order for ASV tables
# meta.wk.small.deep <- meta.wk.small[!rownames(meta.wk.small) %in% samples_to_remove, ]
# kept_idx <- match(rownames(meta.wk.small.deep), rownames(meta.wk.small))
# tab.16S.deep <- tab.16S[kept_idx, ]
# tab.18S.deep <- tab.18S[kept_idx, ]
# 
# # 3. Set rarefaction depth to minimum reads per domain and see how it looks just for richness
# euk_min_reads <- min(meta.wk.small.deep$euk_reads)
# prok_min_reads <- min(meta.wk.small.deep$prok_reads)
# rarecurve(tab.18S.deep, step = 100, col = "gray40", label = FALSE, main = "Eukaryota"); abline(v = euk_min_reads, col = "red")
# rarecurve(tab.16S.deep, step = 100, col = "gray40", label = FALSE, main = "Prokaryota"); abline(v = prok_min_reads, col = "red")
# 
# # 4. Calculate mean Shannon diversity after 1000 rarefactions and take average - takes 1 h to run
# calculate_mean_shannon <- function(count_table, rarefy_depth, n_reps = 1000) {
#   set.seed(161193)
#   shannon_results <- matrix(0, nrow = nrow(count_table), ncol = n_reps)       # initialize results matrix
#   for (i in 1:n_reps) {                                                       # perform rarefaction and calculate Shannon index 1000 times
#     rarefied <- suppressWarnings(rrarefy(count_table, sample = rarefy_depth))
#     shannon_results[, i] <- vegan::diversity(rarefied, index = "shannon")
#   }
#   rowMeans(shannon_results)                                                   # calculate mean across rarefactions
# }
# meta.wk.small.deep$Shannon_Euk <- calculate_mean_shannon(tab.18S.deep, euk_min_reads)
# meta.wk.small.deep$Shannon_Prok <- calculate_mean_shannon(tab.16S.deep, prok_min_reads)
# saveRDS(meta.wk.small.deep %>% rownames_to_column("Sample") %>% select(Sample, Shannon_Euk, Shannon_Prok), "OUTPUT/shannon_diversity.rds") # only Shannon

# attach Shannon: rmv if you have it to avoid errors on rerun; if your meta changes, you can use pre-calc-ed Shannon as far sample names are the same
meta.wk.small <- meta.wk.small %>% select(-contains("Shannon")) %>%
  rownames_to_column("Sample") %>% left_join(readRDS("OUTPUT/shannon_diversity.rds"), by = "Sample") %>% column_to_rownames("Sample")


####################################################################################
### PLOT MOST ABUDANT GROUPS AND SHANNON & STATISTICALLY TEST SHANNON DIFFERENCE ###

## 1. Specify plankton groups and calculate and assign colours
parse_tax <- function(asv) {
  tx <- asv %>% str_remove(";\\s*ASV_\\d+$") %>% stringr::str_remove_all("\\s*\\[\\d+\\]"); parts <- strsplit(tx, ";\\s*")[[1]]
  tibble(Taxonomy = tx, Domain = if(length(parts) >= 1) parts[1] else NA_character_, Phylum = if(length(parts) >= 2) parts[2] else NA_character_,
         Order = if(length(parts) >= 3) parts[3] else NA_character_, Genus = if(length(parts) >= 4) parts[4] else NA_character_)}

all_asvs_tax <- tibble(ASV_full = colnames(proptab.small.taxan.saman)) %>% rownames_to_column() %>% rowwise() %>% mutate(bind_rows(parse_tax(ASV_full))) %>% 
  ungroup() %>% mutate(ASV_id = str_extract(ASV_full, "ASV_\\d+$")) %>%
  left_join(bind_rows(asv16, asv18) %>% select(ASV, Class), by = c("ASV_id" = "ASV")) %>%
  transmute(ASV = ASV_full, Domain2 = if_else(Domain %in% c("Bacteria", "Archaea"), "Prokaryote", "Eukaryote"),
            Group = if_else(Domain == "Archaea", "Archaea",
                             case_when(Class == "Alphaproteobacteria" ~ "α-Proteobacteria", Class == "Gammaproteobacteria" ~ "γ-Proteobacteria",
                                       # in Gyrista, focus on Diatoms and other photosynthetic Gyrista: 
                                       Class %in% c("Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae") ~ "Diatoms", # if all of them, then "Ochrophyta"
                                       Class %in% c("Pelagophyceae", "Chrysophyceae", "Dictyochophyceae", "Bolidophyceae", "Phaeophyceae", "Eustigmatophyceae",
                                                    "Raphidophyceae", "MOCH-2") ~ "non-diatom Ochrophyta", # Gyrista: NO Diatoms/Peronosporomycetes/Pirsoniales/MAST
                                       Order == "Crustacea" ~ "Crustacea",
                                       Class == "Dinophyceae" ~ "Dinophyceae", Class == "Syndiniales" ~ "Syndiniales",
                                       TRUE ~ str_replace(Phylum, "^uncl.\\s+", ""))))

# Colour major groups
euk_colors <- c("Diatoms" = "#008080", "Chlorophyta" = "#00D100", "non-diatom Ochrophyta" = "#29A229", "Cryptophyta" = "#86C435", "Dinophyceae" = "#00CDDA",
                "Ciliophora" = "#0000FF", "Cercozoa" = "#1A4AB2", "Telonemia" = "#3A86FF", "Crustacea" = "#8E2E9F", "Syndiniales" = "#8E74E6", "Fungi" = "#D4BF26")
prok_colors <- c("Cyanobacteriota" = "#02c78f", "α-Proteobacteria" = "#B84040", "γ-Proteobacteria" = "#f62334", "Bacteroidota" = "#FF8C00", 
                 "Verrucomicrobiota" = "#cc8d8d", "Actinomycetota" = "#f330b2", "Planctomycetota" = "#b623f6", "Archaea" = "#786D58")

# Calculate Relab of major groups
cluster_group <- meta.wk.small %>% rownames_to_column("Sample") %>% select(Sample, Cluster) %>%
  inner_join(proptab.small.taxan.saman %>% as_tibble(rownames = "Sample") %>% pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb"), by = "Sample") %>%
  left_join(all_asvs_tax, by = "ASV") %>% group_by(Sample, Cluster, Domain2, Group) %>% summarise(GroupAbundSample = sum(RelAb), .groups = "drop") %>%
  group_by(Cluster, Domain2, Group) %>% summarise(AvgAbund = mean(GroupAbundSample) * 100, .groups = "drop") %>%
  mutate(ClusterName = factor(cl_labels[as.numeric(Cluster)], levels = rev(cl_labels), ordered = TRUE))


## 2. Explore Relab of top prokaryotic and eukaryotic groups as well as of specific ASVs, genera, etc. to see how best to communicate the dynamics of the plankton
# # First of all, create a stacked bar plot to see what top 10 major groups per cluster look like
# # Then choose what to colour in order to present as much relative abundance as possible without highlighting too many groups.
#
# # Calculate abundances
# cluster_top10_groups <- cluster_group %>% group_by(Cluster, Domain2) %>% slice_max(AvgAbund, n = 10, with_ties = FALSE) %>% ungroup()
# 
# # Prepare data with correct ordering
# euk_data <- cluster_top10_groups %>% filter(Domain2 == "Eukaryote") %>% mutate(Group = factor(Group, levels = rev(names(euk_colors)), ordered = TRUE))
# prok_data <- cluster_top10_groups %>% filter(Domain2 == "Prokaryote") %>% mutate(Group = factor(Group, levels = rev(names(prok_colors)), ordered = TRUE))
# 
# # Create plots with position_stack(reverse=TRUE)
# euk_plot <- ggplot(euk_data, aes(x = AvgAbund, y = ClusterName, fill = Group)) +
#   geom_bar(stat = "identity", position = position_stack(reverse=FALSE), width = 0.7, color = "black", linewidth = 0.1) +
#   scale_fill_manual(values = euk_colors, breaks = names(euk_colors), drop = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(0,0.05)), limits = c(0,100)) +
#   labs(x = "Average Relative Abundance (%)", y = "Cluster", title = "Eukaryota") +
#   theme_bw() + theme(axis.text.y = element_text(face = "bold", color = "black", size = 15), panel.grid.major.y = element_blank(),
#                      legend.position = "bottom", legend.title = element_blank()) +
#   guides(fill = guide_legend(ncol = 3, byrow = FALSE, title.position = "top", label.position= "right"))
# 
# prok_plot <- ggplot(prok_data, aes(x = AvgAbund, y = ClusterName, fill = Group)) +
#   geom_bar(stat = "identity", position = position_stack(reverse=FALSE), width = 0.7, color = "black", linewidth = 0.1) +
#   scale_fill_manual(values = prok_colors, breaks = names(prok_colors), drop = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(0,0.05)), limits = c(0,100)) +
#   labs(x="Average Relative Abundance (%)", y = NULL, title = "Prokaryota") +
#   theme_bw() + theme(axis.text.y = element_text(face = "bold", color = "black", size = 15), panel.grid.major.y = element_blank(),
#                      legend.position = "bottom", legend.title = element_blank()) +
#   guides(fill = guide_legend(ncol = 3, byrow=FALSE, title.position = "top", label.position = "right"))
# 
# # Combine and save
# combined_plot <- euk_plot + prok_plot + plot_layout(nrow = 1, widths = c(1, 1))
# cairo_pdf("OUTPUT/cluster_top10_group_stackedbar.pdf", width = 16, height = 11, onefile = FALSE); print(combined_plot); dev.off()
#
#
# # This is for further exploration and to make it easier to report specific numbers in the Results/Discussion
#
# # Mean Relab of ASVs across all samples
# tibble(ASV = str_extract(colnames(proptab.small.taxan.saman), "ASV_\\d+"), Relab = colMeans(proptab.small.taxan.saman) * 100) %>% View
# 
# # Relab of a specific ASV in each samples to see the most abundant ones overall
# proptab.small.taxan.saman %>% select(ASV = contains("ASV_210")) %>% rownames_to_column("Sample") %>% mutate(RelAb = ASV * 100) %>% print() %>% View
# 
# # Mean Relab of a specific ASV within each cluster to see where it peaks
# sapply(split(meta.wk.small, meta.wk.small$Cluster), function(d) { # 1 = oHB, 2 = BI22, 3 = cHB, 4 = dBI, 5 = dJB, 6 = EST, 7 = MR, 8 = nJB, 9 = sJB
#   samps <- rownames(d)
#   asv_col <- grep("ASV_714$", colnames(proptab.small.taxan.saman), value = TRUE)
#   gen_cols <- grep("Eukaryota", colnames(proptab.small.taxan.saman), value = TRUE)
#   mean(proptab.small.taxan.saman[samps, asv_col] / rowSums(proptab.small.taxan.saman[samps, gen_cols], na.rm = TRUE), na.rm = TRUE) * 100})
# 
# # Mean Relab of a Genus within each cluster to see how certain genera behave in different clusters
# sapply(split(meta.wk.small, meta.wk.small$Cluster), function(d) {
#   samps <- rownames(d)
#   genus_cols <- grep("Choricystis", colnames(proptab.small.taxan.saman), value = TRUE) # Genus
#   group_cols <- grep("Chlorophyta", colnames(proptab.small.taxan.saman), value = TRUE) # part of what
#   mean(rowSums(proptab.small.taxan.saman[samps, genus_cols, drop = FALSE], na.rm = TRUE) /
#          rowSums(proptab.small.taxan.saman[samps, group_cols, drop = FALSE], na.rm = TRUE), na.rm = TRUE) * 100})
# 
# # For each cyanobacterial bioindicator ASV, calculate:
# # the mean Relab of the ASV within the Cluster it characterizes; the mean Relab of the overall Genus within the Cluster;
# # and how much the bioindicator ASV contributes to its Genus within the Cluster
# params <- tribble(~name, ~clus, ~asv_pat, ~gen_pat, "ASV_17_cHB", 3, "ASV_17(\\D|$)", "Synechococcus",
#                   "ASV_221_EST", 6, "ASV_221(\\D|$)", "Cyanobium", "ASV_210_sJB", 9, "ASV_210(\\D|$)", "Cyanobium")
# params %>%
#   pmap_dfr(function(name, clus, asv_pat, gen_pat) {
#     d <- meta.wk.small %>% filter(Cluster == clus)
#     asv_col <- grep(asv_pat, colnames(proptab.small.taxan.saman), value = TRUE)
#     gen_cols <- grep(gen_pat, colnames(proptab.small.taxan.saman), value = TRUE)
#     asv <- mean(rowSums(proptab.small.taxan.saman[rownames(d), asv_col, drop = FALSE]) * 100, na.rm = TRUE)
#     genus <- mean(rowSums(proptab.small.taxan.saman[rownames(d), gen_cols, drop = FALSE]) * 100, na.rm = TRUE)
#     tibble(ASV_rel_ab = asv, Genus_rel_ab = genus, ASV_pct_of_genus = asv / genus * 100)
#   }) %>% mutate(name = params$name, .before = 1) %>% column_to_rownames("name")
# 
# # Look-up the Relab of all Genera belonging to a specific group or groups to see which ones are the most abundant
# colSums(proptab.small.taxan.saman[rownames(meta.wk.small)[meta.wk.small$`Depth [m]` >= 1 & meta.wk.small$`Depth [m]` <= 6], ]) %>%
#   enframe(name = "ASV", value = "total") %>% mutate (Genus = str_split_i(ASV, "; ", 4)) %>%
#   inner_join(all_asvs_tax %>% filter(Group %in% "Diatoms"), by = "ASV") %>%
#   group_by(Genus) %>% summarise(total = sum(total) * 100) %>% arrange(desc(total))
# 
# # RelAb of a major group to know the precise percentages behind the stacked bars
# View(cluster_group %>% filter(Group %in% c("Ciliophora", "Telonemia", "Cercozoa")))


## 3. Make a bar plots for major groups (once you have decided what to highlight)

# Filter data to include ONLY the chosen groups in the stacked bar plot
euk_data <- cluster_group %>% filter(Domain2 == "Eukaryote", Group %in% names(euk_colors)) %>%
  mutate(Group = factor(Group, levels = rev(names(euk_colors)), ordered = TRUE))
prok_data <- cluster_group %>% filter(Domain2 == "Prokaryote", Group %in% names(prok_colors)) %>%
  mutate(Group = factor(Group, levels = rev(names(prok_colors)), ordered = TRUE))

# Panels A and B: major eukaryotic and prokaryotic groups
euk_plot <- ggplot(euk_data, aes(x = AvgAbund, y = ClusterName, fill = Group)) +
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE),  width = 0.7, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = euk_colors, breaks = names(euk_colors), drop = FALSE) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "Relative Abundance (%)", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_blank(),
        plot.margin = margin(l = 15, r = 15),
        legend.key.spacing.x = unit(2, "cm")) +                       # push legend apart horizontally
  guides(fill = guide_legend(ncol = 3, byrow = FALSE))

prok_plot <- ggplot(prok_data, aes(x = AvgAbund, y = ClusterName, fill = Group)) +
  geom_bar(stat = "identity", position = position_stack(reverse = FALSE),  width = 0.7, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = prok_colors, breaks = names(prok_colors), drop = FALSE) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = "Relative Abundance (%)", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_blank(),
        plot.margin = margin(l = 20, r = 10),
        legend.key.spacing.x = unit(1.5, "cm")) +
  guides(fill = guide_legend(ncol = 3, byrow = FALSE))


## 4. Make a plot of Shannon diversity per cluster
shannon_long <- meta.wk.small %>% rownames_to_column("Sample") %>%
  pivot_longer(c(Shannon_Euk, Shannon_Prok), names_to = "Domain", values_to = "Shannon", values_drop_na = TRUE) %>% # drop samples w/ no Shannon
  mutate(Domain = if_else(Domain == "Shannon_Euk", "Eukaryota", "Prokaryota"), Domain = factor(Domain, levels = c("Eukaryota","Prokaryota")),
         ClusterName = factor(cl_labels[as.numeric(Cluster)], levels = rev(cl_labels), ordered = TRUE))

shannon_plot <- ggplot(shannon_long, aes(x = Shannon, y = ClusterName, color = Domain)) +
  geom_boxplot(fill = "white", width = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("Eukaryota" = "#D55E00", "Prokaryota" = "#0072B2"),
                     breaks = c("Eukaryota", "Prokaryota")) +
  scale_x_continuous(limits = c(1, 6), breaks = 1:6, expand = c(0, 0)) +
  labs(x = "Shannon Diversity Index", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 15),
        legend.title = element_blank(),
        plot.margin = margin(l = 10, r = 20))

## 5. Combine plots with adjusted spacing and output them
final_plot <- (euk_plot + prok_plot + shannon_plot) + plot_layout(nrow = 1, widths = c(1, 1, 0.75))

# Save and make a table of Shannon medians for yourself
cairo_pdf("OUTPUT/figure_community_diversity.pdf", width = 24, height = 9, onefile = FALSE); print(final_plot); dev.off()

png("OUTPUT/figure_community_diversity.png", width = 24, height = 9, units = "in", res = 300)
print(final_plot)
dev.off()


# 6. Explore Shannon and do Dunn statistical tests
# Shannon per-cluster medians
shannon_long %>% group_by(Domain, ClusterName) %>% 
  summarise(n = sum(!is.na(Shannon)), Median = median(Shannon), Q1 = quantile(Shannon, 0.25, na.rm = TRUE), Q3 = quantile(Shannon, 0.75, na.rm = TRUE),
            IQR = IQR(Shannon, na.rm = TRUE), .groups = "drop") %>% mutate(Domain = factor(Domain, levels = c("Eukaryota", "Prokaryota"))) %>%
  arrange(Domain, desc(ClusterName))

# Environmental Parameters vs Shannon - Spearman correlation
shannon_long %>% pivot_longer(cols = where(is.numeric) & !Shannon, names_to = "Parameter", values_to = "Value") %>% 
  group_by(Domain, `Parameter`) %>% filter(sd(Value, na.rm = TRUE) > 0) %>%  # test cor-n of Shannon indices w/ all numeric values w/ SD > 0 
  summarise(n = sum(complete.cases(Shannon, Value)), tidy(cor.test(Shannon, Value, method = "spearman", exact = F)), .groups = "drop") %>%
  filter(p.value < 0.001) %>% arrange(Domain, p.value)

# see the trends in Shannon index with latitude or other env-l parameters for yourself
# par(mfrow = c(2, 1))
# with(shannon_long %>% filter(Domain == "Eukaryota"), {plot(Shannon ~ `Latitude [deg]`); lines(lowess(`Latitude [deg]`, Shannon), col = "red")})
# with(shannon_long %>% filter(Domain == "Prokaryota"), {plot(Shannon ~ `Latitude [deg]`); lines(lowess(`Latitude [deg]`, Shannon), col = "red")})
# par(mfrow = c(1, 1))

# Do ANOVA and check residuals: if they are non-normally distributed -> do Kruskal instead
fits <- lapply(split(shannon_long, shannon_long$Domain), \(d) aov(Shannon ~ ClusterName, data = d)) # ANOVA
data.frame(Domain = factor(names(fits), levels = c("Eukaryota", "Prokaryota")),
           `Shapiro p-value (residuals)` = sapply(fits, \(m) shapiro.test(residuals(m))$p.value), # p < 0.05 -> bad, your residuals aren't normal
           row.names = NULL, check.names = FALSE)
par(mfrow = c(1, 2)); invisible(lapply(names(fits), \(n) hist(residuals(fits[[n]]), main = n, xlab = NA))); par(mfrow = c(1, 1)) # resid. must be bell-curve-y

# Since residuals are not normally distributed, do Kruskal + Dunn
by_dom <- split(shannon_long, shannon_long$Domain)
do.call(rbind, lapply(names(by_dom), \(dom)                                     # Kruskal-Wallis p-values (per Domain)
                      data.frame(Domain = dom,
                                 `Kruskal p-value` = kruskal.test(Shannon ~ ClusterName, by_dom[[dom]])$p.value,
                                 row.names = NULL, check.names = FALSE))) %>%
  mutate(Domain = factor(Domain, levels = c("Eukaryota","Prokaryota"))) %>% arrange(Domain)

do.call(rbind, lapply(names(by_dom), \(dom) {                                   # Dunn post-hoc (BH-adjusted), w/ direction based on cluster medians
  d   <- by_dom[[dom]]
  med <- tapply(d$Shannon, d$ClusterName, median, na.rm = TRUE)
  dt  <- dunn.test(d$Shannon, d$ClusterName, method = "bh", kw = FALSE, list = TRUE)
  out <- data.frame(comparison = dt$comparisons, `p-value` = dt$P.adjusted, check.names = FALSE)
  out <- out[out$`p-value` < 0.05, , drop = FALSE]; if (!nrow(out)) return(NULL)
  g1 <- trimws(sub("\\s*-.*", "", out$comparison)); g2 <- trimws(sub(".*-\\s*", "", out$comparison))
  data.frame(Domain = dom,
             Result = ifelse(med[g1] > med[g2], paste(g1, ">", g2), paste(g2, ">", g1)),
             `p-value` = out$`p-value`,
             row.names = NULL, check.names = FALSE)})) %>%
  mutate(Domain = factor(Domain, levels = c("Eukaryota","Prokaryota"))) %>% arrange(Domain, `p-value`)


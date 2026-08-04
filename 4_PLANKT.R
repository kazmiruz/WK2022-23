# 4------------------------------------------------------------------------
# Purpose:
#   Assess major plankton taxa and Shannon diversity.
#   Specifically, attach pre‑computed Shannon diversity (can compute them here yourself) to the metadata,
#   build a unified taxonomy lookup, compute average relative abundances per cluster for each major group (usually at the phylum level),
#   produce a combined figure (Shannon boxplots + domain‑wise stacked bar charts), and perform
#   statistical tests on Shannon diversity (Spearman correlations, Kruskal–Wallis, Dunn post‑hoc).
#
# Prerequisite: run both 2_TAX and 3_CLUS first, as they generate various objects used here.
#
# Required external files - you can also reproduce them all yourself using the commented out code:
#   EXT/shannon_diversity.rds                      (pre‑computed Shannon indices; can reproduce using commented code)
#
# Outputs:
#   OUT/4_PLANKT/cluster_top10_group_stackedbar.pdf (top 10 major groups per cluster for exploration)
#   OUT/4_PLANKT/plankton_communities.pdf / .png    (major taxonomic groups and Shannon diversity for both eukaryotes and prokaryotes)
#   OUT/4_PLANKT/Shannon_statistics.xlsx            (Shannon diversity statistics: change btwn clusters and corr. w/ env. parameters)
# -------------------------------------------------------------------------

library(tidyverse)
library(vegan)
library(dunn.test)
library(broom)     # for tidy to summarize stats
library(patchwork)
library(openxlsx)  # for writing complex Excel files

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
# euk_min_reads <- min(meta.wk.small.deep$'Eukaryotic Reads')
# prok_min_reads <- min(meta.wk.small.deep$'Prokaryotic Reads')
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
# saveRDS(meta.wk.small.deep %>% rownames_to_column("Sample") %>% select(Sample, Shannon_Euk, Shannon_Prok), "EXT/shannon_diversity.rds") # only Shannon

# attach Shannon: rmv if you have it to avoid errors on rerun; if your meta changes, you can use pre-calc-ed Shannon as far sample names are the same
meta.wk.small <- meta.wk.small %>% select(-contains("Shannon")) %>%
  rownames_to_column("Sample") %>% left_join(readRDS("EXT/shannon_diversity.rds"), by = "Sample") %>% column_to_rownames("Sample")


####################################################################################
### PLOT MOST ABUDANT GROUPS AND SHANNON & STATISTICALLY TEST SHANNON DIFFERENCE ###

## 1. Provide cluster dendrogram (sample tree from 3_CLUS cut to k branches)
# Cut the rotated sample dendrogram at the height giving k groups, keep only the upper part -> one leaf per cluster.
# Converting that upper tree to hclust cleans up the leaf coordinates (evenly spaced 1..k) and drops every tip to a common
# baseline (0), so branches no longer end at uneven heights. Leaf order (1..k) is cluster 1..k = oHB..sJB,
# matching the top->bottom order of the bar/box panels, so no leaf labels are needed (they're read off the panels to the right).
dendro_upper <- cut(dend_row_order,
                    h = dendextend::heights_per_k.dendrogram(dend_row_order)[as.character(k)])$upper
dendro_hc  <- as.hclust(dendro_upper)                          # even tips at height 0, real merge heights kept
dendro_seg <- dendextend::as.ggdend(as.dendrogram(dendro_hc))$segments

# Optionally colour terminal branches per cluster_palette. Not done b/c the abundance panels are colour-dense; set to TRUE to enable.
dendro_color <- FALSE
dendro_seg$branch_col <- "grey20"
if (dendro_color) {tips <- dendro_seg$yend == 0; dendro_seg$branch_col[tips] <- cluster_palette[dendro_seg$x[tips]]}

dendro_plot <- ggplot(dendro_seg) +
  # flip to horizontal (height on x) and map leaf i -> row (k+1-i) so oHB is on top
  geom_segment(aes(x = y, y = (k + 1) - x, xend = yend, yend = (k + 1) - xend, color = branch_col), linewidth = 1, lineend = "round") +
  scale_color_identity() +
  scale_y_continuous(expand = expansion(add = 0.6)) +            # match the discrete bar rows
  scale_x_reverse(expand = expansion(mult = c(0.05, 0.02))) +    # root on the left, tips facing right
  theme_void() +
  theme(plot.margin = margin(l = 10, r = 5, t = 5, b = 5))


## 2. Specify plankton groups and calculate and assign colours
parse_tax <- function(asv) {
  tx <- asv %>% str_remove(";\\s*ASV_\\d+$") %>% stringr::str_remove_all("\\s*\\[\\d+\\]"); parts <- strsplit(tx, ";\\s*")[[1]]
  tibble(Taxonomy = tx, Domain = if(length(parts) >= 1) parts[1] else NA_character_, Phylum = if(length(parts) >= 2) parts[2] else NA_character_,
         Order = if(length(parts) >= 3) parts[3] else NA_character_, Genus = if(length(parts) >= 4) parts[4] else NA_character_)}

all_asvs_tax <- tibble(ASV_full = colnames(proptab.small.taxan.saman)) %>% rownames_to_column() %>% rowwise() %>%
  mutate(bind_rows(parse_tax(ASV_full))) %>% 
  ungroup() %>% mutate(ASV_id = str_extract(ASV_full, "ASV_\\d+$")) %>%
  left_join(bind_rows(asv16, asv18) %>% select(ASV, Class), by = c("ASV_id" = "ASV")) %>%
  transmute(ASV = ASV_full, Domain2 = if_else(Domain %in% c("Bacteria", "Archaea"), "Prokaryote", "Eukaryote"),
            Group = if_else(Domain == "Archaea", "Archaea",
                            case_when(Class == "Alphaproteobacteria" ~ "α-Proteobacteria", Class == "Gammaproteobacteria" ~ "γ-Proteobacteria",
                                      Class %in% c("Bacillariophyceae", "Coscinodiscophyceae", "Mediophyceae") ~ "Diatoms",
                                      Class %in% c("Pelagophyceae", "Chrysophyceae", "Dictyochophyceae", "Bolidophyceae",
                                                   "Phaeophyceae", "Eustigmatophyceae", "Raphidophyceae", "MOCH-2") ~ "non-diatom Ochrophyta",
                                      Order == "Crustacea" ~ "Crustacea",
                                      Class == "Dinophyceae" ~ "Dinophyceae", Class == "Syndiniales" ~ "Syndiniales",
                                      TRUE ~ str_replace(Phylum, "^uncl.\\s+", ""))))

# Colour major groups
euk_colors <- c("Diatoms" = "#008080", "Chlorophyta" = "#00D100", "non-diatom Ochrophyta" = "#29A229", "Cryptophyta" = "#86C435",
                "Dinophyceae" = "#00CDDA", "Ciliophora" = "#0000FF", "Cercozoa" = "#1A4AB2", "Telonemia" = "#3A86FF",
                "Crustacea" = "#8E2E9F", "Syndiniales" = "#8E74E6", "Fungi" = "#D4BF26")
prok_colors <- c("Cyanobacteriota" = "#02c78f", "α-Proteobacteria" = "#B84040", "γ-Proteobacteria" = "#f62334", "Bacteroidota" = "#FF8C00", 
                 "Verrucomicrobiota" = "#cc8d8d", "Actinomycetota" = "#f330b2", "Planctomycetota" = "#b623f6", "Archaea" = "#786D58")

# Calculate RelAb of major groups
cluster_group <- meta.wk.small %>% rownames_to_column("Sample") %>% select(Sample, Cluster) %>%
  inner_join(proptab.small.taxan.saman %>% as_tibble(rownames = "Sample") %>%
               pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb"), by = "Sample") %>%
  left_join(all_asvs_tax, by = "ASV") %>% group_by(Sample, Cluster, Domain2, Group) %>%
  summarise(GroupAbundSample = sum(RelAb), .groups = "drop") %>%
  group_by(Cluster, Domain2, Group) %>% summarise(AvgAbund = mean(GroupAbundSample) * 100, .groups = "drop") %>%
  mutate(Cluster = factor(cl_labels[as.numeric(Cluster)], levels = rev(cl_labels), ordered = TRUE))


# 3. Explore RelAb of top eukaryotic and prokaryotic groups and specific ASVs, genera, etc. to see how best to communicate the community dynamics

## First of all, create a diagnostic stacked bar plot to see what top 10 major groups per cluster look like.
# Then choose what to colour in order to present as much relative abundance as possible without highlighting too many groups.

# Calculate abundances
top10_groups <- cluster_group %>% group_by(Cluster, Domain2) %>% slice_max(AvgAbund, n = 10, with_ties = FALSE) %>% ungroup()

# Prepare data with correct ordering
euk_data <- top10_groups %>% filter(Domain2 == "Eukaryote") %>% mutate(Group = factor(Group, levels = rev(names(euk_colors)), ordered = TRUE))
prok_data <- top10_groups %>% filter(Domain2 == "Prokaryote") %>% mutate(Group = factor(Group, levels = rev(names(prok_colors)), ordered = TRUE))

# # Create plots with position_stack (reverse=TRUE)
# euk_plot <- ggplot(euk_data, aes(x = AvgAbund, y = Cluster, fill = Group)) +
#   geom_bar(stat = "identity", position = position_stack(reverse = FALSE), width = 0.7, color = "black", linewidth = 0.1) +
#   scale_fill_manual(values = euk_colors, breaks = names(euk_colors), drop = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 100)) +
#   labs(x = "Eukaryotic Relative Abundance (%)", y = "Cluster", title = "Eukaryota") +
#   theme_bw() + theme(axis.text.y = element_text(face = "bold", color = "black", size = 15), panel.grid.major.y = element_blank(),
#                      legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10)) +
#   guides(fill = guide_legend(ncol = 4, byrow = FALSE, title.position = "top", label.position= "right"))
# 
# prok_plot <- ggplot(prok_data, aes(x = AvgAbund, y = Cluster, fill = Group)) +
#   geom_bar(stat = "identity", position = position_stack(reverse = FALSE), width = 0.7, color = "black", linewidth = 0.1) +
#   scale_fill_manual(values = prok_colors, breaks = names(prok_colors), drop = FALSE) +
#   scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 100)) +
#   labs(x="Prokaryotic Relative Abundance (%)", y = NULL, title = "Prokaryota") +
#   theme_bw() + theme(axis.text.y = element_text(face = "bold", color = "black", size = 15), panel.grid.major.y = element_blank(),
#                      legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 10)) +
#   guides(fill = guide_legend(ncol = 4, byrow = FALSE, title.position = "top", label.position = "right"))
# 
# # Combine and save
# combined_plot <- euk_plot + prok_plot + plot_layout(nrow = 1, widths = c(1, 1))
# cairo_pdf("OUT/4_PLANKT/top10_per_cluster.pdf", width = 15, height = 6, onefile = FALSE); print(combined_plot); dev.off()
#
#
# ## This is for further exploration and to make it easier to report specific numbers in the Results/Discussion
#
# # Mean RelAb of ASVs across all samples
# tibble(ASV = str_extract(colnames(proptab.small.taxan.saman), "ASV_\\d+"), RelAb = colMeans(proptab.small.taxan.saman) * 100) %>% View
# 
# # RelAb of a specific ASV in each samples to see the most abundant ones overall
# proptab.small.taxan.saman %>% select(ASV = contains("ASV_210")) %>% rownames_to_column("Sample") %>% mutate(RelAb = ASV * 100) %>% print() %>% View
# 
# # Mean RelAb of a specific ASV within each cluster to see where it peaks
# sapply(split(meta.wk.small, meta.wk.small$Cluster), function(d) { # 1 = oHB, 2 = BI22, 3 = cHB, 4 = dBI, 5 = dJB, 6 = EST, 7 = MR, 8 = nJB, 9 = sJB
#   samps <- rownames(d)
#   asv_col <- grep("ASV_714$", colnames(proptab.small.taxan.saman), value = TRUE)
#   gen_cols <- grep("Eukaryota", colnames(proptab.small.taxan.saman), value = TRUE)
#   mean(proptab.small.taxan.saman[samps, asv_col] / rowSums(proptab.small.taxan.saman[samps, gen_cols], na.rm = TRUE), na.rm = TRUE) * 100})
# 
# # Mean RelAb of a Genus within each cluster to see how certain genera behave in different clusters
# sapply(split(meta.wk.small, meta.wk.small$Cluster), function(d) {
#   samps <- rownames(d)
#   genus_cols <- grep("Choricystis", colnames(proptab.small.taxan.saman), value = TRUE) # Genus
#   group_cols <- grep("Chlorophyta", colnames(proptab.small.taxan.saman), value = TRUE) # part of what
#   mean(rowSums(proptab.small.taxan.saman[samps, genus_cols, drop = FALSE], na.rm = TRUE) /
#          rowSums(proptab.small.taxan.saman[samps, group_cols, drop = FALSE], na.rm = TRUE), na.rm = TRUE) * 100})
# 
# # For each cyanobacterial bioindicator ASV, calculate:
# # the mean RelAb of the ASV within the Cluster it characterizes; the mean RelAb of the overall Genus within the Cluster;
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
# # Look-up the RelAb of all Genera belonging to a specific group or groups to see which ones are the most abundant
# colSums(proptab.small.taxan.saman[rownames(meta.wk.small)[meta.wk.small$`Depth [m]` >= 1 & meta.wk.small$`Depth [m]` <= 6], ]) %>%
#   enframe(name = "ASV", value = "total") %>% mutate (Genus = str_split_i(ASV, "; ", 4)) %>%
#   inner_join(all_asvs_tax %>% filter(Group %in% "Diatoms"), by = "ASV") %>%
#   group_by(Genus) %>% summarise(total = sum(total) * 100) %>% arrange(desc(total))
# 
# # RelAb of a major group to know the precise percentages behind the stacked bars
# View(cluster_group %>% filter(Group %in% c("Ciliophora", "Telonemia", "Cercozoa")))


## 4. Make a bar plots for major groups (once you have decided what to highlight)

# Formatting shared by all panels
community_theme <- theme_bw() +
  theme(axis.text.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 15),
        axis.title = element_text(color = "black", size = 16, face = "bold"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(color = "black", size = 16),
        legend.title = element_blank())

# Filter data to include ONLY the chosen groups in the stacked bar plot
major_group_plot <- function(data, domain, colors, xlab, left, right, legend_space) {
  data %>% filter(Domain2 == domain, Group %in% names(colors)) %>%
    mutate(Group = factor(Group, levels = rev(names(colors)), ordered = TRUE)) %>%
    ggplot(aes(AvgAbund, Cluster, fill = Group)) +
    geom_col(orientation = "y", position = position_stack(reverse = FALSE),
             width = 0.7, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = colors, breaks = names(colors), drop = FALSE) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(x = xlab, y = NULL) +
    community_theme +
    theme(plot.margin = margin(l = left, r = right), legend.key.spacing.x = grid::unit(legend_space, "cm")) +
    guides(fill = guide_legend(ncol = 3, byrow = FALSE))
}

# Plot major groups for eukaryotes and prokaryotes
euk_plot <- major_group_plot(cluster_group, "Eukaryote", euk_colors, "Eukaryotic Relative Abundance (%)", 15, 15, 2)
prok_plot <- major_group_plot(cluster_group, "Prokaryote", prok_colors, "Prokaryotic Relative Abundance (%)", 20, 10, 1.5)


## 5. Plot Shannon diversity for both euks and proks per cluster; explore how diversity changes btwn clusters and env. variables 
domain_order <- c("Eukaryota", "Prokaryota")
shannon_long <- meta.wk.small %>% rownames_to_column("Sample") %>%
  pivot_longer(c(Shannon_Euk, Shannon_Prok), names_to = "Domain", values_to = "Shannon", values_drop_na = TRUE) %>%
  mutate(Domain = recode(Domain, Shannon_Euk = "Eukaryota", Shannon_Prok = "Prokaryota"),
         Domain = factor(Domain, levels = domain_order), Cluster = factor(cl_labels[as.numeric(Cluster)], levels = rev(cl_labels), ordered = TRUE))
by_dom <- split(shannon_long, shannon_long$Domain)

# Plot Shannon diversity for eukaryotes and prokaryotes
shannon_plot <- ggplot(shannon_long, aes(Shannon, Cluster, color = Domain)) +
  geom_boxplot(fill = "white", width = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(Eukaryota = "#D55E00", Prokaryota = "#0072B2")) +
  scale_x_continuous(breaks = 1:6, expand = c(0, 0)) + coord_cartesian(xlim = c(1, 6)) +
  labs(x = "Shannon Diversity Index", y = NULL) +
  community_theme + theme(plot.margin = margin(l = 10, r = 20))

# Statistics for Shannon diversity within each cluster (the same as shown on the boxplot)
shannon_summary <- shannon_long %>% group_by(Domain, Cluster) %>%
  summarise(n = n(), Median = median(Shannon), Q1 = quantile(Shannon, 0.25), Q3 = quantile(Shannon, 0.75), IQR = IQR(Shannon), .groups = "drop") %>%
  arrange(Domain, desc(Cluster))

# Check ANOVA residual normality before choosing ANOVA or Kruskal-Wallis
fits <- lapply(by_dom, \(d) aov(Shannon ~ Cluster, data = d))
imap_dfr(fits, \(m, dom) tibble(Domain = dom, `Shapiro p-value` = shapiro.test(residuals(m))$p.value)) %>%
  mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain) %>% print() # Shapiro p < 0.05 indicates non-normal residuals for proks
walk2(fits, names(fits), ~{ qqnorm(residuals(.x), main = .y); qqline(residuals(.x), col = "red") }) # depart from the red line -> non-normal residuals

# Overall differences in Shannon diversity among clusters (Kriskal-Wallis test)
shannon_kruskal <- imap_dfr(by_dom, \(d, dom) tibble(Domain = dom, `Kruskal-Wallis p-value` = kruskal.test(Shannon ~ Cluster, data = d)$p.value)) %>%
  mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain)

# Significant pairwise differences in Shannon diversity among clusters (post hoc Dunn’s test with Benjamini-Hochberg adjustment)
shannon_dunn <- imap_dfr(by_dom, \(d, dom) {
  med <- c(tapply(d$Shannon, d$Cluster, median))
  dt <- dunn.test(d$Shannon, d$Cluster, method = "bh", kw = FALSE, list = TRUE)
  groups <- str_split_fixed(dt$comparisons, "\\s+-\\s+", 2)
  tibble(Domain = dom, `Direction of difference` = case_when(
    med[groups[, 1]] > med[groups[, 2]] ~ paste(groups[, 1], ">", groups[, 2]),
    med[groups[, 2]] > med[groups[, 1]] ~ paste(groups[, 2], ">", groups[, 1]),
    TRUE ~ paste(groups[, 1], "=", groups[, 2])), `BH-adjusted p-value` = dt$P.adjusted)}) %>%
  filter(`BH-adjusted p-value` < 0.05) %>% mutate(Domain = factor(Domain, levels = domain_order)) %>% arrange(Domain, `BH-adjusted p-value`)

# Significant Spearman correlations with spatiotemporal and environmental parameters
shannon_spearman <- shannon_long %>%
  pivot_longer(where(is.numeric) & !any_of(c("Shannon", "Cluster")), names_to = "Parameter", values_to = "Value") %>%
  filter(complete.cases(Shannon, Value)) %>% group_by(Domain, Parameter) %>%
  filter(n() >= 3, sd(Shannon) > 0, sd(Value) > 0) %>%
  summarise(n = n(), rho = cor(Shannon, Value, method = "spearman"),
            `p-value` = cor.test(Shannon, Value, method = "spearman", exact = FALSE)$p.value, .groups = "drop") %>%
  filter(`p-value` < 0.001) %>% arrange(Domain, `p-value`)

# Export cluster results together and environmental results separately
wb <- createWorkbook()
addWorksheet(wb, "Diversity by cluster")
addWorksheet(wb, "Diversity vs environment")

tables <- list("Shannon diversity by cluster" = shannon_summary, "Kruskal-Wallis tests among clusters" = shannon_kruskal,
               "Significant Dunn's pairwise tests" = shannon_dunn)

scientific_txt <- \(p) {
  if (is.na(p)) return(NA_character_)
  x <- str_split_fixed(formatC(p, format = "e", digits = 1), "e", 2)
  paste0(x[1], " × 10", as.integer(x[2]))}

r <- 1
for (nm in names(tables)) {
  tab <- tables[[nm]]
  p_cols <- which(str_detect(names(tab), "p-value$"))
  if (length(p_cols)) tab[p_cols] <- lapply(tab[p_cols], \(x) map_chr(x, scientific_txt))
  writeData(wb, 1, nm, startRow = r)
  writeDataTable(wb, 1, tab, startRow = r + 1)
  r <- r + nrow(tab) + 4
}

writeData(wb, 2, "Significant Spearman correlations", startRow = 1)
writeDataTable(wb, 2, shannon_spearman %>% mutate(`p-value` = map_chr(`p-value`, scientific_txt)), startRow = 2)
for (sheet in 1:2) setColWidths(wb, sheet, cols = 1:7, widths = "auto")

saveWorkbook(wb, "OUT/4_PLANKT/Shannon_statistics.xlsx", overwrite = TRUE)

## 6. Combine plots with adjusted spacing and output them
final_plot <- (dendro_plot + euk_plot + prok_plot + shannon_plot) + plot_layout(nrow = 1, widths = c(0.25, 1, 1, 0.75))
cairo_pdf("OUT/4_PLANKT/plankton_communities.pdf", width = 24, height = 9.5, onefile = FALSE); print(final_plot); dev.off()
png("OUT/4_PLANKT/plankton_communities.png", width = 24, height = 9.5, units = "in", res = 300); print(final_plot); dev.off()


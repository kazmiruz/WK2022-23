# 5------------------------------------------------------------------------
# Purpose:
#   Pipeline to display the most prominent bioindicators and cosmopolitans.
#   Build the indicator/cosmopolitan dot‑plot with dendrogram, total reads, station labels, similar to the ComplexHeatmap.
#
# Prerequisite: run 2_TAX, 3_CLUS, and 4_PLANKT first, as they generate various objects used here.
#
# Outputs:
#   OUT/5_INDIC/dot_ASV_indic_cosmo2.pdf / .png
# -------------------------------------------------------------------------

# ========================= # Libraries # =========================
library(tidyverse)
library(ggdendro)
library(patchwork)

# ========================= # KNOBS (threshold, placement, spacing, sizes) # =========================

# Make the map bigger/smaller for the paper
S <- (24/52)

# RelAb cutoff
thr_euk_indic <- 0.015         # threshold for indic+cosmo
thr_prok_indic <- 0.015        # threshold for indic+cosmo
thr_euk_not_indic <- 0.01      # threshold for non-indic+cosmo
thr_prok_not_indic <- 0.01     # threshold for non-indic+cosmo

# File resolution
width_inches <- 52 * S # export width (inches)
height_inches <- 34 * S # export height (inches)
res_dpi <- 300 # export resolution (DPI)

# Horizontal widths (plots only)
TREE_SCALE <- 0.5 # tree width
W_TREE <- 0.11    # the border around the tree
W_EUK <- 0.8
W_PROK <- 0.8
W_READS <- 0.085
W_YEAR <- 0.038
W_STATION <- 0.026
# W_COAST <- 0.055
W_DEPTH <- 0.034
# W_SAL <- 0.055
# W_TEMP <- 0.055

# Spacer widths (relative). GAP_SMALL can be negative to tighten.
GAP_MED <- -0.002 # gap between Euk and Prok dotplots (go negative if no reaction)
GAP_SMALL <- -0.007 # small gaps between oceanography columns on the right side (go negative if no reaction)

# ASV label placement (labels drawn *inside* dot-panels at y = 0.5)
LABEL_SIZE <- 30 * S
LABEL_ANGLE <- 57  # the smaller - the more to the left the tail
LABEL_Y_AT <- 0.50 # anchor at panel bottom edge
LABEL_VJUST <- 1.4 # >1 places label just below the anchor
LABEL_XNUDGE <- -0.4 # small left shift (counters rightward drift at angle)
LABEL_DNUDGE <- 0.00 # keep 0 so labels never dip below y = 0.5 (prevents clipping)"

BOTTOM_MARGIN_PTS <- 86 * S # bottom-only padding so rotated tails aren't clipped
LEGEND_GAP_HT <- 3.0 # gap between plots row and legends row

euk_colors_dot <- c("Diatoms" = "#008080", "Dinophyceae" = "#00CDDA", "Chlorophyta" = "#00D100", "non-diatom Ochrophyta" = "#29A229",
                    "Cryptophyta" = "#86C435", "Ciliophora" = "#0000FF",
                    # "Cercozoa" = "#1A4AB2", "Telonemia" = "#3A86FF",
                    "Crustacea" = "#8E2E9F", "Syndiniales" = "#8E74E6", "Fungi" = "#D4BF26")
prok_colors_dot <- c("Cyanobacteriota" = "#02c78f",
                     "α-Proteobacteria" = "#B84040", "γ-Proteobacteria" = "#f62334", "Bacteroidota" = "#FF8C00", 
                     # "Verrucomicrobiota" = "#cc8d8d",
                     "Actinomycetota" = "#f330b2", "Planctomycetota" = "#b623f6", "Archaea" = "#786D58")


# ========================= # Helpers # =========================
normalize_sample <- function(x) x %>% str_replace_all("[\u2212\u2013\u2014]", "-") %>% str_squish()
ensure_other_color <- function(cols) if (!"Other" %in% names(cols)) c(cols, Other = "black") else cols
to_hex <- function(cols) {m <- col2rgb(cols); apply(m, 2, \(v) rgb(v[1], v[2], v[3], maxColorValue = 255))}
label_total_reads_ticks <- function(breaks) { sapply(breaks, \(b) if (b == 0) "0" else if (b == 50000) "5×10⁴" else if (b == 100000) "10⁵" else "") }
get_visual_mapping <- function(RelAb) list(alpha = if_else(RelAb <= 0.25, 1 + 0 * (RelAb / 0.25), 1),          # if you want to change alpha too
                                           # Circle size increases with RelAb up to 25%; triangle, diamond, and square sizes were adjusted separately
                                           # because these shapes appear different at the same numeric size
                                           size = S * if_else(RelAb <= 0.25, 2 + 8 * sqrt(RelAb / 0.25),
                                                              case_when(RelAb <= 0.35 ~ 8.6, RelAb <= 0.50 ~ 12.4, TRUE ~ 11.4)),
                                           shape = case_when(RelAb <= 0.25 ~ 16L, RelAb <= 0.35 ~ 17L, RelAb <= 0.50 ~ 18L, TRUE ~ 15L))
add_bottom_margin <- function(p, b_pts = BOTTOM_MARGIN_PTS) p + theme(plot.margin = margin(0, 0, b_pts, 0))

# Clean indicator parser: "1", "2+", "", NA → numeric sort key (9999 if none)
parse_indicator_key <- function(x) {
  x <- replace_na(as.character(x), "")
  vapply(x, function(s) {
    if (!nzchar(s)) return(9999)
    nums <- suppressWarnings(as.integer(strsplit(s, "\\+")[[1]]))
    nums <- nums[!is.na(nums)]
    if (!length(nums)) return(9999)
    2*min(nums) + ifelse(length(nums) > 1, 1, 0)
  }, numeric(1))
}

# ========================= # Ordering + indicator-aware prioritization # =========================
prepare_domain_data <- function(data_plot, data_all, phylum_order, phylum_colors, asv_mean_info) {
  phylum_map <- data_all %>% distinct(ASV, Group)
  
  indicator_info <- data_plot %>% distinct(ASV) %>% mutate(ASV_number = str_extract(ASV, "ASV_\\d+$")) %>%
    left_join(indic_asv_tbl %>% mutate(ASV_number = str_extract(ASV, "ASV_\\d+$"), Clusters = as.character(Clusters)) %>%
                select(ASV_number, Clusters), by = "ASV_number") %>%
    mutate(indicator_sort_key = parse_indicator_key(Clusters)) %>% transmute(ASV, indicator_sort_key)
  
  # Use pre‑computed max‑cluster info instead of recomputing inside the function
  distinct_asvs <- data_plot %>% distinct(ASV) %>% left_join(phylum_map, by = "ASV")
  
  ### ASV Order: Cluster of the indicator -> Phylum -> Cluster with highest mean abundance -> Highest mean abundance in the cluster ###
  asv_ordering <- distinct_asvs %>% left_join(asv_mean_info, by = "ASV") %>% left_join(indicator_info, by = "ASV") %>%
    mutate(Group = factor(Group, levels = phylum_order), MaxCluster_num = suppressWarnings(as.integer(as.character(MaxCluster)))) %>%
    arrange(indicator_sort_key, Group, MaxCluster_num, desc(MaxRelAb), ASV) %>% mutate(asv_x = row_number())
  
  extended_phylum_colors <- ensure_other_color(phylum_colors)
  
  d <- data_plot %>% select(-any_of("Group")) %>% left_join(asv_ordering %>% select(ASV, asv_x, Group), by = "ASV") %>%
    mutate(Group = as.character(Group),
           Group = if_else(Group %in% names(extended_phylum_colors), Group, "Other"),
           Group = factor(Group, levels = names(extended_phylum_colors)))
  
  gl <- asv_ordering %>% left_join(data_plot %>% distinct(ASV, label), by = "ASV") %>% distinct(asv_x, ASV, label, Group) %>%
    mutate(Group = as.character(Group),
           Group = if_else(Group %in% names(extended_phylum_colors), Group, "Other"),
           Group = factor(Group, levels = names(extended_phylum_colors))) %>%
    arrange(asv_x)
  list(d = d, gl = gl, extended_colors = extended_phylum_colors)
}

add_indicator_info <- function(obj) {
  obj$gl <- obj$gl %>% mutate(IsCosmo = ASV %in% cosmo_asvs,                                    # Determine if ASV is cosmopolitan
                              ClusterNum = as.character(asv_max_cluster_vec[ASV]),              # Get cluster number from asv_max_cluster_vec
                              tag = ifelse(IsCosmo, "COSM", cl_labels[as.integer(ClusterNum)]), # Assign tag: cluster label or COSM for cosmopolitans
                              label = paste0(label, " [", tag, "]")) %>% select(-IsCosmo, -ClusterNum)
  obj
}

build_tick_df <- function(gl_tbl, ph_cols) {                          # community labels
  ph_cols_ext <- ensure_other_color(ph_cols)
  gl_tbl %>% distinct(asv_x, ASV, label, Group) %>%
    mutate(ph = if_else(as.character(Group) %in% names(ph_cols_ext), as.character(Group), "Other"), col_hex = to_hex(unname(ph_cols_ext[ph]))) %>%
    arrange(asv_x)
}


# ========================= # Plots (dot-panels draw their own thin grids) # =========================
create_dot_panel <- function(obj, title, ph_cols, label_angle = LABEL_ANGLE, label_size = LABEL_SIZE, label_y_at = LABEL_Y_AT, label_vjust = LABEL_VJUST) {
  d <- obj$d
  vm <- get_visual_mapping(d$RelAb)
  d <- d %>% mutate(mapped_alpha = vm$alpha, mapped_size = vm$size, shape_var = vm$shape)
  xmx <- max(d$asv_x, na.rm = TRUE) + 0.5
  ticks <- build_tick_df(obj$gl, ph_cols)
  xgrid_df <- tibble(x = sort(unique(ticks$asv_x)))
  ygrid_df <- tibble(y = seq_len(max(d$ix, na.rm = TRUE)))
  y_top <- max(ygrid_df$y) + 0.5
  p <- ggplot(d) +
    geom_segment(data = xgrid_df, aes(x = x, xend = x, y = 0.5, yend = max(ygrid_df$y) + 0.5), inherit.aes = FALSE, color = "grey85", linewidth = 0.25) +
    geom_segment(data = ygrid_df, aes(y = y, yend = y, x = 0.5, xend = xmx), inherit.aes = FALSE, color = "grey85", linewidth = 0.25) +
    geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", linewidth = 0.25, linetype = "solid") +
    geom_point(aes(x = asv_x, y = ix, color = Group, alpha = mapped_alpha, size = mapped_size, shape = shape_var), na.rm = TRUE) +
    geom_text(data = ticks, aes(x = asv_x, y = label_y_at, label = label, colour = I(col_hex)),
              angle = label_angle, hjust = 1, vjust = label_vjust, size = label_size / ggplot2::.pt,
              inherit.aes = FALSE, position = position_nudge(x = LABEL_XNUDGE, y = LABEL_DNUDGE)) +
    # geom_text(data = ticks, aes(x = asv_x, y = y_top, label = Community),
    #           fontface = "bold.italic", size = label_size / ggplot2::.pt, inherit.aes = FALSE, vjust = -0.3) +       # community labels if you want
    scale_shape_identity() +
    scale_color_manual(values = ph_cols, limits = names(ph_cols), drop = FALSE, guide = "none") +
    scale_alpha_identity() + scale_size_identity() +
    scale_x_continuous(limits = c(0.5, xmx), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, y_top), expand = c(0, 0)) +
    labs(title = title) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 10) +
    theme(text = element_text(color = "black"),
          plot.title.position = "panel",
          plot.title = element_text(hjust = 0.5, size = LABEL_SIZE + 8, face = "bold", margin = margin(t = 5, r = 0, b = 7, l = 0)),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          legend.position = "none",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  add_bottom_margin(p)
}

# ========================= # Legends (compact helpers) # =========================
make_col_phylum_legend <- function(phylum_names, color_map, title_text, num_cols) {
  n <- length(phylum_names)
  items_per_col <- ceiling(n / num_cols)
  legend_data <- tibble(Phylum = phylum_names,
                        x_dot   = rep(seq(1, by = 1, length.out = num_cols), each = items_per_col)[1:n] - 0.25,
                        x_label = rep(seq(1, by = 1, length.out = num_cols), each = items_per_col)[1:n] - 0.20,
                        y = rep(items_per_col:1, length.out = n))
  ggplot(legend_data, aes(x = x_dot, y = y)) +
    geom_point(aes(color = Phylum), size = 10 * S, show.legend = FALSE) +
    geom_text(aes(x = x_label, y = y, label = Phylum), hjust = 0, vjust = 0.5, size = 13 * S) +
    scale_color_manual(values = color_map, na.translate = FALSE) +
    labs(title = title_text) +
    theme_void() +
    theme(plot.title = element_text(size = LABEL_SIZE + 6, face = "bold", hjust = 0.5), plot.margin = margin(0, 60 * S, 0, 60 * S)) + #incr. margin to part
    scale_x_continuous(limits = c(0.5, max(legend_data$x_label) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, max(legend_data$y) + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off")
}

create_clean_legend_with_dots <- function() {
  # Only show phyla that actually occur in the dotplot for this figure
  euk_names  <- intersect(c(names(euk_colors_dot), "Other"),  euk_phyla_used)
  prok_names <- intersect(c(names(prok_colors_dot), "Other"), prok_phyla_used)
  euk_legend_plot  <- make_col_phylum_legend(euk_names, euk_colors_dot_extended[euk_names], "Major Eukaryotic Groups", 2)
  prok_legend_plot <- make_col_phylum_legend(prok_names, prok_colors_dot_extended[prok_names], "Major Prokaryotic Groups", 2)
  vm <- get_visual_mapping(c(0.01, 0.05, 0.10, 0.25, 0.30, 0.42, 0.62))
  ab_df <- tibble(RelAb = c("≤1%", "5%", "10%", "25%", "26-35%", "36-50%", "51-75%"),
                  alpha_val = vm$alpha, size_val = vm$size, shape_val = vm$shape,
                  x = c(1.01, 1.26, 1.47, 1.68, 1.99, 2.30, 2.61), y = 2)
  
  abundance_legend_plot <- ggplot(ab_df, aes(x = x, y = y)) +
    geom_point(aes(alpha = alpha_val, size = size_val, shape = shape_val), color = "black") +
    geom_text(aes(label = RelAb), nudge_y = -2.5, size = 12 * S, vjust = 0.5) +
    scale_alpha_identity() + scale_size_identity() + scale_shape_identity() +
    labs(title = "Relative Abundance") +
    theme_void() +
    theme(plot.title.position = "panel",
          plot.title = element_text(size = LABEL_SIZE + 6, face = "bold", hjust = 0.5, vjust = 1, margin = margin(t = 0, b = 0)),
          plot.margin = margin(0, 0, 0, 0)) +
    scale_x_continuous(limits = c(0.8, 2.8), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1.0, 2.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off")
  centers <- c(0.19, 0.63); y_mid <- 0.825; pad <- 0.03; box_sz <- 9 * S
  tr_df <- tibble(kind = c("Eukaryota", "Prokaryota"), fill = c("black", "lightgrey"), cx = centers, cy = y_mid, xtext = centers + pad)
  
  total_reads_plot <- ggplot(tr_df) +
    geom_point(aes(x = cx, y = cy, fill = fill), shape = 22, size = box_sz, stroke = 1, color = "black") +
    geom_text(aes(x = xtext, y = cy, label = kind), hjust = 0, vjust = 0.4, size = 13 * S) +
    labs(title = "Total Reads") +
    scale_fill_identity() +
    theme_void() +
    theme(plot.title = element_text(size = LABEL_SIZE + 6, face = "bold", hjust = 0.5, vjust = 0.5, margin = margin(t = 15)),
          plot.margin = margin(0, 40, 0, 40)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.0, 1.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off")
  right_column <- abundance_legend_plot / total_reads_plot + plot_layout(heights = c(1, 1))
  euk_legend_plot + plot_spacer() + prok_legend_plot + plot_spacer() + right_column +
    plot_layout(nrow = 1, widths = c(1.35, 0, 1.4, 0.25, 1.7)) # legend, spacer, legend, spacer, legend
}

# ========================= # Build ONE dataset row (minimal wrapper) # =========================
build_main_row <- function(proptab, thr_euk, thr_prok) {
  
  # --- ASV inclusion based on mean RelAb within affinity clusters (cosmo = all samples) ---
  aff_map <- indic_asv_tbl %>% transmute(ASV, Aff = as.character(Clusters))
  samp_clus <- meta.wk.small %>% rownames_to_column("Sample") %>% transmute(Sample, Cluster = as.character(Cluster))
  
  # Long table ONLY for computing the mean used for filtering
  long_sel <- proptab %>% as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
    left_join(all_asvs_tax %>% select(ASV, Domain2), by = "ASV") %>%
    left_join(samp_clus, by = "Sample") %>%
    left_join(aff_map, by = "ASV") %>%
    rowwise() %>%
    mutate(in_aff = is.na(Aff) || Cluster %in% str_split(Aff, "\\s*\\+\\s*")[[1]]) %>%
    ungroup()
  
  # Select ASVs using mean RelAb computed on:
  #   - indicators: only their affinity cluster(s)
  #   - cosmopolitans (Aff is NA): all samples
  asv_selection <- long_sel %>%
    filter(in_aff) %>%
    group_by(ASV, Domain2) %>%
    summarise(meanRelAb = mean(RelAb, na.rm = TRUE), .groups = "drop") %>%
    filter(meanRelAb >= if_else(Domain2 == "Eukaryote", thr_euk, thr_prok)) %>%
    select(ASV, Domain2)
  
  # Full data for plotting (NOT restricted to affinity clusters)
  dat_all <- proptab %>% as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
    left_join(all_asvs_tax %>% select(ASV, Domain2, Group), by = "ASV") %>%
    left_join(meta.wk.small %>% rownames_to_column("Sample") %>%
                transmute(Sample, Cluster = as.character(Cluster)), by = "Sample") %>%
    mutate(Genus = str_extract(ASV, "(?<=; )[^;]*(?=; ASV_\\d+$)")) %>%
    semi_join(asv_selection, by = c("ASV", "Domain2"))
  
  dat_plot <- dat_all %>% filter(RelAb > 0)
  
  # Compute per-ASV cluster with maximum mean abundance (avoids dplyr slice issues)
  asv_max_clus <- dat_all %>%
    group_by(ASV, Cluster) %>%
    summarise(mean_rel = mean(RelAb, na.rm = TRUE), .groups = "drop") %>%
    group_by(ASV) %>%
    summarise(
      MaxCluster = Cluster[which.max(mean_rel)],
      MaxRelAb   = max(mean_rel, na.rm = TRUE),
      .groups = "drop")
  
  dd <- dendro_data(dend_row_order, type = "rectangle")
  n <- nrow(dd$labels)
  
  # Build idx_lut: get positions from dendrogram, then join with cluster from metadata
  idx_lut <- tibble(Sample = dd$labels$label, ix = max(dd$labels$x) + 1 - dd$labels$x) %>%
    left_join(meta.wk.small %>% rownames_to_column("Sample") %>% select(Sample, Cluster), by = "Sample")
  
  dat_all <- dat_all %>% left_join(idx_lut, by = "Sample") %>% mutate(label = paste0(Genus, " (", str_extract(ASV, "ASV_\\d+$"), ")"))
  dat_plot <- dat_plot %>% left_join(idx_lut, by = "Sample") %>% mutate(label = paste0(Genus, " (", str_extract(ASV, "ASV_\\d+$"), ")"))
  cluster_ranges <<- idx_lut %>% group_by(Cluster) %>%
    summarise(ymin = min(ix) - 0.5, ymax = max(ix) + 0.5, y_center = (min(ix) + max(ix)) / 2, .groups = "drop")
  
  # Panels
  E <- prepare_domain_data(dat_plot %>% filter(Domain2 == "Eukaryote"), dat_all %>% filter(Domain2 == "Eukaryote"),
                           names(euk_colors_dot), euk_colors_dot, asv_max_clus) %>% add_indicator_info()
  P <- prepare_domain_data(dat_plot %>% filter(Domain2 == "Prokaryote"), dat_all %>% filter(Domain2 == "Prokaryote"),
                           names(prok_colors_dot), prok_colors_dot, asv_max_clus) %>% add_indicator_info()
  
  # Which phyla actually appear in this figure (incl. 'Other' only if used)
  used_euk_phyla  <- sort(unique(as.character(E$d$Group)))
  used_prok_phyla <- sort(unique(as.character(P$d$Group)))
  
  euk_cols_ext <- ensure_other_color(E$extended_colors)
  prok_cols_ext <- ensure_other_color(P$extended_colors)
  p_euk <- create_dot_panel(E, "Eukaryota", ph_cols = euk_cols_ext)
  p_prok <- create_dot_panel(P, "Prokaryota", ph_cols = prok_cols_ext)
  
  # Total Reads (stacked bars)
  bar_wide <- idx_lut %>% mutate(Sample_key = normalize_sample(Sample)) %>%
    left_join(meta.wk.small %>% rownames_to_column("Sample") %>% mutate(Sample_key = normalize_sample(Sample)) %>%
                select(Sample_key, `Eukaryotic Reads`, `Prokaryotic Reads`), by = "Sample_key") %>%
    mutate(`Eukaryotic Reads` = as.numeric(coalesce(`Eukaryotic Reads`, 0)),
           `Prokaryotic Reads` = as.numeric(coalesce(`Prokaryotic Reads`, 0)),
           y0 = ix - 0.4, y1 = ix + 0.4)
  bar_segments <- bind_rows(transmute(bar_wide, ix, xmin = 0, xmax = `Eukaryotic Reads`, ymin = y0, ymax = y1, read_type = "Eukaryotes"),
                            transmute(bar_wide, ix, xmin = `Eukaryotic Reads`, xmax = `Eukaryotic Reads` + `Prokaryotic Reads`,
                                      ymin = y0, ymax = y1, read_type = "Prokaryotes"))
  p_stacked_bar <- ggplot(bar_segments) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = read_type), color = "black") +
    scale_fill_manual(values = c("Eukaryotes" = "black", "Prokaryotes" = "lightgrey")) +
    geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.25, linetype = "solid") +
    scale_x_continuous(breaks = c(0, 50000, 100000), labels = label_total_reads_ticks,
                       expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    labs(x = "Total Reads", y = NULL) +
    guides(fill = "none") +
    theme_minimal() +
    theme(text = element_text(color = "black"),
          axis.line.x = element_line(linewidth = 0.25, color = "black"),
          axis.ticks.x = element_line(linewidth = 0.25, color = "black"),
          panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 27 * S, color = "black"),
          axis.title.x = element_text(size = 29 * S, face = "bold", margin = margin(t = -60, r = 0, b = 0, l = 0), vjust = -0.25))
  p_stacked_bar <- add_bottom_margin(p_stacked_bar)
  
  # Dendrogram
  dd_segs <- dd$segments
  xmax <- max(dd_segs$y, dd_segs$yend, na.rm = TRUE)
  seg_tree <- dd_segs %>% mutate(leaf_x = x, leaf_xend = xend) %>%
    transmute(x = TREE_SCALE * (xmax - y), xend = TREE_SCALE * (xmax - yend), y = n + 1 - leaf_x, yend = n + 1 - leaf_xend)
  p_tree <- ggplot(seg_tree) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, lineend = "round", na.rm = TRUE, color = "black") +
    geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", linewidth = 0.25, linetype = "solid", inherit.aes = FALSE) +
    geom_text(data = cluster_ranges, aes(x = -0.7, y = y_center, label = cl_labels[as.integer(as.character(Cluster))]),
              color = "black", size = 14 * S, fontface = "bold", inherit.aes = FALSE, hjust = 0.5, vjust = 0.5) +
    # annotate("text",  x = xmax - 0.2, y = n + 0.4, label = "Community", hjust = 1, vjust = -0.3, fontface = "bold.italic", size = 10) + # left/top/right-adj
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-1.5, TREE_SCALE * xmax), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() + theme(plot.margin = margin(t = 20 * S, r = 0, b = BOTTOM_MARGIN_PTS, l = 0))
  p_tree <- add_bottom_margin(p_tree)
  
  # Side-panels
  panel_df <- idx_lut %>% left_join(meta.wk.small %>% rownames_to_column("Sample") %>%
                                      select(Sample, Yr = Year, Stn = Station,
                                             # `C (km)` = `Distance to coast [km]`,
                                             `D (m)` = `Depth [m]` 
                                             # S = `Salinity`, `T (°C)` = `Temperature [deg C]`
                                      ), by = "Sample") %>%
    mutate(across(c(`D (m)`
                    # S, `T (°C)`
    ), ~ format(round(as.numeric(.), 1), nsmall = 1, trim = TRUE)))
  
  create_component_plot <- function(component, title_txt, box_right = 0.8, pad_rel = 0.05) {
    pad_data <- box_right * pad_rel
    ggplot(panel_df) +
      geom_rect(data = cluster_ranges, aes(xmin = -Inf, xmax = box_right, ymin = ymin, ymax = ymax),
                fill = NA, color = "black", linewidth = 0.25, linetype = "solid", inherit.aes = FALSE) +
      geom_text(aes(y = ix, label = !!sym(component)), x = box_right - pad_data, hjust = 1, size = 7.8 * S, color = "black") +
      scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
      scale_x_continuous(limits = c(0, box_right), expand = c(0, 0)) +
      coord_cartesian(clip = "on") +
      labs(title = title_txt) +
      theme_void() +
      theme(plot.title = element_text(size = LABEL_SIZE-1, face = "bold", hjust = 0.5, vjust = 0, margin = margin(b = 5)))
  }
  side_specs <- c("Yr", "Stn",
                  # "C (km)",
                  "D (m)"
                  # "S", "T (°C)"
  )
  width_map <- c("Yr" = W_YEAR, "Stn" = W_STATION,
                 # "C (km)" = W_COAST,
                 "D (m)" = W_DEPTH
                 # "S" = W_SAL, "T (°C)" = W_TEMP
  )
  side_plots <- side_specs %>% set_names() %>% map(~ add_bottom_margin(create_component_plot(.x, .x)))
  plots_row <- list(p_tree, p_euk, plot_spacer(), p_prok, plot_spacer(), p_stacked_bar)
  widths_row <- c(W_TREE, W_EUK, GAP_MED, W_PROK, GAP_SMALL, W_READS)
  for (nm in names(side_plots)) {
    plots_row <- append(plots_row, list(plot_spacer(), side_plots[[nm]]))
    widths_row <- c(widths_row, GAP_SMALL, unname(width_map[nm]))
  }
  list(
    plot_row        = wrap_plots(plots_row, nrow = 1) + plot_layout(widths = widths_row),
    euk_cols_ext    = euk_cols_ext,
    prok_cols_ext   = prok_cols_ext,
    euk_phyla_used  = used_euk_phyla,
    prok_phyla_used = used_prok_phyla,
    E               = E,
    P               = P
  )
}

# ========================= # Build & SAVE: separate PNGs for indic and not-indic # =========================

# --- Saving the plots ---
save_plot <- function(p, path, w, h, dpi){
  if (tolower(tools::file_ext(path)) == "pdf")
    cairo_pdf(path, width = w, height = h, onefile = FALSE)
  else
    png(path, width = w, height = h, units = "in", res = dpi, bg = "white", type = "cairo")
  print(p); dev.off()
}

# --- INDIC + COSMO ---
row_indic <- build_main_row(proptab.indic.cosmo, thr_euk = thr_euk_indic, thr_prok = thr_prok_indic)
E <- row_indic$E
P <- row_indic$P
euk_colors_dot_extended <- row_indic$euk_cols_ext
prok_colors_dot_extended <- row_indic$prok_cols_ext
euk_phyla_used       <- row_indic$euk_phyla_used
prok_phyla_used      <- row_indic$prok_phyla_used
legend_plot <- create_clean_legend_with_dots()
final_indic <- row_indic$plot_row / plot_spacer() / legend_plot + plot_layout(heights = c(10, LEGEND_GAP_HT, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
save_plot(final_indic, "OUT/5_INDIC/dot_ASV_indic_cosmo2.pdf", width_inches, height_inches, res_dpi)
save_plot(final_indic, "OUT/5_INDIC/dot_ASV_indic_cosmo2.png", width_inches, height_inches, res_dpi)

# # --- NON-INDIC & NON-COSMO ---
# row_not_indic <- build_main_row(proptab.not.indic.cosmo, thr_euk = thr_euk_not_indic, thr_prok = thr_prok_not_indic)
# euk_colors_dot_extended <- row_not_indic$euk_cols_ext
# prok_colors_dot_extended <- row_not_indic$prok_cols_ext
# euk_phyla_used       <- row_not_indic$euk_phyla_used
# prok_phyla_used      <- row_not_indic$prok_phyla_used
# legend_plot <- create_clean_legend_with_dots()
# final_not_indic <- row_not_indic$plot_row / plot_spacer() / legend_plot + plot_layout(heights = c(10, LEGEND_GAP_HT, 1)) +
#   plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
# save_plot(final_not_indic, "OUT/5_INDIC/dot_ASV_not_indic_or_cosmo2.pdf", width_inches, height_inches, res_dpi)


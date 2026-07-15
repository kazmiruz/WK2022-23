# WK_6---------------------------------------------------------------------
# Purpose:
#   Produce the faceted cluster map, T‑S-D and fChl-a‑fCDOM-Attenuation bubble plots, the cluster environmental summary table,
#   and the maps displaying distribution of PP taxa in surface waters (diatoms, dinoflagellates, non‑diatom algae, cyanobacteria).
#
# Prerequisites (objects must be in the workspace):
#   meta.wk.small, proptab.small.taxan.saman, cl_labels, cluster_palette, k (from WK_3_CLUS_TITAN_FT)
#   ne_10m/ (folder containing Natural Earth shapefiles)
#
# Outputs:
#   OUTPUT/wk_maps.pdf / .png                 (faceted cluster map)
#   OUTPUT/wk_ts.pdf / .png                   (T‑S-D bubble plot)
#   OUTPUT/wk_fc.pdf / .png                   (fChl-a‑fCDOM-Attenuation bubble plot)
#   OUTPUT/cluster_summary.csv                (env‑parameter summary table)
#   OUTPUT/diatoms_chaetoceros...pdf / .png   (surface diatom map)
#   OUTPUT/dinos_gyrodinium...pdf / .png      (surface dinoflagellate map)
#   OUTPUT/algae_picochlorum...pdf / .png     (surface non‑diatom algae map)
#   OUTPUT/cyano_surface.pdf / .png           (surface cyanobacteria map)
# -------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(sf)        # sf::sf_use_s2(FALSE); st_bbox, st_read, st_transform, st_geometry_type, st_make_valid, st_buffer, st_crop
library(ggspatial) # for scale and arrow
library(ggrepel)
library(scales)
library(vegan)     # for statistics (Kruskal)
library(dunn.test) # Dunn's test

# Divide samples by year and depth groups for plotting
set.seed(161193)   # reproducible jitter offsets for overlapping points (runif below); WK_6 may be run standalone
meta_grouped <- meta.wk.small %>% rownames_to_column("Sample") %>%
  mutate(`Longitude [deg]` = as.numeric(`Longitude [deg]`),
         `Latitude [deg]`  = as.numeric(`Latitude [deg]`),
         Depth_group = case_when(`Depth [m]` >= 1  & `Depth [m]` <= 6  ~ "Surface (1-6 m)",
                                 `Depth [m]` >= 9.7 & `Depth [m]` <= 75   ~ "Subsurface (10-75 m)"),
         Depth_group = factor(Depth_group, levels = c("Surface (1-6 m)", "Subsurface (10-75 m)"))) %>%
  # Round Lat/Long to 1 decimal to detect close points
  mutate(lon_round = round(`Longitude [deg]`, 1),
         lat_round = round(`Latitude [deg]`, 1)) %>%
  group_by(Depth_group, Year, lon_round, lat_round) %>%
  mutate(n_close = n()) %>%
  ungroup() %>%
  # Apply jitter only to groups with >1 close point
  mutate(jitter_lon = if_else(n_close > 1, `Longitude [deg]` + runif(n(), -0.1, 0.1), `Longitude [deg]`),
         jitter_lat = if_else(n_close > 1, `Latitude [deg]` + runif(n(), -0.1, 0.1), `Latitude [deg]`)) %>%
  select(-lon_round, -lat_round, -n_close)

# 0) one-time for this session: do all shapefile ops with GEOS
sf::sf_use_s2(FALSE)

# 1) your window as a bbox OBJECT (rectangular crop)s
bb   <- st_bbox(c(xmin = -95.7, xmax = -74.7, ymin = 50.6, ymax = 59), crs = 4326)
xlim <- c(bb["xmin"], bb["xmax"]); ylim <- c(bb["ymin"], bb["ymax"])

# 2) tiny, robust reader: read -> WGS84 -> repair polygons -> rectangular crop
rd <- function(path, bb) {
  if (!file.exists(path)) return(NULL)
  g <- suppressWarnings(st_read(path, quiet = TRUE)) |> st_transform(4326)
  if (any(st_geometry_type(g) %in% c("POLYGON","MULTIPOLYGON"))) {
    g <- suppressWarnings(st_make_valid(g))
    g <- suppressWarnings(st_buffer(g, 0))  # fixes degenerate/duplicate vertices
  }
  suppressWarnings(st_crop(g, bb)) |>
    dplyr::select(geometry, any_of(c("name", "name_en")))
}

# 3) read layers
ne <- "ne_10m"
land   <- rd(file.path(ne, "ne_10m_land.shp"), bb)
coast  <- rd(file.path(ne, "ne_10m_coastline.shp"), bb)
isles  <- rd(file.path(ne, "ne_10m_minor_islands.shp"), bb)
islesc <- rd(file.path(ne, "ne_10m_minor_islands_coastline.shp"), bb)
r1     <- rd(file.path(ne, "ne_10m_rivers_lake_centerlines_scale_rank.shp"), bb)
r2     <- rd(file.path(ne, "ne_10m_rivers_north_america.shp"), bb)
l1     <- rd(file.path(ne, "ne_10m_lakes.shp"), bb)
l2     <- rd(file.path(ne, "ne_10m_lakes_north_america.shp"), bb)
bathy  <- rd(file.path(ne, "ne_10m_bathymetry_L_0.shp"), bb)

# 4) highlight specific rivers
highlight_rivers <- function(river_sf, pattern) {
  if (is.null(river_sf)) return(NULL)
  # Natural Earth uses 'name' or 'name_en' – adapt if needed
  nm_col <- intersect(c("name", "name_en"), names(river_sf))[1]
  if (is.na(nm_col)) {
    # If no name column, can't filter; return empty
    return(river_sf[0, ])
  }
  river_sf[grepl(pattern, river_sf[[nm_col]], ignore.case = TRUE), ]
}
riv_pattern <- "Churchill|Nelson|Hayes|Winisk|Moose|Maquatua|La Grande"
r1_highlight <- highlight_rivers(r1, riv_pattern)
r2_highlight <- highlight_rivers(r2, riv_pattern)
l1_highlight <- highlight_rivers(l1, riv_pattern)
l2_highlight <- highlight_rivers(l2, riv_pattern)

# faceted cluster map
p_maps <- ggplot() +
  geom_sf(data = bathy, fill = "aliceblue", linewidth = 0.4) +
  geom_sf(data = land, fill = "grey100",  linewidth = 0.4) +
  geom_sf(data = isles, linewidth = 0.4) +
  geom_sf(data = coast, color = "black", linewidth = 0.4) +
  geom_sf(data = islesc, color = "black", linewidth = 0.4) +
  geom_sf(data = l1, fill = "cyan3", color = "cyan3", linewidth = 0.4) +
  geom_sf(data = l2, fill = "cyan3", color = "cyan3", linewidth = 0.4) +
  geom_sf(data = r1, color = "cyan3", linewidth = 0.4) +
  geom_sf(data = r2, color = "cyan3", linewidth = 0.4) +
  geom_sf(data = l1_highlight, fill = "navy", color = "navy", linewidth = 1.2) +
  geom_sf(data = l2_highlight, fill = "navy", color = "navy", linewidth = 1.2) +
  geom_sf(data = r1_highlight, color = "navy", linewidth = 1.2) +
  geom_sf(data = r2_highlight, color = "navy", linewidth = 1.2) +
  # Layer 1: thick border in cluster colour (no fill)
  geom_point(data = meta_grouped, aes(x = jitter_lon, y = jitter_lat, color = factor(Cluster)),
             shape = 21, size = 6.0, stroke = 2, fill = NA) +
  # Layer 2: semi‑transparent interior fill (no stroke) – colour is not used for border
  geom_point(data = meta_grouped, aes(x = jitter_lon, y = jitter_lat, fill = factor(Cluster)),
             shape = 21, size = 6.0, stroke = 0, alpha = 1.0) +
  geom_text(aes(x = x, y = y, label = river), data = data.frame(river = c("Churchill R.", "Nelson R.", "Hayes R.", "Winisk R.", "Moose R.", "Maquatua R.", "La Grande R."),
                                                                x = c(-94.0, -94.6, -92.6, -86.0, -82.5, -76.2, -76.2), y = c(57.7, 56.15, 55.5, 54.5, 50.90, 53.2, 54.4)),
            inherit.aes = FALSE, size = 5.5, fontface = "bold", color = "navy") +
  geom_label_repel(data = meta_grouped, aes(x = jitter_lon, y = jitter_lat, label = Station),
                  color = "black", fill = alpha(c("white"), 0.9), size = 5.6, fontface = "bold", box.padding = 0.3, point.padding = 0.3,
                  segment.size = 0.6, max.overlaps = Inf, force = 15, force_pull = 50, max.iter = 10000, min.segment.length = 0.05, seed = 161193) +
  scale_color_manual(values = cluster_palette, guide = "none") +   # border uses same palette, no extra legend
  scale_fill_manual(values = cluster_palette, labels = function(x) setNames(cl_labels, 1:k)[as.character(x)], name = "Cluster") +
  guides(fill = guide_legend(byrow = TRUE, ncol = 1, 
                             keyheight = unit(1.6, "cm"), keywidth = unit(1.6, "cm"),
                             override.aes = list(shape = 21, size = 9, stroke = 2, colour = NA))) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_grid(Depth_group ~ Year, labeller = labeller(Year = \(x) paste("Year", x)), switch = "y") +
  annotation_scale(location = "bl", text_cex = 1.5) +
  annotation_north_arrow(location = "tr", which_north = "true", height = unit(2.3, "cm"), width = unit(2.3, "cm")) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "grey85", linewidth = 0.2), panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"),
        axis.title = element_blank(), axis.text  = element_text(size = 20, colour = "black"), axis.ticks = element_line(color = "black"),
        strip.placement = "outside", strip.text.x = element_text(size = 33, face = "bold"), strip.text.y = element_text(size = 33, face = "bold"),
        legend.title = element_text(size = 33, face = "bold"), legend.text  = element_text(size = 30, face = "bold"), legend.position = "right",
        legend.key.height = unit(1.6, "cm"), legend.key.width  = unit(1.6, "cm"), legend.spacing.y= unit(0.8, "cm"), legend.box.spacing = unit(1.0, "cm"))

cairo_pdf("OUTPUT/wk_maps.pdf", width = 24, height = 16)
print(p_maps)
dev.off()

png("OUTPUT/wk_maps.png", width = 24, height = 16, units = "in", res = 300)
print(p_maps)
dev.off()


################
### TSD plot ###

p_ts <- ggplot(meta_grouped, aes(`Salinity [PSU]`, `Temperature [deg C]`)) +
  scale_x_continuous(limits = c(0, 31), trans = exp_trans(1.15),  # exp transformation
                     breaks = c(0, 10, 15, 20, 25, 30, 31), minor_breaks = NULL) +
  scale_y_continuous(limits = c(-1, 15), breaks = c(-1, 0, 3, 6, 9, 12, 15), minor_breaks = NULL) +
  geom_point(aes(color = factor(Cluster), shape = factor(Year)), fill = NA, size = 7, stroke = 2.5) +
  geom_point(aes(color = factor(Cluster), shape = factor(Year), alpha = `Depth [m]`, fill = after_scale(colour)), size = 7, stroke = 0) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0))) +
  scale_color_manual(values = cluster_palette, labels = \(x) cl_labels[as.integer(x)], name = "Cluster") +
  scale_alpha_continuous(name = "Depth (m)", range = c(0.1, 1), limits = range(meta_grouped$`Depth [m]`), breaks = c(1, 5, 20, 40, 75)) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2, override.aes = list(size = 7, shape = 21, stroke = 2.4),
                              keyheight = unit(1.1, "cm"), keywidth = unit(1.1, "cm"), default.unit = "cm", byrow = TRUE),
         alpha = guide_legend(order = 3, override.aes = list(shape = 21, size = 7, stroke = 0, color = "black", fill = "black"))) +
  # facet_grid(. ~ Year, labeller = labeller(Year = ~ paste("Year", .x))) +
  labs(x = "Salinity", y = "Temperature (°C)") +
  theme_minimal() +
  theme(axis.text  = element_text(size = 23, colour = "black"), axis.title = element_text(size = 28, face = "bold"),
        strip.text.x = element_text(size = 31, face = "bold"), legend.text  = element_text(size = 28, face = "bold"),
        legend.title = element_text(size = 28, face = "bold"), legend.position = "right", legend.spacing.y = unit(0.3, "cm"),
        plot.margin = margin(t = 18, r = 36, b = 18, l = 36),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), panel.spacing = unit(0.6, "lines"))

grDevices::cairo_pdf("OUTPUT/wk_ts.pdf", width = 24, height = 11); print(p_ts); dev.off()

png("OUTPUT/wk_ts.png", width = 24, height = 11, units = "in", res = 300)
print(p_ts)
dev.off()


#################
### Fluo plot ###

meta_grouped_FC <- meta_grouped %>% filter(!is.na(`Beam attenuation [1/m]`))
p_fc <- ggplot(meta_grouped_FC, aes(`CDOM fluorescence [mg/m^3]`, `Chlorophyll-a fluorescence [mg/m^3]`)) +
  scale_x_continuous(limits = c(1, 100), trans = log_trans(),  # log transformation
                     breaks = c(1, 2, 5, 10, 25, 50, 100), minor_breaks = NULL) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 11), minor_breaks = NULL) +
  geom_point(aes(color = factor(Cluster), shape = factor(Year)), fill = NA, size = 7, stroke = 2.5) +
  geom_point(aes(color = factor(Cluster), shape = factor(Year), alpha = `Beam attenuation [1/m]`, fill = after_scale(colour)), size = 7, stroke = 0) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0))) +
  scale_color_manual(values = cluster_palette, labels = \(x) cl_labels[as.integer(x)], name = "Cluster") +
  scale_alpha_continuous(
    name = expression(atop(textstyle(bold("Beam")), atop(textstyle(bold("Attenuation")), textstyle(bold("(m"^{-1}*")"))))), 
    range = c(0.1, 1), limits = c(0.1, 6.5), breaks = c(0.2, 0.5, 1, 3, 6.5)) + 
  guides(shape = guide_legend(order = 1), 
         color = guide_legend(order = 2, override.aes = list(size = 7, shape = 21, stroke = 2.4), keyheight = unit(1.1, "cm"), keywidth = unit(1.1, "cm"), default.unit = "cm", byrow = TRUE), 
         alpha = guide_legend(order = 3, override.aes = list(shape = 21, size = 7, stroke = 0, color = "black", fill = "black"))) + 
  labs(x = expression(bold(CDOM~fluorescence~(mg~m^{-3}))), 
       y = expression(bold(Chlorophyll-a~fluorescence~(mg~m^{-3})))) +
  # facet_grid(. ~ Year, labeller = labeller(Year = ~ paste("Year", .x))) +
  theme_minimal() +
  theme(axis.text  = element_text(size = 23, colour = "black"),
        axis.title = element_text(size = 28, face = "bold"),
        strip.text.x = element_text(size = 31, face = "bold"),
        legend.text  = element_text(size = 28, face = "bold"),
        legend.title = element_text(size = 28, face = "bold"),
        legend.position = "right",
        legend.spacing.y = unit(0.3, "cm"),
        plot.margin = margin(t = 18, r = 36, b = 18, l = 36),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
        panel.spacing = unit(0.6, "lines"))

grDevices::cairo_pdf("OUTPUT/wk_fc.pdf", width = 24, height = 11); print(p_fc); dev.off()

png("OUTPUT/wk_fc.png", width = 24, height = 11, units = "in", res = 300, type = "cairo")
print(p_fc)
dev.off()

# median chl a in surface waters
# median(meta.wk.small %>% filter(`Depth [m]` <= 6) %>% .$`Chlorophyll-a fluorescence [mg/m^3]`, na.rm = TRUE)
# meta.wk.small %>% filter(`Depth [m]` <= 6, !is.na(`Chlorophyll-a fluorescence [mg/m^3]`)) %>% nrow()


#############################################################################################################
### Summary table for the clusters - includes env-l parameters found to be significant in db-RDA and more ###
# What env-l properties characterize each cluster?

### Environmental cluster stats
env_vars <- c("Depth [m]", "Distance to coast [km]", "Temperature [deg C]", "Salinity [PSU]",
              "Chlorophyll-a fluorescence [mg/m^3]", "CDOM fluorescence [mg/m^3]", "Beam attenuation [1/m]")

# Median (range) summary function
med_range_txt <- \(v) {
  n <- length(v); m <- sum(!is.na(v))
  if (m/n < 0.5) return(NA_character_) # at least half of the samples should have readings, otherwise NA
  sprintf("%.1f\n[%.1f-%.1f]%s",
          stats::median(v, na.rm = TRUE),
          min(v, na.rm = TRUE),
          max(v, na.rm = TRUE),
          if (m < n) "*" else "")
}

# Normality check histograms (red - fail b/c low Shapiro p; green - pass b/c high Shapiro p)
par(mfrow = c(2, 2))
walk(env_vars, \(x) {
  v <- meta.wk.small[[x]]
  if (sum(!is.na(v)) >= 3 && n_distinct(meta.wk.small$Cluster[!is.na(v)]) > 1) {
    fit <- aov(v ~ factor(meta.wk.small$Cluster))
    p <- shapiro.test(residuals(fit))$p.value
    hist(residuals(fit), main = paste("Resid:", x), xlab = paste("Shapiro p =", round(p, 4)),
         col = ifelse(p < 0.05, "salmon", "lightgreen"))
  }
})
par(mfrow = c(1, 1))

# Kruskal & Dunn tests - for each env-l parameter, which cluster has significantly higher/lower readings
test_var <- \(var) {
  d <- meta.wk.small %>% filter(!is.na(.data[[var]]))
  if (n_distinct(d$Cluster) < 2) return(list(kw_p = NA_real_, dunn_res = NULL))
  kw_p <- kruskal.test(d[[var]] ~ d$Cluster)$p.value
  dunn_res <- if (kw_p < 0.05) {
    capture.output(dt <- dunn.test(d[[var]], d$Cluster, method = "bh", kw = FALSE, list = TRUE), file = nullfile())
    med <- tapply(d[[var]], d$Cluster, \(z) median(z, na.rm = TRUE))
    comps <- str_split(dt$comparisons, " - ", simplify = TRUE)
    tibble(group1 = comps[,1], group2 = comps[,2], p_val = dt$P.adjusted) %>%
      filter(p_val < 0.05) %>%
      mutate(Result = ifelse(med[group1] > med[group2],
                             paste(cl_labels[as.numeric(group1)], ">", cl_labels[as.numeric(group2)]),
                             paste(cl_labels[as.numeric(group2)], ">", cl_labels[as.numeric(group1)]))) %>%
      select(Result, p_val)
  } else NULL
  list(kw_p = kw_p, dunn_res = dunn_res)
}

test_results <- map(env_vars, test_var) %>% setNames(env_vars)
bind_rows(test_results)

# Create summary table
cluster_summary <- meta.wk.small %>%
  mutate(Cluster = as.character(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(`Most frequent\nsample type` = names(which.max(table(Code))),
            across(all_of(env_vars), med_range_txt, .names = "{.col}"), .groups = "drop") %>%
  arrange(as.numeric(Cluster)) %>%
  mutate(Cluster = cl_labels[as.numeric(Cluster)]) %>%
  bind_rows(tibble(Cluster = "Kruskal-Wallis p") %>%
              bind_cols(map_dfc(env_vars, \(x) {
                v <- meta.wk.small[[x]]
                g <- meta.wk.small$Cluster
                keep <- !is.na(v)
                tibble(!!x := if (n_distinct(g[keep]) < 2) NA_character_
                       else formatC(kruskal.test(v[keep] ~ g[keep])$p.value, format = "e", digits = 1))
              }))) %>%
  rename_with(~ str_replace_all(., "(.*?) \\[(.*?)\\]", "\\1\n(\\2)"))

write.csv(cluster_summary, "OUTPUT/cluster_summary.csv", row.names = FALSE)


###########################################################################################

# Make special maps for the surface waters alone for diatoms, dinos, non-diatom algae, and cyanos.

# ============================================================
# 1. Base theme
# ============================================================
big_black_theme <- theme_minimal() +
  theme(text = element_text(color = "black", size = 18),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title = element_text(color = "black", size = 20, face = "bold"),
        axis.text = element_text(color = "black", size = 16),
        legend.title = element_text(color = "black", size = 20, face = "bold"),
        legend.text = element_text(color = "black", size = 20),
        strip.text = element_text(color = "black", size = 20, face = "bold"),
        panel.grid.major = element_line(color = "grey85", linewidth = 0.2))

# ============================================================
# 2. Conditional jitter for overlapping map points
# ============================================================
add_conditional_jitter <- function(df) {
  set.seed(161193)   # reproducible jitter offsets (runif) on every call
  df %>%
    mutate(lon_round = round(`Longitude [deg]`, 1),
           lat_round = round(`Latitude [deg]`, 1)) %>%
    group_by(lon_round, lat_round) %>%
    mutate(n_close = n()) %>%
    ungroup() %>%
    mutate(jitter_lon = if_else(n_close > 1, `Longitude [deg]` + runif(n(), -0.1, 0.1), `Longitude [deg]`),
           jitter_lat = if_else(n_close > 1, `Latitude [deg]` + runif(n(), -0.1, 0.1), `Latitude [deg]`)) %>%
    select(-lon_round, -lat_round, -n_close)
}

# ============================================================
# 3. Plotting function: 2022 map + 2023 map above full-width T-S plot
# ============================================================
plot_group <- function(data, group_col, colors, title, filename,
                       max_size, upper_limit, breaks, labels, legend_labels = waiver()) {
  
  make_map <- function(df, year_lab) {
    ggplot() +
      geom_sf(data = bathy, fill = "aliceblue", linewidth = 0.2) +
      geom_sf(data = land, fill = "grey90", linewidth = 0.2) +
      geom_sf(data = coast, color = "black", linewidth = 0.4) +
      geom_sf(data = l1, fill = "cyan2", color = "cyan2", linewidth = 0.2) +
      geom_sf(data = l2, fill = "cyan2", color = "cyan2", linewidth = 0.2) +
      geom_sf(data = r1, color = "cyan2", linewidth = 0.2) +
      geom_sf(data = r2, color = "cyan2", linewidth = 0.2) +
      geom_point(data = df,
                 aes(x = jitter_lon, y = jitter_lat, size = RelAb, color = !!sym(group_col)),
                 shape = 21, fill = NA, stroke = 2, alpha = 0.9) +
      scale_size_area(max_size = max_size, limits = c(0, upper_limit),
                      breaks = breaks, labels = labels, name = "Relative abundance") +
      scale_color_manual(values = colors, labels = legend_labels, name = title) +
      coord_sf(xlim = c(-94.5, -76.5), ylim = c(50.8, 59), expand = FALSE) +
      labs(title = year_lab) +
      big_black_theme +
      theme(axis.title = element_blank(),
            plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
            legend.position = "none")
  }
  
  p_map_2022 <- make_map(data %>% filter(Year == 2022), "Year 2022")
  p_map_2023 <- make_map(data %>% filter(Year == 2023), "Year 2023")
  
  p_ts <- ggplot(data, aes(x = `Salinity [PSU]`, y = `Temperature [deg C]`,
                           size = RelAb, color = !!sym(group_col))) +
    geom_point(shape = 21, fill = NA, stroke = 2, alpha = 0.9,
               position = position_jitter(width = 0.1, height = 0.1)) +
    scale_size_area(max_size = max_size, limits = c(0, upper_limit),
                    breaks = breaks, labels = labels, name = "Relative abundance") +
    scale_color_manual(values = colors, labels = legend_labels, name = title) +
    scale_x_continuous(expand = c(0, 0.5), limits = c(-0.5, 30.5)) +
    scale_y_continuous(expand = c(0, 0.5), limits = c(-0.5, 15.5)) +
    labs(x = "Salinity", y = "Temperature (°C)") +
    big_black_theme +
    guides(color = guide_legend(title.position = "top", direction = "vertical", order = 1,
                                override.aes = list(fill = colors, shape = 21, stroke = 1, size = 5)),
           size = guide_legend(title.position = "top", direction = "horizontal",
                               label.position = "bottom", order = 2)) +
    theme(legend.position = "bottom",
          legend.box = "horizontal",
          legend.box.just = "top",
          legend.spacing.x = unit(5, "cm"),
          legend.key.width = unit(2.5, "cm"))
  
  fig <- (p_map_2022 + p_map_2023) / p_ts + plot_layout(widths = c(1.2, 1.2), heights = c(1, 0.5))
  
  ggsave(paste0("OUTPUT/", filename, ".pdf"), fig, width = 22, height = 15, device = cairo_pdf)
  ggsave(paste0("OUTPUT/", filename, ".png"), fig, width = 22, height = 15)
}

# ============================================================
# 4. Eukaryotic groups
# ============================================================
groups <- list(
  list(genera = c("Chaetoceros", "Minidiscus", "Arcocellulus"),
       colors = c("Chaetoceros" = "red2", "Minidiscus" = "green3", "Arcocellulus" = "blue2"),
       title = "Diatoms", file = "diatoms_chaetoceros_minidiscus_arcocellulus",
       max_size = 15, lim = 35, breaks = c(1, 5, 25, 35), labels = c("1%", "5%", "25%", "35%")),
  list(genera = c("Gyrodinium", "Heterocapsa", "Gymnodinium"),
       colors = c("Gyrodinium" = "red2", "Heterocapsa" = "green3", "Gymnodinium" = "blue2"),
       title = "Dinoflagellates", file = "dinos_gyrodinium_heterocapsa_gymnodinium",
       max_size = 15 * sqrt(40 / 35), lim = 40, breaks = c(1, 5, 25, 40), labels = c("1%", "5%", "25%", "40%")),
  list(genera = c("Picochlorum", "Micromonas", "Aureococcus"),
       colors = c("Picochlorum" = "red2", "Micromonas" = "green3", "Aureococcus" = "blue2"),
       title = "Non-diatom algae", file = "algae_picochlorum_micromonas_aureococcus",
       max_size = 15 * sqrt(55 / 35), lim = 55, breaks = c(1, 5, 25, 55), labels = c("1%", "5%", "25%", "55%"))
)

for (g in groups) {
  d <- proptab.small.taxan.saman %>%
    as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
    mutate(Genus = sapply(strsplit(ASV, "; "), `[`, 4)) %>%
    filter(Genus %in% g$genera) %>%
    group_by(Sample, Genus) %>%
    summarise(RelAb = sum(RelAb), .groups = "drop") %>%
    left_join(meta.wk.small %>%
                rownames_to_column("Sample") %>%
                mutate(Depth_group = ifelse(`Depth [m]` >= 1 & `Depth [m]` <= 6, "Surface", NA)),
              by = "Sample") %>%
    filter(!is.na(Depth_group), RelAb > 0) %>%
    mutate(RelAb = RelAb * 100,
           Genus = factor(Genus, levels = g$genera),
           Year = as.numeric(Year)) %>%
    add_conditional_jitter()
  
  plot_group(d, "Genus", g$colors, g$title, g$file,
             g$max_size, g$lim, g$breaks, g$labels)
}

# ============================================================
# 5. Cyanobacteria
# ============================================================
target_asvs <- c("Bacteria; Cyanobacteriota; Synechococcales; Synechococcus CC9902; ASV_17",
                 "Bacteria; Cyanobacteriota; Synechococcales; Cyanobium PCC-6307; ASV_221",
                 "Bacteria; Cyanobacteriota; Synechococcales; Cyanobium PCC-6307; ASV_210")

cols <- setNames(c("red2", "green3", "blue2"), target_asvs)

cyano_labels <- str_replace(target_asvs, "^.*; (.*); (ASV_\\d+)$", "\\1 (\\2)")

d <- proptab.small.taxan.saman[, target_asvs] %>%
  as_tibble(rownames = "Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
  left_join(meta.wk.small %>%
              rownames_to_column("Sample") %>%
              mutate(Depth_group = ifelse(`Depth [m]` >= 1 & `Depth [m]` <= 6, "Surface", NA)),
            by = "Sample") %>%
  filter(!is.na(Depth_group), RelAb > 0) %>%
  mutate(ASV = factor(ASV, levels = target_asvs),
         RelAb = RelAb * 100,
         Year = as.numeric(Year)) %>%
  add_conditional_jitter()

plot_group(d, "ASV", cols, "Cyanobacteria", "cyano_surface",
           max_size = 15, upper_limit = 15,
           breaks = c(0.1, 1, 5, 15),
           labels = c("0.1%", "1%", "5%", "15%"),
           legend_labels = cyano_labels)


# # Explore: at what T or S does a certain ASV reach max Relab? what is the distribution?
# ASV = "ASV_17"
# meta.wk.small %>% rownames_to_column("Sample") %>% filter(`Depth [m]` %in% 1:6) %>%
#   left_join(proptab.small.taxan.saman %>% as.data.frame() %>% rownames_to_column("Sample") %>%
#               select(Sample, asv = matches(paste0(ASV, "(\\D|$)"))), by = "Sample") %>%
#   {
#     # occurrence range (all surface samples with RelAb > 0)
#     cat("Occurrence range (RelAb > 0):\n")
#     filter(., asv > 0) %>% summarise(Temp_range = list(range(`Temperature [deg C]`, na.rm = TRUE)),
#                                      Sal_range  = list(range(`Salinity [PSU]`, na.rm = TRUE))) %>%
#       unlist() %>% print()
#     
#     # at which T and S doe it reach maximum Relab? ----
#     cat("\nAt maximum relative abundance:\n")
#     filter(., asv == max(asv, na.rm = TRUE)) %>% select(`Temperature [deg C]`, `Salinity [PSU]`) %>% print()
#     
#     # plot your ASV vs T and S for Relab > 0
#     pos <- filter(., asv > 0)   # drop zeros
#     par(mfrow = c(1, 2))
#     plot(pos$asv, pos$`Temperature [deg C]`, xlab = ASV, ylab = "Temperature (°C)", log = "x")
#     plot(pos$asv, pos$`Salinity [PSU]`, xlab = ASV, ylab = "Salinity (PSU)", log = "x")
#     par(mfrow = c(1, 1))}


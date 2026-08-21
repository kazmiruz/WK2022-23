# 6------------------------------------------------------------------------
# Purpose:
#   Display stations/samples, including by cluster, and relative abundance of interesting ASVs in surface waters.
#   Produce faceted cluster maps, the map displaying what have already been sampled in HBC and where my survey fits,
#   and the maps displaying distribution of PP taxa in surface waters (diatoms, non‑diatom algae, dinoflagellates, and cyanobacteria).
#
# Prerequisite: run 3W_CLUS first, as it generate various objects used here.
#
# Required external files:
#   EXT/ne_10m/                            (includes various shape and associated files for creating maps)
#   EXT/marmap/                            (bathymetry)
#   EXT/HBC_stations.csv                   (HBC stations sampled in this and previous studies)
#
# Outputs:
#   OUT/6W_MAPS/wk_maps.png                        (faceted cluster map)
#   OUT/6W_MAPS/HBC_station_map.png                (map of stations sampled in HBC)
#   OUT/6W_MAPS/diatoms_surface.png                (surface diatom map)
#   OUT/6W_MAPS/dinos_surface.png                  (surface dinoflagellate map)
#   OUT/6W_MAPS/algae_surface.png                  (surface non‑diatom algae map)
#   OUT/6W_MAPS/cyano_surface.png                  (surface cyanobacteria map)
# -------------------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(cowplot)
library(sf)        # sf::sf_use_s2(FALSE); st_bbox, st_read, st_transform, st_geometry_type, st_make_valid, st_buffer, st_crop
library(ggspatial) # for scale and arrow
library(ggnewscale)
library(ggrepel)
library(marmap)
library(scales)
library(vegan)     # for statistics (Kruskal)
library(dunn.test) # Dunn's test


# Divide samples by year and depth groups for plotting

set.seed(161193)   # reproducible jitter offsets for overlapping points: 0.1 degrees aka ~10 km is necessary to see all overlapping points
meta_grouped <- meta.wk %>% rownames_to_column("Sample") %>%
  mutate(`Longitude [deg]` = as.numeric(`Longitude [deg]`), `Latitude [deg]` = as.numeric(`Latitude [deg]`),
         Depth_group = case_when(`Depth [m]` >= 1 & `Depth [m]` <= 6 ~ "Surface (1-6 m)",
                                 `Depth [m]` >= 9.7 & `Depth [m]` <= 75 ~ "Subsurface (10-75 m)"),
         Depth_group = factor(Depth_group, levels = c("Surface (1-6 m)", "Subsurface (10-75 m)")),
         lon_round = round(`Longitude [deg]`, 1), lat_round = round(`Latitude [deg]`, 1)) %>%    # "similar" locations within roughly 0.1°
  group_by(Depth_group, Year, lon_round, lat_round) %>% mutate(n_close = n()) %>% ungroup() %>%
  mutate(jitter_lon = if_else(n_close > 1, `Longitude [deg]` + runif(n(), -0.1, 0.1), `Longitude [deg]`),
         jitter_lat = if_else(n_close > 1, `Latitude [deg]` + runif(n(), -0.1, 0.1), `Latitude [deg]`)) %>%
  select(-lon_round, -lat_round, -n_close)

# Map extent and shapefiles
sf::sf_use_s2(FALSE)
bb <- st_bbox(c(xmin = -95.7, xmax = -74.7, ymin = 50.6, ymax = 59), crs = 4326)
xlim <- c(bb["xmin"], bb["xmax"]); ylim <- c(bb["ymin"], bb["ymax"])

rd <- function(path, bb) {
  if (!file.exists(path)) return(NULL)
  g <- suppressWarnings(st_read(path, quiet = TRUE)) %>% st_transform(4326)
  if (any(st_geometry_type(g) %in% c("POLYGON", "MULTIPOLYGON")))
    g <- suppressWarnings(st_buffer(st_make_valid(g), 0))
  suppressWarnings(st_crop(g, bb)) %>% dplyr::select(geometry, any_of(c("name", "name_en")))
}

ne <- "EXT/ne_10m"
land <- rd(file.path(ne, "ne_10m_land.shp"), bb)
coast <- rd(file.path(ne, "ne_10m_coastline.shp"), bb)
isles <- rd(file.path(ne, "ne_10m_minor_islands.shp"), bb)
islesc <- rd(file.path(ne, "ne_10m_minor_islands_coastline.shp"), bb)
r1 <- rd(file.path(ne, "ne_10m_rivers_lake_centerlines_scale_rank.shp"), bb)
r2 <- rd(file.path(ne, "ne_10m_rivers_north_america.shp"), bb)
l1 <- rd(file.path(ne, "ne_10m_lakes.shp"), bb)
l2 <- rd(file.path(ne, "ne_10m_lakes_north_america.shp"), bb)

# Download once if the file does not already exist
# marmap::getNOAA.bathy(lon1 = as.numeric(bb["xmin"]), lon2 = as.numeric(bb["xmax"]),
#                       lat1 = as.numeric(bb["ymin"]), lat2 = as.numeric(bb["ymax"]),
#                       resolution = 0.5, keep = TRUE, path = "EXT/marmap")

# Bathymetry classes use exactly the same palette as the larger map
bathy_labs <- c("0–25", "25–50", "50–100", "100–250", "≥250")
bathy_cols <- setNames(c("#D7EAF3", "#B8D3E0", "#95B7C7", "#7198A9", "#4D7180"), bathy_labs)

bathy_file <- file.path("EXT/marmap", paste0("marmap_coord_", as.numeric(bb["xmin"]), ";",
                                             as.numeric(bb["ymin"]), ";",
                                             as.numeric(bb["xmax"]), ";",
                                             as.numeric(bb["ymax"]), "_res_0.5.csv"))

bathy_df <- fortify.bathy(marmap::read.bathy(bathy_file, header = TRUE)) %>% filter(z <= 0) %>%
  mutate(depth = -z, depth_class = cut(depth, breaks = c(0, 25, 50, 100, 250, Inf), labels = bathy_labs, include.lowest = TRUE, right = FALSE))
# range(bathy_df$depth, na.rm = TRUE) # depth range check


# Highlight selected rivers and lakes
highlight_rivers <- function(x, pattern) {
  if (is.null(x)) return(NULL)
  nm_col <- intersect(c("name", "name_en"), names(x))[1]
  if (is.na(nm_col)) return(x[0, ])
  x[grepl(pattern, x[[nm_col]], ignore.case = TRUE), ]
}

riv_pattern <- "Churchill|Nelson|Hayes|Winisk|Moose|Maquatua|La Grande"
r1_highlight <- highlight_rivers(r1, riv_pattern)
r2_highlight <- highlight_rivers(r2, riv_pattern)
l1_highlight <- highlight_rivers(l1, riv_pattern)
l2_highlight <- highlight_rivers(l2, riv_pattern)

# Faceted cluster map
options(ggspatial.scale_bar_warn = FALSE)
meta_labels <- meta_grouped %>% distinct(Year, Depth_group, Station, .keep_all = TRUE)

river_labels <- data.frame(
  river = c("Churchill R.", "Nelson R.", "Hayes R.", "Winisk R.", "Moose R.", "Maquatua R.", "La Grande R."),
  x = c(-94.0, -94.6, -92.6, -86.0, -82.5, -76.2, -76.2),
  y = c(57.7, 56.15, 55.5, 54.5, 50.90, 53.2, 54.4))

p_maps <- ggplot() +
  geom_tile(data = bathy_df, aes(x = x, y = y, fill = depth_class)) +
  scale_fill_manual(values = bathy_cols, breaks = levels(bathy_df$depth_class),
                    drop = TRUE, na.translate = FALSE, name = "Depth (m)",
                    guide = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                                         keyheight = unit(1.6, "cm"), keywidth = unit(1.6, "cm"))) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = land, fill = "grey100", colour = NA, linewidth = 0.4) +
  geom_sf(data = isles, fill = "grey100", colour = NA, linewidth = 0.4) +
  geom_sf(data = coast, colour = "black", linewidth = 0.4) +
  geom_sf(data = islesc, colour = "black", linewidth = 0.4) +
  geom_sf(data = l1, fill = "#9ecae1", colour = "#9ecae1", linewidth = 0.4) +
  geom_sf(data = l2, fill = "#9ecae1", colour = "#9ecae1", linewidth = 0.4) +
  geom_sf(data = r1, colour = "#9ecae1", linewidth = 0.4) +
  geom_sf(data = r2, colour = "#9ecae1", linewidth = 0.4) +
  geom_sf(data = l1_highlight, fill = "#08519c", colour = "#08519c", linewidth = 0.9) +
  geom_sf(data = l2_highlight, fill = "#08519c", colour = "#08519c", linewidth = 0.9) +
  geom_sf(data = r1_highlight, colour = "#08519c", linewidth = 0.9) +
  geom_sf(data = r2_highlight, colour = "#08519c", linewidth = 0.9) +
  geom_point(data = meta_grouped,
             aes(jitter_lon, jitter_lat, colour = factor(Cluster), fill = factor(Cluster)),
             shape = 21, size = 4.5, stroke = 1.2, colour = "black") +
  geom_text(data = river_labels, aes(x, y, label = river), inherit.aes = FALSE,
            size = 5.5, fontface = "bold", colour = "#08299c") +
  geom_label_repel(data = meta_labels, aes(jitter_lon, jitter_lat, label = Station),
                   colour = "black", fill = alpha("white", 0.9), size = 5.5, fontface = "bold",
                   box.padding = 0.3, point.padding = 0.2, segment.size = 0.6,
                   max.overlaps = Inf, force = 20, force_pull = 50, max.iter = 10000,
                   min.segment.length = 0.05, seed = 161193) +
  scale_color_manual(values = cluster_palette, guide = "none") +
  scale_fill_manual(values = cluster_palette,
                    labels = \(x) setNames(cl_labels, 1:k)[as.character(x)], name = "Cluster") +
  guides(fill = guide_legend(order = 2, title.position = "top", title.hjust = 0.5, ncol = 1,
                             keyheight = unit(1.6, "cm"), keywidth = unit(1.6, "cm"),
                             override.aes = list(shape = 21, size = 9, stroke = 1.2, colour = "black"))) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  facet_grid(Depth_group ~ Year, labeller = labeller(Year = \(x) paste("Year", x)), switch = "y") +
  annotation_scale(location = "bl", text_cex = 1.5) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(2.3, "cm"), width = unit(2.3, "cm")) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.2),
        panel.grid.minor = element_blank(), panel.spacing = unit(2, "lines"),
        axis.title = element_blank(), axis.text = element_text(size = 18, colour = "black"),
        axis.ticks = element_line(colour = "black"), strip.placement = "outside",
        strip.text.x = element_text(size = 30, face = "bold"),
        strip.text.y = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 28, face = "bold"),
        legend.position = "right", legend.key.height = unit(1.6, "cm"),
        legend.key.width = unit(1.6, "cm"), legend.box = "vertical",
        legend.box.just = "center", legend.spacing.y = unit(1, "cm"))

# cairo_pdf("OUT/6W_MAPS/wk_maps.pdf", width = 23, height = 16)
# print(p_maps)
# dev.off()

png("OUT/6W_MAPS/wk_maps.png", width = 23, height = 16, units = "in", res = 300)
print(p_maps)
dev.off()


################################################
### Hudson Bay Complex sampling coverage map ###
################################################

sf::sf_use_s2(FALSE)

land_col <- "grey90"
sampling_cols <- c("July 2017" = "#F05ACD", "June–July 2018" = "#9B5DE5", "August 2018" = "#18BE5B",
                   "July 2022" = "#F2CC00", "August 2022; August 2023" = "#FF6040")

eDNA_shapes <- c(Eukaryotes = 24, Prokaryotes = 21, `All organisms` = 22)
eDNA_sizes <- c(Eukaryotes = 2.9, Prokaryotes = 4, `All organisms` = 4)

period <- c("2017-07" = "July 2017",
            "2018-06" = "June–July 2018", "2018-07" = "June–July 2018",
            "2018-08" = "August 2018", "2022-07" = "July 2022",
            "2022-08" = "August 2022; August 2023",
            "2023-08" = "August 2022; August 2023")

# Published and current stations
hb <- bind_rows(read.csv("EXT/HBC_stations.csv") %>% transmute(Date = as.Date(Date), eDNA, Latitude, Longitude),
                meta.wk %>% transmute(Date = as.Date(paste(Year, Month, Day, sep = "-")), eDNA = "All organisms",
                                            Latitude = `Latitude [deg]`, Longitude = `Longitude [deg]`)) %>%
  mutate(Sampling = period[format(Date, "%Y-%m")], lon0 = round(Longitude, 4), lat0 = round(Latitude, 4))

# Combine duplicate coordinates and jitter similar stations
set.seed(161193)   # reproducible jitter offsets for overlapping points: 0.20 degrees aka ~20 km is necessary to see all overlapping points
hb_all <- hb %>% group_by(Sampling, lon0, lat0) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude),
            eDNA = if (n_distinct(eDNA) == 1) (eDNA)[1] else "All organisms", .groups = "drop") %>%
  mutate(lon_round = round(Longitude, 1), lat_round = round(Latitude, 1)) %>%                    # "similar" locations within roughly 0.1°
  group_by(lon_round, lat_round) %>% mutate(n_close = n()) %>% ungroup() %>%
  mutate(plot_lon = if_else(n_close > 1, Longitude + runif(n(), -0.2, 0.2), Longitude),
         plot_lat = if_else(n_close > 1, Latitude + runif(n(), -0.2, 0.2), Latitude)) %>%
  select(-lon_round, -lat_round, -n_close)

# Map extent
bb <- st_bbox(c(xmin = min(hb_all$Longitude) - 1, xmax = max(hb_all$Longitude) + 1,
                ymin = min(hb_all$Latitude) - 1, ymax = max(hb_all$Latitude) + 1), crs = 4326)
xlim <- bb[c("xmin", "xmax")]
ylim <- bb[c("ymin", "ymax")]
x_breaks <- seq(floor(xlim[1] / 5) * 5, ceiling(xlim[2] / 5) * 5, 5)
y_breaks <- seq(floor(ylim[1]), ceiling(ylim[2]), 2)

# Natural Earth layers
ne <- "EXT/ne_10m"
ne_files <- c(land = "ne_10m_land.shp", coast = "ne_10m_coastline.shp",
              isles = "ne_10m_minor_islands.shp", l1 = "ne_10m_lakes.shp",
              l2 = "ne_10m_lakes_north_america.shp",
              r1 = "ne_10m_rivers_lake_centerlines_scale_rank.shp",
              r2 = "ne_10m_rivers_north_america.shp")

read_ne <- \(file, bounds) suppressWarnings(st_read(file.path(ne, file), quiet = TRUE) %>% st_transform(4326) %>% st_make_valid() %>% st_crop(bounds))
L <- lapply(ne_files, read_ne, bounds = bb)

# Bathymetry
bb_bathy <- c(xmin = floor(as.numeric(bb["xmin"])), ymin = floor(as.numeric(bb["ymin"])),
              xmax = ceiling(as.numeric(bb["xmax"])), ymax = ceiling(as.numeric(bb["ymax"]))) # rounded out box

# Download once for the general HBC extent, then leave commented
# marmap::getNOAA.bathy(lon1 = as.numeric(bb_bathy["xmin"]), lon2 = as.numeric(bb_bathy["xmax"]),
#                       lat1 = as.numeric(bb_bathy["ymin"]), lat2 = as.numeric(bb_bathy["ymax"]),
#                       resolution = 0.5, keep = TRUE, path = "EXT/marmap")

bathy_labs <- c("0–25", "25–50", "50–100", "100–250", "≥250")
bathy_cols <- setNames(c("#D7EAF3", "#B8D3E0", "#95B7C7", "#7198A9", "#4D7180"), bathy_labs)

bathy_file <- file.path("EXT/marmap", paste0("marmap_coord_", as.numeric(bb_bathy["xmin"]), ";", as.numeric(bb_bathy["ymin"]), ";",
                                             as.numeric(bb_bathy["xmax"]), ";", as.numeric(bb_bathy["ymax"]), "_res_0.5.csv"))

bathy_df <- fortify.bathy(read.bathy(bathy_file, header = TRUE)) %>% filter(z <= 0) %>%
  mutate(depth = -z, depth_class = cut(depth, c(0, 25, 50, 100, 250, Inf), labels = bathy_labs, include.lowest = TRUE, right = FALSE))
# range(bathy_df$depth, na.rm = TRUE) # depth range check

# Densified bounding-box polygon for the Arctic projection
outline <- \(b, pad = c(0, 0, 0, 0), n = 200) {
  x <- c(b["xmin"] - pad[1], b["xmax"] + pad[2])
  y <- c(b["ymin"] - pad[3], b["ymax"] + pad[4])
  xy <- rbind(cbind(seq(x[1], x[2], length.out = n), y[1]),
              cbind(x[2], seq(y[1], y[2], length.out = n)),
              cbind(seq(x[2], x[1], length.out = n), y[2]),
              cbind(x[1], seq(y[2], y[1], length.out = n)))
  st_sfc(st_polygon(list(rbind(xy, xy[1, ]))), crs = 4326)
}

# Arctic inset
crs_arctic <- 3573
bb_arctic <- st_bbox(c(xmin = -180, xmax = 180, ymin = 35, ymax = 90), crs = 4326)
bb_view <- st_bbox(c(xmin = -179.9, xmax = -10, ymin = 47.5, ymax = 60), crs = 4326)

A <- lapply(ne_files[c("land", "isles", "coast")], read_ne, bounds = bb_arctic)
A <- lapply(A, st_transform, crs = crs_arctic)
view <- st_bbox(st_transform(outline(bb_view, n = 300), crs_arctic))
hb_box <- st_transform(st_sf(geometry = outline(bb, c(0.3, 0.3, 1, 0.2))), crs_arctic)

p_inset <- ggplot() +
  geom_sf(data = A$land, fill = land_col, colour = NA) +
  geom_sf(data = A$isles, fill = land_col, colour = NA) +
  geom_sf(data = A$coast, colour = "black", linewidth = 0.25) +
  geom_sf(data = hb_box, fill = NA, colour = "red2", linewidth = 1) +
  coord_sf(crs = crs_arctic,
           xlim = c(view["xmin"], view["xmax"] - 1000000),
           ylim = view[c("ymin", "ymax")], expand = FALSE, datum = NA) +
  theme_void() +
  theme(panel.background = element_rect(fill = "#dbeaf4", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black", linewidth = 1.2),
        plot.margin = margin(0, 0, 0, 0))

# Custom legend
sampling_key <- tibble(Sampling = names(sampling_cols), y = c(8.75, 7.35, 5.32, 3.92, 2.52),
                       Authors = c("Jacquemot et al. (2021; 2024)", "Jacquemot et al. (2021; 2022; 2024)\nMorency et al. (2022)",
                                   "Blais et al. (2022)", "Loutet et al. (2025)", "This study"))
edna_key <- tibble(eDNA = names(eDNA_shapes), y = c(8.75, 8, 7.25))
depth_key <- tibble(Depth = names(bathy_cols), y = c(4.9, 4.2, 3.5, 2.8, 2.1))

legend_panel <- ggplot() +
  annotate("text", x = 0.2, y = 9.55, label = "Sampling", hjust = 0, size = 6.5, fontface = "bold") +
  geom_point(data = sampling_key, aes(0.36, y, fill = Sampling), shape = 21, size = 4.5, stroke = 1.2, colour = "black") +
  geom_text(data = sampling_key, aes(0.66, y, label = Sampling), hjust = 0, size = 6) +
  geom_text(data = sampling_key, aes(0.66, y - 0.33, label = Authors), hjust = 0, vjust = 1, size = 6, lineheight = 1.02) +
  scale_fill_manual(values = sampling_cols, guide = "none") +
  annotate("text", x = 7.15, y = 9.55, label = "eDNA", hjust = 0, size = 6.5, fontface = "bold") +
  geom_point(data = edna_key, aes(7.35, y, shape = eDNA, size = eDNA), fill = NA, stroke = 1.5, colour = "black") +
  geom_text(data = edna_key, aes(7.58, y, label = eDNA), hjust = 0, size = 6) +
  scale_shape_manual(values = eDNA_shapes, guide = "none") +
  scale_size_manual(values = eDNA_sizes, guide = "none") +
  annotate("text", x = 7.15, y = 5.75, label = "Depth (m)", hjust = 0, size = 6.5, fontface = "bold") +
  ggnewscale::new_scale_fill() +
  geom_rect(data = depth_key, aes(xmin = 7.2, xmax = 7.62, ymin = y - 0.19, ymax = y + 0.19, fill = Depth), colour = "black", linewidth = 0.25) +
  geom_text(data = depth_key, aes(7.87, y, label = Depth), hjust = 0, size = 5.9) +
  scale_fill_manual(values = bathy_cols, guide = "none") +
  coord_cartesian(c(0, 10), c(1.65, 9.85), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.background = element_rect(fill = alpha("white", 0.92), colour = "black", linewidth = 0.9), plot.margin = margin(0, 0, 0, 0))

# Main map
p_hb <- ggplot() +
  geom_tile(data = bathy_df, aes(x, y, fill = depth_class)) +
  scale_fill_manual(values = bathy_cols, drop = FALSE, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = L$land, fill = land_col, colour = NA) +
  geom_sf(data = L$isles, fill = land_col, colour = NA) +
  geom_sf(data = L$l1, fill = "#9ecae1", colour = "#9ecae1", linewidth = 0.2) +
  geom_sf(data = L$l2, fill = "#9ecae1", colour = "#9ecae1", linewidth = 0.2) +
  geom_sf(data = L$r1, colour = "#6baed6", linewidth = 0.2) +
  geom_sf(data = L$r2, colour = "#6baed6", linewidth = 0.2) +
  geom_sf(data = L$coast, colour = "black", linewidth = 0.3) +
  geom_point(data = hb_all %>%                 # establish priority of points showing to display all stations as much as possible
               arrange(match(Sampling, c("August 2018", "July 2017", "June–July 2018", "July 2022", "August 2022; August 2023"))),
             aes(plot_lon, plot_lat, fill = Sampling, shape = eDNA, size = eDNA), stroke = 1.1, colour = "black") +
  scale_fill_manual(values = sampling_cols, guide = "none") +
  scale_shape_manual(values = eDNA_shapes, guide = "none") +
  scale_size_manual(values = eDNA_sizes, guide = "none") +         # smaller on the map than on the legend to see all stations
  scale_x_continuous(breaks = x_breaks, labels = \(x) paste0(abs(x), "°W"),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, labels = \(y) paste0(y, "°N"),
                     expand = c(0, 0)) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  annotation_scale(location = "bl", text_cex = 1.4) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         height = unit(1.8, "cm"), width = unit(1.8, "cm")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text = element_text(size = 15, colour = "black"),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 0, 0, 0)) +
  inset_element(p_inset, left = 0.005, bottom = 0.03, right = 0.305,  top = 0.275, align_to = "panel") +
  inset_element(legend_panel, left = 0.56, bottom = 0.001, right = 0.999, top = 0.36, align_to = "panel")

ggsave("OUT/6W_MAPS/HBC_station_map.png", p_hb, width = 15.2, height = 12, dpi = 300, bg = "white")


####################################################################
### Relative abundance of primary producers in the surface layer ###
####################################################################

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
# 2. Plotting function: 2022 map, 2023 map, and T-S plot
# ============================================================
plot_group <- function(data, group_col, colors, title, filename, max_size, upper_limit,
                       breaks, labels, legend_labels = waiver()) {
  map_xlim <- c(-95.7, -74.7); map_ylim <- c(50.6, 59)
  
  make_map <- function(df, year_lab) ggplot() +
    geom_sf(data = land, fill = "grey100", colour = NA, linewidth = 0.4) +
    geom_sf(data = isles, fill = "grey100", colour = NA, linewidth = 0.4) +
    geom_sf(data = coast, colour = "black", linewidth = 0.4) +
    geom_sf(data = islesc, colour = "black", linewidth = 0.4) +
    geom_sf(data = l1, fill = "#9ecae1", colour = "#9ecae1", linewidth = 0.4) +
    geom_sf(data = l2, fill = "#9ecae1", colour = "#9ecae1", linewidth = 0.4) +
    geom_sf(data = r1, colour = "#9ecae1", linewidth = 0.4) +
    geom_sf(data = r2, colour = "#9ecae1", linewidth = 0.4) +
    geom_sf(data = l1_highlight, fill = "#08519c", colour = "#08519c", linewidth = 0.9) +
    geom_sf(data = l2_highlight, fill = "#08519c", colour = "#08519c", linewidth = 0.9) +
    geom_sf(data = r1_highlight, colour = "#08519c", linewidth = 0.9) +
    geom_sf(data = r2_highlight, colour = "#08519c", linewidth = 0.9) +
    geom_point(data = df %>% arrange(desc(RelAb)), aes(`Longitude [deg]`, `Latitude [deg]`, size = RelAb, colour = .data[[group_col]]),
               shape = 21, fill = NA, stroke = 2, alpha = 0.9) +
    geom_text(data = river_labels, aes(x, y, label = river), inherit.aes = FALSE,
              size = 4.1, fontface = "bold", colour = "#08299c") +
    scale_size_area(max_size = max_size, limits = c(0, upper_limit), breaks = breaks,
                    labels = labels, name = "Relative abundance") +
    scale_color_manual(values = colors, labels = legend_labels, name = title) +
    scale_x_continuous(breaks = seq(-95, -75, 5), labels = \(x) paste0(abs(x), "°W")) +
    scale_y_continuous(breaks = seq(51, 59, 2), labels = \(y) paste0(y, "°N")) +
    coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
    labs(title = year_lab) +
    guides(colour = "none", size = "none") +
    theme_minimal(base_size = 18) +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.line = element_blank(),
          axis.title = element_blank(), axis.text = element_text(size = 13.5, colour = "black"),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.ticks.length = unit(0.14, "cm"),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 3)))
  
  p_map_2022 <- make_map(filter(data, Year == 2022), "Year 2022")
  p_map_2023 <- make_map(filter(data, Year == 2023), "Year 2023")
  
  p_ts <- ggplot(data, aes(Salinity, `Temperature [deg C]`, size = RelAb, colour = .data[[group_col]])) +
    geom_point(data = data %>% arrange(desc(RelAb)), shape = 21, fill = NA, stroke = 2, alpha = 0.9) +
    scale_size_area(max_size = max_size, limits = c(0, upper_limit), breaks = breaks,
                    labels = labels, name = "Relative abundance") +
    scale_color_manual(values = colors, labels = legend_labels, name = title) +
    scale_x_continuous(limits = c(-0.5, 30.5), expand = c(0, 0.5)) +
    scale_y_continuous(limits = c(-0.5, 15.5), expand = c(0, 0.5)) +
    labs(x = "Salinity", y = "Temperature (°C)") +
    big_black_theme +
    guides(colour = guide_legend(title.position = "top", title.hjust = 0,
                                 title.theme = element_text(size = 22.5, face = "bold", hjust = 0, margin = margin(b = 13)),
                                 order = 1, direction = "horizontal", nrow = 1, byrow = TRUE,
                                 keywidth = unit(0.85, "cm"), keyheight = unit(0.75, "cm"),
                                 override.aes = list(shape = 21, fill = NA, stroke = 1.8, size = 7.5)),
           size = guide_legend(title.position = "top", title.hjust = 0,
                               title.theme = element_text(size = 22.5, face = "bold", hjust = 0, margin = margin(b = 2)),
                               order = 2, direction = "horizontal", nrow = 1, byrow = TRUE,
                               keywidth = unit(0.45, "cm"), keyheight = unit(0.75, "cm"))) + theme(aspect.ratio = 0.72)
  
  legend_theme <- theme(legend.position = "bottom", legend.justification = "left",
                        legend.box.just = "left", legend.box = "horizontal",
                        legend.spacing.x = unit(4.5, "cm"), legend.spacing.y = unit(0.05, "cm"),
                        legend.box.spacing = unit(1, "cm"), legend.key.spacing.x = unit(0.6, "cm"),
                        legend.title = element_text(size = 22.5, face = "bold", hjust = 0),
                        legend.text = element_text(size = 19.5), legend.key.height = unit(0.75, "cm"),
                        legend.key = element_blank(), legend.background = element_blank(),
                        legend.box.background = element_blank(), legend.margin = margin(0),
                        legend.box.margin = margin(0), plot.margin = margin(0))
  
  fig <- ((p_map_2022 | plot_spacer() | p_map_2023 | plot_spacer() | p_ts) +
            plot_layout(widths = c(1, 0.05, 1, 0.06, 1.05), guides = "collect")) & legend_theme
  
  # ggsave(paste0("OUT/6W_MAPS/", filename, ".pdf"), fig, width = 24, height = 7.5, device = cairo_pdf, bg = "white")
  ggsave(paste0("OUT/6W_MAPS/", filename, ".png"), fig, width = 24, height = 7.5, dpi = 300, bg = "white")
}

# ============================================================
# 3. Eukaryotic groups
# ============================================================
groups <- list(
  list(genera = c("Chaetoceros", "Minidiscus", "Arcocellulus"),
       colors = c("Chaetoceros" = "red2", "Minidiscus" = "green3", "Arcocellulus" = "#8F00FF"),
       title = "Diatoms", file = "diatoms_surface",
       max_size = 15, lim = 35, breaks = c(1, 5, 25, 35), labels = c("1%", "5%", "25%", "35%")),
  list(genera = c("Gyrodinium", "Heterocapsa", "Gymnodinium"),
       colors = c("Gyrodinium" = "red2", "Heterocapsa" = "green3", "Gymnodinium" = "#8F00FF"),
       title = "Dinoflagellates", file = "dinos_surface",
       max_size = 15 * sqrt(40 / 35), lim = 40, breaks = c(1, 5, 25, 40), labels = c("1%", "5%", "25%", "40%")),
  list(genera = c("Picochlorum", "Micromonas", "Aureococcus"),
       colors = c("Picochlorum" = "red2", "Micromonas" = "green3", "Aureococcus" = "#8F00FF"),
       title = "Non-diatom algae", file = "algae_surface",
       max_size = 15 * sqrt(55 / 35), lim = 55, breaks = c(1, 5, 25, 55), labels = c("1%", "5%", "25%", "55%")))

for (g in groups) {
  d <- proptab.wk %>% as_tibble(rownames = "Sample") %>%
    pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
    mutate(Genus = sapply(strsplit(ASV, "; "), `[`, 4)) %>%
    filter(Genus %in% g$genera) %>% group_by(Sample, Genus) %>%
    summarise(RelAb = sum(RelAb), .groups = "drop") %>%
    left_join(meta.wk %>% rownames_to_column("Sample") %>%
                mutate(Depth_group = ifelse(`Depth [m]` >= 1 & `Depth [m]` <= 6, "Surface", NA)), by = "Sample") %>%
    filter(!is.na(Depth_group), RelAb > 0) %>%
    mutate(RelAb = RelAb * 100, Genus = factor(Genus, levels = g$genera), Year = as.numeric(Year))
  
  plot_group(d, "Genus", g$colors, g$title, g$file, g$max_size, g$lim, g$breaks, g$labels)
}

# ============================================================
# 4. Cyanobacteria
# ============================================================
target_asvs <- c("Bacteria; Cyanobacteriota; Synechococcales; Synechococcus CC9902; ASV_17",
                 "Bacteria; Cyanobacteriota; Synechococcales; Cyanobium PCC-6307; ASV_221",
                 "Bacteria; Cyanobacteriota; Synechococcales; Cyanobium PCC-6307; ASV_210")
cols <- setNames(c("red2", "green3", "#8F00FF"), target_asvs)
cyano_labels <- str_replace(target_asvs, "^.*; (.*); (ASV_\\d+)$", "\\1 (\\2)")

d <- proptab.wk[, target_asvs] %>% as_tibble(rownames = "Sample") %>%
  pivot_longer(-Sample, names_to = "ASV", values_to = "RelAb") %>%
  left_join(meta.wk %>% rownames_to_column("Sample") %>%
              mutate(Depth_group = ifelse(`Depth [m]` >= 1 & `Depth [m]` <= 6, "Surface", NA)), by = "Sample") %>%
  filter(!is.na(Depth_group), RelAb > 0) %>%
  mutate(ASV = factor(ASV, levels = target_asvs), RelAb = RelAb * 100, Year = as.numeric(Year))

plot_group(d, "ASV", cols, "Cyanobacteria", "cyano_surface",
           max_size = 15, upper_limit = 15, breaks = c(0.1, 1, 5, 15), labels = c("0.1%", "1%", "5%", "15%"), legend_labels = cyano_labels)

 
# # Explore: at what T or S does a certain ASV reach max RelAb? what is the distribution?
# ASV = "ASV_17"
# meta.wk %>% rownames_to_column("Sample") %>% filter(`Depth [m]` %in% 1:6) %>%
#   left_join(proptab.wk %>% as.data.frame() %>% rownames_to_column("Sample") %>%
#               select(Sample, asv = matches(paste0(ASV, "(\\D|$)"))), by = "Sample") %>%
#   {
#     # occurrence range (all surface samples with RelAb > 0)
#     cat("Occurrence range (RelAb > 0):\n")
#     filter(., asv > 0) %>% summarise(Temp_range = list(range(`Temperature [deg C]`, na.rm = TRUE)),
#                                      Sal_range  = list(range(`Salinity`, na.rm = TRUE))) %>%
#       unlist() %>% print()
#     
#     # at which T and S doe it reach maximum RelAb? ----
#     cat("\nAt maximum relative abundance:\n")
#     filter(., asv == max(asv, na.rm = TRUE)) %>% select(`Temperature [deg C]`, `Salinity`) %>% print()
#     
#     # plot your ASV vs T and S for RelAb > 0
#     pos <- filter(., asv > 0)   # drop zeros
#     par(mfrow = c(1, 2))
#     plot(pos$asv, pos$`Temperature [deg C]`, xlab = ASV, ylab = "Temperature (°C)", log = "x")
#     plot(pos$asv, pos$`Salinity`, xlab = ASV, ylab = "Salinity", log = "x")
#     par(mfrow = c(1, 1))}

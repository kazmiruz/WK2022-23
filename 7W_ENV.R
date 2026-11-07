# 7---------------------------------------------------------------------
# Purpose:
#   This last script generally assesses how environmental parameters drive the composition and diversity of microbial communities.
#   Perform distance‑based redundancy analysis (db‑RDA) on the fourth‑root transformed Bray‑Curtis dissimilarities
#   to resolve the influence of depth, distance to coast, temperature, and salinity on plankton community structure.
#   Produce a 2‑panel ordination plot (dbRDA1 vs dbRDA2 and dbRDA3 vs dbRDA2), a 3‑D interactive plot, an ANOVA table,
#   and a multi‑page diagnostic PDF with sample points sized by various environmental parameters.
#   Plot T‑S-D and fChl-a‑fCDOM-attenuation bubble plots, and produce tables comparing environments between clusters.
#   Put together MIMARKS and SRA tables for NCBI submission.
#
# Prerequisite: run 2W_TAX, 3W_CLUS, 4W_PLANKT, and 6W_MAPS first, as they generate various objects used here.
#
# Required external files:
#   DNA/WK_filelist.txt                                         (list of sequencing files for the SRA table)
#
# Outputs:
#   OUT/7W_ENV/wk_dbRDA_many.pdf                                (exploration: multi‑page diagnostic, samples sized by each environmental variable)
#   OUT/7W_ENV/wk_dbRDA.pdf / .png                              (Figure 4: 2‑panel ordination; optional .svg for Inkscape in commented code)
#   OUT/7W_ENV/wk_dbRDA_3d.html                                 (exploration: interactive 3‑D plot of dbRDA1-3)
#   OUT/7W_ENV/wk_ts.pdf / .png                                 (Figure 3A: T‑S-D bubble plot)
#   OUT/7W_ENV/wk_fc.pdf / .png                                 (Figure 3B: fChl-a‑fCDOM-attenuation bubble plot)
#   OUT/7W_ENV/Environmental_statistics.xlsx                    (Tables S1-S2: env properties of clusters and statistically sign-t diff-ces btwn them)
#   OUT/7W_ENV/NCBI_submission/WK_MIMARKS.survey.water.6.0.csv  (MIMARKS table for NCBI)
#   OUT/7W_ENV/NCBI_submission/WK_SRA_metadata.csv              (SRA table for NCBI)

# -------------------------------------------------------------------------

library(tidyverse)
library(vegan)      # for dbRDA and standardization: decostand(), rankindex(), dbrda(), vif.cca(), anova.cca(), RsquareAdj(), scores()
library(GGally)     # for ggpairs
library(ggrepel)    # to place labels nicely
library(scales)     # for rescale
library(grid)
library(gridExtra)
library(plotly)     # for 3D dbRDA
library(htmlwidgets)
library(dunn.test)
library(openxlsx)   # for writing complex Excel files

####################################################################################################################
### dbRDA (distance-based Redundancy Analysis; in this case, we will base it on the Bray Curtis Distance matrix) ###
####################################################################################################################

# dbRDA is used to identify how environmental variables influence microbial relative abundance.
# Vectors point in the direction of their greatest +ve influence and away from their greatest -ve influence.
# Many know PCA - an unconstrained form of ordination that searches for ANY variables best explaining the community dynamics.
# In contrast, RDA is constrained as it tries to explain the community dynamics with the explanatory (environmental) variables GIVEN.

# Info that helped me to put together a dbRDA pipeline below:
# https://youtu.be/Wqa7_18jf5A
# https://youtu.be/_bpQ16uUqs4
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# https://wiki.qcbs.ca/r_workshop10
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/rda-and-dbrda/
# https://sites.ualberta.ca/~ahamann/teaching/renr690/seminars/dbRDA.pdf
# https://vegandevs.github.io/vegan/reference/goodness.cca.html

#
## Select and transform environmental/explanatory and microbial/response data 
metadata <- meta.wk        # explanatory variables
microbial_data <- proptab.trans  # response variables (already normalized + transformed)

# Select numeric environmental variables (keep spaces; drop units in names)
env_vars <- metadata %>% select("Distance to coast [km]", "Depth [m]", "Salinity", "Temperature [deg C]") %>%
  rename_with(~ str_replace_all(., "\\s*\\[.*?\\]", "")) # remove units

# Inspect pairwise relationships
ggpairs(env_vars, lower = list(continuous = wrap("smooth", method = "lm", se = TRUE, color = "blue", size = 0.8)),
        diag  = list(continuous = wrap("barDiag", bins = 20, color = "cyan")),
        upper = list(continuous = wrap("cor", size = 3, color = "black")))

# Keep complete cases across env vars
keep <- complete.cases(env_vars)
metadata <- metadata[keep, , drop = FALSE]
microbial_data <- microbial_data[keep, , drop = FALSE]
env_vars <- env_vars[keep, , drop = FALSE]

# Standardize env vars
env_stand <- decostand(env_vars, method = "standardize")
apply(env_stand, 2, mean)  # check: must be ~0
apply(env_stand, 2, sd)    # check: must be ~1

## Distance index ranking (optional)
rankindex(env_stand, microbial_data, indices = c("euc","man","gow","bra","kul"), stepacross = FALSE, method = "pearson")

## Model
rda_model <- dbrda(microbial_data ~ ., data = env_stand, distance = "bray", dfun = vegdist)

## VIFs
vif.cca(rda_model) %>% as.data.frame() %>% rownames_to_column("Variable") %>% rename(VIF = ".") %>% mutate(VIF = round(VIF, 0)) %>% arrange(desc(VIF))

## ANOVAs (model, axes, vars)

# set.seed(161193)
# anova.cca(rda_model, permutations = 99999)                                    # test model sign-ce (*** = very significant)
#          Df  SumOfSqs     F  Pr(>F)    
# Model     4   5.4174 8.8772  1e-05 ***
# Residual 67  10.2218  

# set.seed(161193)
# anova.cca(rda_model, permutations = 99999, by = "axis")                       # test axes sign-ce
#   Df SumOfSqs       F Pr(>F)    
#   dbRDA1    1   2.8090 18.4121   1e-05 ***
#   dbRDA2    1   1.2811  8.5226   1e-05 ***
#   dbRDA3    1   0.9366  6.3225   1e-05 ***
#   dbRDA4    1   0.3906  2.6748 0.00727 ** 
#   Residual 67  10.2218  

# set.seed(161193)
# anova.cca(rda_model, permutations = 99999, by = "term") %>% arrange(desc(F))  # test var. sign-ce (biased towards the vars appearing first in the formula)
# Df SumOfSqs       F Pr(>F)    
#   Salinity             1   2.0717 13.5789  1e-05 ***
#   `Distance to coast`  1   1.2045  7.8947  1e-05 ***
#   Temperature          1   1.1566  7.5811  1e-05 ***
#   Depth                1   0.9847  6.4540  1e-05 ***
#   Residual            67  10.2218   

set.seed(161193)
anova_per_var <- anova.cca(rda_model, permutations = 99999, by = "margin") %>%  # test var. sign-ce if marginal (last in the formula) and store it neatly
  as.data.frame() %>% rownames_to_column("Variable") %>% filter(Variable != "Residual") %>% 
  transmute(Variable = gsub("`", "", Variable), F = round(F, 1), P = `Pr(>F)`,
            Sig = as.character(symnum(P, cutpoints = c(0, 0.00009, 0.0009, 0.0049, 0.009, 1), symbols = c("***", "**", "*", ".", " ")))) %>% arrange(desc(F))
#   Variable   F       P Sig
#   1          Salinity 8.3 0.00001 ***
#   2       Temperature 7.6 0.00001 ***
#   3 Distance to coast 7.1 0.00001 ***
#   4             Depth 4.3 0.00017  **

## Variance info
perc <- round(100 * (summary(rda_model)$cont$importance[2, 1:4]), 1) %>% print() # how much of the variation each axis explains
# dbRDA1 dbRDA2 dbRDA3 dbRDA4 
# 18.0    8.2    6.0    2.5 

RsquareAdj(rda_model)$adj.r.squared                                              # adjusted R2 (total variaiton explained but adjusted for number of vars)
# [1] 0.3075065

## Scores and arrow scaling (remove any backticks in bp labels)
sc_bp <- scores(rda_model, display = "bp", scaling = 2, choices = 1:rda_model$CCA$rank) %>%
  as.data.frame() %>% rownames_to_column("Variable") %>% mutate(Variable = gsub("`", "", Variable))
sc_si <- scores(rda_model, display = "sites", choices = 1:rda_model$CCA$rank) %>%
  as.data.frame() %>% cbind(metadata) %>% mutate(Cluster = as.factor(Cluster))
arrow_scaling_factor <- max(abs(as.matrix(sc_bp[grep("^dbRDA", names(sc_bp))])), na.rm = TRUE) * 1.7

## Exploration: dbRDA plots with one page per numeric variable sized onto points
pdf(file = "OUT/7W_ENV/wk_dbRDA_many.pdf", width = 35, height = 20)
numeric_vars <- meta.wk %>% select(Year, "Bottom depth [m]", "Depth [m]", "Temperature [deg C]", "Salinity",
                                   "Chlorophyll-a fluorescence [mg/m^3]", "CDOM fluorescence [mg/m^3]", "Beam attenuation [1/m]",
                                   "Dissolved oxygen [uM]", "Mixed layer depth [m]", "Distance to coast [km]", "DNA [mg/L]") %>% colnames()
size_range <- c(5, 30)

for (param in numeric_vars) {
  param_vals <- suppressWarnings(as.numeric(sc_si[[param]])); good <- is.finite(param_vals)
  if (!any(good)) next
  rng <- range(param_vals[good])
  size_val <- scales::rescale(param_vals, to = size_range, from = rng)
  pretty_vals <- pretty(param_vals[good], n = 4); pretty_vals <- pretty_vals[pretty_vals >= rng[1] & pretty_vals <= rng[2]]
  size_breaks <- scales::rescale(pretty_vals, to = size_range, from = rng)
  
  param_label <- gsub("\\s*\\[.*?\\]", "", param)
  df <- mutate(sc_si, size_val = size_val)
  
  param_plot <- ggplot(df, aes(dbRDA1, dbRDA2)) +
    geom_point(aes(color = factor(Cluster), size = size_val), alpha = 0.5) +
    geom_segment(data = sc_bp,
                 aes(x = 0, y = 0, xend = dbRDA1 * arrow_scaling_factor, yend = dbRDA2 * arrow_scaling_factor),
                 arrow = arrow(length = unit(0.5, "cm")), linewidth = 1.2, color = "black") +
    geom_text_repel(data = sc_bp,
                    aes(x = dbRDA1 * arrow_scaling_factor, y = dbRDA2 * arrow_scaling_factor, label = Variable),
                    color = "black", size = 10, fontface = "bold",
                    max.overlaps = Inf, box.padding = 1, point.padding = 0.5, segment.color = "grey50", seed = 161193) +
    scale_color_manual("Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)]) +
    scale_size_identity(param_label, guide = "legend", breaks = size_breaks, labels = pretty_vals) +
    guides(color = guide_legend(order = 1, override.aes = list(size = 15)),
           size  = guide_legend(order = 2, override.aes = list(color = "black"))) +
    labs(title = paste0("dbRDA plot for the 2022 and 2023 WK expeditions: ", param_label),
         x = paste0("dbRDA1 (", perc[1], "%)"), y = paste0("dbRDA2 (", perc[2], "%)")) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2), expand = FALSE) +
    theme_minimal() +
    theme(plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
          axis.text = element_text(size = 16), axis.title.x = element_text(size = 35, face = "bold"),
          axis.title.y = element_text(size = 35, face = "bold"),
          plot.margin = margin(t = 40, r = 60, b = 40, l = 60),
          legend.text = element_text(size = 30), legend.title = element_text(size = 35, face = "bold"),
          legend.key.size = unit(2, "cm"), legend.spacing.y = unit(3, "cm"))
  
  print(param_plot)
}

grid.newpage(); grid.draw(tableGrob(anova_per_var, rows = NULL,
                                    theme = ttheme_minimal(base_size = 60, padding = unit(c(10, 10), "mm"),
                                                           core = list(bg_params = list(fill = c("white", "grey95"), col = "black")),
                                                           colhead = list(bg_params = list(fill = "grey80", col = "black")))))
dev.off()

## dbRDA plots for the manuscript
# Shared ranges: dbRDA1 and dbRDA3 are horizontal; dbRDA2 is vertical in both panels
rbind(x = range(c(sc_si$dbRDA1, sc_si$dbRDA3)), y = range(sc_si$dbRDA2))

# Choose shared limits after inspecting the ranges; leave room for arrows and labels
db_xlim <- c(-2, 2)
db_ylim <- c(-1, 2.5)

p1 <- ggplot(sc_si, aes(dbRDA1, dbRDA2)) +
  geom_point(aes(fill = factor(Cluster), shape = factor(Year)), size = 7, alpha = 0.9, color = "black", stroke = 0.5) +
  geom_segment(data = sc_bp, aes(x = 0, y = 0, xend = dbRDA1*arrow_scaling_factor, yend = dbRDA2*arrow_scaling_factor),
               arrow = arrow(length = unit(1, "cm")), linewidth = 1, color = alpha("black", 0.9)) +
  geom_text_repel(data = sc_bp %>% mutate(Variable = gsub("`", "", Variable)),
                  aes(x = dbRDA1*arrow_scaling_factor, y = dbRDA2*arrow_scaling_factor, label = Variable),
                  size = 8, fontface = "bold", force = 15, box.padding = 5, seed = 161193
                  # segment.color = NA                                           # employ for svg b/c Inkscape doesn't see it
  ) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0.5),
                                          keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 1)) +
  scale_fill_manual(name = "Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)],
                    guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 7, alpha = 0.9, stroke = 0.5),
                                         keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 2)) +
  labs(x = sprintf("dbRDA1 (%.1f%%)", perc[1]), y = sprintf("dbRDA2 (%.1f%%)", perc[2])) +
  coord_cartesian(xlim = db_xlim, ylim = db_ylim, expand = FALSE) +
  theme_minimal() +
  theme(axis.title = element_text(size = 25, face = "bold"),
        axis.text  = element_text(size = 22, color = "black"),
        legend.title = element_text(size = 25, face = "bold"),
        legend.text  = element_text(size = 25, face = "bold"),
        legend.spacing.y  = unit(0.3, "cm"),
        plot.margin = margin(20, 30, 20, 30),
        legend.position = "none")

p2 <- ggplot(sc_si, aes(dbRDA3, dbRDA2)) +
  geom_point(aes(fill = factor(Cluster), shape = factor(Year)), size = 7, alpha = 0.9, color = "black", stroke = 0.5) +
  geom_segment(data = sc_bp, aes(x = 0, y = 0, xend = dbRDA3*arrow_scaling_factor, yend = dbRDA2*arrow_scaling_factor),
               arrow = arrow(length = unit(1, "cm")), linewidth = 1, color = alpha("black", 0.9)) +
  geom_text_repel(data = sc_bp %>% mutate(Variable = gsub("`", "", Variable)),
                  aes(x = dbRDA3*arrow_scaling_factor, y = dbRDA2*arrow_scaling_factor, label = Variable),
                  size = 8, fontface = "bold", force = 15, box.padding = 5, seed = 161193
                  # segment.color = NA                                           # employ for svg b/c Inkscape doesn't see it
  ) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0.5),
                                          keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 1)) +
  scale_fill_manual(name = "Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)],
                    guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 7, alpha = 0.9, stroke = 0.5),
                                         keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 2)) +
  labs(x = sprintf("dbRDA3 (%.1f%%)", perc[3]), y = NULL) +
  coord_cartesian(xlim = db_xlim, ylim = db_ylim, expand = FALSE) +
  theme_minimal() +
  theme(axis.title = element_text(size = 25, face = "bold"),
        axis.text  = element_text(size = 22, color = "black"),
        legend.title = element_text(size = 25, face = "bold"),
        legend.text  = element_text(size = 25, face = "bold"),
        legend.spacing.y  = unit(0.3, "cm"),
        plot.margin = margin(20, 30, 20, 30),
        axis.title.y = element_blank())

cairo_pdf("OUT/7W_ENV/wk_dbRDA.pdf", width = 24, height = 13.2)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

png("OUT/7W_ENV/wk_dbRDA.png", width = 24, height = 13.2, units = "in", res = 300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off() 

# # Optional: save side-by-side as .svg for manual editing in Inkscape
# svg("OUT/7W_ENV/wk_dbRDA.svg", width = 24, height = 13.2)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
# dev.off()


# Exploration: interactive 3D dbRDA
p3d <- plot_ly() %>% add_markers(data = sc_si, x = ~dbRDA1, y = ~dbRDA2, z = ~dbRDA3, color = ~Cluster, colors = cluster_palette, # sites
                                 symbol = ~Year, symbols = c(`2022` = "circle", `2023` = "diamond"), mode = "markers",
                                 marker = list(size = 6, line = list(color = "black", width = 1)), hoverinfo = "text",
                                 text = ~paste0(Station, "<br>Cluster ", cl_labels[as.integer(as.character(Cluster))], "<br>Year ", Year)) %>%
  add_trace(data = {                                                                                            # arrows + labels, built inline
    b <- sc_bp %>% left_join(anova_per_var, by = "Variable"); m <- nrow(b)
    df <- data.frame(Variable = rep(b$Variable, each = 2), F = rep(b$F, each = 2),
                     x = c(rbind(rep(0, m), b$dbRDA1 * arrow_scaling_factor)),
                     y = c(rbind(rep(0, m), b$dbRDA2 * arrow_scaling_factor)),
                     z = c(rbind(rep(0, m), b$dbRDA3 * arrow_scaling_factor)))
    df$label <- ""; df$label[seq(2, 2 * m, by = 2)] <- b$Variable; df}, x = ~x, y = ~y, z = ~z, type = "scatter3d",
    mode = "lines+text", split = ~Variable, text = ~label, textposition = "top center", line = list(color = "black", width = 5),
    textfont = list(color = "black", size = 16, family = "Arial Black"), showlegend = FALSE, hoverinfo = "text") %>%
  layout(scene = list(xaxis = list(title = list(text = paste0("dbRDA1 (", perc[1], "%)"), font = list(color="black")),
                                   tickfont = list(color="black"), range = c(-3,2)),
                      yaxis = list(title = list(text = paste0("dbRDA2 (", perc[2], "%)"), font = list(color="black")),
                                   tickfont = list(color="black"), range = c(-3,2)),
                      zaxis = list(title = list(text = paste0("dbRDA3 (", perc[3], "%)"), font = list(color="black")),
                                   tickfont = list(color="black"), range = c(-3,2)),
                      aspectmode = "data"),
         legend = list(orientation = "h"), showlegend = FALSE)
htmlwidgets::saveWidget(p3d, "OUT/7W_ENV/wk_dbRDA_3d.html", selfcontained = TRUE)  # self-contained HTML; requires pandoc (included with RStudio)
unlink("OUT/7W_ENV/wk_dbRDA_3d_files", recursive = TRUE)                              # remove the temporary library folder left by saveWidget

# #PCA inspection
# # if it generally agrees w/ dbRDA, then your env. variables are doing a good job explaining the Bray-Curtis dissimilarity
# 
# t1 <- prcomp(microbial_data) 
# t2 <- scores(t1, display = "sites") %>% as.data.frame() %>% cbind(metadata) %>% mutate(Cluster = as.factor(Cluster))
# 
# ggplot(t2, aes(PC1, PC3)) +
#   geom_point(aes(fill = factor(Cluster), shape = factor(Year)), size = 8.5, alpha = 0.9, color = "black", stroke = 0.5) +
#   scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
#                      guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 8.5, stroke = 0.5),
#                                           keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 1)) +
#   scale_fill_manual(name = "Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)],
#                     guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 8.5, alpha = 0.9, stroke = 0.5),
#                                          keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 2))


################
### TSD plot ###
################

p_ts <- ggplot(meta_grouped, aes(`Salinity`, `Temperature [deg C]`)) +
  scale_x_continuous(limits = c(0, 31), trans = transform_exp(1.15),  # exp transformation
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

grDevices::cairo_pdf("OUT/7W_ENV/wk_ts.pdf", width = 24, height = 11); print(p_ts); dev.off()

png("OUT/7W_ENV/wk_ts.png", width = 24, height = 12, units = "in", res = 300)
print(p_ts)
dev.off()


#################
### Fluo plot ###
#################

meta_grouped_FC <- meta_grouped %>% filter(!is.na(`Beam attenuation [1/m]`))
p_fc <- ggplot(meta_grouped_FC, aes(`CDOM fluorescence [mg/m^3]`, `Chlorophyll-a fluorescence [mg/m^3]`)) +
  scale_x_continuous(limits = c(1, 100), trans = transform_log(),  # log transformation
                     breaks = c(1, 2, 5, 10, 25, 50, 100), minor_breaks = NULL) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 11), minor_breaks = NULL) +
  geom_point(aes(color = factor(Cluster), shape = factor(Year)), fill = NA, size = 7, stroke = 2.5) +
  geom_point(aes(color = factor(Cluster), shape = factor(Year), alpha = `Beam attenuation [1/m]`, fill = after_scale(colour)), size = 7, stroke = 0) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0))) +
  scale_color_manual(values = cluster_palette, labels = \(x) cl_labels[as.integer(x)], name = "Cluster") +
  scale_alpha_continuous(
    name = expression(atop(textstyle(bold("Beam")), atop(textstyle(bold("attenuation")), textstyle(bold("(m"^{-1}*")"))))), 
    range = c(0.1, 1), limits = c(0.1, 6.5), breaks = c(0.2, 0.5, 1, 3, 6.5)) + 
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2, override.aes = list(size = 7, shape = 21, stroke = 2.4), keyheight = unit(1.1, "cm"),
                              keywidth = unit(1.1, "cm"), default.unit = "cm", byrow = TRUE), 
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

grDevices::cairo_pdf("OUT/7W_ENV/wk_fc.pdf", width = 24, height = 11); print(p_fc); dev.off()

png("OUT/7W_ENV/wk_fc.png", width = 24, height = 11, units = "in", res = 300, type = "cairo")
print(p_fc)
dev.off()

# median chl a in surface waters
# median(meta.wk %>% filter(`Depth [m]` <= 6) %>% .$`Chlorophyll-a fluorescence [mg/m^3]`, na.rm = TRUE)
# meta.wk %>% filter(`Depth [m]` <= 6, !is.na(`Chlorophyll-a fluorescence [mg/m^3]`)) %>% nrow()


############################################################
### Environmental properties characterizing each cluster ###
############################################################

# Metadata columns and publication-ready headings
env_headers <- c("Depth\n(m)" = "Depth [m]", "Distance to coast\n(km)" = "Distance to coast [km]",
                 "Temperature\n(°C)" = "Temperature [deg C]", "Salinity" = "Salinity",
                 "fChl-a\n(mg m-3)" = "Chlorophyll-a fluorescence [mg/m^3]", "fCDOM\n(mg m-3)" = "CDOM fluorescence [mg/m^3]",
                 "Beam attenuation\n(m-1)" = "Beam attenuation [1/m]")
env_vars <- unname(env_headers)
env_labels <- setNames(str_replace_all(names(env_headers), "\n", " "), env_vars)
year_col <- if ("Year" %in% names(meta.wk)) "Year" else "Cruise"
cluster_full <- c(oHB = "oHB (offshore Hudson Bay)", BI22 = "BI22 (the Belcher Islands)", cHB = "cHB (coastal Hudson Bay)",
                  dBI = "dBI (the Belcher Islands, subsurface)", dJB = "dJB (James Bay, subsurface)", EST = "EST (the Churchill and Winisk estuaries)",
                  MR = "MR (the Maquatua River)", nJB = "nJB (northern James Bay)", sJB = "sJB (southern James Bay)")

# Report median [range] when at least half the cluster has measurements; add * when some values are missing
med_range_txt <- \(x) {
  n <- length(x); x <- x[!is.na(x)]
  if (length(x)/n < 0.5) return(NA_character_)
  star <- if (length(x) < n) "*" else ""
  if (min(x) == max(x)) sprintf("%.1f%s", median(x), star) else sprintf("%.1f%s\n[%.1f-%.1f]", median(x), star, min(x), max(x))
}

year_txt <- \(x) {
  years <- sort(unique(unlist(str_extract_all(as.character(x), "20\\d{2}"))))
  if (length(years)) paste(years, collapse = "-") else NA_character_}

scientific_txt <- \(p) {
  if (is.na(p)) return(NA_character_)
  x <- str_split_fixed(formatC(p, format = "e", digits = 1), "e", 2)
  paste0(x[1], " × 10", as.integer(x[2]))}

# Check ANOVA residual normality before using non-parametric tests
environment_fits <- set_names(lapply(env_vars, \(var) {
  d <- meta.wk %>% transmute(Cluster, Value = .data[[var]]) %>% drop_na()
  aov(Value ~ factor(Cluster), data = d)}), env_vars)

imap_dfr(environment_fits, \(m, var) tibble(Parameter = env_labels[[var]], `Shapiro p-value` = shapiro.test(residuals(m))$p.value)) %>% print()
walk2(environment_fits, names(environment_fits), ~{qqnorm(residuals(.x), main = env_labels[[.y]]); qqline(residuals(.x), col = "red")})

# Overall differences in environmental properties among clusters
environment_kruskal <- map_dfr(env_vars, \(var) {
  d <- meta.wk %>% transmute(Cluster, Value = .data[[var]]) %>% drop_na()
  tibble(Parameter = var, `Kruskal-Wallis p-value` = if (n_distinct(d$Cluster) < 2) NA_real_ else kruskal.test(Value ~ Cluster, data = d)$p.value)})
kw_p <- setNames(environment_kruskal$`Kruskal-Wallis p-value`, environment_kruskal$Parameter)

# Significant Dunn pairwise differences among clusters
environment_dunn <- map_dfr(env_vars, \(var) {
  if (is.na(kw_p[[var]]) || kw_p[[var]] >= 0.05) return(tibble())
  d <- meta.wk %>% transmute(Cluster = as.integer(as.character(Cluster)), Value = .data[[var]]) %>% drop_na()
  dt <- dunn.test(d$Value, d$Cluster, method = "bh", kw = FALSE, altp = TRUE, table = FALSE, list = FALSE, interpret = FALSE)
  med <- tapply(d$Value, d$Cluster, median)
  med <- setNames(as.numeric(med), names(med))  # remove array dimension
  groups <- str_split_fixed(dt$comparisons, "\\s+-\\s+", 2)
  g1 <- groups[, 1]
  g2 <- groups[, 2]
  group1 <- unname(cl_labels[as.numeric(g1)])
  group2 <- unname(cl_labels[as.numeric(g2)])
  tibble(Parameter = env_labels[[var]],
         `Direction of difference` = case_when(med[g1] > med[g2] ~ paste(group1, ">", group2), med[g2] > med[g1] ~ paste(group2, ">", group1),
                                               TRUE ~ paste(group1, "=", group2)), `BH-adjusted p-value` = dt$altP.adjusted)}) %>%
  filter(`BH-adjusted p-value` < 0.05) %>% arrange(Parameter, `BH-adjusted p-value`) %>%
  mutate(`BH-adjusted p-value` = map_chr(`BH-adjusted p-value`, scientific_txt))

# Publication-ready median [range] table
environment_summary <- meta.wk %>% mutate(Cluster = as.character(Cluster)) %>% group_by(Cluster) %>%
  summarise(Year = year_txt(.data[[year_col]]), across(all_of(env_vars), med_range_txt), .groups = "drop") %>%
  arrange(as.numeric(Cluster)) %>%
  mutate(Cluster = unname(cluster_full[cl_labels[as.numeric(Cluster)]])) %>%
  bind_rows(tibble(Cluster = "Kruskal-Wallis p", Year = NA_character_) %>%
              bind_cols(as_tibble_row(set_names(map_chr(env_vars, \(var) scientific_txt(kw_p[[var]])), env_vars)))) %>%
  rename(all_of(env_headers)) %>% mutate(across(everything(), \(x) replace(x, is.na(x), "NA")))

# Export summary and significant Dunn comparisons
wb <- createWorkbook()
addWorksheet(wb, "Environment by cluster")
writeData(wb, 1, "Environmental properties by cluster", startRow = 1)
writeDataTable(wb, 1, environment_summary, startRow = 2, tableName = "EnvironmentalSummary")
summary_end <- nrow(environment_summary) + 2
note_row <- summary_end + 2
dunn_row <- note_row + 3
writeData(wb, 1, "Values are medians [ranges]. *Some samples in the cluster (no more than half) had no measurements.", startRow = note_row)
mergeCells(wb, 1, cols = 1:ncol(environment_summary), rows = note_row)
writeData(wb, 1, "Significant Dunn's pairwise tests", startRow = dunn_row)
writeDataTable(wb, 1, environment_dunn, startRow = dunn_row + 1, tableName = "EnvironmentalDunn")
wrap <- createStyle(wrapText = TRUE, valign = "top")
addStyle(wb, 1, wrap, rows = 2:summary_end, cols = 1:ncol(environment_summary), gridExpand = TRUE, stack = TRUE)
addStyle(wb, 1, wrap, rows = note_row, cols = 1, stack = TRUE)
setRowHeights(wb, 1, rows = 2, heights = 32)
setRowHeights(wb, 1, rows = 3:(summary_end - 1), heights = 30)
setColWidths(wb, 1, cols = 1:9, widths = c(42, 14, 16, 22, 18, 14, 18, 18, 24))
saveWorkbook(wb, "OUT/7W_ENV/Environmental_statistics.xlsx", overwrite = TRUE)


##########################################
### Prepare MIMARKS BioSample metadata ###
##########################################

bioproject_accession <- "PRJNA421293"
times_are_utc <- TRUE
filelist_path <- "DNA/WK_filelist.txt"

fmt_num <- function(x, digits = 3) if_else(is.na(x), NA_character_, sub("\\.?0+$", "", sprintf(paste0("%.", digits, "f"), x)))

wk <- meta.wk %>%
  mutate(depth_id = fmt_num(`Depth [m]`),
         coastal = case_when(is.na(`Distance to coast [km]`) ~ NA, `Distance to coast [km]` <= 10 ~ TRUE, TRUE ~ FALSE),
         # River samples (salinity < 5, i.e. fresh to oligohaline in the Venice system; the Maquatua River samples, MR) get freshwater terms
         riverine = !is.na(Salinity) & Salinity < 5,
         sample_name = if_else(is.na(Cruise) | is.na(Station) | is.na(depth_id), NA_character_,
                               str_replace_all(paste(Cruise, Station, paste0(depth_id, "m"), sep = "_"), "[^A-Za-z0-9._-]", "-")),
         sample_title = if_else(is.na(sample_name), NA_character_,
                                paste0("Hudson Bay microbial communities: ",
                                       Cruise, ", station ", Station, ", ", depth_id, " m")),
         collection_date = case_when(
           times_are_utc & if_all(c(Year, Month, Day, Hour, Minute, Second), ~ !is.na(.x)) ~
             sprintf("%04d-%02d-%02dT%02d:%02d:%02dZ", Year, Month, Day, Hour, Minute, Second),
           !times_are_utc & if_all(c(Year, Month, Day), ~ !is.na(.x)) ~ sprintf("%04d-%02d-%02d", Year, Month, Day),
           TRUE ~ NA_character_))

mimarks_wk <- wk %>%
  transmute(`*sample_name` = sample_name, sample_title, bioproject_accession = .env$bioproject_accession,
            `*organism` = if_else(riverine, "freshwater metagenome", "marine metagenome"), `*collection_date` = collection_date, `*depth` = `Depth [m]`,
            `*env_broad_scale` = if_else(riverine, "freshwater river biome [ENVO:01000253]", "marine biome [ENVO:00000447]"),
            `*env_local_scale` = case_when(riverine ~ "river [ENVO:00000022]", coastal ~ "coastal water body [ENVO:02000049]",
                                           TRUE ~ "marine pelagic zone [ENVO:00000208]"),
            `*env_medium` = case_when(riverine ~ "river water [ENVO:01000599]", coastal ~ "coastal sea water [ENVO:00002150]",
                                      TRUE ~ "sea water [ENVO:00002149]"),
            `*geo_loc_name` = geo_loc_name,
            `*lat_lon` = if_else(is.na(`Latitude [deg]`) | is.na(`Longitude [deg]`), NA_character_,
                                 sprintf("%.5f %s %.5f %s", abs(`Latitude [deg]`), if_else(`Latitude [deg]` >= 0, "N", "S"),
                                         abs(`Longitude [deg]`), if_else(`Longitude [deg]` >= 0, "E", "W"))),
            tot_depth_water_col = `Bottom depth [m]`, size_frac_up = 0.22, samp_store_temp = -80,
            temp = `Temperature [deg C]`, salinity = Salinity, fluor = `Chlorophyll-a fluorescence [mg/m^3]`,
            cdom = `CDOM fluorescence [mg/m^3]`, beam_attenuation = `Beam attenuation [1/m]`,
            mixed_layer_depth = `Mixed layer depth [m]`, distance_to_coast = `Distance to coast [km]`,
            dna_conc = `DNA [mg/L]`,
            `parameter (units)` = paste(
              if (times_are_utc) "collection_date (UTC);" else "collection_date;",
              "depth and tot_depth_water_col (m); size_frac_up (um); samp_store_temp and temp (degrees C); salinity;",
              "fluor and cdom (mg m-3); beam_attenuation (m-1); mixed_layer_depth (m);",
              "distance_to_coast (km); dna_conc (mg L-1)"))

required <- c("*sample_name", "sample_title", "bioproject_accession", "*organism", "*collection_date", "*depth",
              "*env_broad_scale", "*env_local_scale", "*env_medium", "*geo_loc_name", "*lat_lon")
bad <- vapply(mimarks_wk[required], \(x) any(is.na(x) | trimws(as.character(x)) == ""), logical(1))
if (any(bad)) stop("Missing required MIMARKS values: ", paste(names(bad)[bad], collapse = ", "))
if (anyDuplicated(mimarks_wk$`*sample_name`)) stop("Duplicate MIMARKS sample_name values.")
write.csv(mimarks_wk, "OUT/7W_ENV/NCBI_submission/WK_MIMARKS.survey.water.6.0.csv", row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

build_lookup <- function(x, marker) {
  p <- str_split_fixed(x, "\\|", 4)
  tibble(plate.geno = toupper(paste(str_split_fixed(p[, 3], "-", 6)[, 2], p[, 4], sep = "_")),
         library_ID = paste0(p[, 3], "_", p[, 4]), marker)
}

file_map <- tibble(file = basename(str_trim(readLines(filelist_path)))) %>%
  filter(str_detect(file, "_R[12]\\.fastq\\.gz$")) %>%
  transmute(library_ID = str_remove(file, "_R[12]\\.fastq\\.gz$"),
            read = str_extract(file, "R[12](?=\\.fastq\\.gz$)"), file) %>% distinct()
if (anyDuplicated(file_map[c("library_ID", "read")])) stop("Duplicate library/read entries in ", filelist_path)
file_map <- file_map %>% pivot_wider(names_from = read, values_from = file)

sra_wk <- bind_rows(build_lookup(rownames(seqtab.16S), "16S"), build_lookup(rownames(seqtab.18S), "18S")) %>%
  inner_join(wk %>% select(plate.geno, sample_name, sample_title), by = "plate.geno") %>%
  left_join(file_map, by = "library_ID") %>%
  arrange(match(plate.geno, wk$plate.geno), match(marker, c("16S", "18S"))) %>%
  transmute(sample_name, library_ID, title = paste0(marker, " rRNA gene metabarcoding of marine microbial communities: ", sub("^Hudson Bay microbial communities: ", "", sample_title)),  # as submitted to the SRA
            library_strategy = "AMPLICON", library_source = "METAGENOMIC", library_selection = "PCR",
            library_layout = "paired", platform = "ILLUMINA", instrument_model = "Illumina NovaSeq 6000",
            design_description = paste0("as per CD Genomics; Earth Microbiome Project ", marker, " rRNA gene primers"),
            filetype = "fastq", filename = R1, filename2 = R2)

if (nrow(sra_wk) != 2 * nrow(wk) || anyNA(sra_wk[c("sample_name", "library_ID", "filename", "filename2")]) ||
    anyDuplicated(sra_wk$library_ID)) stop("SRA matching failed: each sample must have one unique 16S and 18S library with R1/R2 files.")
write.csv(sra_wk, "OUT/7W_ENV/NCBI_submission/WK_SRA_metadata.csv",  row.names = FALSE, na = "NA", fileEncoding = "UTF-8")


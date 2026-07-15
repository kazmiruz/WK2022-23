# WK_7---------------------------------------------------------------------

# Purpose:
#   Perform distance‑based redundancy analysis (db‑RDA) on the fourth‑root transformed Bray‑Curtis dissimilarities
#   to resolve the influence of depth, distance to coast, temperature, and salinity on plankton community structure.
#   Produce a 2‑panel ordination plot (dbRDA1‑2 and dbRDA3‑2), a 3‑D interactive plot, an ANOVA table,
#   and a multi‑page diagnostic PDF with sample points sized by various environmental parameters.
#
# Prerequisites (objects must be in the workspace, created by earlier scripts):
#   proptab.trans, cluster_palette, cl_labels, meta.wk.small (WK_3_CLUS_TITAN_FT)
#
# Outputs:
#   OUTPUT/wk_db-RDA.pdf / .png         (2‑panel ordination)
#   OUTPUT/wk_dbRDA_3d.html             (interactive 3‑D plot)
#   OUTPUT/wk_db-RDA_many.pdf           (multi‑page diagnostic, each env variable)
#   Console: ANOVA table, adjusted R², axis percentages
# -------------------------------------------------------------------------

library(tidyverse)
library(vegan)    # for RDA, standardization: decostand(), rankindex(), dbrda(), vif.cca(), anova.cca(), RsquareAdj(), scores()
library(GGally)   # for ggpairs
library(ggrepel)  # to place labels nicely
library(scales)   # for rescale
library(grid)
library(gridExtra)
library(plotly)   # for 3D
library(htmlwidgets)

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
metadata <- meta.wk.small        # explanatory variables
microbial_data <- proptab.trans  # response variables (already normalized + transformed)

# Select numeric environmental variables (keep spaces; drop units in names)
env_vars <- metadata %>% select("Distance to coast [km]", "Depth [m]", "Salinity [PSU]", "Temperature [deg C]") %>%
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

# Distance index ranking (optional)
rankindex(env_stand, microbial_data, indices = c("euc","man","gow","bra","kul"), stepacross = FALSE, method = "pearson")

# Model
rda_model <- dbrda(microbial_data ~ ., data = env_stand, distance = "bray", dfun = vegdist)

# VIFs
vif.cca(rda_model) %>% as.data.frame() %>% rownames_to_column("Variable") %>% rename(VIF = ".") %>% mutate(VIF = round(VIF, 0)) %>% arrange(desc(VIF))

# ANOVAs (model, axes, vars)
# anova.cca(rda_model, permutations = 99999)                                      # test model sign-ce

# anova.cca(rda_model, permutations = 99999, by = "axis")                         # test axes sign-ce
# Model: dbrda(formula = microbial_data ~ `Distance to coast` + Depth + Salinity + Temperature, data = env_stand, distance = "bray", dfun = vegdist)
# Df SumOfSqs       F  Pr(>F)    
# dbRDA1    1   2.7425 17.7794   1e-05 ***
#   dbRDA2    1   1.2604  8.1711   1e-05 ***
#   dbRDA3    1   0.9215  5.9742   2e-05 ***
#   dbRDA4    1   0.3906  2.5321 0.00951 ** 
#   Residual 67  10.3347

# anova.cca(rda_model, permutations = 99,999, by = "term") %>% arrange(desc(F))  # test var. sign-ce (biased towards the vars appearing first in the formula)
anova_per_var <- anova.cca(rda_model, permutations = 99999, by = "margin") %>%  # test var. sign-ce if marginal (last in the formula) and store it neatly
  as.data.frame() %>% rownames_to_column("Variable") %>% filter(Variable != "Residual") %>% 
  transmute(Variable = gsub("`", "", Variable), F = round(F, 1), P = `Pr(>F)`,
            Sig = as.character(symnum(P, cutpoints = c(0, 0.00009, 0.0009, 0.0049, 0.009, 1), symbols = c("***", "**", "*", ".", " ")))) %>% arrange(desc(F))

# Variance info for titles
perc <- round(100 * (summary(rda_model)$cont$importance[2, 1:4]), 1) %>% print()  # Bray-Curtis dis-ty explained by your env var-s: but only only dbRA1 and dbRDA2
RsquareAdj(rda_model)$adj.r.squared                                               # adjusted R2

# Scores and arrow scaling (remove any backticks in bp labels)
sc_bp <- scores(rda_model, display = "bp", scaling = 2, choices = 1:rda_model$CCA$rank) %>%
  as.data.frame() %>% rownames_to_column("Variable") %>% mutate(Variable = gsub("`", "", Variable))
sc_si <- scores(rda_model, display = "sites", choices = 1:rda_model$CCA$rank) %>%
  as.data.frame() %>% cbind(metadata) %>% mutate(Cluster = as.factor(Cluster))
arrow_scaling_factor <- max(abs(as.matrix(sc_bp[grep("^dbRDA", names(sc_bp))])), na.rm = TRUE) * 1.7

# Multipage plot: one page per numeric variable sized onto points
pdf(file = "OUTPUT/wk_db-RDA_many.pdf", width = 35, height = 20)
numeric_vars <- meta.wk.small %>% select(Year, "Bot. Depth [m]", "Depth [m]", "Distance to coast [km]", "Temperature [deg C]", "Salinity [PSU]",
                                         "PAR [umol photons/m^2/sec]", "Chlorophyll-a fluorescence [mg/m^3]", "CDOM fluorescence [mg/m^3]", "Oxygen [umol/kg]", dna_conc) %>% colnames()
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
                    max.overlaps = Inf, box.padding = 1, point.padding = 0.5, segment.color = "grey50") +
    scale_color_manual("Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)]) +
    scale_size_identity(param_label, guide = "legend", breaks = size_breaks, labels = pretty_vals) +
    guides(color = guide_legend(order = 1, override.aes = list(size = 15)),
           size  = guide_legend(order = 2, override.aes = list(color = "black"))) +
    labs(title = paste0("db-RDA plot for the 2022 and 2023 WK expeditions: ", param_label),
         x = paste0("dbRDA1 (", perc[1], "%)"), y = paste0("dbRDA2 (", perc[2], "%)")) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2), expand = FALSE) +
    theme_minimal() +
    theme(plot.title = element_text(size = 50, face = "bold", hjust = 0.5, margin = margin(b = 20)),
          axis.text = element_text(size = 16), axis.title.x = element_text(size = 35, face = "bold"), axis.title.y = element_text(size = 35, face = "bold"),
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

################################################################################

# --- db-RDA plot ---

p1 <- ggplot(sc_si, aes(dbRDA1, dbRDA2)) +
  geom_point(aes(fill = factor(Cluster), shape = factor(Year)), size = 7, alpha = 0.9, color = "black", stroke = 0.5) +
  geom_segment(data = sc_bp, aes(x = 0, y = 0, xend = dbRDA1*arrow_scaling_factor, yend = dbRDA2*arrow_scaling_factor),
               arrow = arrow(length = unit(1, "cm")), linewidth = 1, color = alpha("black", 0.9)) +
  geom_text_repel(data = sc_bp %>% mutate(Variable = gsub("`", "", Variable)),
                  aes(x = dbRDA1*arrow_scaling_factor, y = dbRDA2*arrow_scaling_factor, label = Variable),
                  size = 8, fontface = "bold", force = 15, box.padding = 5,
                  # segment.color = NA # employ for svg b/c Inkscape doesn't see it
                  ) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0.5),
                                          keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 1)) +
  scale_fill_manual(name = "Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)],
                    guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 7, alpha = 0.9, stroke = 0.5),
                                         keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 2)) +
  labs(x = sprintf("dbRDA1 (%.1f%%)", perc[1]), y = sprintf("dbRDA2 (%.1f%%)", perc[2])) +
  coord_cartesian(xlim = c(-2, 1.5), ylim = c(-3, 1), expand = FALSE) +
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
                  size = 8, fontface = "bold", force = 15, box.padding = 5,
                  # segment.color = NA # employ for svg b/c Inkscape doesn't see it
                  ) +
  scale_shape_manual(name = "Year", values = c(`2022` = 21, `2023` = 23),
                     guide = guide_legend(override.aes = list(fill = alpha("black", 0.9), color = "black", size = 7, stroke = 0.5),
                                          keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 1)) +
  scale_fill_manual(name = "Cluster", values = cluster_palette, labels = \(x) cl_labels[as.integer(x)],
                    guide = guide_legend(override.aes = list(shape = 21, color = "black", size = 7, alpha = 0.9, stroke = 0.5),
                                         keywidth = unit(1.1, "cm"), keyheight = unit(1.1, "cm"), order = 2)) +
  labs(x = sprintf("dbRDA3 (%.1f%%)", perc[3]), y = NULL) +
  coord_cartesian(xlim = c(-2, 1.5), ylim = c(-3, 1), expand = FALSE) +
  theme_minimal() +
  theme(axis.title = element_text(size = 25, face = "bold"),
        axis.text  = element_text(size = 22, color = "black"),
        legend.title = element_text(size = 25, face = "bold"),
        legend.text  = element_text(size = 25, face = "bold"),
        legend.spacing.y  = unit(0.3, "cm"),
        plot.margin = margin(20, 30, 20, 30),
        axis.title.y = element_blank())

cairo_pdf("OUTPUT/wk_db-RDA.pdf", width = 24, height = 13.2)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

png("OUTPUT/wk_db-RDA.png", width = 24, height = 13.2, units = "in", res = 300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off() 

# # save side-by-side (for Inkscape)
# svg("OUTPUT/wk_db-RDA.svg", width = 30, height = 15)
# gridExtra::grid.arrange(p1, p2, ncol = 2)
# dev.off()

# dbRDA axis ranges
# min(sc_si$dbRDA1); min(sc_si$dbRDA2); min(sc_si$dbRDA3)
# max(sc_si$dbRDA1); max(sc_si$dbRDA2); max(sc_si$dbRDA3)


# 3D dbRDA
p3d <- plot_ly() %>% add_markers(data = sc_si, x = ~dbRDA1, y = ~dbRDA2, z = ~dbRDA3, color = ~Cluster, colors = cluster_palette, # sites
                                 symbol = ~Year, symbols = c(`2022` = "circle", `2023` = "diamond"), mode = "markers",
                                 marker = list(size = 6, line = list(color = "black", width = 1)), hoverinfo = "text",
                                 text = ~paste0(Station, "<br>Cluster ", cl_labels[as.integer(as.character(Cluster))], "<br>Year ", Year)) %>%
  add_trace(data = {                                                                                                              # arrows + labels, built inline
    b <- sc_bp %>% left_join(anova_per_var, by = "Variable"); m <- nrow(b)
    df <- data.frame(Variable = rep(b$Variable, each = 2), F = rep(b$F, each = 2),
                     x = c(rbind(rep(0, m), b$dbRDA1 * arrow_scaling_factor)),
                     y = c(rbind(rep(0, m), b$dbRDA2 * arrow_scaling_factor)),
                     z = c(rbind(rep(0, m), b$dbRDA3 * arrow_scaling_factor)))
    df$label <- ""; df$label[seq(2, 2 * m, by = 2)] <- b$Variable; df},
    x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "lines+text", split = ~Variable, text = ~label, textposition = "top center",
    line = list(color = "black", width = 5), textfont = list(color = "black", size = 16, family = "Arial Black"), showlegend = FALSE, hoverinfo = "text") %>%
  layout(scene = list(xaxis = list(title = list(text = paste0("dbRDA1 (", perc[1], "%)"), font = list(color="black")),
                                   tickfont = list(color="black"), range = c(-3,2)),
                      yaxis = list(title = list(text = paste0("dbRDA2 (", perc[2], "%)"), font = list(color="black")),
                                   tickfont = list(color="black"), range = c(-3,2)),
                      zaxis = list(title = list(text = paste0("dbRDA3 (", perc[3], "%)"), font = list(color="black")),
                                   tickfont = list(color="black"), range = c(-3,2)),
                      aspectmode = "data"),
         legend = list(orientation = "h"), showlegend = FALSE)
htmlwidgets::saveWidget(p3d, "OUTPUT/wk_dbRDA_3d.html", selfcontained = TRUE)


# ######################
# ### PCA inspection ###
# # if it generally agrees w/ db-RDA, then your env. variables are doing a good job explaining the Bray-Curtis dissimilarity
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

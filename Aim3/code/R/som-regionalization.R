###########################################
# Self Organizing Maps (SOM) implementation
# for management priority landscapes 
# "archetypes" of aspen management (?)

library(tidyverse)
library(sf)
library(kohonen)
library(reshape2)
library(scales)
library(ggspatial)
library(gridExtra)
library(ggradar)
library(cluster)
library(ggcorrplot)
library(patchwork)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"

# load and prep the firesheds data
firesheds <- st_read(paste0(projdir,"data/spatial/mod/srm_firesheds_model_data_wBivar.gpkg")) %>%
 mutate(bi_class = as.factor(bi_class))
glimpse(firesheds)

####################################
# select numeric variables for SOM #
X <- firesheds %>%
 # drop the geometry
 st_drop_geometry() %>%
 # keep only the subset we want to use for SOM
 select(
  trend_count, delta585, # future future / future aspen
  aspen_ha, patch_size, patch_den, big_patch, # contemporary aspen presence/patches
  combust, pop_den, wui, wui_dist, msbf_count, # built environment metrics
  burned_pct, disturbed, # proportion of fireshed that has experienced disturbance
  whp_p90, hui_p90, # wildfire risk to communities
  forest_pct, lf_canopy, lf_height, cfl_p90 # forested area and mean canopy cover
 ) %>%
 # Fill NAs with zero and scale variables
 mutate(
  across(everything(), ~ replace_na(.x, 0))
 ) %>%
 mutate(
  across(everything(), ~ as.numeric(scale(., center=T)))
 )
# Check structure
str(X)

#########################
# correlation structure #
# Compute correlation matrix
cor_mat <- cor(X, use = "complete.obs", method = "spearman")
# Plot correlation matrix
ggcorrplot(cor_mat, method = "circle",
 type = "lower", lab = TRUE, lab_size = 3, tl.cex = 10,
 colors = c("blue", "white", "red")
)
rm(cor_mat) # tidy up



#==============VARIABLE WEIGHTING================#

##########################
# WUI codes for reference:
# 1: 'Forest/Shrubland/Wetland-dominated Intermix WU'
# 2: 'Forest/Shrubland/Wetland-dominated Interface WUI'
# 3: 'Grassland-dominated Intermix WUI'
# 4: 'Grassland -dominated Interface WUI'

# # weight by aspen canopy and patch size
# # Get un-scaled variable to use as a weight
# aspen_wt <- scale(firesheds$aspen10_ha, center = FALSE)  # avoid zero centering
# # Apply to entire matrix (row-wise)
# X_wt <- X * aspen_wt

# assign "domains" to variables
# use "domains" to assign weighting schemes
domains <- list(
 aspen = X[, c("delta585", "aspen_ha", "patch_size", "patch_den", "big_patch")],
 fire = X[, c("trend_count", "burned_pct", "disturbed")],
 socio = X[, c("combust", "pop_den", "wui", "wui_dist", "msbf_count")],
 hazard = X[, c("hui_p90", "whp_p90")],
 forest = X[, c("forest_pct", "lf_canopy", "lf_height", "cfl_p90")]
)

# Define domains and domain-level weights
domain_weights <- c(
 aspen = 1.0,
 fire = 0.6,
 socio = 0.8,
 hazard = 1.0,
 forest = 0.4
)

# weight the input variables by domain
X_domains <- mapply(function(df, w) df * w, domains, domain_weights, SIMPLIFY = FALSE)
X_domains <- do.call(cbind, X_domains)


#==============SOM SETUP================#
set.seed(123)

# Define SOM grid size (adjust xdim & ydim as needed)
som.grid <- kohonen::somgrid(
 xdim = 12, 
 ydim = 12, 
 topo = "hexagonal"
)

# Train the SOM
som.model <- som(
 as.matrix(X_domains),
 grid = som.grid,
 rlen = 5001,   # Number of training iterations
 alpha = c(0.1, 0.02), # Learning rate (start, end)
 keep.data = TRUE
)

# View summary of the trained model
summary(som.model)
# Plot the nodes
plot(som.model)
plot(som.model, type = "changes")

# look at SOM node SD
# high = variable helps differentiate SOM nodes
apply(som.model$codes[[1]], 2, sd)



#==============MAP GRIDS================#

# define the cluster assignments
n_cluster = 8 # number of clusters to compute
# Compute distance matrix on SOM codebook vectors
dist_mat <- dist(som.model$codes[[1]])
# hierarchical clustering on SOM node weights
hc <- hclust(dist_mat, method = "ward.D2")
som.cluster <- cutree(hc, k = n_cluster)  

# Assign clusters to each fireshed
firesheds$som.cluster <- som.cluster[som.model$unit.classif]
# Check distribution of clusters
table(firesheds$som.cluster)

# Create a color mapping for the clusters
cluster_colors <- c(
 "#1b9e77",  # Teal Green
 "#d95f02",  # Burnt Orange
 "#7570b3",  # Slate Blue
 "#e7298a",  # Magenta
 "#66a61e",  # Olive Green
 "#e6ab02",  # Mustard
 "#a6761d",  # Earth Brown
 "#666666"   # Charcoal Gray
)
names(cluster_colors) <- as.character(1:8)

# make the spatial map
(cl.map <- ggplot(firesheds) +
 geom_sf(aes(fill = as.factor(som.cluster)), color = NA) +
 # scale_fill_brewer(palette = "Accent", name = "SOM Cluster") +
 scale_fill_manual(values = cluster_colors, name = "SOM Cluster") +
 theme_void() +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(
  location = "br", width_hint = 0.1,
  pad_x = unit(0.15,"in"), pad_y = unit(0.35,"in"),
  line_width = 0.5, text_pad = unit(0.15,"cm"),
  height = unit(0.15,"cm")
 ) +
 ggspatial::annotation_north_arrow(
  location = "br", which_north = "true",
  pad_x = unit(0.15,"in"), pad_y= unit(0.55,"in"),
  width = unit(0.8,"cm"), height = unit(0.8,"cm"),
  style = north_arrow_fancy_orienteering
 ) +
 theme(
  legend.title = element_text(angle = 90, size = 11, vjust = 1.2, hjust = 0.5),
  legend.text = element_text(size = 10),
  legend.key.size = unit(0.4, "cm"),  
  legend.position = c(0.24, 0.80)
 ) +
 guides(fill = guide_legend(title.position = "left", title.vjust = 1)))

# save the plot.
out_png <- paste0(projdir,'figures/SRM_Firesheds_SOMclusters.png')
ggsave(out_png, plot = cl.map, dpi = 500, width = 8, height = 6, bg = 'white')



#==============RADAR PLOTS================#

# Scale the SOM input variables (already done)
# Add cluster labels to each observation (in X)
X_clustered <- X %>%
 mutate(cluster = as.factor(firesheds$som.cluster))

# Compute mean for each variable per cluster
cl_means <- X_clustered %>%
 group_by(cluster) %>%
 summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Rescale for radar plot aesthetics 
cl_means_r <- cl_means %>%
 mutate(across(-cluster, scales::rescale)) %>%
 mutate(group = paste0("Cluster ", cluster)) %>%
 select(group, everything(), -cluster)

###########################################################
# Function to generate ggplot-based radar chart per cluster
make_ggradar <- function(df_row, cluster_id, fill_color) {
 # Plot with ggradar
 ggradar(df_row,
         values.radar = c("0", "0.5", "1"),
         grid.min = 0, grid.mid = 0.5, grid.max = 1,
         group.line.width = 1.0,
         group.point.size = 1.6,
         group.colours = fill_color,
         fill = TRUE, fill.alpha = 0.3,
         axis.label.size = 1.8,
         grid.label.size = 4,
         background.circle.colour = "grey90",
         gridline.mid.colour = "grey80",
         font.radar = "sans",
         legend.position = "none") +
  theme(
   plot.margin = margin(t = 15, r = 80, b = 15, l = 80)  
  )
}

# List of radar plots
(radar_plots <- lapply(1:n_cluster, function(i) {
 df_i <- cl_means_r %>% filter(group == paste0("Cluster ", i))
 make_ggradar(df_i, paste0("Cluster ", i), cluster_colors[[i]])
}))
# panel the radar plots
(radar_grid <- patchwork::wrap_plots(radar_plots, ncol = 2) + 
  plot_layout(
   widths = c(1, 1), 
   heights = rep(1, ceiling(n_cluster / 2)), 
   guides = "collect"
  ) & 
  theme(plot.margin = margin(1,1,1,1)))
# combine with the map
(map_radar <- cl.map + radar_grid + 
 patchwork::plot_layout(
  ncol = 2, widths = c(1.2, 1.2)))

# save the plot.
out_png <- paste0(projdir,'figures/SRM_Firesheds_SOMclusters_wRadar.png')
ggsave(out_png, plot = map_radar, dpi = 300, 
       width = 10, height = 7, bg = 'white')



#==============Cluster Characteristics================#

###########
# Boxplot #
df <- X %>%
 mutate(
  cluster = factor(firesheds$som.cluster),
  cluster_label = paste0("Cluster ", firesheds$som.cluster) 
 ) %>%
 pivot_longer(cols = -c(cluster, cluster_label), 
              names_to = "variable", values_to = "value")
# assign clusters
df$cluster <- factor(df$cluster, levels = as.character(1:n_cluster))
# Plot variable contributions by cluster with custom colors
(p.box <- ggplot(df, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(outlier.size = 0.4, lwd = 0.2) +
  facet_wrap(~ variable, scales = "free") +
  theme_classic(base_size = 10) +
  labs(x = "Cluster", y = "Value") +
  scale_fill_manual(values = cluster_colors, name = "SOM Cluster") +
  theme(
   axis.text.x = element_text(angle = 0, hjust = 0.5),
   legend.position = c(0.92, 0.08),  # manual position
   legend.direction = "vertical",
   legend.title = element_text(angle = 0, vjust = 0.5, size = 11),
   legend.text = element_text(size = 10),
   legend.key.size = unit(0.6, "cm")
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)))
# save the plot.
out_png <- paste0(projdir,'figures/SRM_Firesheds_SOMclusters_VarMetrics.png')
ggsave(out_png, plot = p.box, dpi = 500, 
       width = 8, height = 6, bg = 'white')


#==============LINK TO BIVARIATE CLASSES================#
# Create a table: rows = bi_class, columns = SOM clusters
cluster_biv <- firesheds %>%
 st_drop_geometry() %>%
 count(label, som.cluster) %>%
 pivot_wider(names_from = som.cluster, values_from = n, values_fill = 0) %>%
 arrange(label)
# View the table
print(cluster_biv)
# save CSV
write_csv(cluster_biv, paste0(projdir, "data/tabular/bivariate_cluster_table.csv"))


#==============SUMMARY TABLE WITHIN CLUSTERS================#
# create a summary table for clusters with all variables:
sum.tbl <- firesheds %>%
 st_drop_geometry() %>%
 select(-c(fs_id, sfs_id)) %>%
 group_by(som.cluster) %>%
 summarise(
  n_firesheds = n(),
  sfs_area_km2 = sum(sfs_area_ha) * 0.01,
  aspen_ha = sum(aspen_ha, na.rm=T),
  aspen_km2 = sum(aspen_ha, na.rm=T) * 0.01,
  msbf_count_sum = sum(msbf_count, na.rm=T),
  pop_n = sum(pop_n, na.rm=T),
  across(c(
   whp_p90, hui_p90, cfl_p90, exposure, patch_size,
   delta585, trend_count, aspen10_pct
  ), \(x) mean(x, na.rm = TRUE)),
  .groups = "drop"
 ) %>%
 mutate(som.cluster = as.factor(som.cluster))
sum.tbl %>%
 select(som.cluster, n_firesheds, sfs_area_km2, 
        exposure, aspen_ha, aspen_km2, patch_size,
        msbf_count_sum, pop_n)
# save the table.
write_csv(sum.tbl, paste0(projdir, "data/tabular/som-cluster_summary_table.csv"))




# #==============CLUSTER VIZ================#
# 

# ###########
# # Heatmap #
# # Long format for ggplot heatmap
# cl_melt <- melt(cl_means, id.vars = "cluster")
# # plot it.
# ggplot(cl_melt, aes(x = variable, y = cluster, fill = value)) +
#  geom_tile(color = "white") +
#  scale_fill_viridis_c(option = "C", name = "Z-Score") +
#  theme_minimal(base_size = 11) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1),
#        panel.grid = element_blank()) +
#  labs(title = "SOM Cluster Characterization (Standardized)",
#       x = "Variable", y = "Cluster")

# # Extract grid info
# som_codes <- som.model$codes[[1]]
# unit_coords <- som.model$grid$pts
# unit_classif <- som.model$unit.classif
# # Get SOM unit positions and codebook vectors
# xdim <- som.model$grid$xdim
# ydim <- som.model$grid$ydim
# 
# # Function to get neighbor indices
# get_neighbors <- function(i, coords) {
#  this_coord <- coords[i, ]
#  diffs <- sweep(coords, 2, this_coord)
#  dists <- sqrt(rowSums(diffs^2))
#  neighbors <- which(dists > 0 & dists <= 1.1)  # typical threshold for adjacent hexes
#  return(neighbors)
# }
# 
# # Calculate average distance to neighbors for each SOM unit
# u_matrix_vals <- map_dbl(1:nrow(som_codes), function(i) {
#  neighbors <- get_neighbors(i, unit_coords)
#  if (length(neighbors) == 0) return(NA)
#  dists <- apply(som_codes[neighbors, , drop = FALSE], 1, function(n) {
#   sqrt(sum((som_codes[i, ] - n)^2))
#  })
#  mean(dists)
# })
# 
# # Hierarchical clustering
# dist_mat <- dist(som_codes)
# hc <- hclust(dist_mat, method = "ward.D2")
# n_clusters <- 6
# som_cluster <- cutree(hc, k = n_clusters)
# 
# # Make a SOM dataframe
# som_df <- as_tibble(unit_coords) %>%
#  rename(x = x, y = y) %>%
#  mutate(node = row_number(),
#         uval = u_matrix_vals,
#         cluster = as.factor(som_cluster))
# 
# # Set hex size
# hex_size <- 0.5
# # Create hexagon coordinates
# hex_coords <- som_df %>%
#  rowwise() %>%
#  mutate(hex = list(tibble(
#   hx = x + hex_size * cos(seq(0, 2 * pi, length.out = 7)),
#   hy = y + hex_size * sin(seq(0, 2 * pi, length.out = 7))
#  ))) %>%
#  unnest(hex)
# 
# # make the plot:
# uplot <- ggplot(hex_coords, aes(x = hx, y = hy, group = node)) +
#  geom_polygon(aes(fill = uval), color = "white", size = 0.1) +
#  scale_fill_viridis_c(name = "U-Matrix Distance", option = "C", na.value = "grey90") +
#  geom_path(aes(group = cluster), color = "black", linewidth = 0.4) +
#  coord_equal() +
#  theme_void() +
#  theme(
#   legend.position = "right",
#   plot.title = element_text(hjust = 0.5)
#  ) +
#  labs(title = "SOM U-Matrix with Cluster Outlines")
# uplot


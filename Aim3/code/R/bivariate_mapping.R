#############################################
# Bivariate map of future fire / future aspen
# for management priority landscapes 

# Load the libraries
library(tidyverse)
library(sf)
library(biscale)
library(cowplot)
library(viridis)
library(ggspatial)
library(gridExtra)
library(scales)
library(patchwork)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/"

# Load the data (firesheds)
firesheds <- st_read(paste0(projdir, 'Aim3/data/spatial/mod/srm_firesheds_model_data.gpkg')) %>%
 rename(
  # tidy some of the column names
  patch_size = mean_patch_ha,
  n_patches = number_of_patches,
  patch_den = patch_density,
  big_patch = largest_patch_index,
  pop_den = pop_density_max,
  pop_n = pop_count_sum,
  wui_dist = wui_dist_mean,
  combust = combust_sum,
  msbf_count = msbf_count_sum,
  lf_canopy = forest_cc_mean,
  lf_height = forest_ch_mean,
  exposure = sfs_exposure,
  disturbed = pct_disturbed_y,
  burned_pct = burned_pct_c
 ) %>%
 mutate(
  # calculate contemporary aspen hectares
  aspen_ha = aspen10_pixn * 0.01,
  forest_ha = forest_pixels * 0.09,
  # percent of forested area that is aspen?
  aspen_forest = aspen_ha / forest_ha,
  # combine the WUI classes
  wui = wui1 + wui2 + wui3 + wui4
 )
glimpse(firesheds)


############################################
# Bivariate map: Future Fire -> Future Aspen
# SSP585

# scale the variables
firesheds.biv <- firesheds %>%
 mutate(across(c(trend_count, delta585), ~scale(.)))

# create the bivariate classes
firesheds.biv <- bi_class(
 firesheds, x = trend_count, y = delta585, 
 style = "quantile", dim = 3)

# create the bivariate chloropleth map
p1 <- ggplot() +
 geom_sf(data = firesheds.biv, 
         aes(fill = bi_class), 
         color = NA, linewidth = 0,
         show.legend = F) +
 bi_scale_fill(pal = "GrPink", dim = 3) +
 bi_theme() +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(
  location = "br", width_hint = 0.1,
  pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
  line_width = 0.5, text_pad = unit(0.15,"cm"),
  height = unit(0.15,"cm")
 ) +
 ggspatial::annotation_north_arrow(
  location = "br", which_north = "true",
  pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
  width = unit(0.8,"cm"), height = unit(0.8,"cm"),
  style = north_arrow_fancy_orienteering)

# create the legend
legend <- bi_legend(
 pal = "GrPink",
 dim = 3,
 xlab = "Future Fire ",
 ylab = "Future Aspen ",
 pad_width = 1.4,
 size = 7)

# merge the plot and legend
(p1.1 <- ggdraw() +
 draw_plot(p1, 0, 0, 1, 1) +
 draw_plot(legend, 0.20, .72, 0.25, 0.25))

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_Bivar_FutureFireAspen_SSP585.png')
ggsave(out_png, plot = p1.1, dpi = 500, width = 6, height = 6, bg = 'white')



#####################################################
# Add inset plot showing the total area in each class
# first, define meaningful labels for classes
firesheds.biv <- firesheds.biv %>%
 mutate(
  label = recode(
   bi_class,
   "1-1" = "Low \u2191 Fire,\nDecreasing Aspen",
   "1-2" = "Low \u2191 Fire,\nLow Change Aspen",
   "1-3" = "Low \u2191 Fire,\nIncreasing Aspen",
   "2-1" = "Moderate \u2191 Fire,\nDecreasing Aspen",
   "2-2" = "Moderate \u2191 Fire,\nLow Change Aspen",
   "2-3" = "Moderate \u2191 Fire,\nIncreasing Aspen",
   "3-1" = "High \u2191 Fire,\nDecreasing Aspen",
   "3-2" = "High \u2191 Fire,\nLow Change Aspen",
   "3-3" = "High \u2191 Fire,\nIncreasing Aspen"
 ))

# Compute total area per bivariate class
(biv_area <- firesheds.biv %>%
 as_tibble() %>%
 mutate(sfs_area_km2 = sfs_area_ha * 0.01) %>%
 group_by(label, bi_class) %>%
 summarize(n_firesheds = n(),
           area_ha = sum(sfs_area_ha),
           area_km2 = sum(sfs_area_km2)))

# Plot
(p_area <- ggplot(biv_area, 
                 aes(x = reorder(label, area_ha), 
                     y = area_ha, fill = bi_class)) +
 geom_col() +
 geom_text(aes(label = paste0("N = ",n_firesheds)),
           hjust = -0.1, size = 3) +
 bi_scale_fill(pal = "GrPink", dim = 3) +
 theme_classic() +
  theme(
   legend.position = "none",
   plot.margin = margin(t = 5, r = 30, b = 5, l = 5)
  ) +
 labs(x = "Fire-Aspen Class", y = "Total Fireshed Area (ha)") +
 scale_y_continuous(
  breaks = c(1e5, 1e6, 2e6, 3e6), 
  labels = label_number(scale_cut = cut_short_scale())  
 ) +
 coord_flip(clip = "off"))

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_area.png')
ggsave(out_png, plot = p_area, height = 4, width = 4, bg = 'white')

# combine
(p3 <- ggdraw() +
 draw_plot(p1.1, -0.2, 0, 1, 1) +  # Main map takes full space
 draw_plot(p_area, 0.5, 0.15, 0.48, 0.75))

# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_Panelarea.png')
ggsave(out_png, plot = p3, height = 7, width = 9, bg = 'white')


##########################################
# Join firesheds to get unscaled variables
firesheds <- firesheds %>%
 left_join(firesheds.biv%>%as_tibble()%>%select(c(sfs_id, bi_class, label)),
           by = "sfs_id")
glimpse(firesheds)


##########################################
# Create a summary table for all variables
firesheds <- firesheds %>%
 mutate(
  label = recode(
   bi_class,
   "1-1" = "Low Increasing Fire, Decreasing Aspen",
   "1-2" = "Low Increasing Fire, Low Change Aspen",
   "1-3" = "Low Increasing Fire, Increasing Aspen",
   "2-1" = "Moderate Increasing Fire, Decreasing Aspen",
   "2-2" = "Moderate Increasing Fire, Low Change Aspen",
   "2-3" = "Moderate Increasing Fire, Increasing Aspen",
   "3-1" = "High Increasing Fire, Decreasing Aspen",
   "3-2" = "High Increasing Fire, Low Change Aspen",
   "3-3" = "High Increasing Fire, Increasing Aspen"
  ))
sum.tbl <- firesheds %>%
 st_drop_geometry() %>%
 select(-c(fs_id, sfs_id)) %>%
 group_by(label, bi_class) %>%
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
 mutate(bi_class = as.factor(bi_class))
sum.tbl %>%
 select(label, bi_class, n_firesheds, sfs_area_km2, 
        exposure, aspen_ha, aspen_km2, patch_size,
        msbf_count_sum, pop_n)
# save the table.
write_csv(sum.tbl, paste0(projdir, "Aim3/data/tabular/bivariate_summary_table.csv"))


##############################
# Write out the spatial data #
out_fp = paste0(projdir,"Aim3/data/spatial/mod/srm_firesheds_model_data_wBivar.gpkg")
st_write(firesheds, out_fp, append=F)



# #########################################
# # Summarize current aspen cover per class
# biv_aspen <- firesheds.biv %>%
#  group_by(label, bi_class) %>%
#  summarise(aspen_ha = sum(aspen10_pixn) * 0.01,
#            area_ha = sum(sfs_area_ha))
# # Plot
# p_aspen <- ggplot(biv_aspen, aes(x = reorder(label, area_ha), 
#                                  y = aspen_ha, fill = bi_class)) +
#  geom_col() +
#  bi_scale_fill(pal = "GrPink", dim = 3) +
#  theme_minimal() +
#  labs(x = "Bivariate Class", y = "Current Aspen Cover (ha)") +
#  scale_y_continuous(
#   breaks = c(1e4, 5e4, 1e5, 5e6), 
#   labels = label_number(scale_cut = cut_short_scale())  
#  ) +
#  theme(legend.position = "none") +
#  coord_flip()
# p_aspen
# 
# # combine
# p_area / p_aspen  # Stacked layout
# 
# #####
# # try a combined version
# biv_combined <- firesheds.biv %>%
#  group_by(label) %>%
#  summarise(
#   total_area_ha = sum(sfs_area_ha),
#   aspen_area_ha = sum(aspen10_pixn) * 0.01  # Convert pixels to ha
#  ) %>%
#  pivot_longer(cols = c(total_area_ha, aspen_area_ha), 
#               names_to = "Metric", values_to = "Value")
# 
# ggplot(biv_combined, aes(x = reorder(label, Value), 
#                          y = Value, fill = Metric)) +
#  geom_col(position = "dodge") +  # Side-by-side bars
#  scale_fill_manual(values = c("total_area_ha" = "gray70", "aspen_area_ha" = "forestgreen"),
#                    labels = c("Total Area", "Current Aspen Area")) +
#  theme_minimal() +
#  labs(x = "Fire-Aspen Class", y = "Area (ha)", fill = "Metric") +
#  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
#  coord_flip()
# 
# # proportion line chart
# biv_prop <- firesheds.biv %>%
#  group_by(label) %>%
#  summarise(
#   total_area_ha = sum(sfs_area_ha),
#   aspen_area_ha = sum(aspen10_pixn) * 0.01
#  ) %>%
#  mutate(aspen_ratio = aspen_area_ha / total_area_ha)  # Proportion of area that is aspen
# 
# p_total <- ggplot(biv_prop, aes(x = reorder(label, total_area_ha), 
#                                 y = total_area_ha)) +
#  geom_line(group = 1, color = "gray70", size = 1.2) +
#  geom_point(color = "gray40", size = 3) +
#  theme_minimal() +
#  labs(x = NULL, y = "Total Area (ha)") +
#  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
#  coord_flip()
# 
# p_aspen <- ggplot(biv_prop, aes(x = reorder(label, aspen_ratio), 
#                                 y = aspen_ratio)) +
#  geom_line(group = 1, color = "forestgreen", size = 1.2) +
#  geom_point(color = "darkgreen", size = 3) +
#  theme_minimal() +
#  labs(x = "Fire-Aspen Class", y = "Aspen Cover (%)") +
#  scale_y_continuous(labels = percent_format()) +
#  coord_flip()
# 
# # Combine them in a multipanel layout
# p_total / p_aspen  # Stacked layout
# 
# # create a merged plot
# p3 <- p2.1 | (p_total / p_aspen) + 
#  plot_layout(widths = c(3, 0.5))  # Map takes 2x space compared to graphs
# p3
# 
# # save the plot.
# out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_PanelMap.png')
# ggsave(out_png, plot = p3, height = 8, width = 8, bg = 'white')

# #############
# # WUI summary
# wui_summary <- future.fire.bi %>%
#  as_tibble() %>%
#  group_by(bi_class) %>%
#  summarize(
#   WUI_Interface = mean(wui3, na.rm = TRUE),
#   WUI_Intermix = mean(wui4, na.rm = TRUE),
#   Distance_to_WUI = mean(wui_dist_mean, na.rm = TRUE)
#  ) %>%
#  pivot_longer(cols = c(WUI_Interface, WUI_Intermix), 
#               names_to = "WUI_Class", values_to = "Percent_Cover")
# glimpse(wui_summary)
# 
# # plot it
# wui_barplot <- ggplot(wui_summary, aes(x = bi_class, y = Percent_Cover, fill = WUI_Class)) +
#  geom_col(position = "dodge", color = "black") +  # Grouped bars with outlines
#  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")) +  # Distinct colors
#  theme_minimal(base_size = 12) +
#  labs(x = "Bivariate Class", y = "Mean % Cover", fill = "WUI Class") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# wui_barplot
# 
# # distance to wui
# wui_boxplot <- ggplot(future.fire.bi, aes(x = bi_class, y = wui_dist_mean, fill = bi_class)) +
#  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxes
#  geom_jitter(alpha = 0.2, width = 0.2, size = 0.8) +  # Light scatter points
#  bi_scale_fill(pal = "DkViolet", dim = 3) +  # Keep consistent with bivariate map
#  theme_minimal(base_size = 12) +
#  labs(x = "Bivariate Class", y = "Distance to WUI (m)") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# wui_boxplot
# 
# 
# # create a combined plot
# # Arrange the two summary plots in a vertical stack (right side)
# # Merge with the bivariate map (left side)
# final_plot <- plot_grid(
#  p1.1,        # Bivariate map
#  counts_plot,  # Stacked plots
#  ncol = 2,    # Two-column layout
#  labels = c("A", "B"), 
#  rel_widths = c(2, 1.3)
# )
# final_plot
# 
# # save the plot.
# out_png <- paste0(projdir,'Aim3/figures/SouthernRockies_Bivar_FutureFireAspen_wInsets.png')
# ggsave(out_png, plot = final_plot, dpi = 500, bg = 'white')



# Load the libraries
library(tidyverse)
library(sf)
library(biscale)
library(cowplot)
library(viridis)
library(ggspatial)
library(gridExtra)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/"

# Load the data (firesheds)
firesheds <- st_read(paste0(projdir, 'Aim3/data/spatial/mod/srm_firesheds_model_data.gpkg')) %>%
 rename(
  # tidy some of the column names
  n_patches = number_of_patches,
  patch_den = patch_density,
  prop_ls = proportion_of_landscape,
  big_patch = largest_patch_index,
  pop_den = pop_density_max,
  pop_n = pop_count_sum,
  wui_dist = wui_dist_mean,
  combust = combust_sum,
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

###################
# plotting function
plot_metric <- function(data, var, legend_name) {
 ggplot(data) +
  geom_sf(aes_string(fill = var), color = NA) + 
  scale_fill_distiller(palette = "Greens", direction = 1, 
                       name = legend_name, na.value = "white") +
  guides(fill = guide_colourbar(direction = "vertical", 
                                barwidth = 0.6, 
                                barheight = 8, 
                                ticks.colour = NA, 
                                title.position = "left")) +
  theme_void() +
  theme(
   legend.title = element_text(angle = 90, size = 10),
   legend.text = element_text(size = 9),
   legend.position = c(0.35, 0.95),  
   legend.justification = c(1, 1)
  ) +
  coord_sf(expand = FALSE) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                              pad_x = unit(0.01, "in"), pad_y = unit(0.05, "in"),
                              line_width = 0.5, text_pad = unit(0.15, "cm"),
                              height = unit(0.15, "cm")) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                    pad_x = unit(0.01, "in"), pad_y = unit(0.2, "in"),
                                    width = unit(0.8, "cm"), height = unit(0.8, "cm"),
                                    style = north_arrow_fancy_orienteering)
}

######################################################
#============= aspen suitability map ================#

(p1 <- plot_metric(firesheds, "historic", "Historic aspen suitability"))
# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_FutureAspen_Historic.png')
ggsave(out_png, plot = p1, dpi = 500, bg = 'white')

################################
# change in future suitability #
(p2 <- plot_metric(firesheds, "delta585", "Aspen suitability change") +
  scale_fill_distiller(palette = "RdBu", direction = 1, 
                       name = "Aspen Suitability Change"))
# save the plot.
out_png <- paste0(projdir,'Aim3/figures/SRM_Firesheds_FutureAspen_SSP585.png')
ggsave(out_png, plot = p2, dpi = 500, bg = 'white')

######################
# merge the two maps #
p_merged <- grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ncol = 2)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_FutureAspen_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 8, height = 6, bg = "white")


################################################
#============= future fire map ================#
# fire count
(p1 <- plot_metric(firesheds, "trend_count", "Fire occurrence trend") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       name = "Fire occurrence trend"))
# area burned
(p2 <- plot_metric(firesheds, "trend_area", "Burned area trend") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       name = "Burned area trend"))
######################
# merge the two maps #
p_merged <- grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ncol = 2)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_FutureFire_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 8, height = 6, bg = "white")


##################################################################
#============= Fireshed exposure and disturbance ================#
# fire count
(p1 <- plot_metric(firesheds, "exposure", "Wildfire exposure") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       name = "Wildfire exposure"))
# area burned
(p2 <- plot_metric(firesheds, "disturbed", "Disturbed (%)") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       name = "Disturbed (%)"))
######################
# merge the two maps #
p_merged <- grid.arrange(ggplotGrob(p1), ggplotGrob(p2), ncol = 2)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_ExposureDisturbed_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 8, height = 6, bg = "white")


######################################################
#============= built environment maps================#

####################
# WUI designations #
# calculate the sum of WUI classes
firesheds <- firesheds %>%
 mutate(wui_sum = wui1 + wui2 + wui3 + wui4)

p1 <- ggplot(firesheds) +
 geom_sf(data = firesheds %>% filter(wui_sum == 0), fill = "gray90", color = NA) + 
 geom_sf(data = firesheds %>% filter(wui_sum > 0), aes(fill = wui_sum), color = NA) +  
 scale_fill_viridis(option = "rocket", name = "Total WUI %", direction=-1) +
 guides(fill = guide_colourbar(
  direction = "vertical", 
  barwidth = 0.4, 
  barheight = 6, 
  ticks.colour = NA, 
  title.position = "left")
 ) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.35, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p1

#################
# COMBUST (sum) #
p2 <- ggplot(firesheds) +
 geom_sf(data = firesheds %>% filter(combust_sum == 0), fill = "gray90", color = NA) + 
 geom_sf(data = firesheds %>% filter(combust_sum > 0), aes(fill = combust_sum), color = NA) + 
 scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                    name = "COMBUST (sum)", trans="sqrt",
                    labels = scales::label_number(scale = 1e-4, suffix = "k")) +
 guides(fill = guide_colourbar(
  direction = "vertical", 
  barwidth = 0.4, 
  barheight = 6, 
  ticks.colour = NA, 
  title.position = "left")
 ) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.35, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p2

###################
# Distance to WUI #
p3 <- ggplot(firesheds) +
 geom_sf(data = firesheds, aes(fill = wui_dist_mean), color = NA) + 
 scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                    name = "Distance to WUI") +
 # scale_fill_distiller(palette = "Greens", direction = 1, 
 #                      name = "Historic aspen suitability") +
 guides(fill = guide_colourbar(direction = "vertical", 
                               barwidth = 0.4, 
                               barheight = 6, 
                               ticks.colour = NA, 
                               title.position = "left")) +
 theme_void() +
 theme(legend.title = element_text(angle = 90, size = 10),
       legend.text = element_text(size = 9),
       legend.position = c(0.35, 0.95),  
       legend.justification = c(1, 1)) +
 coord_sf(expand = FALSE) +
 ggspatial::annotation_scale(location = "br", width_hint = 0.1,
                             pad_x = unit(0.01,"in"), pad_y = unit(0.05,"in"),
                             line_width = 0.5, text_pad = unit(0.15,"cm"),
                             height = unit(0.15,"cm")) +
 ggspatial::annotation_north_arrow(location = "br", which_north = "true",
                                   pad_x = unit(0.01,"in"), pad_y= unit(0.2,"in"),
                                   width = unit(0.8,"cm"), height = unit(0.8,"cm"),
                                   style = north_arrow_fancy_orienteering)
p3
######################
# merge the two maps #
p_merged <- grid.arrange(
 ggplotGrob(p1), 
 ggplotGrob(p2), 
 ggplotGrob(p3), 
 ncol = 3
)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_BuiltEnv_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 9, height = 6, bg = "white")


################################################
#============= vegetation maps ================#

# Reshape data to long format
firesheds_l <- firesheds %>%
 pivot_longer(cols = c(lodgepole, aspen, ponderosa, pinon_juniper, 
                       sagebrush, spruce_fir, douglas_fir, white_fir, gambel_oak),
              names_to = "forest", values_to = "pct_cover")

# create a facet map
p1 <- ggplot(firesheds_l) +
 geom_sf(aes(fill = pct_cover), color = NA) +
 scale_fill_distiller(palette = "Greens", direction = 1,
                      name = "Percent cover", na.value = "gray80") +
 facet_wrap(~forest, nrow = 2) +  # Creates small multiples by forest type
 theme_void() +
 theme(legend.position = c(0.95, 0.26),
       legend.title = element_text(angle = 90, size = 10),
       legend.key.width = unit(4, "cm"),
       legend.key.height = unit(0.4, "cm"),
       strip.text = element_text(size = 9, face = "italic",
                                 margin = margin(b = 2)),
       panel.spacing = unit(0.5, "lines")) +
 guides(fill = guide_colorbar(
  ticks.colour = NA, 
  title.position = "left",
  title.hjust = 0.5,  
  barwidth = unit(0.65, "cm"), 
  barheight = unit(6.5, "cm")  
 ))
p1

# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_EVTSAF_Panel.png")
ggsave(out_png, plot = p1, dpi = 500, width = 9, height = 6, bg = "white")


####################################################
#============= risk to communities ================#
(p1 <- plot_metric(firesheds, "hui_p90", "Housing Unit Impact") +
  scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                     name = "Housing Unit Impact"))
(p2 <- plot_metric(firesheds, "whp_p90", "Wildfire Hazard Potential") +
  scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                     name = "Wildfire Hazard Potential"))
(p3 <- plot_metric(firesheds, "cfl_p90", "Cond. Flame Length") +
  scale_fill_viridis(option = "rocket", na.value = "gray80", direction=-1,
                     name = "Cond. Flame Length"))
#########################
p_merged <- grid.arrange(
 ggplotGrob(p3), 
 ggplotGrob(p2), 
 ggplotGrob(p1), 
 ncol = 3
)
p_merged
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_RiskToCommunities_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 500, width = 9, height = 6, bg = "white")


####################################################
#============= aspen patch metrics ================#
p1 <- plot_metric(firesheds, "number_of_patches", "Number of patches")
p2 <- plot_metric(firesheds, "patch_density", "Patch Density")
p3 <- plot_metric(firesheds, "largest_patch_index", "Largest patch")
p4 <- plot_metric(firesheds, "mean_patch_ha", "Average patch size")
#########################
(p_merged <- grid.arrange(
 ggplotGrob(p1), # number of patches
 ggplotGrob(p4), # mean patch size
 ggplotGrob(p2), # patch density
 ggplotGrob(p3), # largest patch index
 ncol = 4, nrow = 1
))
# Save the merged plot
out_png <- paste0(projdir, "Aim3/figures/SRM_Firesheds_AspenPatches_Panel.png")
ggsave(out_png, plot = p_merged, dpi = 300, width = 10, height = 6, bg = "white")



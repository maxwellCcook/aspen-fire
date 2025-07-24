
################################
# Plotting INLA model results. #

library(tidyverse)
library(ggridges)
library(patchwork)

maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# Load the data
frp.df <- read_csv(paste0(maindir, "data/tabular/mod/results/INLA_tidy_effects_FRP.csv"))
cbi.df <- read_csv(paste0(maindir, "data/tabular/mod/results/INLA_tidy_effects_CBI.csv"))

# merge into one dataframe
frp.df$response <- 'FRPc'
cbi.df$response <- 'CBIbc'
combined.df <- bind_rows(frp.df, cbi.df)

glimpse(combined.df)

rm(frp.df, cbi.df)

##########################
# Set some figure parameters

# response color map
response_cmap <- c(
 "CBIbc" = "#800026",
 "FRPc" = "#FEB24C"
)

# order for structure metrics
param_order <- c(
 "Forest type proportion",
 "Live basal area",
 "Diameter/height ratio",
 "Quadratic mean diameter",
 "Stand density index"
)

# alpha mapping (transparency)
alpha_map <- c(
 "Quaking aspen" = 0.8,
 "Lodgepole pine" = 0.6,
 "Douglas-fir" = 0.6,
 "White fir" = 0.6,
 "Gambel oak" = 0.6,
 "Piñon-juniper" = 0.6,
 "Ponderosa pine" = 0.6,
 "Spruce-fir" = 0.6
)


#========Plot 1: Forest type effects relative to aspen========#

###################################
# Extract the forest type effects #

fortyp_effects <- combined.df %>%
 filter(
  str_detect(parameter, "fortypnm_gp"),
  !str_detect(parameter, ":fortyp_pct"),
 ) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(
  fill_species = fct_reorder(fill_species, -mean_effect)
 )

# order the forest types from low->high on the main effect
fortyp_order <- fortyp_effects %>%
 filter(response == "FRPc") %>%
 group_by(fill_species) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(mean_effect)

# apply the ordering
fortyp_effects$fill_species <- factor(
 fortyp_effects$fill_species, 
 levels = fortyp_order$fill_species
)

# glimpse the resulting table
glimpse(fortyp_effects)


############
# ridge plot

(p1 <- fortyp_effects %>%
  ggplot(., aes(x = x, y = fill_species, height = y, fill = response)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.85) +
  # geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect relative to aspen",
   y = "Forest Type",
   # fill = "Response"
   fill = ""
  ) +
  scale_fill_manual(
   values = response_cmap,
   labels = c("FRPc" = "Cumulative FRP",
              "CBIbc" = "90th Percentile CBI")
  ) +
  scale_color_manual(values = rep(scales::alpha("black", 0.6), 2)) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size = 9),
   axis.text.x = element_text(angle = 0, hjust = 0, size = 9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_text(size = 10, margin = margin(t = 12)),
   legend.position = c(0.24, 0.94),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.4),
    color = NA, linewidth = 0.8
   ),
   legend.title = element_text(size = 9),
   legend.text = element_text(size = 8)
  ))

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP_species.png')
ggsave(out_png, plot = p1, dpi = 500, width = 7, height = 4, bg = 'white')


#========Plot 1b: Effect of increasing proportional cover of forest types========#

###################################
# Extract the forest type effects #

fortyp_pct_effects <- combined.df %>%
 filter(
  str_detect(parameter, ":fortyp_pct"),
 ) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(
  fill_species = fct_reorder(fill_species, -mean_effect)
 )

# make the ridge plot
(p1b <- fortyp_pct_effects %>%
  ggplot(., aes(x = x, y = fill_species, height = y, fill = response)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect size",
   y = "Forest Type",
   fill = ""
  ) +
  scale_fill_manual(
   values = response_cmap,
   labels = c("FRPc" = "Cumulative FRP",
              "CBIbc" = "90th Percentile CBI")
  ) +
  scale_color_manual(values = rep(scales::alpha("black", 0.6), 2)) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size = 9),
   axis.text.x = element_text(angle = 0, hjust = 0, size = 9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_text(size = 10, margin = margin(t = 12)),
   legend.position = c(0.82, 0.92),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.4),
    color = NA, linewidth = 0.8
   ),
   legend.title = element_text(size = 9),
   legend.text = element_text(size = 8)
  ))


#========Plot 2: Forest structure effects========#

# filter and prepare data
spp.effects <- combined.df %>%
 filter(
  # str_detect(parameter, ":fortyp_pct"),
  effect %in% param_order
 ) %>%  # keep species-specific effects
 mutate(
  effect = factor(effect, levels = param_order),
  fill_species = factor(fill_species, levels = sort(unique(fill_species), decreasing = TRUE))
 )
glimpse(spp.effects)

# order species by proportion forest area effect

species_order <- combined.df %>%
 filter(
  str_detect(parameter, ":fortyp_pct"),
  response == "FRPc"
 ) %>%
 group_by(fill_species) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(mean_effect) %>%  # lowest effect first (top of plot)
 pull(fill_species)

# order the dataframe:

spp.effects <- spp.effects %>%
 mutate(
  fill_species = factor(fill_species, levels = species_order)
 )

# ridge plot
(p2 <- ggplot(spp.effects, aes(x = x, y = fill_species, height = y, fill=response)) +
  geom_density_ridges(
   stat = "identity", 
   scale = 1.3,
   color = scales::alpha("black", 0.4),
   aes(alpha = fill_species),
   linewidth = 0.4
  ) +
  scale_fill_manual(
   values = response_cmap,
   labels = c("FRPc" = "Cumulative FRP",
              "CBIbc" = "90th Percentile CBI")
  ) +
  scale_alpha_manual(values = alpha_map, guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~effect, ncol = 4, scales = "free_y") +
  labs(
   x = "Effect size",
   y = NULL  # remove y-label from each subplot
  ) +
  theme_classic() +
  theme(
   strip.background = element_blank(),
   strip.text = element_text(size = 10, face = "bold"),
   axis.text.x = element_text(size = 9),
   axis.text.y = element_text(size = 8),
   axis.title.x = element_text(size = 11, margin = margin(t = 10)),
   axis.title.y = element_blank(),
   panel.spacing = unit(1.4, "lines")
  ))


#========Plot 3: All other fixed effects========#

#####################################################
# create the ridge plot for all other fixed effects #
fixed_effects.frp <- tidy.effects.frp %>%
 filter(!effect %in% c("aspen1","overlap","log_fire_size"),
        !str_detect(parameter,"fortypnm"),
        !str_detect(parameter,"species"))

# plot the ridges
(frp.p3 <- ggplot(fixed_effects.frp, 
                  aes(x = x, y = effect, height = y)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect Size",
   y = "Fixed Effect",
  ) +
  theme_minimal() +
  theme(
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 8),
   axis.text.y = element_text(size = 10),
   axis.title.y = element_text(size = 12)
  ))




#========90th Percentile Composite Burn Index (CBIbc90)========#

###################################
# Extract the forest type effects #
fortyp_effects.cbi <- tidy.effects.cbi %>%
 filter(str_detect(parameter, "fortypnm_gp"),
        !str_detect(parameter, ":fortyp_pct")) %>%
 mutate(
  group = case_when(
   str_detect(parameter, "vpd:") ~ "VPD-mediated",
   TRUE ~ "Forest type"  # Default for all other fixed effects
  )) %>%
 group_by(fill_species) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() %>%
 mutate(fill_species = fct_reorder(fill_species, -mean_effect))
glimpse(fortyp_effects.cbi)

# make the ridge plot
(cbi.p1 <- ggplot(fortyp_effects.cbi,
                  aes(x = x, y = fill_species, height = y,
                      fill = group, alpha=group)) +
  geom_density_ridges(stat = "identity", scale = 1.5, aes(alpha = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect relative to aspen",
   y = "Forest Type"
  ) +
  scale_fill_manual(
   values = c(
    "VPD-mediated" = "#800026",
    "Forest type" = "#800026"
   )
  ) +
  scale_alpha_manual(
   values = c(
    "VPD-mediated" = 0.2,
    "Forest type" = 0.98
   )
  ) +
  theme_minimal() +
  theme(
   legend.position = "none",
   axis.text.y = element_text(size = 10),
   axis.title.y = element_text(size = 11)
  ))


#######################################################
# create the ridge plot for species structure metrics #
# filter to structure effects
sp_effects.cbi <- tidy.effects.cbi %>%
 filter(effect %in% param_order)

# plot it
(cbi.p2 <- sp_effects.cbi %>%
  mutate(effect = factor(effect, levels = param_order)) %>%
  ggplot(., aes(
   x = x, y = effect, height = y, 
   fill = fill_species, 
   alpha = fill_species,
   color = fill_species
  )) +
  geom_density_ridges(
   stat = "identity", scale = 1.5, show.legend=T,
   color = scales::alpha("black", 0.2),
   linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "",
   y = "Fixed Effect",
   fill = "Species"
  ) +
  # add a subplot label
  annotate("text", x = -0.20, y = 5.2,
           label = expression(bold("(B)") ~ "CBIbc"),
           size = 4, hjust = 0.2) +
  scale_fill_manual(
   values = color_map,
   breaks = names(alpha_map),
   guide = guide_legend(
    title = "Forest type",
    override.aes = list(
     fill = scales::alpha(unname(color_map[names(alpha_map)]), 
                          unname(alpha_map)),
     alpha = unname(alpha_map)
    )
   )
  ) +
  scale_alpha_manual(
   values = alpha_map,
   guide = "none"  # Still hide alpha scale
  ) +
  xlim(-0.20, 0.20) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size=10),
   axis.text.x = element_text(angle = 0, hjust = 0, size=10),
   axis.title.y = element_text(size = 11, margin = margin(r = 8)),
   axis.title.x = element_text(size = 11, margin = margin(t = 8)),
   legend.position = c(0.94, 0.52),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.6), 
    color = NA, linewidth = 1),
   legend.title = element_text(size = 9),
   legend.text = element_text(size = 8)
  ))
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_CBI_species.png')
ggsave(out_png, plot = cbi.p2, dpi = 500, width = 7, height = 5, bg = 'white')

##############################
# Option #2: faceted by metric

# filter and prepare data
sp_effects.cbi_ <- tidy.effects.cbi %>%
 filter(effect %in% param_order,
        !str_detect(parameter, "vpd"),
        !is.na(species)) %>%  # keep species-specific effects
 mutate(
  effect = factor(effect, levels = param_order),
  fill_species = factor(fill_species, levels = rev(levels(fill_species)))  # reverse species order
 )

# plot
(frp.p2_f <- ggplot(sp_effects.cbi_, aes(
 x = x, y = fill_species, height = y
)) +
  geom_density_ridges(
   stat = "identity", scale = 1.3,
   fill = "#B2182B",  # deep red
   color = scales::alpha("black", 0.4),
   alpha = 0.8, linewidth = 0.4
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~effect, ncol = 4, scales = "free_y") +
  labs(
   x = "Effect size",
   y = NULL  # remove y-label from each subplot
  ) +
  theme_classic() +
  theme(
   strip.background = element_blank(),
   strip.text = element_text(size = 10, face = "bold"),
   axis.text.x = element_text(size = 9),
   axis.text.y = element_text(size = 8),
   axis.title.x = element_text(size = 11, margin = margin(t = 10)),
   axis.title.y = element_blank(),
   panel.spacing = unit(1.4, "lines")
  ))


#####################################################
# create the ridge plot for all other fixed effects #
fixed_effects.cbi <- tidy.effects.cbi %>%
 filter(!effect %in% c("aspen1","overlap"),
        !str_detect(parameter,"fortypnm"),
        !str_detect(parameter,"species"),
        !fill_species %in% c("Gambel oak","White fir"))

# plot the ridges
(cbi.p3 <- ggplot(fixed_effects.cbi, 
                  aes(x = x, y = effect, height = y)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect Size",
   y = "Fixed Effect",
  ) +
  theme_minimal() +
  theme(
   legend.title = element_text(size = 10),
   legend.text = element_text(size = 8),
   axis.text.y = element_text(size = 10),
   axis.title.y = element_text(size = 12)
  ))


#######################################
#######################################
#######################################
#######################################
# Make combined figures (FRP and CBI) #

# Create the combined forest type relative to aspen
fortyp_effects.frp$response <- "FRP"
fortyp_effects.cbi$response <- "CBI"
fortyp_effects <- bind_rows(fortyp_effects.frp, fortyp_effects.cbi) %>%
 # reorder the effects so main effect is drawn on top
 mutate(
  group = factor(group, levels = c("VPD-mediated", "Forest type")),
  fill_cat = factor(
   paste(response, group),
   levels = c(
    "CBI VPD-mediated", "FRP VPD-mediated", 
    "CBI Forest type", "FRP Forest type"
   ))
 )
# order the forest types from low->high on the main effect
forest_order <- fortyp_effects %>%
 filter(group == "Forest type", response == "FRP") %>%
 group_by(fill_species) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(mean_effect)

# apply the ordering
fortyp_effects$fill_species <- factor(
 fortyp_effects$fill_species, 
 levels = rev(forest_order$fill_species)
)

############################
# plot without VPD-mediation
# Filter to just forest type effects
(p4 <- fortyp_effects %>%
  filter(group == "Forest type") %>%
  ggplot(., aes(
   x = x, y = fill_species, height = y,
   fill = fill_cat
  )) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.85) +
  # geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect relative to aspen",
   y = "Predominant Forest Type",
   fill = "Effect Type"
  ) +
  scale_fill_manual(
   values = cmap,
   labels = c("FRP Forest type" = "Cumulative FRP",
              "CBI Forest type" = "90th percentile CBIbc")
  ) +
  scale_color_manual(values = rep(scales::alpha("black", 0.6), 2)) +
  # scale_color_manual(values = rep(scales::alpha("black", 0.3), 2)) +
  # coord_cartesian(xlim = c(-0.30, 0.36)) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size = 9),
   axis.text.x = element_text(angle = 0, hjust = 0, size = 9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_text(size = 10, margin = margin(t = 12)),
   legend.position = c(0.18, 0.18),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.4),
    color = NA, linewidth = 0.8
   ),
   legend.title = element_text(size = 9),
   legend.text = element_text(size = 8)
  ))
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP-CBI_species.png')
ggsave(out_png, dpi = 500, width = 7, height = 4, bg = 'white')

############################
# plot without VPD-mediation
# Filter to just forest type effects
(p4 <- fortyp_effects %>%
  filter(group == "Forest type") %>%
  ggplot(., aes(
   x = x, y = fill_species, height = y,
   fill = fill_cat
  )) +
  # geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.95) +
  geom_density_ridges(stat = "identity", scale = 1.5, color = NA, alpha = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect relative to aspen",
   y = "Predominant Forest Type",
   fill = "Effect Type"
  ) +
  scale_fill_manual(
   values = cmap,
   labels = c("FRP Forest type" = "Cumulative FRP",
              "CBI Forest type" = "90th percentile CBIbc")
  ) +
  # scale_color_manual(values = rep(scales::alpha("black", 0.6), 2)) +
  scale_color_manual(values = rep(scales::alpha("black", 0.3), 2)) +
  # coord_cartesian(xlim = c(-0.30, 0.36)) +
  theme_classic() +
  theme(
   axis.text.y = element_text(angle = 0, hjust = 1, size = 9),
   axis.text.x = element_text(angle = 0, hjust = 0, size = 9),
   axis.title.y = element_text(size = 10, margin = margin(r = 12)),
   axis.title.x = element_text(size = 10, margin = margin(t = 12)),
   legend.position = c(0.18, 0.18),
   legend.background = element_rect(
    fill = scales::alpha("white", 0.4),
    color = NA, linewidth = 0.8
   ),
   legend.title = element_text(size = 9),
   legend.text = element_text(size = 8)
  ))
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP-CBI_BLANK.png')
ggsave(out_png, dpi = 500, width = 7, height = 4, bg = 'white')


# ##########################
# # Plot with VPD-mediation
# (p4.1 <- ggplot(fortyp_effects, aes(x = x, y = fill_species, height = y, 
#                                  fill = fill_cat, 
#                                  alpha = fill_cat,
#                                  color = fill_cat)) +
#  geom_density_ridges(stat = "identity", scale = 1.5) +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#  labs(
#   x = "Effect relative to aspen",
#   y = "Predominant Forest Type",
#   fill = "Effect Type",
#  ) +
#  scale_fill_manual(
#   values = c(
#    "FRP Forest type" = "#FEB24C",  
#    "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
#    "CBI Forest type" = "#800026",  
#    "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
#   ),
#   labels = c(
#    "Cumulative FRP", "VPD-mediated (FRPc)",
#    "90th percentile CBIbc", "VPD-mediated (CBIbc)"
#   ),
#   guide = guide_legend(
#    override.aes = list(
#     fill = c("#FEB24C", scales::alpha("#FEB24C", 0.3), 
#              "#800026", scales::alpha("#800026", 0.3)),
#     alpha = c(1, 0.3, 1, 0.3),  # Match transparency
#     color = c(
#      scales::alpha("black", 0.7),
#      scales::alpha("#FEB24C", 0.3), 
#      scales::alpha("black", 0.7), 
#      scales::alpha("#800026", 0.3))
#    ),
#    order = 1  # Ensure this legend appears first
#   )
#  ) +
#  scale_alpha_manual(
#   values = c(
#    "FRP Forest type" = 1,
#    "FRP VPD-mediated" = 0.3,
#    "CBI Forest type" = 1,
#    "CBI VPD-mediated" = 0.3
#   ),
#   guide = "none"  
#  ) +
#  scale_color_manual(
#   values = c(
#    "FRP Forest type" = scales::alpha("black", 0.6),
#    "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
#    "CBI Forest type" = scales::alpha("black", 0.6),  
#    "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
#   ),
#   guide = "none"
#  ) +
#  coord_cartesian(xlim=c(-0.38,0.44)) +
#  theme_classic() +
#  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
#        axis.text.x = element_text(angle = 0, hjust = 0, size=9),
#        axis.title.y = element_text(size = 10, margin = margin(r = 12)),
#        axis.title.x = element_text(size = 10, margin = margin(t = 12)),
#        legend.position = c(0.18, 0.30),
#        legend.background = element_rect(
#         fill = scales::alpha("white", 0.4), 
#         color = NA, size = 0.8),
#        legend.title = element_text(size = 9),
#        legend.text = element_text(size = 8)))
# # Save the plot
# out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP-CBI_species_VPD.png')
# ggsave(out_png, dpi = 500, width = 7, height = 4, bg = 'white')


# #############################################
# # version 2: without Gambel oak and white fir
# p4.1 <- fortyp_effects %>%
#  filter(!fill_species %in% c("Gambel oak")) %>%
#  ggplot(., aes(x = x, y = fill_species, height = y, 
#                fill = fill_cat, 
#                alpha = fill_cat,
#                color = fill_cat)) +
#  geom_density_ridges(stat = "identity", scale = 1.5) +
#  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#  labs(
#   x = "Effect relative to aspen",
#   y = "Predominant Forest Type",
#   fill = "Effect Type",
#  ) +
#  scale_fill_manual(
#   values = c(
#    "FRP Forest type" = "#FEB24C",  
#    "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
#    "CBI Forest type" = "#800026",  
#    "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
#   ),
#   labels = c(
#    "Cumulative FRP", "VPD-mediated (FRPc)",
#    "90th percentile CBIbc", "VPD-mediated (CBIbc)"
#   ),
#   guide = guide_legend(
#    override.aes = list(
#     fill = c("#FEB24C", scales::alpha("#FEB24C", 0.3), 
#              "#800026", scales::alpha("#800026", 0.3)),
#     alpha = c(1, 0.3, 1, 0.3),  # Match transparency
#     color = c(
#      scales::alpha("black", 0.7),
#      scales::alpha("#FEB24C", 0.3), 
#      scales::alpha("black", 0.7), 
#      scales::alpha("#800026", 0.3))
#    ),
#    order = 1  # Ensure this legend appears first
#   )
#  ) +
#  scale_alpha_manual(
#   values = c(
#    "FRP Forest type" = 1,
#    "FRP VPD-mediated" = 0.3,
#    "CBI Forest type" = 1,
#    "CBI VPD-mediated" = 0.3
#   ),
#   guide = "none"  
#  ) +
#  scale_color_manual(
#   values = c(
#    "FRP Forest type" = scales::alpha("black", 0.6),
#    "FRP VPD-mediated" = scales::alpha("#FEB24C", 0.3),
#    "CBI Forest type" = scales::alpha("black", 0.6),  
#    "CBI VPD-mediated" = scales::alpha("#800026", 0.3)
#   ),
#   guide = "none"
#  ) +
#  coord_cartesian(xlim=c(-0.38,0.44)) +
#  theme_classic() +
#  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
#        axis.text.x = element_text(angle = 0, hjust = 0, size=9),
#        axis.title.y = element_text(size = 10, margin = margin(r = 12)),
#        axis.title.x = element_text(size = 10, margin = margin(t = 12)),
#        legend.position = c(0.18, 0.26),
#        legend.background = element_rect(
#         fill = scales::alpha("white", 0.4), 
#         color = NA, size = 0.8),
#        legend.title = element_text(size = 9),
#        legend.text = element_text(size = 8))
# p4.1
# # Save the plot
# out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FORTYPCD_FRP-CBI_species.png')
# ggsave(out_png, dpi = 500, width = 8, height = 4, bg = 'white')


#=============COMBINED PLOTS=============#

###################################
# Species composition / structure #
p5 <- frp.p2 / cbi.p2 # "/" stacks them, "|" would place them side-by-side
p5
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_species.png')
ggsave(out_png, plot = p5, dpi = 500, width = 9, height = 7, bg = 'white')

# # version 2
# p5.1 <- frp.p2.1 / cbi.p2.1 # "/" stacks them, "|" would place them side-by-side
# p5.1
# # Save the plot
# out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_species_v2.png')
# ggsave(out_png, plot = p5.1, dpi = 500, width = 8, height = 6, bg = 'white')

#######################
# All other effects ...
# Create the combined forest type relative to aspen
fixed_effects.frp$response <- "FRP"
fixed_effects.cbi$response <- "CBI"
fixed_effects <- bind_rows(fixed_effects.frp, fixed_effects.cbi)

# plot it
(p6 <- fixed_effects %>%
  filter(effect != "Re-burn") %>%
  ggplot(., aes(x = x, y = effect, height = y, fill = response)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect size",
   y = "Parameter",
   fill = "Response",
  ) +
  scale_fill_manual(values = c("FRP" = "#FEB24C", "CBI" = "#800026"),
                    labels = c(
                     "FRP" = "Cumulative FRP",
                     "CBI" = expression("90"^"th" ~ "Percentile CBI"))) +
  # coord_cartesian(xlim=c(-0.16,0.20)) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
        axis.text.x = element_text(angle = 0, hjust = 0, size=9),
        axis.title.y = element_text(size = 10, margin = margin(r = 12)),
        axis.title.x = element_text(size = 10, margin = margin(t = 12)),
        legend.position = c(0.82, 0.88),
        legend.background = element_rect(
         fill = scales::alpha("white", 0.4), 
         color = NA, size = 0.8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))) 
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_fixedeffects.png')
ggsave(out_png, plot = p6, dpi = 500, width = 8, height = 6, bg = 'white')

gc()
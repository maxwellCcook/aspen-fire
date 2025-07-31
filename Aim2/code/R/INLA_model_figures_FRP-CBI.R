
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
fixed_effects <- combined.df %>%
 filter(
  !effect %in% c("overlap"),
  !str_detect(parameter, "fortypnm"),
  !str_detect(parameter, "species")
 ) 

# reorder the elements
# order the forest types from low->high on the main effect
fixed_effects_order <- fixed_effects %>%
 filter(response == "FRPc") %>%
 group_by(effect) %>%
 summarize(mean_effect = mean(x, na.rm = TRUE)) %>%
 arrange(mean_effect)

# apply the ordering
fixed_effects$effect <- factor(
 fixed_effects$effect, 
 levels = fixed_effects_order$effect
)

# make the plot
(p3 <- fixed_effects %>%
  ggplot(., aes(x = x, y = effect, height = y, fill = response)) +
  geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
   x = "Effect size",
   y = "Parameter",
   fill = "Response",
  ) +
  scale_fill_manual(values = c("FRPc" = "#FEB24C", "CBIbc" = "#800026"),
                    labels = c(
                     "FRPc" = "Cumulative FRP",
                     "CBIbc" = expression("90"^"th" ~ "Percentile CBI"))) +
  theme_classic() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size=9),
        axis.text.x = element_text(angle = 0, hjust = 0, size=9),
        axis.title.y = element_text(size = 10, margin = margin(r = 12)),
        axis.title.x = element_text(size = 10, margin = margin(t = 12)),
        legend.position = c(0.16, 0.88),
        legend.background = element_rect(
         fill = scales::alpha("white", 0.4), 
         color = NA, size = 0.8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8))) 

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_fixedeffects.png')
ggsave(out_png, plot = p3, dpi = 500, width = 8, height = 6, bg = 'white')

gc()

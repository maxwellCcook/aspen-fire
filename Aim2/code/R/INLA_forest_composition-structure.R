
# libraries
library(tidyverse)
library(sf) 
library(INLA) 
library(ggcorrplot)
library(ggridges)
library(reshape2)
library(spdep)
library(patchwork)
library(forcats)
library(gstat)


#=========Load gridcell data=========#

# see INLA_gridcell-prep.R
fp <- "data/tabular/mod/model_data_cleaned.csv"
grid_tm <- read_csv(fp) %>%
 # filter for response variables
 filter(
  log_frp_csum > 0,
  CBIbc_p90 > 0
 ) %>%
 # remove NA on BALIVE
 filter(
  !is.na(BALIVE), 
  !is.infinite(BALIVE)
 ) %>%
 # create a forest proportion
 group_by(grid_idx) %>%
 mutate(
  forest_prop = sum(proportion, na.rm=T),
  balive_total = sum(BALIVE, na.rm=T)) %>%
 ungroup() %>%
 mutate(
  balive_prop = BALIVE / balive_total,
  dom_fortyp = dom_fortyp %>%
   tolower() %>%
   gsub("–", "_", .) %>%  # en-dash to underscore
   gsub("-", "_", .) %>%
   gsub(" ", "_", .),
  fortypnm_gp = fortypnm_gp %>%
   tolower() %>%
   gsub("–", "_", .) %>%  # en-dash to underscore
   gsub("-", "_", .) %>%
   gsub(" ", "_", .)
 ) %>%
 # select the model attributes ...
 select(
  fire_id, grid_idx, # ID fields
  log_frp_csum, CBIbc_p90, # response variables
  fortypnm_gp, dom_fortyp, proportion, 
  forest_prop, balive_prop, balive_total, # forest type information
  first_obs_date, day_max_frp, # temporal attributes
  BALIVE, SDI, QMD, # TreeMap attributes
  CC, CH, CBD, CBH, # Fuel attributes
  erc, erc_dv, vpd, vpd_dv, vs, # fire weather
  elev, slope, northness, eastness, tpi, # topography
  day_prop, overlap, # VIIRS aggregation
  dist_to_perim, log_fire_size, # gridcell position and fire size
  fire_aspen, grid_aspen, # aspen presence
  x, y # gridcell center coordinates
 ) %>%
 distinct(grid_idx, fortypnm_gp, .keep_all=T)

glimpse(grid_tm)

print(length(unique(grid_tm$grid_idx))) # unique gridcells

# Check NAs
grid_tm %>%
 summarize(across(
  .cols = c(BALIVE, SDI, QMD, CC, vs, balive_prop),
  .fns = list(
   na = ~sum(is.na(.)),
   inf = ~sum(is.infinite(.))
  ),
  .names = "{.col}_{.fn}"
 ))

# check on BALIVE
grid_tm %>%
 summarize(
  min = min(BALIVE, na.rm = TRUE),
  q10 = quantile(BALIVE, 0.10, na.rm = TRUE),
  q50 = quantile(BALIVE, 0.50, na.rm = TRUE),
  q90 = quantile(BALIVE, 0.90, na.rm = TRUE),
  q99 = quantile(BALIVE, 0.99, na.rm = TRUE),
  max = max(BALIVE, na.rm = TRUE),
  n_extreme = sum(BALIVE > 1e4, na.rm = TRUE),
  total = sum(!is.na(BALIVE))
 )

# filter extreme values of BALIVE
# remove SDI and QMD
grid_tm <- grid_tm %>%
 filter(BALIVE < 1e30) %>%
 select(-c(SDI, QMD))
 
# check the distribution
ggplot(grid_tm, aes(x = BALIVE)) +
 geom_histogram(bins = 100) +
 scale_x_continuous(trans = "log10") +
 theme_minimal() +
 labs(title = "Distribution of BALIVE", x = "BALIVE (log10)", y = "Count")

# get the BALIVE of the dominant forest type
balive_dom <- grid_tm %>%
 filter(fortypnm_gp == dom_fortyp) %>%
 select(grid_idx, balive_dom = BALIVE, balive_dom_pr = balive_prop)

# join back
grid_tm <- grid_tm %>%
 left_join(balive_dom, by = "grid_idx") %>%
 filter(!is.na(balive_dom))

rm(balive_dom)

###################################################
# check how many grids are majority forested, etc #
# percent of gridcells majority forested
# Check how many grids are > 50% forested
print(paste0(
 "Percent forested gridcells (>50% gridcell area): ",
 round(dim(grid_tm %>% filter(forest_prop >= 0.50) %>% distinct(grid_idx))[1]/
        dim(grid_tm %>% distinct(grid_idx))[1], 3) * 100, "%"
))

# # retain predominantly forested gridcells ...
# grid_tm <- grid_tm %>%
#  filter(forest_prop >= 0.50)
# # Check how many gridcells and fires remain
# print(length(unique(grid_tm$fire_id))) # unique fires
# print(length(unique(grid_tm$grid_idx))) # unique gridcells


##################################################
# Assess the proportion of live basal area and TPP
# filter potential noise in the data ...
(qt.ba <- tibble::enframe(
 quantile(grid_tm$balive_prop, probs = seq(0.1, 0.9, by = .1)),
 name = "qt", value = "val"
))
(qt10.ba = qt.ba[1,]$val) # 10%

# plot the distribution and 10th percentile line
ggplot(grid_tm, aes(x = balive_prop)) +
 geom_histogram(
  bins = 30, 
  fill = "grey", 
  color="grey40", 
  alpha = 0.4
 ) +
 geom_vline(
  xintercept = 0.05, color = "darkred", linewidth=1, linetype = "dashed"
 ) +
 labs(x = "Proportion live basal area", y = "Count") +
 theme_classic()

# check percentage of rows below the threshold
1 - dim(grid_tm %>% filter(
 balive_prop >= qt10.ba
))[1] / dim(grid_tm)[1]

# 5% of BA live
dim(grid_tm %>% filter(balive_prop >= 0.05))[1] / dim(grid_tm)[1]

# filter species contributing less than the Nth percentile
grid_tm <- grid_tm %>%
 # filter noise in species data
 filter(
  # remove noise from small species proportions
  # e.g., less than 5% of live basal area
  balive_prop >= 0.05
 )

# Tidy up !
rm(qt.ba, qt10.ba)


##########################################
# filter fires with not enough gridcells
# should have at least 10 (?) for spatial model

# check on the grid cell counts
gridcell_counts <- grid_tm %>%
 distinct(fire_id, grid_idx) %>% # keep only distinct rows
 group_by(fire_id) %>%
 summarise(n = n())
# check the distribution
summary(gridcell_counts$n)

# calculate the quantile distribution
(qt <- tibble::enframe(
 round(quantile(gridcell_counts$n, probs = seq(.1, .9, by = .1))),
 name = "qt", value = "val"
))
(qt10 = qt[1,]$val) # 10th percentile

# filter fires below the 10th percentile
fires_keep <- gridcell_counts %>%
 filter(n >= qt10) %>%
 pull(fire_id)

# filter to retain fires with enough gridcells
grid_tm <- grid_tm %>%
 filter(fire_id %in% fires_keep)

# Check how many gridcells and fires remain
print(length(unique(grid_tm$fire_id))) # unique fires
print(length(unique(grid_tm$grid_idx))) # unique gridcells

# tidy up!
rm(gridcell_counts, fires_keep, qt, qt10)

###############################
# check on the species counts #
grid_tm %>% 
 group_by(fortypnm_gp) %>%
 summarise(n = length(fortypnm_gp)) %>%
 arrange(desc(n))

grid_tm %>%
 distinct(grid_idx, .keep_all = T) %>%
 count(dom_fortyp, sort = TRUE)


#===========MODEL SETUP==============#

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # set factors
  fire_id = as.factor(fire_id),
  grid_idx = as.factor(grid_idx),
  first_obs_date = as.factor(first_obs_date),
  day_max_frp = as.factor(day_max_frp),
  grid_aspen = as.factor(grid_aspen),
  fire_aspen = as.factor(fire_aspen),
  # force aspen to be the baseline factor #
  fortypnm_gp = fct_relevel(fortypnm_gp, "aspen"),
  dom_fortyp = fct_relevel(dom_fortyp, "aspen"),
  # center/scale metrics (fixed effects)
  across(
   # only scale numeric columns not in exclude list
   .cols = where(is.numeric) & !all_of(c("log_frp_csum","CBIbc_p90",
                                         "log_fire_size","x","y")),  
   .fns = ~ as.numeric(scale(.))
  )
  ) %>%
 # let's rename the outcome vars
 rename(
  frpc = log_frp_csum,
  cbibc90 = CBIbc_p90
 ) %>%
 arrange(grid_idx) # arrange by grid index

# check the factor levels
# make sure aspen is first
levels(da$dom_fortyp)
levels(da$fortypnm_gp)


#########################################################
# save the mean and standard deviation for back-transform
(sc <- da %>%
 select(where(is.numeric)) %>%
 select(-c(frpc, cbibc90, log_fire_size, x, y)) %>%
 names())

# Create scaling lookup table (mean and SD from unscaled data)
sc_lookup <- grid_tm %>%
 summarise(across(all_of(sc), list(
  mean = ~mean(.x, na.rm = TRUE), 
  sd = ~sd(.x, na.rm = TRUE))
 ))

# Save as CSV
write_csv(sc_lookup, "data/tabular/mod/var_scaling_lookup.csv")

# tidy up
rm(grid_tm, sc_lookup)
gc()


#===========CORRELATIONS==============#

########################################
# correlation matrix for fixed effects #
cor_da <- da %>%
 select(where(is.numeric)) %>%
 select(-c("x","y"))
# Compute correlation matrix
cor_mat <- cor(cor_da, method = "spearman")
# Plot correlation matrix
ggcorrplot(cor_mat, method = "circle",
           type = "lower", lab = TRUE, lab_size = 3, tl.cex = 10,
           colors = c("blue", "white", "red")
)

# save the plot.
out_png <- 'figures/INLA_CorrelationMatrix_FixedEffects_Model.png'
ggsave(out_png, dpi=500, width=12, height=12, bg = 'white')

rm(cor_da, cor_mat) # tidy up


#######################################
# species-specific correlation matrix #
sp_cor <- da %>%
 select(grid_idx, fortypnm_gp, BALIVE) %>%
 group_by(grid_idx, fortypnm_gp) %>%
 summarize(BALIVE = sum(BALIVE, na.rm = TRUE), .groups = "drop") %>%
 pivot_wider(
  names_from = fortypnm_gp,
  values_from = BALIVE,
  values_fill = 0
 )
# compute the correlation matrix
sp_cormat <- cor(sp_cor[,-1], use = "complete.obs", method = "spearman")
ggcorrplot(sp_cormat, method = "circle", type = "lower",
           lab = TRUE, lab_size = 3, colors = c("blue", "white", "red"))
# save the plot.
out_png <- 'figures/INLA_CorrelationMatrix_SpeciesTPP_sc.png'
ggsave(out_png, dpi=500, bg = 'white')

rm(sp_cor,sp_cormat) # tidy up


###################################################
# Save a fire perimeter dataset after the filtering
fp <- "data/spatial/mod/srm_fire_census_2017_to_2023_ics_perims.gpkg"
fires <- st_read(fp)
# subset and save out
fires <- fires %>% 
 filter(Fire_ID %in% unique(da$fire_id)) %>%
 st_as_sf() %>%
 st_transform(st_crs(5070))
st_write(fires, "data/spatial/mod/srm_fire_census_model_data.gpkg",
         append=F)
# tidy up
rm(fires)


###########################################
# create the summary table for sample sizes
# of the dominant forest type by grid cell
dom_occurrence <- da %>%
 distinct(grid_idx, .keep_all = T) %>%
 group_by(dom_fortyp) %>%
 summarise(
  n_cells = n(),
  n_fires = n_distinct(fire_id),
  pct_fires = n_fires / length(unique(da$fire_id)),
  .groups = "drop"
 ) %>%
 arrange(desc(n_cells))  
print(dom_occurrence)

# do this for all occurrences
fortyp_occurrence <- da %>%
 group_by(fortypnm_gp) %>%
 summarise(
  n_occurrences = n(),
  n_gridcells = n_distinct(grid_idx),
  pct_gridcells = n_gridcells / length(unique(da$grid_idx)),
  .groups = "drop"
 ) %>%
 arrange(desc(n_occurrences))
print(fortyp_occurrence)

# combine them 
sample_sizes <- full_join(
 dom_occurrence,
 fortyp_occurrence,
 by = c("dom_fortyp" = "fortypnm_gp")
) %>%
 arrange(desc(n_occurrences))

# Save this file out.
write_csv(sample_sizes, paste0("data/tabular/mod/results/model_sample_size_fortyp.csv"))

rm(sample_sizes, dom_occurrence, fortyp_occurrence)


#========CREATE THE SPATIAL MESH GRID=========#

# extract gridcell coordinates
grid_sf <- da %>%
 arrange(grid_idx) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
# extract coordinates
coords <- grid_sf %>% st_coordinates(.)
# convert coordinates to a matrix for INLA
coords_mat <- as.matrix(coords)

# Define the spatial mesh
mesh <- inla.mesh.2d(
 loc = coords_mat, # Locations (grid centroids)
 max.edge = c(1, 20), # Maximum edge lengths (inner and outer)
 cutoff = 0.01, # Minimum distance between points (0.01deg = ~1.1km)
 offset = c(0.5, 0.1) # Boundary buffer
)
# # Plot mesh to check
# plot(mesh, main = "SPDE Mesh for FRP Model")
# points(coords, col = "red", pch = 20)
rm(coords, grid_sf)

# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.3, 0.5), # 50% certainty that range is below ~5km
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.05) # variance
)

# Compute the projector matrix (A)
A <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

# Assign the spatial index
field.idx <- inla.spde.make.index(
 name = "mesh.idx",
 n.spde = spde.ml$n.spde
)
str(field.idx)



#===========MODEL FITTING==============#
set.seed(456)

####################
# define some priors
# fixed effects
pr.fixed <- list(
 prec.intercept = 1e-6, 
 prec = 0.001  # shrinkage prior (very weak)
)
# for the Gaussian precision
pr.family <- list(
 hyper = list(
  prec = list(
   prior = "pc.prec", param = c(1, 0.01))
  )
)
# pc.prec for random effects, etc.
pr.pc_prec <- list(
 prec = list(
  prior = "pc.prec", param = c(1, 0.05))
)


#========Cumulative Fire Radiative Power (FRPc)========#


# Create the INLA data stack for FRPc
X <- da %>% select(-c(frpc,cbibc90))
stack.frp <- inla.stack(
 data = list(frpc = da$frpc),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)

rm(X)

dim(inla.stack.A(stack.frp))


#################################################################
# Model with spatial process (SPDE) and fire-level random effects

# setup the model formula
mf.frp <- frpc ~ 1 +
 # tree-level structure X dominance
 fortypnm_gp:proportion + # proportion of gridcell
 fortypnm_gp:BALIVE + # species live basal area
 fortypnm_gp:CC + # species canopy cover percent
 fortypnm_gp:CBH + # species canopy base height
 # gridcell totals
 # other fixed effects
 forest_prop + # gridcell forest proportion
 erc_dv + vpd + vs + # day-of-burn fire weather
 elev + slope + northness + eastness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop +  # gridcell proportion daytime detections
 dist_to_perim + # gridcell distance to perimeter
 # random/latent effects
 f(fire_id, model = "iid", hyper = pr.pc_prec) + # fire-level random effects
 f(mesh.idx, model = spde.ml) # spatial process model

# fit the model
ml.frp <- inla(
 mf.frp, 
 data = inla.stack.data(stack.frp),
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute=T),
 control.compute = list(
  dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE,
  return.marginals = TRUE
 ),
 control.fixed = pr.fixed, # regularization on fixed effects
 control.family = pr.family # on the Gaussian observation
)

# print model outputs
summary(ml.frp) # model summary table
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)

# extract the hyperparameter summary
(hyperpar.frp <- ml.frp$summary.hyperpar)
write_csv(hyperpar.frp,paste0("data/tabular/mod/results/frp_hyperpar.csv"))

# Save the model.
saveRDS(ml.frp, file = "code/R/models/ml_frp_re_sp.rds")
saveRDS(da, "code/R/models/frp_model_da.rds")

gc()


#==========EXTRACTING THE SPATIAL EFFECT===========#

# Compute spatial effect for each grid location
spat.eff <- A %*% ml.frp$summary.random$mesh.idx$mean  # FRP spatial field
# Compute 2.5% (lower bound) credible interval
lower <- A %*% ml.frp$summary.random$mesh.idx$`0.025quant`
# Compute 97.5% (upper bound) credible interval
upper <- A %*% ml.frp$summary.random$mesh.idx$`0.975quant`
# Compute uncertainty as the credible interval width (upper - lower)
uncertainty <- upper - lower

# convert to a spatial data frame
spat.eff.frp <- data.frame(
 x = coords_mat[, 1],  # Longitude
 y = coords_mat[, 2],  # Latitude
 spat_effect = as.vector(spat.eff),  # Convert matrix to vector
 lower = as.vector(lower),  # 2.5% quantile
 upper = as.vector(upper),  # 97.5% quantile
 uncertainty = as.vector(uncertainty),  # CI width
 response = "FRPc"
) %>%
 distinct(x, y, .keep_all=T) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.frp)

# write this a spatial file
out_fp = "data/spatial/mod/spatial_effect_FRP_spp.gpkg"
st_write(spat.eff.frp, out_fp, append=F)

# Tidy up!
rm(spat.eff, spat.eff.frp, lower, upper, uncertainty, stack.frp)


#=================MODEL STATEMENTS=================#

# compute the exponentiated fixed effects
exp.frp <- ml.frp$summary.fixed %>%
 rownames_to_column(var = "parameter") %>%
 mutate(
  exp_mean = exp(mean) - 1,  # Convert log(FRP) effect to % difference
  lower_ci = exp(`0.025quant`) - 1,  # 2.5% CI bound
  upper_ci = exp(`0.975quant`) - 1   # 97.5% CI bound
 )

# check results
exp.frp %>% select(parameter, exp_mean, lower_ci, upper_ci)

# save this table
write_csv(exp.frp, "data/tabular/mod/results/INLA_exp_FRP.csv")


#===========POSTERIOR EFFECTS===========#

#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
frp_marginals <- ml.frp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.frp <- tibble::tibble(
 parameter = names(frp_marginals),
 data = purrr::map(frp_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 # Exclude the intercept
 filter(!parameter %in% c(
  "(Intercept)","overlap","day_prop"
  )
 ) %>%  
 mutate(
  effect = case_when(
   str_detect(parameter, "BALIVE") ~ "Live basal area",
   str_detect(parameter, "balive_prop") ~ "Proportion of live basal area",
   str_detect(parameter, "proportion") ~ "Forest type proportion",
   str_detect(parameter, "lf_canopy") ~ " Canopy cover percent (average)",
   str_detect(parameter, "lf_height") ~ " Canopy height (average)",
   str_detect(parameter, "BALIVE_total") ~ "Live basal area (total)",
   str_detect(parameter, "ba_dead_total") ~ "Dead basal area (total)",
   str_detect(parameter, "tpp_live_total") ~ "Live trees/acre (total)",
   str_detect(parameter, "vs") ~ "Wind speed",
   str_detect(parameter, "elev") ~ "Elevation",
   str_detect(parameter, "northness") ~ "Northness",
   str_detect(parameter, "eastness") ~ "Eastness",
   str_detect(parameter, "slope") ~ "Slope",
   str_detect(parameter, "tpi") ~ "Topographic position",
   str_detect(parameter, "H_tpp") ~ "Shannon diversity index (H-TPP)",
   str_detect(parameter, "erc_dv") ~ "Energy release component\n(15-year deviation)",
   str_detect(parameter, "vpd") & !str_detect(parameter, ":") ~ "Vapor pressure deficit",
   str_detect(parameter, "dist_to_perim") ~ "Distance to fire edge",
   TRUE ~ parameter  # Default for all other fixed effects
  ),
  # Extract species names from parameters
  species = case_when(
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 

# Order effects: Species-specific effects first, 
# global effects by mean effect size
tidy.effects.frp <- tidy.effects.frp %>%
 mutate(
  species = case_when(
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   parameter == "fortyp_pct" ~ "quaking_aspen",  # baseline
   TRUE ~ NA_character_
  )
 ) %>%
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.frp %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global effect", species)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole_pine" = "Lodgepole pine",
  "douglas_fir" = "Douglas-fir",
  "white_fir" = "White fir",
  "gambel_oak" = "Gambel oak",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa_pine" = "Ponderosa pine",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(
   fill_species,
   levels = c("Lodgepole pine", "Douglas-fir", "White fir",
              "Gambel oak", "Piñon-juniper", "Ponderosa pine",
              "Spruce-fir", "Quaking aspen"))) %>%
 mutate(exp_effect = exp(x))

# check on the species name extraction
unique(tidy.effects.frp$fill_species)
spps_breaks <- unique(tidy.effects.frp$fill_species)

# Save the tidy model effects to a CSV for plotting.
write_csv(tidy.effects.frp, paste0(maindir, "data/tabular/mod/results/INLA_tidy_effects_FRP.csv"))



#========90th Percentile Composite Burn Index (CBIbc90)========#

# # prep the data frame
# da <- da %>%
#  # either add a constant or remove zeros
#  # gamma requires strictly positive
#  # CBI = 0 is likely "unburned" but had some FRP detection (?)
#  # filter zeros
#  # filter(cbibc90 > 0)
#  # OR, add a very small constant
#  mutate(cbibc90 = cbibc90 + 1e-4)

# Create the INLA data stack

X <- da %>% select(-c(frpc,cbibc90)) # isolate vars

stack.cbi <- inla.stack(
 data = list(cbibc90 = da$cbibc90),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)

rm(X)

dim(inla.stack.A(stack.cbi))

#######################################
# 1. Baseline model (no latent effects)

# setup the model formula
mf.cbi <- cbibc90 ~ 1 + 
 # forest composition and structure
 fortypnm_gp + # majority forest type (factor, aspen baseline)
 fortypnm_gp:fortyp_pct + # by proportional area
 # tree-level structure X dominance
 fortypnm_gp:BALIVE + # species live basal area
 fortypnm_gp:qmd_live + # species quadratic mean diameter
 fortypnm_gp:hdr_live + # species height-diameter ratio
 # gridcell totals
 BALIVE_total + tpp_live_total + # total live BA and TPP
 # other fixed effects 
 lf_canopy + # gridcell forest canopy cover percent
 erc_dv + vpd + vs + # day-of-burn fire weather
 elev + slope + northness + eastness + tpi + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop +  # gridcell proportion daytime detections
 dist_to_perim + # gridcell distance to perimeter
 # random/latent effects
 f(fire_id, model = "iid", hyper = pr.pc_prec) + # fire-level random effects
 f(mesh.idx, model = spde.ml) # spatial process model

# fit the model
ml.cbi <- inla(
 mf.cbi, # the model formula
 data = inla.stack.data(stack.cbi),
 family = "gamma",
 control.predictor = list(A = inla.stack.A(stack.cbi), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = pr.fixed, # regularization on fixed effects
 control.family = pr.family # on the gaussian observation
)

# print the model results
summary(ml.cbi)
mean(ml.cbi$cpo$cpo, na.rm = TRUE)
# extract the hyperparameter summary
(hyperpar.cbi <- ml.cbi$summary.hyperpar)
write_csv(hyperpar.cbi,paste0("data/tabular/mod/results/cbi_hyperpar.csv"))

# Save the model.
saveRDS(ml.cbi, file = "code/R/models/ml_cbi_re_sp.rds")
saveRDS(da, "code/R/models/cbi_model_da.rds")


#==========EXTRACTING THE SPATIAL EFFECT===========#

# Compute spatial effect for each grid location
spat.eff <- A %*% ml.cbi$summary.random$mesh.idx$mean  # FRP spatial field
# Compute 2.5% (lower bound) credible interval
lower <- A %*% ml.cbi$summary.random$mesh.idx$`0.025quant`
# Compute 97.5% (upper bound) credible interval
upper <- A %*% ml.cbi$summary.random$mesh.idx$`0.975quant`
# Compute uncertainty as the credible interval width (upper - lower)
uncertainty <- upper - lower

# convert to a spatial data frame
spat.eff.cbi <- data.frame(
 x = coords_mat[, 1],  # Longitude
 y = coords_mat[, 2],  # Latitude
 spat_effect = as.vector(spat.eff),  # Convert matrix to vector
 lower = as.vector(lower),  # 2.5% quantile
 upper = as.vector(upper),  # 97.5% quantile
 uncertainty = as.vector(uncertainty),  # CI width
 response = "CBIbc"
) %>%
 distinct(x, y, .keep_all=T) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.cbi)

# write this a spatial file
out_fp = paste0(maindir,"data/spatial/mod/spatial_effect_CBI_spp.gpkg")
st_write(spat.eff.cbi, out_fp, append=F)

# Tidy up!
rm(spat.eff, spat.eff.cbi, lower, upper, uncertainty,
   A, coords_mat, field.idx, mesh, spde.ml, stack.cbi)


#=================MODEL STATEMENTS=================#

fixed_effects <- ml.cbi$summary.fixed

# Compute percentage change relative to aspen
exp.cbi <- fixed_effects %>%
 rownames_to_column(var = "parameter") %>%
 mutate(
  exp_mean = exp(mean) - 1,  # Convert log(FRP) effect to % difference
  lower_ci = exp(`0.025quant`) - 1,  # 2.5% CI bound
  upper_ci = exp(`0.975quant`) - 1   # 97.5% CI bound
 )
# save this table
write_csv(exp.cbi, paste0(maindir,"data/tabular/mod/results/INLA_exp_CBI.csv"))
# check results
exp.cbi%>%select(parameter,exp_mean,lower_ci,upper_ci)


#===========POSTERIOR EFFECTS===========#


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
cbi_marginals <- ml.cbi$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.cbi <- tibble::tibble(
 parameter = names(cbi_marginals),
 data = purrr::map(cbi_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 # Exclude the intercept
 filter(!parameter %in% c(
  "(Intercept)","day_prop","overlap"
 )
 ) %>%  
 mutate(
  effect = case_when(
   str_detect(parameter, "BALIVE") ~ "Live basal area",
   str_detect(parameter, "BALIVE_pr") ~ "Proportion of live basal area",
   str_detect(parameter, "hdr_live") ~ "Diameter/height ratio",
   str_detect(parameter, "sdi_live") ~ "Stand density index",
   str_detect(parameter, "qmd_live") ~ "Quadratic mean diameter",
   str_detect(parameter, "fortyp_pct") ~ "Forest type proportion",
   str_detect(parameter, "lf_canopy") ~ " Canopy cover percent (average)",
   str_detect(parameter, "BALIVE_total") ~ "Live basal area (total)",
   str_detect(parameter, "ba_dead_total") ~ "Dead basal area (total)",
   str_detect(parameter, "tpp_live_total") ~ "Live trees/acre (total)",
   str_detect(parameter, "vs") ~ "Wind speed",
   str_detect(parameter, "elev") ~ "Elevation",
   str_detect(parameter, "northness") ~ "Northness",
   str_detect(parameter, "eastness") ~ "Eastness",
   str_detect(parameter, "slope") ~ "Slope",
   str_detect(parameter, "tpi") ~ "Topographic position",
   str_detect(parameter, "H_tpp") ~ "Shannon diversity index (H-TPP)",
   str_detect(parameter, "erc_dv") ~ "Energy release component\n(15-year deviation)",
   str_detect(parameter, "vpd") & !str_detect(parameter, ":") ~ "Vapor pressure deficit",
   str_detect(parameter, "dist_to_perim") ~ "Distance to fire edge",
   TRUE ~ parameter  # Default for all other fixed effects
  ),
  # Extract species names from parameters
  species = case_when(
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 

# Order effects: Species-specific effects first, 
# global effects by mean effect size
tidy.effects.cbi <- tidy.effects.cbi %>%
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.cbi %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global effect", species)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole_pine" = "Lodgepole pine",
  "douglas_fir" = "Douglas-fir",
  "white_fir" = "White fir",
  "gambel_oak" = "Gambel oak",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa_pine" = "Ponderosa pine",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(
   fill_species,
   levels = c("Lodgepole pine", "Douglas-fir", "White fir",
              "Gambel oak", "Piñon-juniper", "Ponderosa pine",
              "Spruce-fir", "Quaking aspen")))

# check on the species name extraction
unique(tidy.effects.cbi$fill_species)
spps_breaks <- unique(tidy.effects.cbi$fill_species)

# Save the tidy model effects to a CSV for plotting.
write_csv(tidy.effects.cbi, paste0(maindir, "data/tabular/mod/results/INLA_tidy_effects_CBI.csv"))



#====================RESULTS TABLES=======================#

alpha <- ml.cbi$marginals.fixed[[1]]
ggplot(data.frame(inla.smarginal(alpha)), aes(x, y)) +
 geom_line() +
 theme_bw()

# set the results directory
results_dir = paste0(maindir,"data/tabular/mod/results/")

# exponentiated effects
exp.cbi$response <- "CBIbc"
exp.frp$response <- "FRPc"
(exp.fixed <- bind_rows(exp.frp, exp.cbi))
write_csv(exp.fixed, paste0(results_dir,"fixed_effects_FRP-CBI.csv"))


######################################################
# Extract fixed effects summary for both FRP and CBI #

# FRPc
summary.frp <- as.data.frame(ml.frp$summary.fixed) %>%
 rownames_to_column("Parameter") %>%
 rename(
  Mean = mean,
  SD = sd,
  `2.5% CI` = `0.025quant`,
  `50% CI` = `0.5quant`,
  `97.5% CI` = `0.975quant`
 ) %>%
 mutate(Response = "FRP")

# CBIbc p90
summary.cbi <- as.data.frame(ml.cbi$summary.fixed) %>%
 rownames_to_column("Parameter") %>%
 rename(
  Mean = mean,
  SD = sd,
  `2.5% CI` = `0.025quant`,
  `50% CI` = `0.5quant`,
  `97.5% CI` = `0.975quant`
 ) %>%
 mutate(Response = "CBI")

# combine them
ml_summary.fixed <- bind_rows(summary.frp, summary.cbi)
head(ml_summary.fixed)

# save to CSV
write_csv(ml_summary.fixed, paste0(results_dir,"model_summaries_FRP-CBI.csv"))

gc()
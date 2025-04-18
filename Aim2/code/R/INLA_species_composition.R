
# Load the required libraries
library(tidyverse)
library(sf) # spatial
library(INLA) # for spatial Bayes model
library(ggcorrplot)
library(lubridate)
library(ggridges)
library(reshape2)
library(spdep)
library(patchwork)
library(forcats)

# Environment variables
maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load the fire data
fires <- paste0(maindir,"data/spatial/mod/srm_fire_census_2017_to_2023_ics_perims.gpkg")
fires <- st_read(fires) %>%
 as_tibble() %>%
 rename(
  fire_ig_dt = DISCOVERY_DATE,
  fire_acres = ICS_ACRES,
  fire_name = Fire_Name,
  fire_aspenpct = Aspen_Pct
 ) %>%
 select(Fire_ID, fire_name, fire_ig_dt, fire_acres, fire_aspenpct) %>%
 mutate(
  Fire_ID = as.numeric(Fire_ID),
  fire_acres = as.numeric(fire_acres),
  fire_aspenpct = as.numeric(fire_aspenpct))
head(fires)


#=========Prep the grid data=========#

# Format the species composition data frame 
# (long-format) each grid has rows for species co-occurring with the dominant type
# climate and topography are summarized at the grid level

# load the aggregated FRP grid with TreeMap and climate/topography
fp <- paste0(maindir,'data/tabular/mod/gridstats_fortypnm_gp_tm_ct_frp-cbi.csv')
grid_tm <-  read_csv(fp)  %>% # read in the file
 select(-'...1') %>%
 # join in some of the fire information
 left_join(fires, by="Fire_ID") %>%
 # remove missing FRP, prep columns
 # Ensure daytime FRP and CBIbc is greater than 0
 filter(
  frp_csum > 0, # make sure cumulative FRP > 0
  CBIbc_p90 > 0 # make sure CBI > 0
 ) %>% 
 # create a numeric fire ID
 mutate(
  # tidy the temporal fields
  fire_ig_dt = as.Date(fire_ig_dt),  
  fire_year = year(fire_ig_dt),              
  fire_month = month(fire_ig_dt),
  first_obs_date = as.Date(first_obs_date), # the day of first observation
  first_obs_date = as.factor(first_obs_date),
  first_obs_doy = as.numeric(lubridate::yday(first_obs_date)),
  # Create a "season" variable based on the month
  season = case_when(
   fire_month %in% c(3, 4, 5) ~ "spring",
   fire_month %in% c(6, 7, 8) ~ "summer",
   fire_month %in% c(9, 10, 11) ~ "fall"
  ),
  # Year/season interaction
  year_season = interaction(fire_year, season, drop = TRUE),
  # Factorize the temporal attributes
  season = as.factor(season),
  fire_year = as.factor(fire_year),
  year_season = as.factor(year_season),
  # Format species names consistently
  fortypnm_gp = str_to_lower(fortypnm_gp),
  fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
  fortypnm_gp = str_replace_all(fortypnm_gp, " ", "_"),
  # Format species names consistently
  species_gp_n = str_to_lower(species_gp_n),
  species_gp_n = str_replace_all(species_gp_n, "-", "_"),
  species_gp_n = str_replace_all(species_gp_n, " ", "_"),
  # make sure factor variables are set for species names
  fortypnm_gp = factor(fortypnm_gp),
  species_gp_n = factor(species_gp_n),
  Fire_ID = factor(Fire_ID),
  # proportion of contributing AFD during daytime observations
  day_prop = day_count / afd_count,
  # log-scaled response variables
  log_frp_max = log(frp_max),
  log_frp_csum = log(frp_csum),
  log_frp_max_day = log(frp_max_day),
  log_frp_max_night = log(frp_max_night),
  log_frp_csum_day = log(frp_csum_day),
  log_frp_csum_night = log(frp_csum_night),
  # scale the percentages
  forest_pct = forest_pct / 100,
  fire_aspenpct = fire_aspenpct / 100
 ) %>%
 # filter noise in species data
 filter(
  # forest_pct > 0.50, # majority forested gridcells
  # remove noise from  small species proportions
  (ba_live_pr > 0.10 | tpp_live_pr > 0.10) # at least 10% of BA or TPP
  # ba_live_pr > 0.10 # at least 10% of live basal area
 ) %>%
 group_by(Fire_ID) %>%
 mutate(aspen = ifelse(any(fortypnm_gp == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates

# load and filter the duplicate (reburn) grids
dup_grids <- st_read(paste0(maindir,"data/spatial/mod/grid_model_data_duplicates.gpkg"))
dup_idx <- unique(dup_grids$grid_index)
grid_tm <- grid_tm %>%
 filter(!grid_index %in% dup_idx)
rm(dup_grids, dup_idx)
gc()

glimpse(grid_tm)

rm(fires)


#####################################
# get the grid-level aspen proportion
aspen_pr <- grid_tm %>%
 group_by(grid_idx) %>%
 reframe(
  # Aspen-specific live BA and TPP
  aspen_ba_live = sum(ba_live[species_gp_n == "quaking_aspen"], na.rm = TRUE),
  aspen_tpp_live = sum(tpp_live[species_gp_n == "quaking_aspen"], na.rm = TRUE),
  # Calculate proportions (avoid division by zero)
  aspen_ba_pr = if_else(ba_live_total > 0, aspen_ba_live / ba_live_total, 0),
  aspen_tpp_pr = if_else(tpp_live_total > 0, aspen_tpp_live / tpp_live_total, 0)
 ) %>%
 select(grid_idx, aspen_ba_pr, aspen_tpp_pr)
head(aspen_pr)


########################################
# merge attributes back to the grid data
grid_tm <- grid_tm %>%
 left_join(aspen_pr, by = "grid_idx") %>%
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm) # check the results

# Check on the grid cell counts for daytime observations
grid_counts <- grid_tm %>%
 distinct(Fire_ID, grid_idx) %>% # keep only distinct rows
 group_by(Fire_ID) %>%
 summarise(n_grids = n())
summary(grid_counts$n_grids)
qt <- quantile(grid_counts$n_grids, probs = seq(.1, .9, by = .1))
qt

# how many fires have >= 10 grids?
dim(grid_counts %>% filter(n_grids >= 10))[1]
# how many fires have >= 50 grids?
dim(grid_counts %>% filter(n_grids >= 50))[1]
# how many fires have >= 100 grids?
dim(grid_counts %>% filter(n_grids >= 100))[1]

# Identify fires with n_grids below the Nth percentile
idx <- grid_counts %>%
 filter(n_grids < 10) %>%
 pull(Fire_ID)
length(idx)

# filter the data frame to remove these fires
grid_tm <- grid_tm %>%
 filter(!Fire_ID %in% idx)
# tidy up!
rm(grid_counts, aspen_pr, idx)
gc()


# get the species counts
grid_tm %>% 
 group_by(species_gp_n) %>%
 summarise(n = length(species_gp_n)) %>%
 arrange(desc(n))

# check how many grids and fires
length(unique(grid_tm$grid_idx))
length(unique(grid_tm$Fire_ID))


#==========EXPLORE THE DATA==========#

######################################
# correlation matrix for fixed effects
# Select only numeric columns
cor_da <- grid_tm %>%
 select(
  tpp_live, tpp_live_pr, ba_live, ba_live_pr,  # species-level proportion live metrics
  aspen_ba_pr, aspen_tpp_pr, forest_pct, # grid-level aspen proportion and forest percent
  tpp_live_total, ba_live_total, qmd_live_mean, # grid-level live biomass
  qmd_live, tree_ht_live, tree_dia_live, # species tree height/diameter
  tree_ht_live_mean, tree_dia_live_mean, # grid-level mean tree height and diameter
  H_ba, H_tpp, # species diversity
  erc, erc_dv, vpd, vpd_dv, # climate
  fm1000, fm1000_dv, rmin, rmin_dv, # climate
  tmmx, tmmx_dv, vs, vs_dv, #climate
  elev, slope, tpi, chili, aspect,  # topography
  forest_pct, canopypct_mean, balive_sum, # grid-level forest characteristics
  lf_forest_cc_mean, lf_forest_ch_mean # landfire forest canopy cover/height
 ) %>%
 mutate_if(is.factor, as.numeric)  # Convert factors to numeric (if needed)

# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")

# Plot correlation matrix
ggcorrplot(
 cor_mat,
 method = "circle",  # Circle or square for visualization
 type = "lower",  # Lower triangle of the correlation matrix
 lab = TRUE,  # Show correlation values
 lab_size = 3,
 tl.cex = 10,  # Text label size
 colors = c("blue", "white", "red")  # Color gradient
)

rm(cor_da, cor_mat)
gc()

# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_SppComp.png')
ggsave(out_png, dpi=500, bg = 'white')



#===========MODEL SETUP==============#

# force aspen to be the baseline
grid_tm <- grid_tm %>%
 mutate(
  species_gp_n = fct_relevel(species_gp_n, "quaking_aspen"),
  fortypnm_gp = fct_relevel(fortypnm_gp, "quaking_aspen")
 )

# check the factor levels
# make sure aspen is first
levels(grid_tm$species_gp_n)
levels(grid_tm$fortypnm_gp)

# prep the model data frame
# center and scale fixed effects
da <- grid_tm %>%
 mutate(
  # create a gridcell-level live/dead proportion
  live_dead_pr_ba = ba_dead_total / (ba_live_total + ba_dead_total),
  live_dead_pr_tpp = tpp_dead_total / (tpp_live_total + tpp_dead_total),
  # center/scale metrics / fixed effects
  across(
   c(ba_live, ba_live_pr, # species-specific live basal area (dominance)
     tpp_live, tpp_live_pr, # species-specific tree/pixel (abundance)
     qmd_live, # species-specific quadratic mean diameter (stand structure)
     ba_live_total, tpp_live_total, qmd_live_mean, # gridcell live biomass
     live_dead_pr_ba, live_dead_pr_tpp, # proportion dead BA/TPP
     tree_ht_live, tree_ht_dead, tree_ht_live_mean,  # species- and gridcell-level tree height
     tree_dia_live, tree_dia_dead, tree_dia_live_mean,  # tree diameter
     H_ba, H_tpp, # gridcell diversity metrics
     aspen_ba_pr, aspen_tpp_pr, # gridcell aspen proportion
     forest_pct, fortyp_pct, canopypct_mean, balive_sum, # gridcell forest characteristics
     lf_forest_cc_mean, lf_forest_ch_mean, # landfire forest canopy cover/height
     erc, erc_dv, vpd, vpd_dv, # climate
     fm1000, fm1000_dv, rmin, rmin_dv, # climate
     tmmx, tmmx_dv, vs, vs_dv, #climate
     elev, slope, tpi, chili, aspect, # topography
     day_prop, overlap, afd_count, # VIIRS detection characteristics
    ), ~ as.numeric(scale(.))
  ),
  log_fire_size = log(fire_acres) # log-scale fire size
 ) %>%
 # rename the response variables
 rename(
  fire_id = Fire_ID,
  tm_canopypct = canopypct_mean,
  tm_balive = balive_sum,
  lf_canopypct = lf_forest_cc_mean,
  lf_canopyht = lf_forest_ch_mean
 ) %>%
 # calculate the predominant species by raw BA and TPP
 group_by(grid_idx) %>%
 mutate(
  dom_spp_ba = as.factor(species_gp_n[which.max(ba_live)]),
  dom_spp_tpp = as.factor(species_gp_n[which.max(tpp_live)]),
  dom_spp_qmd = as.factor(species_gp_n[which.max(qmd_live)]),
 ) %>%
 ungroup() %>%
 # keep just one row per grid cell by predominant type
 distinct(grid_idx, species_gp_n, .keep_all = TRUE) %>%
 arrange(grid_idx)

# tidy up
rm(grid_tm) 
gc()


#####################################################
# get value counts for dominance/abundance by species
da %>%
 group_by(fortypnm_gp) %>%
 distinct(grid_idx) %>%
 summarize(n = n())# fortypcd
da %>%
 group_by(dom_spp_ba) %>%
 distinct(grid_idx) %>%
 summarize(n = n())# dominance
da %>%
 group_by(dom_spp_tpp) %>%
 distinct(grid_idx) %>%
 summarize(n = n()) # abundance

# (optionally) remove rare species
rm.spps <- c("gambel_oak")
da <- da %>%
 filter(!species_gp_n %in% rm.spps,
        !fortypnm_gp %in% rm.spps) %>%
 mutate(species_gp_n = droplevels(species_gp_n),
        fortypnm_gp = droplevels(fortypnm_gp),
        dom_spp_ba = droplevels(dom_spp_ba),
        dom_spp_tpp = droplevels(dom_spp_tpp)) %>%
 # re-calculate the predominant species by raw BA and TPP
 group_by(grid_idx) %>%
 mutate(
  dom_spp_ba = as.factor(species_gp_n[which.max(ba_live)]),
  dom_spp_tpp = as.factor(species_gp_n[which.max(tpp_live)]),
  dom_spp_qmd = as.factor(species_gp_n[which.max(qmd_live)]),
 ) %>%
 ungroup()
# levels(da$species_gp_n)


#####################################################




# check on factor levels
levels(da$dom_spp_ba)
levels(da$dom_spp_tpp)
levels(da$dom_spp_qmd)

# check how many time dominant species match fortypcd
dim(da)[1]
dim(da%>%filter(fortypnm_gp == dom_spp_ba))[1]
dim(da%>%filter(fortypnm_gp == dom_spp_tpp))[1]
dim(da%>%filter(fortypnm_gp == dom_spp_qmd))[1]
dim(da%>%filter(dom_spp_tpp == dom_spp_ba))[1]



#===========MODEL FITTING==============#

set.seed(456)

# make aspen presence a factor
da$aspen <- as.factor(da$aspen)
levels(da$aspen)

########################
# FIRE RADIATIVE POWER #

#################################################
# 1. Baseline model (no random or latent effects)

# setup the model formula
mf.frp <- log_frp_csum ~ 1 + 
 fortypnm_gp:fortyp_pct +
 species_gp_n:ba_live + # species-specific live basal area proportion
 species_gp_n:tree_ht_live + # species-specific live quadratic mean diameter
 species_gp_n:tree_dia_live +
 H_tpp +
 lf_canopypct + # grid-level canopy percent
 erc_dv + vpd + vs + # day-of-burn climate/weather
 elev + slope + aspect + tpi + # grid-level topography 
 overlap + # grid-level VIIRS overlap (sum)
 day_prop + # grid-level proportion of daytime detections 
 aspen # fire-level aspen presence

# fit the model
ml.frp <- inla(
 mf.frp, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp)
# check on predictive power of the random effects models
mean(ml.frp$cpo$cpo, na.rm = TRUE)


#####################################
# 2. Baseline + Fire-ID Random Effect

# update the model formula
mf.frp.re <- update(
 mf.frp, . ~ 1 + . + 
  f(fire_id, model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.5))
    )) 
)

# fit the model
ml.frp.re <- inla(
 mf.frp.re, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp.re)
# check on predictive power of the random effects models
mean(ml.frp.re$cpo$cpo, na.rm = TRUE)


###################################################
# 3. Baseline + Temporal + Fire-level Random Effect

# update the model formula
mf.frp.re2 <- update(
 mf.frp.re, . ~ 1 + . + 
  f(first_obs_date, # Fire-ID + day-of-burn
    model = "iid", 
    hyper = list(
     prec = list(prior = "pc.prec", param = c(1, 0.5))
    )) 
)

# fit the model
ml.frp.re2 <- inla(
 mf.frp.re2, data = da,
 family = "gaussian",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.family = list(
  # relax the variance assumptions
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 # controls for posterior approximation and hyperparameter search
 control.inla = list(strategy = "laplace", int.strategy = "grid")
)
summary(ml.frp.re2)
# check on predictive power of the random effects models
mean(ml.frp.re2$cpo$cpo, na.rm = TRUE)

gc()


#######################
# 4. Spatial SPDE model

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
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
# Plot mesh to check
plot(mesh, main = "SPDE Mesh for FRP Model")
points(coords, col = "red", pch = 20)

rm(coords)
gc()

######################
# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.1, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(0.5, 0.01) # variance
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

# Extract the predictor variables
X <- da %>% 
 select(fire_id, first_obs_date, grid_index, 
        dom_spp_tpp, dom_spp_ba, dom_spp_qmd,
        fortypnm_gp, species_gp_n, forest_pct,  
        tpp_live, tpp_live_total, tpp_live_pr, ba_live, ba_live_total, ba_live_pr,    
        canopypct, H_tpp, H_ba, 
        tree_ht_live, tree_dia_live, qmd_live, 
        live_dead_pr_tpp, live_dead_pr_ba, 
        erc_dv, vpd, vs, slope, tpi, chili, elev,
        aspen, day_prop, afd_count, overlap)
head(X)

# Create the INLA data stack
stack.frp <- inla.stack(
 data = list(log_frp_csum_day = da$log_frp_csum_day),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.frp))

rm(grid_sf, X)
gc()

##########################
# update the model formula
mf.frp.re.sp <- update(
 mf.frp.re2, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.frp.re.sp <- inla(
 mf.frp.re.sp, 
 data = inla.stack.data(stack.frp), ,
 family = "gaussian",
 control.predictor = list(A = inla.stack.A(stack.frp), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = list(
  prec = list(prior = "pc.prec", param = c(1, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.5)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.frp.re.sp)
mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)



#==========EXTRACTING THE SPATIAL EFFECT===========#

# Compute spatial effect for each grid location
# mean
spat.eff <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.frp.re.sp$summary.random$mesh.idx$mean  # FRP spatial field
# credible intervals
# Compute 2.5% (lower bound) credible interval
lower <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.frp.re.sp$summary.random$mesh.idx$`0.025quant`
# Compute 97.5% (upper bound) credible interval
upper <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.frp.re.sp$summary.random$mesh.idx$`0.975quant`
# Compute uncertainty as the credible interval width (upper - lower)
uncertainty <- upper - lower

# convert to a spatial dataframe
spat.eff.frp <- data.frame(
 x = coords_mat[, 1],  # Longitude
 y = coords_mat[, 2],  # Latitude
 spat_effect = as.vector(spat.eff),  # Convert matrix to vector
 lower = as.vector(lower),  # 2.5% quantile
 upper = as.vector(upper),  # 97.5% quantile
 uncertainty = as.vector(uncertainty),  # CI width
 response = "FRP"
) %>%
 distinct(x, y, .keep_all=T) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.frp)
# write this a spatial file
st_write(spat.eff.frp, paste0(maindir,"data/spatial/mod/spatial_effect_FRP_spp.gpkg"), append=F)

# Tidy up!
rm(spat.eff, lower, upper, uncertainty)
gc()



#=================MODEL COMPARISON=================#

# Create a model comparison table
ml_comparison.frp <- tibble(
 Model = c("Baseline", 
           "W/Fire Random Effect", 
           "W/Fire + Temporal Effect",
           "W/Fire + Temporal + Spatial Effect"),
 Response = "FRP",
 DIC = c(
  ml.frp$dic$dic,
  ml.frp.re$dic$dic,
  ml.frp.re2$dic$dic,
  ml.frp.re.sp$dic$dic
 ),
 WAIC = c(
  ml.frp$waic$waic,
  ml.frp.re$waic$waic,
  ml.frp.re2$waic$waic,
  ml.frp.re.sp$waic$waic
 ),
 Marginal_LogLikelihood = c(
  ml.frp$mlik[1],
  ml.frp.re$mlik[1],
  ml.frp.re2$mlik[1],
  ml.frp.re.sp$mlik[1]
 ),
 Effective_Params = c(
  ml.frp$dic$p.eff,
  ml.frp.re$dic$p.eff,
  ml.frp.re2$dic$p.eff,
  ml.frp.re.sp$dic$p.eff
 ),
 Mean_CPO = c(
  mean(ml.frp$cpo$cpo, na.rm = TRUE),
  mean(ml.frp.re$cpo$cpo, na.rm = TRUE),
  mean(ml.frp.re2$cpo$cpo, na.rm = TRUE),
  mean(ml.frp.re.sp$cpo$cpo, na.rm = TRUE)
 )
) %>% arrange(DIC)
# Print the comparison table
print(ml_comparison.frp)


# tidy up !
rm(A, coords_mat, field.idx, mesh, 
   ml.frp, ml.frp.re, ml.frp.re2, spde.ml, stack.frp)
gc()



#===========POSTERIOR EFFECTS===========#


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
frp_marginals <- ml.frp.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.frp <- tibble::tibble(
 parameter = names(frp_marginals),
 data = purrr::map(frp_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 filter(!parameter %in% c("(Intercept)","aspen","day_prop")) %>%  # Exclude the intercept
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live_pr") ~ "Proportion of \nLive Basal Area",
   str_detect(parameter, "qmd_live") ~ "Quadratic\nMean Diameter",
   str_detect(parameter, "H_tpp") ~ "Shannon Index\n(Abundance-based)",
   # str_detect(parameter, "tree_ht_live") ~ "Tree height (m)",
   TRUE ~ parameter  # Default for all other fixed effects
  ),
  # Extract species names from parameters
  species = case_when(
   str_detect(parameter, "species_gp_n") ~ str_extract(parameter, "(?<=species_gp_n)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 

tidy.effects.frp <- tidy.effects.frp %>%
 # Order effects: Species-specific effects first, global effects by mean effect size
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.frp %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global Effect", species)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole" = "Lodgepole pine",
  "douglas_fir" = "Douglas-fir",
  "white_fir" = "White fir",
  "gambel_oak" = "Gambel oak",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa" = "Ponderosa pine",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(fill_species, 
                        levels = c("Lodgepole pine", "Douglas-fir", "White fir",
                                   "Gambel oak", "Piñon-juniper", "Ponderosa pine", 
                                   "Spruce-fir", "Quaking aspen")))

# check on the species name extraction
unique(tidy.effects.frp$fill_species)
spps_breaks <- unique(tidy.effects.frp$fill_species)


#######################
# create the ridge plot
frp.p1 <- tidy.effects.frp %>%
 filter(!parameter %in% c("(Intercept)", "aspen1", "overlap", "day_prop", "afd_count"),
        !str_starts(parameter, "fortyp")) %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = fill_species)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 scale_fill_manual(
  values = c(
   "Global Effect" = "gray",  # Neutral color for global effects
   "quaking_aspen" = "#e6ab02",
   "lodgepole_pine" = "#d95f02",
   "douglas_fir" = "#7570b3",
   "white_fir" = "#a6cee3",
   "piñon_juniper" = "#e7298a",
   "ponderosa_pine" = "#66a61e",
   "spruce_fir" = "#1b9e77"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )
frp.p1

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_full.png')
ggsave(out_png, dpi = 500, bg = 'white')


####################################
# plot with just the species metrics
frp.p2 <- tidy.effects.frp %>%
 filter(fill_species != "Global Effect") %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = fill_species, alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5, show.legend=F) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # add a subplot label
 annotate("text", x = -0.13, y = 5, 
          label = "(A)", fontface = "bold", 
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = c(
   "Quaking aspen" = "#1b9e77",
   "Lodgepole" = "#d95f02",
   "Mixed-conifer" = "#7570b3",
   "Piñon-juniper" = "#e7298a",
   "Ponderosa" = "#66a61e",
   "Spruce-fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 scale_alpha_manual(
  values = c(
   "Quaking aspen" = 0.9, 
   "Lodgepole" = 0.7,  
   "Mixed-conifer" = 0.7,
   "Piñon-juniper" = 0.7,
   "Ponderosa" = 0.7,
   "Spruce-fir" = 0.7
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 xlim(-0.13, 0.15) +
 theme_minimal() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 8)),
  axis.title.x = element_text(size = 10, margin = margin(t = 8)),
  legend.position = c(0.85, 0.70),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.6), 
   color = NA, size = 1),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8)
 )
frp.p2

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP_species.png')
ggsave(out_png, dpi = 500, bg = 'white')
# 
# # Tidy up!
# rm(frp_marginals, tidy.effects)
# gc()




########################
# COMPOSITE BURN INDEX #


#######################################
# 1. Baseline model (no latent effects)

# setup the model formula
mf.cbi <- CBIbc_p90 ~ 1 + 
 fortypnm_gp:H_tpp + # species diversity (abundance) by predominant forest type
 species_gp_n:ba_live_pr + # species proportion of live basal area
 species_gp_n:qmd_live + # species-level live quadratic mean diameter
 species_gp_n:tree_ht_live + # species-level live tree height
 # live_dead_pr_tpp + # grid-level proportion dead TPP
 canopypct + # grid-level canopy percent
 erc_dv + vpd + vs + # day-of-burn climate/weather
 elev + slope + tpi + chili + # grid-level topography 
 aspen # fire-level aspen presence
# fit the model
ml.cbi <- inla(
 mf.cbi, data = da,
 family = "gamma",
 control.predictor = list(compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi)
# check on predictive power of the random effects models
mean(ml.cbi$cpo$cpo, na.rm = TRUE)


##############################################
# 2. Baseline model + temporal random effect

# update the model formula
mf.cbi.re <- update(
 mf.cbi, . ~ 1 + . + 
  f(first_obs_date, model = "iid", 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  ) # temporal random effect
)
# fit the model                     
ml.cbi.re <- inla(
 mf.cbi.re, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re)
# check on predictive power of the random effects models
mean(ml.cbi.re$cpo$cpo, na.rm = TRUE)


#######################################################################
# 3. Baseline model + temporal random effect + fire-level random effect 

# update the model formula
mf.cbi.re2 <- update(
 mf.cbi.re, . ~ 1 + . + 
  f(fire_id, model = "iid", 
    hyper = list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  ) # fire-level random effect
)
# fit the model                     
ml.cbi.re2 <- inla(
 mf.cbi.re2, data = da, 
 family = "gamma",
 control.predictor = list(compute = TRUE),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re2)
# check on predictive power of the random effects models
mean(ml.cbi.re2$cpo$cpo, na.rm = TRUE)


#######################
# 4. Spatial SPDE model

# extracting spatial coordinates for grid centroids
grid_sf <- da %>%
 arrange(grid_index) %>%
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
# Plot mesh to check
plot(mesh)
points(coords, col = "red", pch = 20)

rm(coords)
gc()

######################
# Build the SPDE model
spde.ml <- inla.spde2.pcmatern(
 # Mesh and smoothness parameter
 mesh = mesh, 
 alpha = 2,
 # P(practic.range < 0.3) = 0.5
 prior.range = c(0.1, 0.5),
 # P(sigma > 1) = 0.01
 prior.sigma = c(1, 0.01) # variance
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

# Extract the predictor variables
X <- da %>% 
 select(fire_id, first_obs_date, grid_index, 
        fortypnm_gp, species_gp_n, total_tpp, forest_pct,      
        canopypct, H_tpp, H_ba, ba_live_pr, tpp_live_pr, qmd_live, 
        tree_ht_live, live_dead_pr_tpp, live_dead_pr_ba,
        erc_dv, vpd, vs, slope, tpi, chili, elev,
        aspen, day_prop)
head(X)

# Create the INLA data stack
stack.cbi <- inla.stack(
 data = list(CBIbc_p90 = da$CBIbc_p90),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)
dim(inla.stack.A(stack.cbi))

rm(grid_sf, X)
gc()

##########################
# update the model formula
mf.cbi.re.sp <- update(
 mf.cbi.re2, . ~ 1 + . + 
  f(mesh.idx, model = spde.ml) # spatial process model
)
# fit the model
ml.cbi.re.sp <- inla(
 mf.cbi.re.sp, 
 data = inla.stack.data(stack.cbi), ,
 family = "gamma",
 control.predictor = list(A = inla.stack.A(stack.cbi), compute=T),
 control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 control.fixed = list(
  prec = list(prior = "loggamma", param = c(0.5, 0.1))
 ),
 control.family = list(
  hyper = list(prec = list(prior = "loggamma", param = c(0.5, 0.1)))
 ),
 control.inla = list(strategy = "adaptive", int.strategy = "grid")
)
summary(ml.cbi.re.sp)

# check on predictive power of the random effects models
mean(ml.cbi.re.sp$cpo$cpo, na.rm = TRUE)



#==========EXTRACTING THE SPATIAL EFFECT===========#

# Compute spatial effect for each grid location
# mean
spat.eff <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.cbi.re.sp$summary.random$mesh.idx$mean  # FRP spatial field
# credible intervals
# Compute 2.5% (lower bound) credible interval
lower <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.cbi.re.sp$summary.random$mesh.idx$`0.025quant`
# Compute 97.5% (upper bound) credible interval
upper <- inla.spde.make.A(mesh, coords_mat) %*% 
 ml.cbi.re.sp$summary.random$mesh.idx$`0.975quant`
# Compute uncertainty as the credible interval width (upper - lower)
uncertainty <- upper - lower

# convert to a spatial dataframe
spat.eff.cbi <- data.frame(
 x = coords_mat[, 1],  # Longitude
 y = coords_mat[, 2],  # Latitude
 spat_effect = as.vector(spat.eff),  # Convert matrix to vector
 lower = as.vector(lower),  # 2.5% quantile
 upper = as.vector(upper),  # 97.5% quantile
 uncertainty = as.vector(uncertainty),  # CI width
 response = "CBI"
) %>%
 distinct(x, y, .keep_all=T) %>%
 st_as_sf(., coords = c("x", "y"), crs = 4326)
glimpse(spat.eff.cbi)

# write this a spatial file
st_write(spat.eff.cbi, paste0(maindir,"data/spatial/mod/spatial_effect_CBI_spp.gpkg"), append=F)

# Tidy up!
rm(spat.eff, lower, upper, uncertainty)
gc()



#=================MODEL COMPARISON=================#

# Create a model comparison table
ml_comparison.cbi <- tibble(
 Model = c("Baseline", 
           "W/Fire Random Effect", 
           "W/Fire + Temporal Effect",
           "W/Fire + Temporal + Spatial Effect"),
 Response = "CBI",
 DIC = c(
  ml.cbi$dic$dic,
  ml.cbi.re$dic$dic,
  ml.cbi.re2$dic$dic,
  ml.cbi.re.sp$dic$dic
 ),
 WAIC = c(
  ml.cbi$waic$waic,
  ml.cbi.re$waic$waic,
  ml.cbi.re2$waic$waic,
  ml.cbi.re.sp$waic$waic
 ),
 Marginal_LogLikelihood = c(
  ml.cbi$mlik[1],
  ml.cbi.re$mlik[1],
  ml.cbi.re2$mlik[1],
  ml.cbi.re.sp$mlik[1]
 ),
 Effective_Params = c(
  ml.cbi$dic$p.eff,
  ml.cbi.re$dic$p.eff,
  ml.cbi.re2$dic$p.eff,
  ml.cbi.re.sp$dic$p.eff
 ),
 Mean_CPO = c(
  mean(ml.cbi$cpo$cpo, na.rm = TRUE),
  mean(ml.cbi.re$cpo$cpo, na.rm = TRUE),
  mean(ml.cbi.re2$cpo$cpo, na.rm = TRUE),
  mean(ml.cbi.re.sp$cpo$cpo, na.rm = TRUE)
 )
) %>% arrange(DIC)
# Print the comparison table
print(ml_comparison.cbi)


# tidy up !
rm(A, coords_mat, field.idx, mesh, 
   ml.cbi, ml.cbi.re, ml.cbi.re2, spde.ml, stack.cbi)
gc()



#===========POSTERIOR EFFECTS===========#


#########################################
# Plot all of the posterior fixed effects
# Extract fixed effect marginals
cbi_marginals <- ml.cbi.re.sp$marginals.fixed
# Tidy marginals for all fixed effects
tidy.effects.cbi <- tibble::tibble(
 parameter = names(cbi_marginals),
 data = purrr::map(cbi_marginals, ~ as.data.frame(.x))
) %>%
 unnest(data) %>%
 filter(!parameter %in% c("(Intercept)","aspen","day_prop")) %>%  # Exclude the intercept
 mutate(
  effect = case_when(
   str_detect(parameter, "ba_live_pr") ~ "Proportion of \nLive Basal Area",
   str_detect(parameter, "qmd_live") ~ "Quadratic\nMean Diameter",
   str_detect(parameter, "H_tpp") ~ "Shannon Index\n(Abundance-based)",
   str_detect(parameter, "tree_ht_live") ~ "Tree height (m)",
   TRUE ~ parameter  # Default for all other fixed effects
  ),
  # Extract species names from parameters
  species = case_when(
   str_detect(parameter, "species_gp_n") ~ str_extract(parameter, "(?<=species_gp_n)[a-zA-Z_ñ]+"),
   str_detect(parameter, "fortypnm_gp") ~ str_extract(parameter, "(?<=fortypnm_gp)[a-zA-Z_ñ]+"),
   TRUE ~ NA_character_  # For non-species effects
  )
 ) %>%
 # Calculate mean effect size for ordering
 group_by(effect) %>%
 mutate(mean_effect = mean(x, na.rm = TRUE)) %>%
 ungroup() 
tidy.effects.cbi <- tidy.effects.cbi %>%
 # Order effects: Species-specific effects first, global effects by mean effect size
 mutate(
  effect_order = case_when(
   !is.na(species) ~ 2,  # Species-specific effects
   TRUE ~ 2              # Global effects
  ),
  effect = factor(effect, levels = tidy.effects.cbi %>%
                   arrange(effect_order, desc(mean_effect)) %>%
                   pull(effect) %>%
                   unique()),
  fill_species = ifelse(is.na(species), "Global Effect", species)
 ) %>%
 mutate(fill_species = recode(
  fill_species,
  "quaking_aspen" = "Quaking aspen",
  "lodgepole" = "Lodgepole",
  "mixed_conifer" = "Mixed-conifer",
  "piñon_juniper" = "Piñon-juniper",
  "ponderosa" = "Ponderosa",
  "spruce_fir" = "Spruce-fir"
 )) %>%
 mutate(
  fill_species = factor(fill_species, 
                        levels = c("Lodgepole", "Mixed-conifer", 
                                   "Piñon-juniper", "Ponderosa", 
                                   "Spruce-fir", "Quaking aspen")))


# check on the species name extraction
unique(tidy.effects.cbi$fill_species)
spps_breaks <- unique(tidy.effects.cbi$fill_species)


#######################
# create the ridge plot
cbi.p1 <- tidy.effects.cbi %>%
 filter(!parameter %in% c("aspen1"),
        !str_starts(parameter, "fortyp")) %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = fill_species)) +
 geom_density_ridges(stat = "identity", scale = 1.5, alpha = 0.7) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 scale_fill_manual(
  values = c(
   "Global Effect" = "gray",  # Neutral color for global effects
   "Quaking aspen" = "#1b9e77",
   "Lodgepole" = "#d95f02",
   "Mixed-conifer" = "#7570b3",
   "Piñon-juniper" = "#e7298a",
   "Ponderosa" = "#66a61e",
   "Spruce-fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species")
 ) +
 theme_minimal() +
 theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12)
 )
cbi.p1
# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_CBI_full.png')
ggsave(out_png, dpi = 500, bg = 'white')


####################################
# plot with just the species metrics
cbi.p2 <- tidy.effects.cbi %>%
 filter(fill_species != "Global Effect") %>%
 ggplot(., aes(x = x, y = effect, height = y, fill = fill_species, alpha = fill_species)) +
 geom_density_ridges(
  stat = "identity", scale = 1.5,
  # alpha = 0.7,
  show.legend=T
 ) +
 geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
 labs(
  x = "Effect Size",
  y = "Fixed Effect",
  fill = "Species"
 ) +
 # add a subplot label
 annotate("text", x = -0.13, y = 5, 
          label = "(B)", fontface = "bold", 
          size = 4, hjust = 0) +
 scale_fill_manual(
  values = c(
   "Quaking aspen" = "#1b9e77",
   "Lodgepole" = "#d95f02",
   "Mixed-conifer" = "#7570b3",
   "Piñon-juniper" = "#e7298a",
   "Ponderosa" = "#66a61e",
   "Spruce-fir" = "#e6ab02"
  ),
  # Exclude "Global Effect" from the legend
  breaks = spps_breaks,  
  guide = guide_legend(title = "Species", override.aes = list(size = 2))
 ) +
 scale_alpha_manual(
  values = c(
   "Quaking aspen" = 0.9, 
   "Lodgepole" = 0.7,  
   "Mixed-conifer" = 0.7,
   "Piñon-juniper" = 0.7,
   "Ponderosa" = 0.7,
   "Spruce-fir" = 0.7
  ),
  guide = "none"  # Hide alpha from legend
 ) +
 coord_cartesian(xlim=c(-0.13,0.15)) +
 theme_minimal() +
 theme(
  axis.text.y = element_text(angle = 0, hjust = 1, size=9),
  axis.text.x = element_text(angle = 0, hjust = 0, size=9),
  axis.title.y = element_text(size = 10, margin = margin(r = 8)),
  axis.title.x = element_text(size = 10, margin = margin(t = 10)),
  legend.position = c(0.91, 0.74),
  legend.background = element_rect(
   fill = scales::alpha("white", 0.7), 
   color = NA, size = 1),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  legend.key.height = unit(0.4, "cm"), 
  legend.key.width = unit(0.5, "cm"),
  legend.spacing.y = unit(0.8, "cm"),    
  legend.spacing.x = unit(0.6, "cm")
 )
cbi.p2

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_CBI_species.png')
ggsave(out_png, dpi = 500, bg = 'white')
# 
# # Tidy up!
# rm(cbi_marginals, tidy.effects)
# gc()


# Make a combined figure (FRP and CBI)
# Load the two ridge plots (assuming they are `frp_plot` and `cbi_plot`)
frp.cbi.p3 <- frp.p2 / cbi.p2 # "/" stacks them, "|" would place them side-by-side
frp.cbi.p3

# Save the plot
out_png <- paste0(maindir, 'figures/INLA_ForestComposition_FixedEffects_FRP-CBI_species.png')
ggsave(out_png, dpi = 500, width = 8, height = 6, bg = 'white')


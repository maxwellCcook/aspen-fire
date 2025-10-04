#---------
## INLA model of dominant forest type
#---------

library(tidyverse)
library(sf)
library(INLA)
inla.setOption(num.threads = "1:1")
getwd() # working directory


## ------------------------------
## 1. Read in the model data frame
## ------------------------------
# load the model data frame
da <- read_csv('data/tabular/mod/gridstats_model_data_c.csv') %>%
 # make sure we do not have any duplicate rows
 distinct(grid_idx, species_gp, .keep_all=TRUE)

## -------------------------------------
## 2. Set the species factor levels once
## -------------------------------------
unique(da$species_gp)
species_levels <- c("quaking_aspen","lodgepole_pine","ponderosa_pine",
                    "spruce_fir","douglas_fir","pinon_juniper","white_fir")
da <- da %>%
 mutate(
  fortypnm_gp = factor(fortypnm_gp, levels = species_levels),
  species_gp = factor(species_gp, levels = species_levels),
  dom_sp_ba = factor(dom_sp_ba, levels = species_levels),
  dom_sp_tpa = factor(dom_sp_tpa, levels = species_levels))
levels(da$species_gp)

## -------------------------------------------------------------
## 3. For dominant forest type model, keep one row per grid cell
## -------------------------------------------------------------
da_fortyp = da %>%
 distinct(grid_idx, fortypnm_gp, .keep_all=TRUE) %>%
 select(
  fire_id, grid_idx, # grid and fire ID
  cfrp, cbibc, # response variables
  first_obs_date, dt_max_frp, # date fields
  fortypnm_gp, fortypcd_pct, forest_prop, dom_sp_ba, dom_sp_tpa, # forest type
  cc, ch, cbd, cbh, # LF attributes
  vpd, vpd_dv, erc, erc_dv, vs, # fire weather
  elev, slope, tpi, eastness, northness, # topography
  ba_total, tpa_total, H_ba, H_tpa,
  log_fire_size, fire_aspen, grid_aspen,
  dist_to_perim, overlap, afd_count, day_prop,
  x, y
 )

glimpse(da_fortyp)

rm(da)

## -------------------------------------
## . Filtering
## -------------------------------------

###################################################
# check how many grids are majority forested, etc #
# percent of gridcells majority forested
# Check how many grids are > 50% forested
print(paste0(
 "Percent forested gridcells (>50% gridcell area): ",
 round(dim(da_fortyp %>% filter(forest_prop >= 0.50) %>% distinct(grid_idx))[1]/
        dim(da_fortyp %>% distinct(grid_idx))[1], 3) * 100, "%"
))

# retain predominantly forested gridcells ...
da_fortyp <- da_fortyp %>%
 filter(forest_prop >= 0.50)

# ###################################################
# # Forest type dominance
# print(paste0(
#  "Percent majority gridcells (>50% gridcell area): ",
#  round(dim(da_fortyp %>% filter(fortypcd_pct > 50) %>% distinct(grid_idx))[1]/
#         dim(da_fortyp %>% distinct(grid_idx))[1], 3) * 100, "%"
# ))
# 
# # retain predominantly forested gridcells ...
# da_fortyp <- da_fortyp %>%
#  filter(fortypcd_pct > 50)
# # Check how many gridcells and fires remain
# print(length(unique(da_fortyp$fire_id))) # unique fires
# print(length(unique(da_fortyp$grid_idx))) # unique gridcells

##########################################
# filter fires with not enough gridcells
# should have at least 10 (?) for spatial model

# check on the grid cell counts
gridcell_counts <- da_fortyp %>%
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
da_fortyp <- da_fortyp %>%
 filter(fire_id %in% fires_keep)

# Check how many gridcells and fires remain
print(length(unique(da_fortyp$fire_id))) # unique fires
print(length(unique(da_fortyp$grid_idx))) # unique gridcells

# tidy up!
rm(gridcell_counts, fires_keep, qt, qt10)


## -------------------------------------
## . Scaling
## -------------------------------------
da_fortyp <- da_fortyp %>%
 mutate(
  # set factors
  fire_id = as.factor(fire_id),
  grid_idx = as.factor(grid_idx),
  first_obs_date = as.factor(first_obs_date),
  dt_max_frp = as.factor(dt_max_frp),
  grid_aspen = as.factor(grid_aspen),
  fire_aspen = as.factor(fire_aspen),
  # scale percentages
  fortypcd_pct = fortypcd_pct / 100,
  # force aspen to be the baseline factor #
  fortypnm_gp = fct_relevel(fortypnm_gp, "quaking_aspen"),
  dom_sp_ba = fct_relevel(dom_sp_ba, "quaking_aspen"),
  dom_sp_tpa = fct_relevel(dom_sp_tpa, "quaking_aspen"),
  # center/scale metrics (fixed effects)
  across(
   # only scale numeric columns not in exclude list
   .cols = where(is.numeric) & !all_of(c("grid_idx","cfrp","cbibc","x","y")),  
   .fns = ~ as.numeric(scale(.))
  )) %>%
 # order canonicaly
 arrange(grid_idx)


## -------------------------------------
## . Correlations
## -------------------------------------
########################################
# correlation matrix for fixed effects #
cor_da <- da_fortyp %>%
 select(c("fortypnm_gp", where(is.numeric))) %>%
 select(-c("x","y")) %>%
 pivot_wider(
  names_from = fortypnm_gp,
  values_from = fortypcd_pct,
  values_fill = 0)
# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")
# Plot correlation matrix
ggcorrplot::ggcorrplot(cor_mat, method = "circle",
           type = "lower", lab = TRUE, lab_size = 3, tl.cex = 10,
           colors = c("blue", "white", "red")
)

# # save the plot.
# out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_FixedEffects_Model-FORTYP.png')
# ggsave(out_png, dpi=500, width=12, height=12, bg = 'white')

rm(cor_da, cor_mat) # tidy up


## ---------------------------------------------
## 4. Set up the spatial mesh grid (SPDE model)
## ---------------------------------------------

# extract gridcell coordinates
grid_sf <- da_fortyp %>%
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
# Plot mesh to check
plot(mesh, main = "SPDE Mesh for FRP Model")
points(coords, col = "red", pch = 20)
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

# sanity check
str(field.idx)


## -------------------------
## 5. Set up the INLA model
## -------------------------

set.seed(456)

#========Cumulative Fire Radiative Power (cfrp)========#

# Create the INLA data stack for cfrp
X <- da_fortyp %>% select(-c(cfrp,cbibc))
stack.frp <- inla.stack(
 data = list(cfrp = da_fortyp$cfrp),
 A = list(A, 1),  
 tag = 'est',
 effects = list(
  c(field.idx),
  list(as.data.frame(X))
 )
)

rm(X)

# dim(inla.stack.A(stack.frp))
# saveRDS(stack.frp, "code/R/models/stack_cfrp_fortyp.rds")


#################################################################
# Model with spatial process (SPDE) and fire-level random effects

####################
# define some priors
# fixed effects
pr.fixed <- list(
 prec.intercept = 1e-6, 
 prec = 0.01  # shrinkage prior (very weak)
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

# setup the model formula
mf.frp <- cfrp ~ 0 + 
 # forest composition and structure
 fortypnm_gp + # mean effect relative to aspen
 fortypnm_gp:fortypcd_pct + # majority forest type (factor, aspen baseline)
 # gridcell totals
 cbd + # landfire CBD/CBH
 H_tpa + # grid cell diversity in trees/acre
 ba_total + 
 # other fixed effects
 erc_dv + vpd + vs + # day-of-burn fire weather
 slope + northness + eastness + # topography 
 overlap + # gridcell VIIRS overlap (cumulative)
 day_prop + # daytime observation proportion
 # day_prop +  # gridcell proportion daytime detections
 dist_to_perim + # gridcell distance to perimeter
 # random/latent effects
 f(mesh.idx, model = spde.ml) + # spatial model
 f(fire_id, model = "iid", hyper = pr.pc_prec) # fire-level random effects

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
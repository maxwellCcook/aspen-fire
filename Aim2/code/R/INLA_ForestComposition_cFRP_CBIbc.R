#####################
# Modeling the effects of forest composition and structure on satellite-derived FRP and CBI
# Author: Maxwell C. Cook; 
# 
#####################

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
glimpse(da)

## ----------------------------------------------------
## 2. Set the factor levels for species (aspen as baseline)
## ----------------------------------------------------
unique(da$species_gp)
species_levels <- c("quaking_aspen","lodgepole_pine","ponderosa_pine",
                    "spruce_fir","douglas_fir","pinon_juniper","white_fir")
S <- length(species_levels)
da <- da %>%
 mutate(species_gp = factor(species_gp, levels = species_levels))
levels(da$species_gp)

## -------------------------------------------------------------
## 3. Build a grid-level table (one row per grid) + canonical order
## -------------------------------------------------------------
grid_covars <- c("cc","ch","cbd","cbh",
                 "erc_dv","vpd","vs",
                 "elev","slope","northness","eastness",
                 "overlap","day_prop","dist_to_perim")
grid_df <- da %>%
 group_by(grid_idx) %>%
 summarise(
  cfrp = first(cfrp),
  cbibc = first(cbibc),
  fire_id = first(fire_id),
  x = first(x),
  y = first(y),
  across(all_of(grid_covars), ~ first(.x)),
  .groups = "drop"
 ) %>%
 arrange(grid_idx)
head(grid_df)

## ---------------------------------------------------------
## 4. Build A-matrices that map species effects to each grid
## ---------------------------------------------------------

# helper for wide matrices: rows are grids, cols are species
make_wide <- function(long_df, grid_ids, species_levels,
                      value_col, weight_with = NULL, center_by_species = FALSE) {
 tmp <- long_df %>%
  filter(species_gp %in% species_levels) %>%
  select(grid_idx, species_gp, ba_live_pr, value = all_of(value_col))
 
 if (center_by_species) {
  tmp <- tmp %>%
   group_by(species_gp) %>%
   mutate(value = value - mean(value, na.rm=TRUE)) %>%
   ungroup()
 }
 if (!is.null(weight_with) && weight_with == "share") {
  tmp <- tmp %>% mutate(value = ba_live_pr * value)
 }
 
 tmp %>%
  group_by(grid_idx, species_gp) %>%
  summarise(val = sum(value, na.rm=TRUE), .groups="drop") %>%
  tidyr::pivot_wider(id_cols = grid_idx,
                     names_from = species_gp,
                     values_from = val,
                     values_fill = 0) %>%
  right_join(tibble(grid_idx = grid_ids), by = "grid_idx") %>%
  arrange(match(grid_idx, grid_ids)) %>%
  {
   for (sp in species_levels) if (!sp %in% names(.)) .[[sp]] <- 0
   as.matrix(.[, species_levels, drop=FALSE])
  }
}

# make the A-matrices
grid_ids <- grid_df$grid_idx

# baseline: species shares
A_baseline <- da %>%
 filter(species_gp %in% species_levels) %>%
 select(grid_idx, species_gp, val = ba_live_pr) %>%
 tidyr::pivot_wider(id_cols = grid_idx,
                    names_from = species_gp,
                    values_from = val,
                    values_fill = 0) %>%
 right_join(tibble(grid_idx = grid_ids), by = "grid_idx") %>%
 arrange(match(grid_idx, grid_ids)) %>%
 { as.matrix(.[, species_levels, drop=FALSE]) }

# traits
A_qmd <- make_wide(da, grid_ids, species_levels,
                   value_col="qmd_live", weight_with="share", center_by_species=TRUE)
A_hdr <- make_wide(da, grid_ids, species_levels,
                   value_col="hdr_live", weight_with="share", center_by_species=TRUE)
A_tpa <- make_wide(da, grid_ids, species_levels,
                   value_col="tpa_live", weight_with="share", center_by_species=TRUE)

# aspen co-occurrence
W_share <- A_baseline
w_aspen <- W_share[, "quaking_aspen", drop=TRUE]
partners <- setdiff(species_levels, "quaking_aspen")
X_aspen_pairs <- sapply(partners, function(k) w_aspen * W_share[, k])
colnames(X_aspen_pairs) <- paste0("cooc_aspen__with__", partners)

#####################
# quick sanity checks
stopifnot(nrow(A_baseline) == nrow(grid_df))
dim(A_baseline)  # should be (nrow(grid_df), 7)
dim(A_qmd)       # same nrow, 7 cols
dim(A_hdr)
dim(A_tpa)
dim(X_aspen_pairs)  # same nrow, 6 cols

colnames(A_baseline)
colnames(A_qmd)
colnames(X_aspen_pairs)


## -------------------
## 5. Build SPDE model
## -------------------
#--- 1. Coordinates from grid_df (already one row per grid) ---
grid_sf <- grid_df %>%
 st_as_sf(coords = c("x", "y"), crs = 4326)  # lon/lat WGS84
coords <- st_coordinates(grid_sf)
coords_mat <- as.matrix(coords)

#--- 2. Build the mesh ---
mesh <- inla.mesh.2d(
 loc = coords_mat,
 max.edge = c(1, 20),   # adjust resolution if needed
 cutoff = 0.01,         # min distance between points (0.01 deg ~ 1.1 km)
 offset = c(0.5, 0.1)   # buffer around domain
)

# check mesh visually
plot(mesh, main = "SPDE Mesh for FRP/CBI Model")
points(coords, col = "red", pch = 20)

#--- 3. SPDE model ---
spde.ml <- inla.spde2.pcmatern(
 mesh = mesh,
 alpha = 2,
 prior.range = c(0.3, 0.5),   # 50% prob that range < 0.3 degrees (~30 km at equator, ~5 km at mid-lat)
 prior.sigma = c(1, 0.05)     # 5% prob sigma > 1
)

#--- 4. Projector matrix ---
A_spde <- inla.spde.make.A(
 mesh = mesh,
 loc = coords_mat
)

#--- 5. Spatial field index ---
spde_idx <- inla.spde.make.index(
 name = "s",                   # name of spatial field
 n.spde = spde.ml$n.spde
)

# sanity check
str(spde_idx)

## -----------------------------
## 5. Build the INLA data stack
## -----------------------------
# -- Aspen co-occurrence from A_baseline (already aligned to grid_df)
W_share  <- A_baseline
w_aspen  <- W_share[, "quaking_aspen", drop=TRUE]
partners <- setdiff(species_levels, "quaking_aspen")

X_aspen_pairs <- sapply(partners, function(k) w_aspen * W_share[, k])
colnames(X_aspen_pairs) <- paste0("cooc_aspen__with__", partners)

# (Optional) center the co-occurrence columns to reduce collinearity
X_aspen_pairs <- scale(X_aspen_pairs, center = TRUE, scale = FALSE)

# -- Grid-level fixed effects + co-occurrence (row order = grid_df)
X_grid <- grid_df %>%
 transmute(
  Intercept = 1,
  fire_id = factor(fire_id),
  cc, ch, cbd, cbh,
  erc_dv, vpd, vs,
  elev, slope, northness, eastness,
  overlap, day_prop, dist_to_perim
 ) %>%
 bind_cols(as.data.frame(X_aspen_pairs))

# -- Stack for FRP
stk_frp <- inla.stack(
 tag  = "frp",
 data = list(y = grid_df$cfrp),
 A    = list(
  A_spde,      # spatial field
  1,           # fixed effects (X_grid)
  A_baseline,  # β0_s
  A_qmd,       # βq_s
  A_hdr,       # βh_s
  A_tpa        # βt_s
 ),
 effects = list(
  c(spde_idx),                       # s
  X_grid,                            # fixed effects + fire_id factor
  list(b0 = 1:length(species_levels)),
  list(bq = 1:length(species_levels)),
  list(bh = 1:length(species_levels)),
  list(bt = 1:length(species_levels))
 )
)

# Extract the data that INLA will actually use (already includes co-occurrence)
dat_frp <- inla.stack.data(stk_frp)

# Quick alignment checks
stopifnot(
 nrow(dat_frp) == nrow(grid_df),
 all(paste0("cooc_aspen__with__", partners) %in% names(dat_frp))
)


## -----------------
## 6. Fit the model
## -----------------
frp_form <- y ~ 0 + Intercept +
 cc + ch + cbd + cbh +
 erc_dv + vpd + vs +
 elev + slope + northness + eastness +
 overlap + day_prop + dist_to_perim +
 # spatial + fire RE
 f(s, model = spde.ml) +
 f(fire_id, model = "iid",
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05)))) +
 # species random-effect blocks (sum-to-zero for identifiability)
 f(b0, model = "iid", constr = TRUE,
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05)))) +
 f(bq, model = "iid", constr = TRUE,
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05)))) +
 f(bh, model = "iid", constr = TRUE,
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05)))) +
 f(bt, model = "iid", constr = TRUE,
   hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.05))))

fit_frp <- inla(
 frp_form,
 data = dat_frp,
 family = "lognormal",
 control.fixed = list(prec = 1),   # mild ridge on fixed effects; adjust/omit as desired
 control.predictor = list(A = inla.stack.A(stk_frp), compute = TRUE),
 control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
 verbose = T
)

summary(fit_frp)


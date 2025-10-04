#####################################################################
# Data prep for forest composition and structure modeling in R-INLA #

library(tidyverse)
library(sf)
library(ggcorrplot)
library(lubridate)
library(geosphere)

maindir <- '/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim2/'

# load and tidy the fire data
# calculate the fire perimeter and fire centroid
fires_fp <- paste0(maindir,"data/spatial/mod/srm_fire_census_2017_to_2023_ics_perims.gpkg")
fires <- st_read(fires_fp) %>%
 st_transform(st_crs(4326)) %>% # make sure its in lat/lon
 mutate(
  fire_perim = st_cast(geom, "MULTILINESTRING")
 ) %>%
 as_tibble() %>%
 rename(
  fire_id = Fire_ID,
  fire_ig_dt = DISCOVERY_DATE,
  fire_acres = ICS_ACRES,
  fire_aspenpct = Aspen_Pct
 ) %>%
 select(fire_id, fire_ig_dt, fire_acres, fire_aspenpct, fire_perim, geom) %>%
 mutate(
  fire_id = as.numeric(fire_id),
  fire_acres = as.numeric(fire_acres),
  log_fire_size = log(fire_acres),
  fire_aspenpct = as.numeric(fire_aspenpct))


#=========Prep the grid data=========#

# Format the forest composition and structure data frame 
# (long-format) each gridcell has rows for all species present
# each gridcell is assigned a majority forest type (FORTYPCD)
# climate and topography are summarized at the gridcell level
fp <- paste0(maindir,'data/tabular/mod/gridstats_model_data.csv')
grid_tm <-  read_csv(fp)  %>% # read in the file
 select(-c('...1')) %>%
 
 # join in some of the fire information
 left_join(fires%>%select(-geom), by="fire_id") %>%
 
 # tidy and create attributes ...
 mutate(
  fire_id = factor(fire_id),
  # tidy the temporal fields
  fire_ig_dt = as.Date(fire_ig_dt),  
  fire_year = year(fire_ig_dt),              
  fire_month = month(fire_ig_dt),
  first_obs_date = as.Date(first_obs_date), # the day of first observation
  first_obs_month = month(first_obs_date), # month of first observation
  first_obs_doy = as.numeric(lubridate::yday(first_obs_date)),
  # Create a "season" variable based on the month of first observation
  season = case_when(
   first_obs_month %in% c(3, 4, 5) ~ "spring",
   first_obs_month %in% c(6, 7, 8) ~ "summer",
   first_obs_month %in% c(9, 10, 11) ~ "fall"
  ),
  # Year/season interaction
  year_season = interaction(fire_year, season, drop = TRUE),
  
  # Factorize the temporal attributes
  season = as.factor(season),
  fire_year = as.factor(fire_year),
  year_season = as.factor(year_season),
  first_obs_date = as.factor(first_obs_date),
  
  # Format species names consistently
  species_gp = str_to_lower(species_gp),
  species_gp = str_replace_all(species_gp, "-", "_"),
  species_gp = str_replace_all(species_gp, " ", "_"),
  fortypnm_gp = str_to_lower(fortypnm_gp),
  fortypnm_gp = str_replace_all(fortypnm_gp, "-", "_"),
  fortypnm_gp = str_replace_all(fortypnm_gp, " ", "_"),
  dom_sp_ba = str_to_lower(dom_sp_ba),
  dom_sp_ba = str_replace_all(dom_sp_ba, "-", "_"),
  dom_sp_ba = str_replace_all(dom_sp_ba, " ", "_"),
  dom_sp_tpa = str_to_lower(dom_sp_tpa),
  dom_sp_tpa = str_replace_all(dom_sp_tpa, "-", "_"),
  dom_sp_tpa = str_replace_all(dom_sp_tpa, " ", "_"),
  
  # make sure factor variables are set for species names
  species_gp = factor(species_gp),
  fortypnm_gp = factor(fortypnm_gp),
  dom_sp_ba = factor(dom_sp_ba),
  dom_sp_tpa = factor(dom_sp_tpa),
  
  # log-scaled response variable (cFRP)
  cfrp = log(frp_csum),
  # beta-scaled response variable (CBIbc)
  cbibc = (cbibc_mean + 1e-4) / (3 + 2 * 1e-4),
  
  # fire size (log-scaled)
  log_fire_size = log(fire_acres), # log-scale fire size

  # scale the percentages
  forest_prop = forest_prop / 100,
  fire_aspenpct = fire_aspenpct / 100,
  
  # create 'northness' and 'eastness' variables from aspect
  # -1 is south-facing, 1 is north-facing
  aspect_rad = aspect * pi / 180, # convert to radians
  northness = cos(aspect_rad), # northness
  eastness = sin(aspect_rad), # eastness

 ) %>%
 
 # fire-level aspen presence
 group_by(fire_id) %>%
 mutate(fire_aspen = ifelse(any(species_gp == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 # gridcell aspen presence
 group_by(grid_idx) %>%
 mutate(grid_aspen = ifelse(any(species_gp == "quaking_aspen"), 1, 0)) %>%
 ungroup() %>%
 mutate(
  fire_aspen = as.factor(fire_aspen),
  grid_aspen = as.factor(grid_aspen),
 ) %>%
 
 # aspen-specific metrics
 group_by(grid_idx) %>%
 mutate(
  ba_total = sum(ba_live), # total live basal area
  tpa_total = sum(tpa_live),
  # summary stats of aspen-specific metrics
  aspen_ba = sum(if_else(species_gp == "quaking_aspen", ba_live, 0.0), na.rm = TRUE),
  aspen_tpa = sum(if_else(species_gp == "quaking_aspen", tpa_live, 0.0), na.rm = TRUE),
  aspen_hdr = mean(if_else(species_gp == "quaking_aspen", hdr_live, NA_real_), na.rm = TRUE),
  aspen_qmd = mean(if_else(species_gp == "quaking_aspen", qmd_live, NA_real_), na.rm = TRUE),
  # proportions 
  aspen_ba_pr = if_else(ba_total > 0, aspen_ba / ba_total, 0.0),
  aspen_tpa_pr = if_else(tpa_total > 0, aspen_tpa / tpa_total, 0.0)
 ) %>%
 ungroup() %>%
 
 # flag re-burn gridcells
 group_by(grid_idx) %>%
 mutate(
  # Get earliest fire date for this gridcell
  first_fire_dt = min(fire_ig_dt, na.rm = TRUE),
  # Flag reburn: if current fire date is after earliest
  reburn = as.factor(fire_ig_dt > first_fire_dt)
 ) %>%
 ungroup() %>%
 
 # be sure there are no duplicate rows for grid/fire/species
 distinct(grid_idx, species_gp, .keep_all = TRUE) # remove duplicates
glimpse(grid_tm)


############################################
# calculate the distance to fire perimeter #
grid_sf <- grid_tm %>%
 distinct(grid_idx, .keep_all=T) %>%
 select(grid_idx, x, y, fire_perim, fire_id) %>%
 st_as_sf(coords = c("x", "y"), crs = 4326) %>%
 # Compute minimum distance from gridcell center to fire perimeter
 rowwise() %>%
 mutate(
  dist_to_perim = as.numeric(min(st_distance(geometry, fire_perim)))
 ) %>%
 ungroup() %>%
 # Normalize distance within each fire
 group_by(fire_id) %>%
 mutate(dist_to_perim = dist_to_perim / max(dist_to_perim, na.rm = TRUE)) %>%
 ungroup() %>%
 # Convert back to tibble for easy manipulation
 st_drop_geometry() %>%
 select(-fire_perim) %>% # remove this large geometry field
 select(grid_idx, dist_to_perim)

# join back to the gridcell data
grid_tm <- grid_tm %>%
 select(-fire_perim) %>% # remove this large geometry field
 left_join(., grid_sf, by = "grid_idx")
head(grid_tm%>%select(grid_idx, dist_to_perim))
rm(grid_sf) # tidy up
gc()


###################################################
# check how many grids are majority forested, etc #
print(length(unique(grid_tm$grid_idx))) # unique gridcells
length(unique(grid_tm$fire_id))
# percent of gridcells majority forested
# Check how many grids are > 50% forested
print(paste0(
 "Proportion forested gridcells (>50% forest pixels): ",
 round(dim(grid_tm %>% filter(forest_prop > 0.50) %>% distinct(grid_idx))[1]/
        dim(grid_tm %>% distinct(grid_idx))[1], 3) * 100, "%"
))


###############################
# check on the species counts #
grid_tm %>% 
 group_by(species_gp) %>%
 summarise(n = length(species_gp)) %>%
 arrange(desc(n))


###################################################
# save the model data frame (tabular and spatial) #
# tabular data (model data)
out_fp = paste0(maindir,"data/tabular/mod/gridstats_model_data_c.csv")
write_csv(grid_tm, out_fp)

# spatial data (keep one row per gridcell)
grid_tm.w <- grid_tm %>%
 pivot_wider(
  names_from = species_gp,
  values_from = ba_live,
  values_fill = 0
 ) %>%
 left_join(., grid_tm %>% distinct(grid_idx, geometry), by="grid_idx") %>%
 distinct(grid_idx, .keep_all = T) %>%
 mutate(geometry = st_as_sfc(geometry)) %>%
 st_as_sf(.) %>% 
 st_set_crs(st_crs(5070))
glimpse(grid_tm.w)
# save it out.
out_fp = paste0(maindir,"data/spatial/mod/gridcell_model_data_cleaned.gpkg")
st_write(grid_tm.w, out_fp, append=F)
rm(grid_tm.w)


#===============Explore Distributions, etc.================#

####################################
# distribution of response variables
resp_plot <- grid_tm %>%
 # pivot longer to facet plot
 pivot_longer(cols = c(log_frp_csum,
                       CBIbc_p90),
              names_to = "variable",
              values_to = "value") %>%
 # Plot with facets
 ggplot(aes(x = value)) +
 geom_histogram(bins = 30, fill = "orange", alpha = 0.7) +
 facet_wrap(
  ~ variable,
  scales = "free",
  labeller = as_labeller(c(log_frp_csum = "log(Cumulative FRP)",
                           CBIbc_p90 = "90th percentile CBIbc"))) +
 labs(x = "value",
      y = "Frequency") +
 theme_classic()
resp_plot
# save the plot.
out_png <- paste0(maindir,'figures/INLA_ResponseDistribution_FRP-CBI.png')
ggsave(out_png, plot = resp_plot, dpi=500, bg = 'white')
rm(resp_plot)

#==========EXPLORE CORRELATIONS==========#

#######################################
# species-specific correlation matrix #
sp_cor <- grid_tm %>%
 select(grid_idx, species_gp_n, ba_live) %>%
 spread(species_gp_n, ba_live, fill = 0)  # Reshape to wide format
# compute the correlation matrix
sp_cormat <- cor(sp_cor[,-1], use = "complete.obs", method = "spearman")
ggcorrplot(sp_cormat, method = "circle", type = "lower",
           lab = TRUE, lab_size = 3, colors = c("blue", "white", "red"))
rm(sp_cor, sp_cormat)
# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_SpeciesBA.png')
ggsave(out_png, dpi=500, bg = 'white')
rm(sp_cor,sp_cormat)

########################################
# correlation matrix for fixed effects #
cor_da <- grid_tm %>%
 mutate(
  log_fire_size = log(fire_acres)
 ) %>%
 select(
  fortypnm_gp,
  # species structure metrics
  tpa_live, ba_live, qmd_live, ht_live, dia_live, hdr_live, ba_per,
  tpa_live_total, ba_live_total, qmd_live_mean,
  tpa_dead_total, ba_dead_total, # dead BA and tpa
  tpa_live_pr, ba_live_pr, # proportions of tpa and BA
  forest_pct, fortyp_pct, # gridcell forest percent and majority forest type percent
  H_ba, H_tpa, # gridcell species diversity (abundance- and dominance-based)
  erc, erc_dv, vpd, vpd_dv, # day-of-burn climate
  fm1000, rmin, tmmx, vs, # day-of-burn climate
  elev, slope, tpi, northness,  # topography
  lf_canopy, lf_height, # gridcell canopy percent/height and BA sum
  day_prop, overlap, # VIIRS detection characteristics
  log_fire_size, dist_to_perim # fire size and distance to perimeter
 ) %>%
 pivot_wider(
  names_from = fortypnm_gp,
  values_from = fortyp_pct,
  values_fill = 0) %>%
 mutate(across(everything(), ~ scale(.) %>% as.numeric()))  # Standardize variables

# Compute correlation matrix
cor_mat <- cor(cor_da, use = "complete.obs", method = "spearman")
# Plot correlation matrix
ggcorrplot(cor_mat, method = "circle",
 type = "lower", lab = TRUE, lab_size = 3, tl.cex = 10,
 colors = c("blue", "white", "red")
)
rm(cor_da, cor_mat) # tidy up
# save the plot.
out_png <- paste0(maindir,'figures/INLA_CorrelationMatrix_FixedEffects_Full.png')
ggsave(out_png, dpi=500, width=12, height=12, bg = 'white')


##################################################
# check on the dominance/abundance distributions #
# get value counts for dominance/abundance by species
(sp_dom_summary <- grid_tm %>%
  distinct(grid_idx, fortypnm_gp, 
           dom_sp_ba, dom_sp_tpa,
           dom_sp_ht, dom_sp_dia,
           dom_sp_qmd) %>%
  group_by(fortypnm_gp) %>%
  summarise(
   n_fortyp = n(),
   n_ba = sum(dom_sp_ba == fortypnm_gp, na.rm = TRUE),
   n_tpa = sum(dom_sp_tpa == fortypnm_gp, na.rm = TRUE),
   n_ht = sum(dom_sp_ht == fortypnm_gp, na.rm = TRUE),
   n_dia = sum(dom_sp_dia == fortypnm_gp, na.rm = TRUE),
   n_qmd = sum(dom_sp_qmd == fortypnm_gp, na.rm = TRUE)
  ) %>%
  arrange(desc(n_fortyp)))
# save this table out
write_csv(sp_dom_summary,
          paste0(maindir,"data/tabular/mod/results/species_dominance_counts.csv"))
rm(sp_dom_summary)










#============================#
# Threshold Analysis Script  #
#============================#

# Load libraries
library(tidyverse)
library(INLA)


#-----------------------------#
# 1. Load model and data
#-----------------------------#

# Load fitted INLA model
ml.cbi.re.sp <- readRDS("code/R/models/ml_cbi_re_sp.rds")
# Load model data frame
da <- readRDS("code/R/models/cbi_model_da.rds")
# load the stack
stack.cbi <- readRDS("code/R/models/stack_cbi.rds")
# Load the mean/sd for fixed effects (pre-scaling)
sc_lookup <- read_csv("data/tabular/mod/var_scaling_lookup.csv")

# ----------------------------------#
# 1b. Extract the fixed effect coefficients
fixed_eff <- ml.cbi.re.sp$summary.fixed
# gather the standard deviation for forest type proportion
sd_fortyp <- sc_lookup$fortyp_pct_sd
# gather the aspen effect
aspen_effect <- fixed_eff["fortypnm_gpquaking_aspen:fortyp_pct", ]

# How much benefit for every 10% increase in aspen cover
thresh <- 0.10  # 10%
(rsc_mean <- aspen_effect$mean * (thresh / sd_fortyp))
(rsc_lwr  <- aspen_effect$`0.025quant` * (thresh / sd_fortyp))
(rsc_upr  <- aspen_effect$`0.975quant` * (thresh / sd_fortyp))
(exp(rsc_mean) - 1) * 100

#####
# FRP

# Load fitted INLA model
ml.frp.re.sp <- readRDS("code/R/models/ml_frp_re_sp.rds")
# Load model data frame
da <- readRDS("code/R/models/frp_model_da.rds")
# load the stack
stack.frp <- readRDS("code/R/models/stack_frp.rds")
# Load the mean/sd for fixed effects (pre-scaling)
sc_lookup <- read_csv("data/tabular/mod/var_scaling_lookup.csv")

# Extract the fixed effect coefficients
fixed_eff <- ml.frp.re.sp$summary.fixed
# gather the standard deviation for forest type proportion
sd_fortyp <- sc_lookup$fortyp_pct_sd
# gather the aspen effect
aspen_effect <- fixed_eff["fortypnm_gpquaking_aspen:fortyp_pct", ]

# How much benefit for every 10% increase in aspen cover
thresh <- 0.10  # 10%
(rsc_mean <- aspen_effect$mean * (thresh / sd_fortyp))
(rsc_lwr  <- aspen_effect$`0.025quant` * (thresh / sd_fortyp))
(rsc_upr  <- aspen_effect$`0.975quant` * (thresh / sd_fortyp))
(exp(rsc_mean) - 1) * 100


#----------------------------------------#
# 2. Threshold function for observed data
#----------------------------------------#


# ----------------------------------#

# Find indices for observed data (not mesh nodes)
idx_est <- inla.stack.index(stack.cbi, tag = "est")$data

# Attach ONLY the correct fitted values to da
da$CBI_fitted <- ml.cbi.re.sp$summary.fitted.values$mean[idx_est]





#----------------------------------------#
# 2. Threshold function for observed data
#----------------------------------------#

thresholding <- function(model, data, forest_type, scaling_lookup,
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         labels = c("0–25%", "25–50%", "50–75%", "75–100%")) {
 # 1. Filter to forest type of interest
 df <- data %>%
  filter(fortypnm_gp == forest_type)
 
 # No need to re-add fitted values here!
 
 # 2. Automatically grab scaling parameters from lookup
 mean_fortyp <- scaling_lookup$fortyp_pct_mean
 sd_fortyp   <- scaling_lookup$fortyp_pct_sd
 
 # 3. Back-transform scaled `fortyp_pct`
 df <- df %>%
  mutate(
   fortyp_pct_raw = fortyp_pct * sd_fortyp + mean_fortyp,
   bin = cut(fortyp_pct_raw, breaks = breaks, labels = labels, include.lowest = TRUE)
  )
 
 # 4. Summarize fitted CBI by bin
 summary_df <- df %>%
  group_by(bin) %>%
  summarise(
   mean_CBI = mean(CBI_fitted, na.rm = TRUE),
   lwr = quantile(CBI_fitted, 0.025, na.rm = TRUE),
   upr = quantile(CBI_fitted, 0.975, na.rm = TRUE),
   n = n(),
   .groups = "drop"
  ) %>%
  filter(!is.na(bin))  # optional: remove any empty bin rows
 
 # 5. Calculate percent reduction from first bin (reference)
 baseline <- summary_df$mean_CBI[1]
 summary_df <- summary_df %>%
  mutate(
   pct_reduction = 100 * (baseline - mean_CBI) / baseline
  )
 
 return(summary_df)
}


#-----------------------------#
# 3. Run Simulation for Aspen
#-----------------------------#

results_aspen <- thresholding(
 model = ml.cbi.re.sp,
 data = da,
 forest_type = "quaking_aspen",
 scaling_lookup = sc_lookup
)

print(results_aspen)



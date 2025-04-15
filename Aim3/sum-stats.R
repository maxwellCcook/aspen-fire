
library(tidyverse)
library(sf)

projdir <- "/Users/max/Library/CloudStorage/OneDrive-Personal/mcook/aspen-fire/Aim3/"

# load and prep the firesheds data
firesheds <- st_read(paste0(projdir,"data/spatial/mod/srm_firesheds_model_data_wBivar.gpkg")) %>%
 mutate(bi_class = as.factor(bi_class))
glimpse(firesheds)

summary(firesheds$patch_size)

# Water Content Data Analysis
# Heili Lowman
# 3/15/22

# The following script will analyze %H2O data collected
# for the sediment nutrient recycling manuscript.

# Load packages
library(tidyverse) 
library(lubridate)
library(plotly)
library(viridis)
library(patchwork)

# Load necessary dataset.
water_dat <- read_csv("data_raw/Sediment_water_020922.csv")

# Calculate additional per tin values.
water_dat <- water_dat %>%
  mutate(WetSed = TinWetSed - `Tin...7`,
         DrySed = TinDrySed - `Tin...7`) %>%
  mutate(Water = WetSed - DrySed) %>%
  mutate(perc_Water = Water/WetSed)

# Summarize by site.
water_sites <- water_dat %>%
  group_by(Site) %>%
  summarize(meanW = mean(perc_Water, na.rm = TRUE),
            sdW = sd(perc_Water, na.rm = TRUE))

# End of script.

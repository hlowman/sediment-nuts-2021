# Additional Data for Reviews
# Heili Lowman
# 9/2/22

# The following script will analyze CTD data collected
# by the SBC LTER during June and August 2016
# for the sediment nutrient recycling manuscript.

# This data was collected by Libe Washburn et al. and can be downloaded
# from the SBC LTER's data repository at the following URL:
# https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.10

# Load packages
library(tidyverse) 
library(lubridate)
library(viridis)
library(patchwork)

# Load necessary dataset.
ctd <- read_csv("data_raw/LTER_monthly_downcasts_registered_stations_20220606.csv")

# Filter for sites of interest.
june <- ctd %>%
  filter(`yyyy-mm-dd` == "6/8/16" & Station == "AB")

august <- ctd %>%
  filter(`yyyy-mm-dd` == "8/3/16" & Station == "MK")

# Plot ctd data.
(june_AB <- ggplot(june, aes(x = ctd_sigmatheta00, 
                             y = (ctd_depth*-1))) +
  geom_line() +
  theme_bw() +
  labs(x = "Sigma (Density)",
       y = "Depth",
       title = "June 8 2016, Arroyo Burro"))

(aug_MK <- ggplot(august, aes(x = ctd_sigmatheta00, 
                              y = (ctd_depth*-1))) +
    geom_line() +
    theme_bw() +
    labs(x = "Sigma (Density)",
         y = "Depth",
         title = "August 3 2016, Mohawk"))

(combined_sigma <- june_AB + aug_MK)

ggsave("density.jpg",
         path = "figures",
         width = 20,
         height = 10,
         units = "cm")

# End of script.

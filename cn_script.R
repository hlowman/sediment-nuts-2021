# C:N Data Analysis
# Heili Lowman
# 2/18/22

# The following script will analyze C:N data collected
# by the SBC LTER during the NSF RAPID study from 2015-2017
# for the sediment nutrient recycling manuscript.

# This data was collected by Mark Page et al. and can be downloaded
# from the SBC LTER's data repository at the following URL:
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.118.1

# Load packages
library(tidyverse) 
library(lubridate)
library(plotly)
library(viridis)
library(patchwork)

# Load necessary dataset.
rapid <- read_csv("data_tidy/RAPID_Master_Datasheet_add_region_streamcode.csv")

# Filter for sites of interest.
marine <- rapid %>%
  filter(Type %in% c("Marine", "MarineRunoff")) %>%
  filter(Site %in% c("REFU", "GOLB", "ABUR", "MOHK", "MICR")) %>%
  filter(Depth == 20) %>%
  mutate(`C/N` = as.numeric(`C/N`))

# Generate summary statistics.
summary <- marine %>%
  group_by(Site) %>%
  summarize(meanOM = mean(perc_OM, na.rm = TRUE),
            sdOM = sd(perc_OM, na.rm = TRUE),
            meanCN = mean(`C/N`, na.rm = TRUE),
            sdCN = sd(`C/N`, na.rm = TRUE),
            medianGS = mean(`Diameter_50.00%`, na.rm = TRUE),
            sdGS = sd(`Diameter_50.00%`, na.rm = TRUE),
            meanClay = mean(perc_Clay, na.rm = TRUE),
            sdClay = sd(perc_Clay, na.rm = TRUE),
            meanSilt = mean(perc_Silt, na.rm = TRUE),
            sdSilt = sd(perc_Silt, na.rm = TRUE),
            meanSand = mean(perc_Sand, na.rm = TRUE),
            sdSand = sd(perc_Sand, na.rm = TRUE))

# Export data for inclusion in manuscript.
write_csv(summary, "data_tidy/sediment_content_stats.csv")

# Overall summary statistics for results section.
summary2 <- marine %>%
  summarize(meanOM = mean(perc_OM, na.rm = TRUE),
            sdOM = sd(perc_OM, na.rm = TRUE),
            meanCN = mean(`C/N`, na.rm = TRUE),
            sdCN = sd(`C/N`, na.rm = TRUE),
            medianGS = mean(`Diameter_50.00%`, na.rm = TRUE),
            sdGS = sd(`Diameter_50.00%`, na.rm = TRUE),
            meanClay = mean(perc_Clay, na.rm = TRUE),
            sdClay = sd(perc_Clay, na.rm = TRUE),
            meanSilt = mean(perc_Silt, na.rm = TRUE),
            sdSilt = sd(perc_Silt, na.rm = TRUE),
            meanSand = mean(perc_Sand, na.rm = TRUE),
            sdSand = sd(perc_Sand, na.rm = TRUE))

# End of script.

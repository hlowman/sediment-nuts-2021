# Table 2 Creation
# Heili Lowman
# 6/4/21

# The following script will calculate summary values to be used in the 
# sediment nutrient recycling manuscript.

# This data was collected by Jason Smith et al. and can be downloaded
# from the SBC LTER's data repository at the following URL:
# https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.117.4

# Load packages
library(tidyverse) 

# Load necessary dataset
watercolumn <- read_csv("data_raw/WaterColumn_NH4_Urea_20200331.csv")

# Summarize ammonium data since some duplicates were analyzed
WCstats <- watercolumn %>% # Takes the original dataset and then
  group_by(STATION, WATER_DEPTH) %>% # Groups the data and then
  summarize(meanConc = mean(NH4_water), sdConc = sd(NH4_water)) %>% # Calculates means and standard deviations and then
  ungroup() # Removes groupings.

# Create new dataset with max water depth included.
WC_w_sites <- watercolumn %>% # using the original dataset
  mutate(TOTAL_DEPTH = c(15,15,15,15,15,15,15,15,15,5,5,5,7,7,7,10,10,10,15,15,15,20,20,20)) # add total water depth

# Export the dataset
write_csv(WC_w_sites, "data_tidy/WC_w_sites.csv")

# Additional summary statistics for the results section:

# minimum concentration
min(WCstats$meanConc)

# maximum concentration
max(WCstats$meanConc)

# mean concentration
mean(WCstats$meanConc)

# s.d. of concentration
sd(WCstats$meanConc)

# quick visualization of overlying water trends
wplot <- WC_w_sites %>%
  mutate(depthf = factor(TOTAL_DEPTH)) %>%
  ggplot(aes(x = WATER_DEPTH, 
             y = NH4_water, 
             color = depthf)) +
  geom_line() +
  facet_wrap(~STATION) +
  theme_bw()

wplot  

# End of script.

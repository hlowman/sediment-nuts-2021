# Table 1 Creation
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
porewater <- read_csv("data_raw/Porewater_NH4_Urea_20180319.csv")

# Summarize ammonium data since duplicates or triplicates were analyzed
PWstats <- porewater %>% # Takes the original dataset and then
  group_by(STATION, WATER_DEPTH, SED_DEPTH) %>% # Groups the data and then
  # not grouping by date because these correspond with sampling site
  summarize(meanConc = mean(NH4_pore), # no NAs to worry about
            sdConc = sd(NH4_pore)) %>% # Calculates means/st deviations
  ungroup() # Don't forget it!!

# Exports summary dataset
write_csv(PWstats, "data_tidy/PWstats.csv")

# Additional summary statistics for the results section:

# minimum concentration
min(PWstats$meanConc)

# maximum concentration
max(PWstats$meanConc)

# mean concentration
mean(PWstats$meanConc)

# s.d. of concentration
sd(PWstats$meanConc)

# quick visualization of porewater trends
pplot <- PWstats %>%
  mutate(depthf = factor(WATER_DEPTH)) %>%
  ggplot(aes(x = SED_DEPTH, 
             y = meanConc, 
             color = depthf)) +
  geom_line() +
  facet_wrap(~STATION) +
  theme_bw()

pplot

# mean concentration in the 0-2cm horizon
top_horiz <- PWstats %>%
  filter(SED_DEPTH == 1) %>%
  summarize(meanNH4 = mean(meanConc))

# mean concentration in all horizons across all sites/depths sampled
presentation_stats <- porewater %>% # Takes the original dataset and then
  mutate(depth = case_when(SED_DEPTH == 8.0 ~ 9.0,
                           TRUE ~ SED_DEPTH)) %>%
  group_by(depth) %>% # Groups the data and then
  summarize(meanConc = mean(NH4_pore), # no NAs to worry about
            sdConc = sd(NH4_pore)) %>% # Calculates means/st deviations
  ungroup() # Don't forget it!!

# End of script.

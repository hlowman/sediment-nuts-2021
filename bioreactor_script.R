# Bioreactor Data Analysis
# Heili Lowman
# 6/18/21

# The following script will analyze nutrient flux data collected
# using the sediment bioreactors run in 2017, 2018, and 2019 for
# the sediment nutrient recycling manuscript.


# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse) 
library(lubridate)
library(viridis)
library(patchwork)

# Load necessary datasets
nh4dat <- read_csv("data_raw/Bioreactor_NH4_2017_2019.csv")
no3dat <- read_csv("data_raw/Bioreactor_NO3_2017_2019.csv")
# Note - using cor_NO3 values below, since these have been corrected
# to report nothing below the limit of detection (0.5uM)
tdndat <- read_csv("data_raw/Bioreactor_DOC_TDN_2017_2019.csv")

# Aggregate replicates 
nh4dat_ed <- nh4dat %>% # Takes the original dataset
  group_by(Year, Site, Treatment, Sample_ID_1) %>% # groups the data
  summarize(meanNH4 = mean(Conc, na.rm = TRUE), 
            sdNH4 = sd(Conc, na.rm = TRUE)) %>% # calculations
  ungroup() %>% # don't forget it!
  mutate(Analyte = "NH4") %>% # add new column for identification
  rename(mean = meanNH4, 
         sd = sdNH4,
         Sample_ID = Sample_ID_1) %>% # rename columns for combining
  filter(Site != "GOSL") %>% # remove estuarine results
  drop_na(Treatment) # remove results from other runs

no3dat_ed <- no3dat %>% 
  group_by(Year, Site, Treatment, Sample_ID) %>%
  summarize(meanNO3 = mean(cor_NO3, na.rm = TRUE), 
            sdNO3 = sd(cor_NO3, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(Analyte = "NO3") %>% 
  rename(mean = meanNO3, 
         sd = sdNO3) %>% 
  filter(Site != "GOSL") %>% 
  drop_na(Treatment)

tdndat_ed <- tdndat %>% 
  group_by(Year, Site, Treatment, Sample_ID) %>%
  summarize(meanTDN = mean(raw_TDN, na.rm = TRUE), 
            sdTDN = sd(raw_TDN, na.rm = TRUE)) %>% # only calculating TDN
  ungroup() %>% 
  mutate(Analyte = "TDN") %>% 
  rename(mean = meanTDN, 
         sd = sdTDN) %>% 
  filter(Site != "GOSL") %>% 
  drop_na(Treatment)

# Stack three datasets on one another
nutdat <- bind_rows(nh4dat_ed, no3dat_ed, tdndat_ed)

# And finally, remove the cores that were pulled due to leakage/other issues
# These filters should cover any and all analytes
# Note, some of these will be missing to begin with, since the analytes were never measured
nutdat_clean <- nutdat %>%
  mutate(mean_clean = case_when(
    Year == 2017 & Site == "MICR" & Sample_ID == "Paul" ~ NA_real_,
    Year == 2018 & Site == "MICR" & Sample_ID == "Orange" ~ NA_real_,
    Year == 2018 & Site == "MICR" & Sample_ID == "Pink" ~ NA_real_,
    Year == 2019 & Site == "ABUR" & Sample_ID == "Blue" ~ NA_real_,
    Year == 2019 & Site == "GOLB" & Sample_ID == "Orange" ~ NA_real_,
    Year == 2019 & Site == "MICR" & Sample_ID == "Red" ~ NA_real_,
    TRUE ~ mean)) # making a new column because filtering wasn't working appropriately

# Flux Calculations -------------------------------------------------------

# Creating a newly revised for loop to calculate changes in concentration

site_list <- c("ABUR", "MICR", "GOLB") # create sortable site list
analyte_list <- c("NH4", "NO3", "TDN") # same for analytes

for(year in 2017:2019){ # iterate over year
  for(site in site_list){ # followed by site
    for(analyte in analyte_list){ # followed by analyte
  
  # filter by iteration
  newdat <- nutdat_clean %>%
    filter(Year == year & Site == site & Analyte == analyte)
  
  # create variable for starting concentration of seawater
  start <- newdat %>%
    filter(Treatment == "Before") %>%
    pull(mean_clean)
  
  # calculate changes and rates
  changes <- newdat %>%
    mutate(change = mean_clean - start) %>%
    mutate(change_hr = change / 3) # all bioreactors run for 3 hours
    
  # export data
  file.name <- paste0("data_analyses/", site, "_", year, "_", 
                      analyte, ".rds") # create file name
  saveRDS(changes, file=file.name)
  
    }
  }
}




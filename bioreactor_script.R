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
library(here)

# Load necessary datasets
nh4dat <- read_csv("data_raw/Bioreactor_NH4_2017_2019.csv")
no3dat <- read_csv("data_raw/Bioreactor_NO3_2017_2019.csv")
# Note - using cor_NO3 values below, since these have been corrected
# to report nothing below the limit of detection (0.5uM in 2017, and 0.2uM in 2018/2019)
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
  file.name <- paste0(site, "_", year, "_", 
                      analyte, ".rds") # create file name
  saveRDS(changes, file=file.name)
  
    }
  }
}

# load in and combine all the created datasets
# I chose to export the files so that I could see they mapped correctly individually
# also for loops need empty receiving dataframes in the environment,
# so it's often easier to export the files to a working data directory instead
nutdat_changes <- here("data_analyses/") %>%
  list.files(pattern = ".rds") %>%
  map(readRDS) %>%
  bind_rows()

# Check to be sure it all looks alright
View(nutdat_changes)
# woohoo!!

# Next, I'll trim down the dataset to help with flux calculations
nutdat_trim <- nutdat_changes %>%
  filter(Treatment %in% c("Control", "Experimental")) %>% # filter out the "before" samples
  # since they have already served their purpose
  select(Year, Site, Treatment, Analyte, change_hr) # trim columns

# Need to aggregate by control/experimental treatment for each run to calculate net flux
nutdat_net <- nutdat_trim %>%
  group_by(Year, Site, Treatment, Analyte) %>%
  summarize(mean_change_hr = mean(change_hr, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Treatment, values_from = mean_change_hr) %>%
  rename(Experimental_change_hr = Experimental, Control_change_hr = Control) %>%
  mutate(Net_change_hr = Experimental_change_hr - Control_change_hr) %>%
  # convert net flux to umol/m^2 * hour
  mutate(Net_change_hr_m2 = (Net_change_hr * 0.25 * 10000) / 19.6)


# Summary Stats -----------------------------------------------------------

# Running some additional calculations for inclusion in the manuscript.
summary <- nutdat_net %>%
  group_by(Analyte) %>%
  summarize(mean_all = mean(Net_change_hr_m2, na.rm = TRUE),
            min_all = min(Net_change_hr_m2, na.rm = TRUE),
            max_all = max(Net_change_hr_m2, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(daily_mean = mean_all * 24) # convert mean flux to umol/m^2 * day

# Creating another dataframe for table 3 in the manuscript.
table3 <- nutdat_net %>%
  select(Year, Site, Analyte, Net_change_hr_m2) %>%
  pivot_wider(names_from = Analyte, values_from = Net_change_hr_m2)

write_csv(table3, path = "data_analyses/table3_output.csv")

# End of script.

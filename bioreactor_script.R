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
nh4dat <- read.csv("data_raw/Bioreactor_NH4_2017_2019.csv")
no3dat <- read.csv("data_raw/Bioreactor_NO3_2017_2019.csv")
# Note - using cor_NO3 values below, since these have been corrected
# to report nothing below the limit of detection (0.5uM in 2017, and 0.2uM in 2018/2019)
tdndat <- read.csv("data_raw/Bioreactor_DOC_TDN_2017_2019.csv")

# All units are currently in uM (micromolar).

# Aggregate replicates 
nh4dat_ed <- nh4dat %>% # Takes the original dataset
  group_by(Year, Site, Treatment, Sample_ID.1) %>% # groups the data
  summarize(meanNH4 = mean(Conc, na.rm = TRUE), 
            sdNH4 = sd(Conc, na.rm = TRUE)) %>% # calculations
  ungroup() %>% # don't forget it!
  mutate(Analyte = "NH4") %>% # add new column for identification
  rename(mean = meanNH4, 
         sd = sdNH4,
         Sample_ID = Sample_ID.1) %>% # rename columns for combining
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
  
  # calculate changes in concentrations (uM) and rates
  changes <- newdat %>%
    mutate(change = mean_clean - start) %>%
    mutate(change_hr = change / 3) # all bioreactors run for 3 hours
    
  # export data
  file.name <- paste0("data_analyses/fluxes/", site, "_", year, "_", 
                      analyte, ".rds") # create file name
  saveRDS(changes, file=file.name)
  
    }
  }
}

# load in and combine all the created datasets
# I chose to export the files so that I could see they mapped correctly individually
# also for loops need empty receiving dataframes in the environment,
# so it's often easier to export the files to a working data directory instead
nutdat_changes <- list.files(path = "data_analyses/fluxes", full.names = TRUE) %>%
  map(readRDS) %>%
  bind_rows()

# Check to be sure it all looks alright
View(nutdat_changes)
# woohoo!!

# Next, I'll trim down the dataset to help with flux calculations
nutdat_trim <- nutdat_changes %>%
  filter(Treatment %in% c("Control", "Experimental")) %>% # filter out "before" samples
  # since they have already served their purpose
  # select function being tough here for some reason
  dplyr::select(Year, Site, Treatment, Analyte, change_hr) # trim columns

# Need to aggregate by control/experimental treatment for each run to calculate net flux
nutdat_net <- nutdat_trim %>%
  group_by(Year, Site, Treatment, Analyte) %>%
  summarize(mean_change_hr = mean(change_hr, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Treatment, values_from = mean_change_hr) %>%
  rename(Experimental_change_hr = Experimental, Control_change_hr = Control) %>%
  # create new column with net change in uM/hour
  mutate(Net_change_uM_hr = Experimental_change_hr - Control_change_hr) %>%
  # convert net flux to umol/m^2 * hour
  mutate(Net_flux_hr_m2 = (Net_change_uM_hr * 0.25 * 10000) / 19.6)

# Size of reservoir: 250mL = 0.25 L
# Surface area of sediment core: 19.6 cm^2
# Conversion from cm^2 to m^2: 10,000 cm^2 / m^2

# calculate turnover times in NH4 concentrations
# first, select only NH4 data
NH4_only <- nutdat_net %>%
  filter(Analyte == "NH4")

# mean change in NH4 concentrations
flux_uM <- mean(NH4_only$Net_change_uM_hr, na.rm = TRUE)
# for manuscript sig figs' purposes, using 0.2uM

# mean umol in the top 2 cm of a core (14 uM)
core_umol_existing <- 14 * 0.011 # With 29% porosity, each core has 11.368mL water

# mean umol produced per hour (250mL size reservoirs)
core_umol_flux <- 0.2 * 0.25

# turnover time a.k.a. time required to replenish existing concentrations
tt_hr <- core_umol_existing/core_umol_flux

# below calcs are for discussion section calculation:
# mean umol in top 15cm of a core (23.4 uM)
# see "calc_horiz" in "table1_script.R" for existing NH4 conc
# volume of water from core volume (39.2 cm^2) * 29% porosity
# umol = uM * L
core_umol_existing15 <- 23.4 * 0.011 * 7.5

# mean umol produced per hour 
# uM * 250mL size reservoirs * 7.5 to reach 15cm of sediment
core_umol_flux15 <- 0.2 * 0.25 * 7.5

# turnover time a.k.a. time required to replenish existing concentrations
tt_hr15 <- core_umol_existing15/core_umol_flux15

tt_min15 <- tt_hr15 * 60

# REDO - only to 2cm
# mean umol in top 2cm of a core (14.0 uM)
# see "calc_horiz" in "table1_script.R" for existing NH4 conc
# volume of water from core volume (39.2 cm^3) * 29% porosity
# umol = uM * L
core_umol_existing2 <- 14.0 * 0.011

# mean umol produced per hour 
# uM * 250mL size reservoirs
core_umol_flux2 <- 0.2 * 0.25

# turnover time a.k.a. time required to replenish existing concentrations
tt_hr2 <- core_umol_existing2/core_umol_flux2

tt_min2 <- tt_hr2 * 60

# Adding water column calculations in here as well

# in 10 m water column, volume of water (m3)
v_m <- 10 * 1 * 1
# volume of water (L)
v_L <- v_m * 1000

# umoles of NH4 using mean diel sampling conc (0.7)
umol_NH4 <- v_L * 0.7
mmol_NH4 <- umol_NH4 / 1000

# contribution with only 2 cm of sediment flux
# using flux of 0.52 mmol m^-2 day^-1
contribution_2cm <- 0.52 / mmol_NH4
# 7%

# contribution with 30 cm of sediment flux
# using flux of 7.8 mmol m^-2 day^-1
contribution_30cm <- 7.8 / mmol_NH4
# 111%

# now, into a 20m water column to mimic turnover
# time in sediments calculations

# in 20 m water column, volume of water (m3)
v_m_20 <- 20 * 1 * 1
# volume of water (L)
v_L_20 <- v_m_20 * 1000

# umoles of NH4 using mean diel sampling conc (0.7)
umol_NH4_20 <- v_L_20 * 0.7
mmol_NH4_20 <- umol_NH4_20 / 1000

# contribution with only 2 cm of sediment flux
# using flux of 0.52 mmol m^-2 day^-1
contribution_2cm_20 <- 0.52 / mmol_NH4_20
# 4%

# contribution with 30 cm of sediment flux
# using flux of 7.8 mmol m^-2 day^-1
contribution_30cm_20 <- 7.8 / mmol_NH4_20
# 56%

# Summary Stats -----------------------------------------------------------

# Running some additional calculations for inclusion in the manuscript.
before <- nutdat_clean %>%
  filter(Treatment == "Before" & Analyte == "NH4")

beforeno3 <- nutdat_clean %>%
  filter(Treatment == "Before" & Analyte == "NO3")

before_tdn <- nutdat_clean %>%
  filter(Treatment == "Before" & Analyte == "TDN")

after <- nutdat_clean %>%
  filter(Treatment == "Experimental" & Analyte == "NH4")

afterno3 <- nutdat_clean %>%
  filter(Treatment == "Experimental" & Analyte == "NO3")

after_tdn <- nutdat_clean %>%
  filter(Treatment == "Experimental" & Analyte == "TDN")

summary <- nutdat_net %>%
  group_by(Analyte) %>%
  summarize(mean_all = mean(Net_flux_hr_m2, na.rm = TRUE),
            min_all = min(Net_flux_hr_m2, na.rm = TRUE),
            max_all = max(Net_flux_hr_m2, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(daily_mean = mean_all * 24, # convert mean flux to umol/m^2 * day
         daily_mean_30cm = mean_all * 24 * 15) # convert mean flux to umol/m^2 * day into 30 cm sediment depth

# Note - the results of the "summary" dataframe creation above are the ones
# cited in the results and the abstract sections in the manuscript.

# Creating another dataframe for table 3 in the manuscript.
table3 <- nutdat_net %>%
  select(Year, Site, Analyte, Net_flux_hr_m2) %>%
  pivot_wider(names_from = Analyte, values_from = Net_flux_hr_m2)

# additional calcs for results/abstract

# first, DON = TDN - NH4 - NO3
don_hr <- 79.111158 - 21.643563 - 8.487654
don_day <- 1898.6678 - 519.4455 - 203.7037

table3_plus <- table3 %>%
  mutate(NH4_Net_flux_day_m2 = NH4 * 24)

write_csv(table3, path = "data_analyses/table3_output.csv")

# And exporting dataset for use in running the linear mixed effects model
# using rate values since they are comparable 
# (uM/surface area of core * hour)
saveRDS(nutdat_trim, file="data_analyses/nutdat_bioreactors.rds")

# Also exporting net flux dataset for publication.
saveRDS(nutdat_net, file="data_analyses/nutdat_net_bioreactors.rds")

# End of script.

# Bioreactor Data Analysis
# Heili Lowman
# 10/28/22

# The following script will prepare data for publication with
# the sediment nutrient recycling manuscript.

# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse) 
library(lubridate)
library(here)

# Load necessary datasets
conc_hrly <- readRDS("data_analyses/nutdat_bioreactors.rds")
net_fluxes <- readRDS("data_analyses/nutdat_net_bioreactors.rds")

# Format -------------------------------------------------------------------

View(conc_hrly)

# All I really need to do is update a column name to better describe the units of change.
conc_hrly <- conc_hrly %>%
  rename("Change_uM_hr" = "change_hr")

# Also round all values to two decimal places.
conc_hrly$Change_uM_hr <- round(conc_hrly$Change_uM_hr, digits = 2)

View(net_fluxes)

# Going to remove the mean control and experimental change columns and rename the remaining net 
# columns that describe the net efflux from sediments.
net_fluxes <- net_fluxes %>%
  select(-Control_change_hr) %>%
  select(-Experimental_change_hr) %>%
  rename("Net_Exp_Change_uM_hr" = "Net_change_uM_hr",
         "Net_Exp_Flux_umol_m2_hr" = "Net_flux_hr_m2")

# Also round all values to two decimal places.
net_fluxes$Net_Exp_Change_uM_hr <- round(net_fluxes$Net_Exp_Change_uM_hr, digits = 2)
net_fluxes$Net_Exp_Flux_umol_m2_hr <- round(net_fluxes$Net_Exp_Flux_umol_m2_hr, digits = 2)

# Export -------------------------------------------------------------------

write_csv(conc_hrly, "data_tidy/bioreactor_nutrient_conc.csv")
write_csv(net_fluxes, "data_tidy/bioreactor_nutrient_fluxes.csv")
  
# End of script.
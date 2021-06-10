# Tidbit Data Processing & Analysis
# Heili Lowman
# 6/8/21

# The following script will analyze temperature data collected
# by the array of tidbit temperature sensors deployed at Mohawk
# Reef for the sediment nutrient recycling manuscript.

# Load packages
library(tidyverse) 
library(lubridate)
library(plotly)

# Load necessary dataset
tidbit <- read_csv("data_raw/mohawk_sedimentTemp_20180319.csv")

# Create less cluttered dataset
tidbit_tidy <- tidbit %>% 
  select(Date_Time_B5, Temp_A50, Temp_B5, Temp_B15, Temp_B30) %>% 
  mutate(DATEtime = mdy_hm(Date_Time_B5))
# Removed B_45 depth from analysis because the noise in the data,
# based on previous analyses, suggests there may have been a mislabeling/
# misentry of this tidbit's data.

# And let's examine where we need to lop off data to clean out the deployment and removal data.

fig_test <- plot_ly(tidbit_tidy, x = tidbit_tidy$DATEtime, y = tidbit_tidy$Temp_B30, type = 'scatter', mode = 'lines')

fig_test

# It appears it took ~48 hours for things to equilibrate at depth and 
# then were disturbed in the afternoon of September 1 for retrieval, 
# so I'm going to remove everything prior to AUGUST 10 0:00 and after 
# SEPTEMBER 1 12:00.

# Remove desired rows
tidbit_tidier <- tidbit_tidy[-c(1:312, 6794:6883),]

# Check to be sure equilibration time is out now
fig_tidier <- plot_ly(tidbit_tidier, x = tidbit_tidier$DATEtime, y = tidbit_tidier$Temp_B30, type = 'scatter', mode = 'lines')

fig_tidier

# 
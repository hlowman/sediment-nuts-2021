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
library(viridis)
library(patchwork)

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

# Perform cross-covariance analysis to determine flushing rates.
# Using the cross-covariance technique from Figure 10 of Fram et al. 
# 2014, "Miniature thermistor chain for determining surficial sediment
# porewater advection," but using the ccf() base R function instead.

# Additional resources referenced:
# https://nwfsc-timeseries.github.io/atsa-labs/sec-tslab-correlation-within-and-among-time-series.html
# http://r-statistics.co/Time-Series-Analysis-With-R.html

# Prior to any cross-covariance calculations, MUST ensure the time 
# series stamps match up, which is not a problem, since all data
# were collected at 5 minute increments here. Keep in mind, this
# will also be the timestep used in analysis calculations.

# Sediment porewaterwater temperature (y) depends on surface water
# temperature (x). So, in each case x = Temp_A50 or "surface".

# Create all horizon vector variables
surface <- tidbit_tidier$Temp_A50
sed5 <- tidbit_tidier$Temp_B5
sed15 <- tidbit_tidier$Temp_B15
sed30 <- tidbit_tidier$Temp_B30

# Calculate cross-covariance at all 3 sediment depths
cc5 <- ccf(surface, sed5, lag.max = NULL, type = c("covariance"))
cc15 <- ccf(surface, sed15, lag.max = NULL, type = c("covariance"))
cc30 <- ccf(surface, sed30, lag.max = 500, type = c("covariance"))
# Changed max lag of cc30 plot to display maximum

# Create new dataframes to examine maximum acf values/timestep lags
# for each sediment horizon.
cc5dat <- data.frame(lag5=cc5$lag, acf5=cc5$acf) 
# max @ lag = -2 == 10 minutes
cc15dat <- data.frame(lag15=cc15$lag, acf15=cc15$acf) 
# max @ lag = -20 == 100 minutes (~1.5 hours)
cc30dat <- data.frame(lag30=cc30$lag, acf30=cc30$acf) 
# max @ lag = -369 == 1845 minutes (~1.25 days)

# Tidy dataset for easier plotting.
tidbit_stacked <- tidbit_tidier %>%
  select(!Date_Time_B5) %>%
  pivot_longer(!DATEtime, names_to = "Horizon", values_to = "Temperature") %>%
  mutate(horizon_f = factor(case_when(Horizon == "Temp_A50" ~ "50cm Water",
                                      Horizon == "Temp_B5" ~ "5cm Sediment",
                                      Horizon == "Temp_B15" ~ "15cm Sediment",
                                      Horizon == "Temp_B30" ~ "30cm Sediment"),
                            levels = c("50cm Water", "5cm Sediment", "15cm Sediment", "30cm Sediment"))) %>% # added leveled/factored horizon column
  mutate(date_ed = ymd_hms(DATEtime)) # and edited date format

# Create main timeseries figure for manuscript.
dev.off() # Had to run to clear out previous figures.
fig3a <- ggplot(tidbit_stacked, aes(DATEtime, Temperature)) +
  geom_line(aes(color = horizon_f)) +
  labs(y = 'Temperature (ºC)', 
       x = 'Date',
       color = 'Horizon') +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14)) 

fig3a

# Check underlying hexadecimal colors so these can be replicated for
# cross-covariance plots.
scales::viridis_pal()(4)

# Create cross-covariance inset plots.
# 5 cm sediment depth
fig3b <- ggplot(cc5dat, aes(x = lag5, y = acf5)) +
  geom_bar(stat = "identity", fill = "#31688EFF") +
  labs(x = "Lag",
       y = "Cross-covariance") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14)) 

fig3b

# 15 cm sediment depth
fig3c <- ggplot(cc15dat, aes(x = lag15, y = acf15)) +
  geom_bar(stat = "identity", fill = "#35B779FF") +
  labs(x = "Lag",
       y = "Cross-covariance") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14)) 

fig3c

# 30 cm sediment depth
fig3d <- ggplot(cc30dat, aes(x = lag30, y = acf30)) +
  geom_bar(stat = "identity", fill = "#FDE725FF") +
  xlim(c(-400, 400)) +
  labs(x = "Lag",
       y = "Cross-covariance") +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14)) 

fig3d

# Combine into full figure.
fig3_full <- fig3a /
  (fig3b | fig3c | fig3d)

fig3_full + plot_annotation(tag_levels = 'A')

# Export figure for manuscript.
# ggsave(("Figure_3.tiff"),
#        path = "/Users/heililowman/Desktop/R_figures/Sediment_N",
#        width = 25,
#        height = 12.5,
#        units = "cm"
#        )

# End of script.

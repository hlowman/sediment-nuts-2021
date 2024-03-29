# Diel Data Analysis
# Heili Lowman
# 6/28/21

# The following script will analyze diel sampling data collected
# in 2018 for the sediment nutrient recycling manuscript.

# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse) 
library(lubridate)
library(viridis)
library(patchwork)
library(here)

# Load necessary datasets

# This data was collected by Heili Lowman et al. and can be downloaded
# from the SBC LTER's data repository at the following URL:
# https://doi.org/10.6073/pasta/7a41a41c52b425c564de7aa079839049
nuts_dat <- read_csv("data_raw/Diel_nutrients.csv")
ctd_dat <- read_csv("data_raw/Diel_CTD.csv")

# Tide data published by NOAA/NOS/CO-OPS
# Disclaimer: These data are based upon the latest information available as of the date of your request, and may differ from the published tide tables.
# Daily Tide Predictions
# StationName: Santa Barbara; State: CA; Stationid: 9411340
# Prediction Type: Harmonic
# From: 20180803 00:00 - 20180812 23:00
# Units: Metric; Time Zone: LST_LDT; Datum: MLLW
# Interval Type: Hourly
# As published at: https://tidesandcurrents.noaa.gov/waterlevels.html?id=9411340&units=metric&bdate=20180803&edate=20180812&timezone=LST/LDT&datum=MLLW&interval=h&action=data

tides_dat <- read_csv("data_raw/bound_tides_pub.csv")

# Additional Summary Statistics -------------------------------------------

# Calculations for inclusion in the manuscript

# mean temperature
mean(ctd_dat$Temperature)

# minimum temperature
min(ctd_dat$Temperature)

# maximum temperature
max(ctd_dat$Temperature)

# mean NH4
mean(nuts_dat$NH4)

# standard deviation NH4
sd(nuts_dat$NH4)

# minimum NH4
min(nuts_dat$NH4)

# maximum NH4
max(nuts_dat$NH4)

# mean NO3
mean(nuts_dat$NO3_NO2)

# standard deviation NO3
sd(nuts_dat$NO3_NO2)

# minimum NO3
min(nuts_dat$NO3_NO2)

# maximum NO3
max(nuts_dat$NO3_NO2)

# mean chla
mean(nuts_dat$Chla, na.rm=TRUE)

# standard deviation chla
sd(nuts_dat$Chla, na.rm=TRUE)

# minimum chla
min(nuts_dat$Chla, na.rm=TRUE)

# maximum chla
max(nuts_dat$Chla, na.rm=TRUE)

# Figures -----------------------------------------------------------------

# Tidal paneled figure

trial_names <- list('A'="August 3-4", 'B'="August 7-8", 'C'="August 10-11")
trial_labeller <- function(variable,value){
  return(trial_names[value])
}

dailytides <- ggplot(data = tides_dat, aes(x = Hour, y = Verified_m)) +
  geom_point(size = 0.5, alpha = 0.7) +
  geom_line() +
  labs(y = 'Tidal Height (m)') +
  facet_wrap(~SampleDate, labeller = trial_labeller) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        #strip.text.x = element_blank(),
        legend.position = "none")

dailytides

# Ammonium paneled figure

fig_nh4 <- nuts_dat %>%
  mutate(Depth_f = factor(Depth),
         Trial_f = factor(case_when(DateTime < "2018-08-05 00:00:00" ~ 1,
                                    DateTime > "2018-08-05 00:00:00" & 
                                      DateTime < "2018-08-09 00:00:00" ~ 2,
                                    DateTime > "2018-08-09 00:00:00" ~ 3))) %>%
  ggplot(aes(DateTime, NH4, color=Depth_f)) +
  scale_color_viridis(discrete=TRUE) +
  geom_point(size = 3, stroke = 0.75, alpha = 0.9) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "8 hours") +
  ylab(expression(paste(NH[4]^{"+"}," (μM)"))) +
  labs(color = "Water Depth (m)") +  
  theme_bw() + 
  facet_grid(.~Trial_f, scales="free") +
  theme(text=element_text(family="Times New Roman", size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "right")

fig_nh4

# Nitrate paneled figure

fig_no3 <- nuts_dat %>%
  mutate(Depth_f = factor(Depth),
         Trial_f = factor(case_when(DateTime < "2018-08-05 00:00:00" ~ 1,
                                    DateTime > "2018-08-05 00:00:00" & 
                                      DateTime < "2018-08-09 00:00:00" ~ 2,
                                    DateTime > "2018-08-09 00:00:00" ~ 3))) %>%
  ggplot(aes(DateTime, NO3_NO2, color=Depth_f)) +
  scale_color_viridis(discrete=TRUE) +
  geom_point(size = 3, stroke = 0.75, alpha = 0.9) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "8 hours") +
  ylab(expression(paste(NO[3]^{"-"}, "+", NO[2]^{"-"}, " (μM)"))) +
  labs(color = "Water Depth (m)") +  
  theme_bw() + 
  facet_grid(.~Trial_f, scales="free") +
  theme(text=element_text(family="Times New Roman", size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")

fig_no3

# Chlorophyll paneled figure

fig_chla <- nuts_dat %>%
  mutate(Depth_f = factor(Depth),
         Trial_f = factor(case_when(DateTime < "2018-08-05 00:00:00" ~ 1,
                                    DateTime > "2018-08-05 00:00:00" & 
                                      DateTime < "2018-08-09 00:00:00" ~ 2,
                                    DateTime > "2018-08-09 00:00:00" ~ 3))) %>%
  ggplot(aes(DateTime, Chla, color=Depth_f)) +
  scale_color_viridis(discrete=TRUE) +
  geom_point(size = 3, stroke = 0.75, alpha = 0.9) +
  scale_x_datetime(date_labels = "%H:%M", date_breaks = "8 hours") +
  ylab(expression(paste(italic("chl a"), " (μg/L)"))) +
  labs(color = "Water Depth (m)",
       x = "Sampling Time") +  
  theme_bw() + 
  facet_grid(.~Trial_f, scales="free") +
  theme(text=element_text(family="Times New Roman", size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        legend.position = "none")

fig_chla

# Temperature paneled figure

dailytemp <- ctd_dat %>%
  mutate(Trial_f = factor(case_when(DateTime < "2018-08-05 00:00:00" ~ 1,
                                    DateTime > "2018-08-05 00:00:00" & 
                                      DateTime < "2018-08-09 00:00:00" ~ 2,
                                    DateTime > "2018-08-09 00:00:00" ~ 3))) %>%
  # group by trial, date, and bin data by 0.1 m intervals
  group_by(Trial_f, DateTime, cut(Depth, breaks=seq(0, 7.3, by = 0.1) )) %>%
  mutate(Temp_binned = mean(Temperature, na.rm = TRUE)) %>%
  ungroup() %>%
  rename(DepthRange = `cut(Depth, breaks = seq(0, 7.3, by = 0.1))`) %>%
  # need to pull out the depth value from the range to make contouring work
  separate(col = DepthRange, into = c("start", "end"), sep = ",") %>%
  mutate(DepthRange_start = parse_number(start),
         DepthRange_end = parse_number(end)) %>%
  # and finally need to group so there's only one time-depth measurement
  group_by(Trial_f, DateTime, DepthRange_start) %>%
  summarize(Temp_binned_grouped = mean(Temp_binned)) %>%
  ungroup() %>%
  ggplot(aes(x = DateTime, y = DepthRange_start, z = Temp_binned_grouped)) + # specify aesthetics
  scale_y_reverse() + # flip the y axis so values increase with depth
  geom_contour_filled() + # create contour plot
  guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) + # edit the legend not to display bins
  geom_point(shape = 20, color = "white", size=0.1, stroke=0.2) + # adds vertical lines for casts
  labs(y = "Water Depth (m)",
       fill = "Temperature (ºC)") +
  facet_grid(~Trial_f, scales = "free") + # scales free in order to only display a given day's data
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

dailytemp

# Combine and export figure -----------------------------------------------

full_fig_4 <- dailytides /
  dailytemp /
  fig_nh4 /
  fig_no3 /
  fig_chla

full_fig_4

# Export map to desktop.
# ggsave(("Figure_4.tiff"),
#        path = "/Users/heililowman/Desktop/R_figures/Sediment_N",
#        width = 25,
#        height = 25,
#        units = "cm"
#        )

# End of script

# Diel Data Analysis, Pt. II
# Heili Lowman
# 6/29/21

# The following script will analyze diel sampling data collected
# in 2018 for the sediment nutrient recycling manuscript using a
# linear mixed effects model.

# I am created a linear mixed effects model (to predict ammonium concentrations) since repeated sampling took place. Model creation begins with fixed effects and random effects using a random intercept structure. Then, model selection follows the protocol outlined by Zuur et al. (2009, Chapter 5), beginning with a linear model, accounting for variance structure, optimizing the fixed structure, and validating best model fit using distribution of residuals and AIC values.

# Based on past analyses, this is the model structure I will investigate:
# NH4 ~ tide*temp + chl a + 1|depth

# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse) 
library(lubridate)
library(here)
library(nlme)
library(multcomp)
library(GGally)

# Load data
nh4_v_lmem <- read_csv("data_raw/Diel_NH4_v_Temp_and_chla.csv") 

# Tidy/rename data
model_nh4_data <- nh4_v_lmem %>%
  rename(NH4 = meanNH4,
         tide = tide_pub,
         temperature = Temp,
         chlorophyll = chl_a) %>%
  mutate(depthf = factor(Depth, levels = c(1, 3, 5))) %>%
  dplyr::select(NH4, depthf, tide, temperature, chlorophyll) %>%
  na.omit() # Removes NAs which cannot be processed by models below.

# Verify structure of data - all numerical except depth
str(model_nh4_data)

# I am not removing any outliers of the dataset.

# Data Exploration --------------------------------------------------------

# Creating boxplots and pairs plots for basic data exploration.
ggpairs(model_nh4_data)

# Trying to determine whether or not the tide and temperature predictors should be combined due to collinearity. Given the regional literature, and my initially proposed structure above, I am going to include tide*temperature as a fixed effect to account for this collinearity. Chlorophyll will remain on its own as a fixed effect.

# Creating histograms for additional data exploration.
hist(model_nh4_data$NH4) # looks skewed left (towards zero)
hist(model_nh4_data$tide) # looks evenly distributed about ~1
hist(model_nh4_data$temperature) # looks evenly distributed about ~20
hist(model_nh4_data$chlorophyll) # looks skewed left (towards zero)

# Transforming NH4 and chl a since they appear less than normally distributed.
model_nh4_data <- model_nh4_data %>%
  mutate(logNH4 = log10(NH4)) %>%
  mutate(logchl = log10(chlorophyll))

hist(model_nh4_data$logNH4) # better-ish
hist(model_nh4_data$logchl) # better

# Per Chapter 19 of Zuur et al. 2009, I'm going to log-transform because otherwise (using different variance structures) you are giving up lots of degrees of freedom and make GLS estimation unstable.

op<- par(mfrow = c(1, 2), mar = c(3, 4, 1, 1))
dotchart(model_nh4_data$NH4, groups = model_nh4_data$depthf) # untransformed on left
dotchart(model_nh4_data$logNH4, groups = model_nh4_data$depthf) # transformed on right
par(op) # looking better!!

# Going into the model workflow, I am using the following structure to start:
# log(NH4) ~ tide*temp + log(chla) + 1|depth

# LME Model Workflow ------------------------------------------------------

# Step 1: Create a linear regression and check residuals.

m1 <- lm(logNH4 ~ tide*temperature + logchl , data = model_nh4_data) # linear model
rm1 <- rstandard(m1) # Assigns standardized residuals to rm1.

m1_data <- model_nh4_data %>%
  mutate(residuals = rm1) # Adds to original dataset.

# Plotting to be sure residuals are evenly distributed about zero.
ggplot(m1_data, aes(x = tide, y = residuals)) +
  geom_point()
ggplot(m1_data, aes(x = temperature, y = residuals)) +
  geom_point()
ggplot(m1_data, aes(x = logchl, y = residuals)) +
  geom_point()
# no patterns emerge

# Step 2: Fit the lm() with GLS and compare to lme().

m2 <- gls(logNH4 ~ tide*temperature + logchl, data = model_nh4_data) # This is effectively a linear regression since it has no additional calls.
m3 <- lme(logNH4 ~ tide*temperature + logchl,
          random = ~1 | depthf,
          data = model_nh4_data) # Creates first LMEM with random effect term.
anova(m2, m3) # Compares the two models. m2 preferred with AIC value of 98.78326, but I choose to proceed with m3 knowing there needs to be random term to account for repeated sampling.

# Step 3: Decide on a variance structure.

plot(m3, col=1) # Plots the residuals of the first lme model. Faint slimming pattern?
qqnorm(m3) # This looks good.

# Based on scatterplots above, let's see what a variance term by temperature does.
m4 <- lme(logNH4 ~ tide*temperature + logchl, 
          random = ~1 | depthf, 
          data = model_nh4_data,
          weights = varIdent(form = ~1 | temperature))

anova(m3, m4) # Compares the two models. m3 still preferred with AIC value of 100.6514. So, removing the variance term.

# Step 4: Fit the lme().

# m3 <- lme(logNH4 ~ tide*temperature + chlorophyll,
#          random = ~1 | depthf,
#          data = model_nh4_data)

# Step 5: Compare the lm() and lme().

# See Step 2.

# Step 6: Check residuals.

# See Step 3.

# Step 7/8: Step-wise Optimal Fixed Structure.

m5_ml <- lme(logNH4 ~ tide*temperature + logchl,
             random = ~1 | depthf,
             method = "ML", 
             data = model_nh4_data) # Full model structure
m5_sub1 <- lme(logNH4 ~ tide*temperature,
               random = ~1 | depthf,
               method = "ML", 
               data = model_nh4_data) # Removes chlorophyll
m5_sub2 <- lme(logNH4 ~ temperature + logchl,
               random = ~1 | depthf,
               method = "ML", 
               data = model_nh4_data) # Removes tide
m5_sub3 <- lme(logNH4 ~ tide + logchl,
               random = ~1 | depthf,
               method = "ML", 
               data = model_nh4_data) # Removes temperature
anova(m5_ml, m5_sub1, m5_sub2, m5_sub3) # Compares all models. m5_ml preferred with AIC value of 82.19151, so keep remaining terms.

# Step 9: Refit with REML

mfinal <- lme(logNH4 ~ tide*temperature + logchl,
              random = ~1 | depthf,
              method = "REML", 
              data = model_nh4_data)

mfinal # Coefficients.
summary(mfinal) # Output of the full model.
plot(mfinal, col=1) # Checking residuals.
qqnorm(mfinal) # Quantile-quantile plot.
anova(mfinal) # Output for manuscript.

# Step 10: What this means in words.

# I applied a linear mixed effects modeling approach because of the repeated sampling (and lack of independence between samples). My model suggests there is a significant effect of tide, temperature, and chlorophyll on summer NH4 concentrations in Goleta Bay. A random intercept by depth was included.

# The final model took the form: log(NH4) ~ tide * temperature + chlorophyll + 1|depth
# This translated to a formula of logNH4 = 0.10[tide*temperature] - 2.43[tide] - 0.24[temperature] - 0.79[chlorophyll] + random

# End of script.

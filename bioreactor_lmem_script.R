# Bioreactor Data Analysis, Pt. II
# Heili Lowman
# 6/22/21

# The following script will run a series of linear mixed effects
# models for the bioreactor data in preparation for the sediment
# nutrient efflux manuscript.

# Models will follow the format:
# Nutrient ~ Treatment + Site + Year + 1/(Site)


# Setup -------------------------------------------------------------------

# Load packages
library(tidyverse) 
library(nlme)
library(multcomp)
library(viridis)
library(patchwork)
library(here)
library(GGally)

# Load necessary dataset
dat <- readRDS("data_analyses/nutdat_bioreactors.rds") %>%
  mutate(Site = factor(Site),
         Treatment = factor(Treatment))

# Make separate datasets for three analytes
dat_nh4 <- dat %>%
  filter(Analyte == "NH4") %>%
  filter(!is.na(change_hr)) # remove all NAs
dat_no3 <- dat %>%
  filter(Analyte == "NO3") %>%
  filter(!is.na(change_hr))
dat_tdn <- dat %>%
  filter(Analyte == "TDN") %>%
  filter(!is.na(change_hr)) %>%
  filter(Year != 2017) # since no samples were run for TDN in 2017

# And briefly examine all three
hist(dat_nh4$change_hr)
hist(dat_no3$change_hr)
hist(dat_tdn$change_hr)

# Going to standardize using means and sd
mean_nh4 <- mean(dat_nh4$change_hr, na.rm=TRUE)
sd_nh4 <- sd(dat_nh4$change_hr, na.rm=TRUE)
dat_nh4$std_change_hr <- (dat_nh4$change_hr - mean_nh4)/sd_nh4

mean_no3 <- mean(dat_no3$change_hr, na.rm=TRUE)
sd_no3 <- sd(dat_no3$change_hr, na.rm=TRUE)
dat_no3$std_change_hr <- (dat_no3$change_hr - mean_no3)/sd_no3

mean_tdn <- mean(dat_tdn$change_hr, na.rm=TRUE)
sd_tdn <- sd(dat_tdn$change_hr, na.rm=TRUE)
dat_tdn$std_change_hr <- (dat_tdn$change_hr - mean_tdn)/sd_tdn

# Examine the standardized values
hist(dat_nh4$std_change_hr)
hist(dat_no3$std_change_hr)
hist(dat_tdn$std_change_hr)

# NH4 Model ---------------------------------------------------------------

### Data Exploration

# I'm examining Treatment, Year, and Site as factors influencing nutrient results. Here are the boxplots and pairplots necessary for initial inspection of the data.
boxplot(std_change_hr ~ Treatment, data = dat_nh4) # Boxplot by treatment - appears to be an effect, as we expect.
boxplot(std_change_hr ~ Site, data = dat_nh4) # Boxplot by site - as expected, little effect.
boxplot(std_change_hr ~ Year, data = dat_nh4) # Boxplot by year - appears to be most variable in 2017 (with prior bioreactor build).

#### STEP 1: Create a linear regression and check residuals.

# Response variable - change in NH4 concentrations
# Explanatory (fixed) variables - Treatment, Site, Year

# Create the initial linear model
a1 <- lm(std_change_hr ~ Treatment + Year + Site, data = dat_nh4) 
# Pull out and examine standardized residuals
ra1 <- rstandard(a1)

dat_nh4 <- dat_nh4 %>%
  mutate(ra1 = ra1)

# Plots said residuals.
ggplot(data = dat_nh4, aes(x = Treatment, y = ra1)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Treatment", y = "Standardised residuals") 
# very different

ggplot(data = dat_nh4, aes(x = Site, y = ra1)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Site", y = "Standardised residuals")

ggplot(data = dat_nh4, aes(x = Year, y = ra1, group = Year)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Year", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().

a2 <- gls(std_change_hr ~ Treatment + Year + Site, data = dat_nh4)
a3 <- lme(std_change_hr ~ Treatment + Year + Site,
          random =~1 | Site, # repeated sampling at each site
          data = dat_nh4)
anova(a2, a3) # Compares the two models.

#### STEP 3: Decide on a variance structure (aka random terms).

plot(a3, col=1) # Check the residuals.
qqnorm(a3) # This looks terrible even with transformed data.

# Adding in variance structure based on residuals by treatment
a4 <- lme(std_change_hr ~ Treatment + Year + Site,
          random =~1 | Site,
          data = dat_nh4,
          weights = varIdent(form =~1 | Treatment))

plot(a4, col=1) # Check the residuals.
qqnorm(a4) # Looking much better.
anova(a3, a4) # Compares the two most recent models.

#### STEP 4: Fit the lme().

# Using a4 above.

#### STEP 5: Compare the lm() and lme().

anova(a2, a4) # Compares the initial lm and the recent lmem, with additional variance term. a4 outperforms a2 using AIC.

#### STEP 6: Everything ok? Check residuals.

# See Step 3.

#### STEP 7/8: Step-wise Optimal Fixed Structure

# I want to retain all effects to generate outputs for each.

#### STEP 9: Refit with REML
afinal <- lme(std_change_hr ~ Treatment + Year + Site,
              random =~1 | Site,
              weights = varIdent(form =~1 | Treatment), 
              method = "REML", 
              data = dat_nh4)

# Output of the model.
summary(afinal)

# Final results.
anova(afinal)

#### STEP 10: What does this mean in WORDS?
# I applied a linear mixed effect modeling approach because the data are nested with samples taken from each site on multiple occasions. My model suggests there is a significant effect of treatment type (experimental/control) on NH4. Random intercepts by site and a variance term by treatment were added.

# Equation: 
# std(NH4) = 8.03 + 0.83[Exp] - 0.004[Year] 
# - 0.07[GOLB] + 0.31[MICR] + random + variance

# BONUS POST HOC:
aHSD <- glht(afinal, linfct=mcp(Treatment="Tukey")) # Run a Tukey's post hoc analysis.
summary(aHSD) # Exp & Control significantly different from one another. 
# Note : Use the results of the first run on this function.

# NO3 Model ---------------------------------------------------------------

### Data Exploration

# Boxplots and pairplots for initial inspection of the data.
boxplot(std_change_hr ~ Treatment, data = dat_no3) # Boxplot by treatment - appears to be an effect, as we expect.
boxplot(std_change_hr ~ Site, data = dat_no3) # Boxplot by site - some effect?
boxplot(std_change_hr ~ Year, data = dat_no3) # Boxplot by year - appears to be most variable in 2019.

#### STEP 1: Create a linear regression and check residuals.

# Response variable - change in NO3 concentrations
# Explanatory (fixed) variables - Treatment, Site, Year

# Create the initial linear model
b1 <- lm(std_change_hr ~ Treatment + Year + Site, data = dat_no3) 
# Pull out and examine standardized residuals
rb1 <- rstandard(b1)

dat_no3 <- dat_no3 %>%
  mutate(rb1 = rb1)

# Plots said residuals.
ggplot(data = dat_no3, aes(x = Treatment, y = rb1)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Treatment", y = "Standardised residuals") 
# very different again

ggplot(data = dat_no3, aes(x = Site, y = rb1)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Site", y = "Standardised residuals")

ggplot(data = dat_no3, aes(x = Year, y = rb1, group = Year)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Year", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().

b2 <- gls(std_change_hr ~ Treatment + Year + Site, data = dat_no3)
b3 <- lme(std_change_hr ~ Treatment + Year + Site,
          random =~1 | Site, # repeated sampling at each site
          data = dat_no3)
anova(b2, b3) # Compares the two models.

#### STEP 3: Decide on a variance structure (aka random terms).

plot(b3, col=1) # Check the residuals.
qqnorm(b3) # This looks better than the last, but let's see what
# an added variance term might do.

# Adding in variance structure based on residuals by treatment
b4 <- lme(std_change_hr ~ Treatment + Year + Site,
          random =~1 | Site,
          data = dat_no3,
          weights = varIdent(form =~1 | Treatment))

plot(b4, col=1) # Looks slightly worse.
qqnorm(b4) # Looks similar.
anova(b3, b4) # Compares the two most recent models.

#### STEP 4: Fit the lme().

# Using b3 above.

#### STEP 5: Compare the lm() and lme().

# Using b3 above.

#### STEP 6: Everything ok? Check residuals.

# See Step 3.

#### STEP 7/8: Step-wise Optimal Fixed Structure

# I want to retain all fixed effects.

#### STEP 9: Refit with REML
bfinal <- lme(std_change_hr ~ Treatment + Year + Site,
              random =~1 | Site,
              method = "REML", 
              data = dat_no3)

# Output of the model.
summary(bfinal)

# Final results.
anova(bfinal)

#### STEP 10: What does this mean in WORDS?
# My model suggests there is a significant effect of year on NO3. Random intercepts by sitewere added.

# Equation: 
# std(NO3) = -821.57 + 0.43[Exp] + 0.41[Year] 
# + 0.28[GOLB] + 0.71[MICR] + random

# TDN Model ---------------------------------------------------------------

### Data Exploration

# Boxplots and pairplots for initial inspection of the data.
boxplot(std_change_hr ~ Treatment, data = dat_tdn) # Boxplot by treatment - appears to be an effect, as we expect.
boxplot(std_change_hr ~ Site, data = dat_tdn) # Boxplot by site - no effect.
boxplot(std_change_hr ~ Year, data = dat_tdn) # Boxplot by year - no effect.

#### STEP 1: Create a linear regression and check residuals.

# Response variable - change in TDN concentrations
# Explanatory (fixed) variables - Treatment, Site, Year

# Create the initial linear model
c1 <- lm(std_change_hr ~ Treatment + Year + Site, data = dat_tdn) 
# Pull out and examine standardized residuals
rc1 <- rstandard(c1)

dat_tdn <- dat_tdn %>%
  mutate(rc1 = rc1)

# Plots said residuals.
ggplot(data = dat_tdn, aes(x = Treatment, y = rc1)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Treatment", y = "Standardised residuals") 
# very different again

ggplot(data = dat_tdn, aes(x = Site, y = rc1)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Site", y = "Standardised residuals")

ggplot(data = dat_tdn, aes(x = Year, y = rc1, group = Year)) +
  geom_boxplot() + 
  geom_point() +
  labs(x = "Year", y = "Standardised residuals")

#### STEP 2: Fit the lm() with GLS and compare to lme().

c2 <- gls(std_change_hr ~ Treatment + Year + Site, data = dat_tdn)
c3 <- lme(std_change_hr ~ Treatment + Year + Site,
          random =~1 | Site, # repeated sampling at each site
          data = dat_tdn)
anova(c2, c3) # Compares the two models.

#### STEP 3: Decide on a variance structure (aka random terms).

plot(c3, col=1) # Check the residuals. Eek!
qqnorm(c3) # This looks better than the last, but let's see what
# an added variance term might do.

# Adding in variance structure based on residuals by treatment
c4 <- lme(std_change_hr ~ Treatment + Year + Site,
          random =~1 | Site,
          data = dat_tdn,
          weights = varIdent(form =~1 | Treatment))

plot(c4, col=1) # Looks slightly better, although gap in the middle persists.
qqnorm(c4) # Looks similar.
anova(c3, c4) # Compares the two most recent models.

#### STEP 4: Fit the lme().

# Using c4 above.

#### STEP 5: Compare the lm() and lme().

# Using c4 above.

#### STEP 6: Everything ok? Check residuals.

# See Step 3.

#### STEP 7/8: Step-wise Optimal Fixed Structure

# I want to retain all fixed effects.

#### STEP 9: Refit with REML
cfinal <- lme(std_change_hr ~ Treatment + Year + Site,
              random =~1 | Site,
              weights = varIdent(form =~1 | Treatment),
              method = "REML", 
              data = dat_tdn)

# Output of the model.
summary(cfinal)

# Final results.
anova(cfinal)

#### STEP 10: What does this mean in WORDS?
# My model suggests there is a significant effect of treatment on TDN. Random intercepts by site and a variance term by treatment were added.

# Equation: 
# std(TDN) = -356.05 + 1.45[Exp] + 0.18[Year] 
# - 0.22[GOLB] - 0.16[MICR] + random + variance

# BONUS POST HOC:
cHSD <- glht(cfinal, linfct=mcp(Treatment="Tukey")) # Tukey's post hoc
summary(cHSD) # Exp & Control significantly different from one another. 

# End of script.

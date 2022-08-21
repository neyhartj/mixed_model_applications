## Analysis of a series of variety trials in perennial crops
##
## from Piepho and Eckl 2013 doi:10.1111/gfs.12054
##

# Load packages
library(sommer)
library(tidyverse)
library(nlme)
library(glmmTMB)
library(lme4)
library(asreml)


# Set directories
setwd("MixedModels/Piepho2013_SeriesTrialAnalysis/")

# Read in the data
pheno_dat <- read_csv(file = "piepho2013_lolium_data.csv") %>%
  # Rename columns like the notation in the paper
  select(yld, GEN, LOC, trl, blk, PLT = PLOT, yr, har) %>%
  rename_all(toupper) %>%
  mutate_at(vars(-YLD), as.factor)

# Set asreml options
asreml.options(workspace = "500mb", ai.sing = TRUE, maxit = 30)


# One location, one trial -------------------------------------------------

pheno_dat_use_142 <- pheno_dat %>%
  filter(LOC == "142") %>%
  droplevels()


pheno_dat_use1 <- pheno_dat_use_142 %>%
  filter(TRL == "2003") %>%
  droplevels()

# Compound symmetry

# First model, fixed gen, random plot and block
fit1_lmer <- lmer(YLD ~ GEN + HAR + GEN:HAR + (1|BLK) + (1|BLK:HAR),
                  data = pheno_dat_use1)

# This should yield the same variance estimates
fit1_asreml <- asreml(fixed = YLD ~ GEN + HAR + GEN:HAR,
                      # random = ~ BLK:cor(HAR) + BLK:PLT:cor(HAR),
                      random = ~ BLK:cor(HAR),
                      residual = ~ units,
                      data = pheno_dat_use1)

# Estimate varcor of BLK across harvest years
summary(fit1_asreml)$varcomp

as.data.frame(VarCorr(fit1_lmer))

# Asreml:
# Variance within blks
# 34.6946902
# Correlation between harvest years
# 0.8733901
# Covariance between harvest years
# 0.8733901 = cov / 34.6946902
# cov = 0.8733901 * 34.6946902
# cov = 30.302
#
# lmer
# variance within blks
# 30.30173 + 4.39266 = 34.69439
# covariance between harvest years
# 30.30173


# AR1

# This should yield the same variance estimates
fit2_asreml <- asreml(fixed = YLD ~ GEN + HAR + GEN:HAR,
                      # random = ~ BLK:cor(HAR) + BLK:PLT:cor(HAR),
                      random = ~ BLK:ar1(HAR),
                      residual = ~ units,
                      data = pheno_dat_use1)

summary(fit2_asreml)$varcomp


# Several locations, one trial --------------------------------------------

pheno_dat_use2 <- pheno_dat %>%
  filter(TRL == "2003") %>%
  droplevels() %>%
  # Create block by plot interaction variable
  mutate(BLK_PLT = interaction(BLK, PLT, sep = ":"))


# First model, fixed gen, random plot and block
fit2 <- lmer(YLD ~ GEN + HAR + GEN:HAR + (1|LOC) + (1|LOC:HAR) + (1|LOC:GEN) + (1|LOC:GEN:HAR) +
               (1|LOC:BLK) + (1|LOC:BLK:HAR) + (1|LOC:BLK:PLT) + (1|LOC:BLK:PLT:HAR),
             data = pheno_dat_use2,
             control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))

as.data.frame(VarCorr(fit2))

# Compare to Table 3 in Piepho 2013
# LOC:HAR variance = 165.64728
# LOC:HAR covariance = 105.05504
#
# LOC:GEN:HAR variance = 22.60766
# LOC:GEN:HAR covariance = 12.62059
#


# ASREML
#
# CS
fit2_asreml <- asreml(fixed = YLD ~ GEN + HAR + GEN:HAR,
                      random = ~ LOC + LOC:HAR + LOC:GEN + LOC:GEN:HAR + LOC:BLK + LOC:BLK:HAR + LOC:BLK:PLT + LOC:BLK:PLT:HAR,
                      residual = ~ units,
                      data = pheno_dat_use2)

summary(fit2_asreml)$varcomp

# # Compare to Table 3 in Piepho 2013
# LOC:HAR variance = 5.896039909
# LOC:HAR covariance = 105.05504
#
# LOC:GEN:HAR variance = 22.60766
# LOC:GEN:HAR covariance = 12.62059
#


# One location, several trials --------------------------------------------


pheno_dat_use3 <- pheno_dat_use_142 %>%
  droplevels()

fit3 <- lmer(yld ~ gen + har + gen:har + (1|year) + (1|year:gen) + (1|year:gen:har) +
               (1|trl) + (1|trl:har) + (1|trl:blk) + (1|trl:blk:har) + (1|trl:blk:plot) + (1|trl:blk:plot:har),
             data = pheno_dat_use3,
             control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))

as.data.frame(VarCorr(fit3))



# Two-stage analysis ------------------------------------------------------



# # One location, one trial -------------------------------------------

# Calculate genotype means per harvest year and location

pheno_dat_use4 <- pheno_dat_use_142 %>%
  filter(trl == "2003") %>%
  droplevels()

# Loop over harvest years and locations
stage_one_means <- pheno_dat_use4 %>%
  group_by(har) %>%
  do({

    df <- .
    df1 <- droplevels(df)

    # Fit a model for this harvest year - location
    fit_i <- lmer(yld ~ 0 + gen + (1|blk) + (1|blk:plot), data = df1,
                  control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))



    stage_one_fit1 <- lmer(yld ~ gen + har + gen:har + (1|loc) + (1|loc:har) + (1|loc:gen) + (1|loc:gen:har) +
                             (1|loc:blk) + (1|loc:blk:har) + (1|loc:blk:plot) + (1|loc:blk:plot:har),
                           data = pheno_dat_use4,
                           control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))






# # Several locations, one trial -------------------------------------------

# Calculate genotype means per harvest year and location

pheno_dat_use4 <- pheno_dat %>%
  filter(trl == "2003") %>%
  droplevels()

# Loop over harvest years and locations
stage_one_means <- pheno_dat_use4 %>%
  group_by(har, loc) %>%
  do({

    df <- .
    df1 <- droplevels(df)

    # Fit a model for this harvest year - location
    fit_i <- lmer(yld ~ 0 + gen + (1|blk) + (1|blk:plot), data = df1,
                  control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))



stage_one_fit1 <- lmer(yld ~ gen + har + gen:har + (1|loc) + (1|loc:har) + (1|loc:gen) + (1|loc:gen:har) +
                         (1|loc:blk) + (1|loc:blk:har) + (1|loc:blk:plot) + (1|loc:blk:plot:har),
                       data = pheno_dat_use4,
                       control = lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore"))




## Analysis of a series of variety trials in perennial crops
##
## from Piepho and Eckl 2013 doi:10.1111/gfs.12054
##

# Load packages
# library(sommer)
library(tidyverse)
# library(nlme)
library(glmmTMB)
# library(lme4)
library(asreml)
library(StageWise)
library(mmrm)


# Read in the data
pheno_dat <- read_csv(file = "data/piepho2013_lolium_data.csv") %>%
  # Rename columns like the notation in the paper
  select(yld, GEN, LOC, trl, blk, PLT = PLOT, yr, har) %>%
  rename_all(toupper) %>%
  mutate_at(vars(-YLD), as.factor)

# Set asreml options
asreml.options(workspace = "500mb", ai.sing = TRUE, maxit = 30)



# Single-stage analysis ---------------------------------------------------

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

# Two-stage analysis will involve calculating genotype-harvest year means and their
# variances per trial and then using this information in a stage-2 analysis

# One location, one trial -------------------------------------------



# One location, several trials --------------------------------------------

# For each trial in the location, fit a model for the single-location, single-trial
pheno_dat_use_142 <- pheno_dat %>%
  subset(LOC == "142") %>%
  droplevels()


# Group by trl
stage1_loc142 <- pheno_dat_use_142 %>%
  group_by(TRL) %>%
  do({
    df <- .

    # # TMB
    # fit0 <- glmmTMB(formula = YLD ~ GEN + HAR + GEN:HAR + (1|BLK:PLT) + (1|BLK),
    #                 data = df, control = glmmTMBControl(rank_check = "adjust"))
    # fit1 <- update(fit0, formula = YLD ~ GEN + HAR + GEN:HAR + (1|BLK:PLT) + (1|BLK:PLT:HAR) +
    #                  (1|BLK) + (1|BLK:HAR))
    # fit2 <- update(fit0, formula = YLD ~ GEN + HAR + GEN:HAR + cs(HAR + 0|BLK:PLT) + cs(HAR + 0|BLK))
    # fit3 <- update(fit0, formula = YLD ~ GEN + HAR + GEN:HAR + ar1(HAR + 0|BLK:PLT) + ar1(HAR + 0|BLK))
    #
    #
    # nd <- distinct(df, GEN, HAR) %>%
    #   mutate(BLK = NA, PLT = NA)
    # pred_fit0 <- predict(fit0, newdata = nd, se.fit = TRUE, cov.fit = TRUE, re.form = ~ 0, allow.new.levels = TRUE)
    # pred_fit1 <- predict(fit1, newdata = nd, se.fit = TRUE, cov.fit = TRUE, re.form = ~ 0, allow.new.levels = TRUE)
    #
    # # MMRM
    # fit0_lmer <- lmer(YLD ~ GEN + HAR + GEN:HAR + (1|BLK:PLT) + (1|BLK),
    #                   data = df)
    # fit1_mmrm <- mmrm(formula = YLD ~ GEN + HAR + GEN:HAR + cs(HAR|BLK/PLT),
    #                   data = df)
    # fit2_mmrm <- mmrm(formula = YLD ~ GEN + HAR + GEN:HAR + csh(HAR|BLK/PLT),
    #                   data = df)
    # fit3_mmrm <- mmrm(formula = YLD ~ GEN + HAR + GEN:HAR + ar1(HAR|BLK/PLT),
    #                   data = df)
    # fit4_mmrm <- mmrm(formula = YLD ~ GEN + HAR + GEN:HAR + ar1h(HAR|PLT),
    #                   data = df)

    # Number of harvest years
    nHAR <- length(unique(df$HAR))

    # Try different single location, single-trial models
    fit0 <- asreml(fixed = YLD ~ GEN + HAR + GEN:HAR,
                   random = ~ BLK:PLT:id(HAR) + BLK:id(HAR),
                   residual = ~ units,
                   data = df)
    fit1_cs <- update(fit0, random = ~ BLK:PLT:corv(HAR) + BLK:corv(HAR))
    fit2_csh <- update(fit0, random = ~ BLK:PLT:corh(HAR) + BLK:corh(HAR))
    # Fit AR1 only if nHAR >= 3
    if (nHAR >= 3) {
      fit3_ar1 <- update(fit0, random = ~ BLK:PLT:ar1(HAR) + BLK:ar1(HAR))
      fit4_ar1h <- update(fit0, random = ~ BLK:PLT:ar1h(HAR) + BLK:ar1h(HAR))
      fit5_us <- update(fit0, random = ~ BLK:PLT:us(HAR) + BLK:us(HAR))
    } else {
      fit3_ar1 <- fit4_ar1h <- fit5_us <- NULL
    }

    # List of models
    fit_list <- list(us = fit5_us, ar1h = fit4_ar1h, ar1 = fit3_ar1, csh = fit2_csh, cs = fit1_cs, id = fit0)
    fit_list <- subset(fit_list, !sapply(fit_list, is.null))

    # Pick the fit with the lowest AIC
    aics <- sapply(fit_list, function(x) summary(x)$aic)
    fit_use <- fit_list[[which.min(aics)]]

    # Predict means
    preds <- predict.asreml(fit_use, classify = "GEN:HAR", vcov = TRUE)

    # Return the means and VCOV
    tibble(covstr = names(which.min(aics)), varcomp = list(summary(fit_use)$varcomp),
           blues = list(as.data.frame(preds$pvals)), vcov = list(preds$vcov))

  }) %>% ungroup()


# Try using Stage2 from the StageWise package

# Prepare the data and VCOV matrix
#
# Set gen x har as id
# Set year as env
# Set trl as loc
#
stage1_loc142_blues <- stage1_loc142 %>%
  select(TRL, blues) %>%
  unnest(blues) %>%
  rename(BLUE = predicted.value) %>%
  mutate(id = paste0(GEN, ":", HAR),
         env = TRL)
stage1_loc142_vcov <- stage1_loc142 %>%
  mutate(vcov = map2(vcov, blues, ~{
    id <- paste0(.y$GEN, ":", .y$HAR)
    `dimnames<-`(.x, list(id, id))
  })) %>%
  pull(vcov)
names(stage1_loc142_vcov) <- levels(stage1_loc142_blues$env)

stage2_out <- Stage2(data = stage1_loc142_blues, vcov = stage1_loc142_vcov,
                     max.iter = 25)

# Try another approach
# id = GEN
# loc = HAR
# env = TRLxHAR

stage1_loc142_blues <- stage1_loc142 %>%
  select(TRL, blues) %>%
  unnest(blues) %>%
  rename(BLUE = predicted.value) %>%
  left_join(., distinct(pheno_dat_use_142, TRL, HAR, YR)) %>%
  mutate(id = GEN, loc = HAR, env = YR)
stage1_loc142_vcov <- stage1_loc142 %>%
  mutate(vcov = map2(vcov, blues, ~{
    id <- paste0(.y$GEN, ":", .y$HAR)
    `dimnames<-`(.x, list(id, id))
  })) %>%
  pull(vcov)
names(stage1_loc142_vcov) <- levels(stage1_loc142_blues$env)

stage2_out <- Stage2(data = stage1_loc142_blues, max.iter = 100)


# BLUP prep
prep <- blup_prep(data = stage1_loc142_blues, vars = stage2_out$vars)
# Calculate blups
har_index <- c("1" = 0.2, "2" = 0.3, "3" = 0.5)
blups <- blup(data = prep, what = "GV", index.coeff = har_index)


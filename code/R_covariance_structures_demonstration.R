## 
## Covariance structures for mixed models in R
## 

library(tidyverse)
library(agridat)
library(lme4)
library(nlme)

# Dataset

dat <- as_tibble(vargas.wheat1.traits) %>%
  select(year, rep, gen, yield) %>%
  mutate(year = as.factor(year),
         dummy = factor(1))

# dat <- dat %>%
#   filter(year %in% c(1990, 1991)) %>%
#   droplevels()


# Empty VCOV matrix for genotypes in environments
G_blank <- matrix(0, nrow = nlevels(dat$year), ncol = nlevels(dat$year),
                  dimnames = list(levels(dat$year), levels(dat$year)))



# Identity - Homogenous variance; no covariance between environments

fit1 <- lmer(yield ~ 1 + year + (1|gen), dat)

var_df <- as.data.frame(VarCorr(fit1))

# What does the G matrix look like for environments?
G_fit1 <- G_blank
diag(G_fit1) <- var_df$vcov[var_df$grp == "gen"]
G_fit1




# Compound symmetry - homogenous variance and homogenous covariance

fitCS <- lmer(yield ~ 1 + year + (1|gen) + (1|gen:year), dat)

(var_df <- as.data.frame(VarCorr(fitCS)))

# Variance in each environment
varG <- var_df$vcov[var_df$grp == "gen"] + var_df$vcov[var_df$grp == "gen:year"]
# Covariance across environments
covGE <- var_df$vcov[var_df$grp == "gen"]

# What does the G matrix look like for environments?
G_fitCS <- G_blank
G_fitCS[] <- covGE
diag(G_fitCS) <- varG
G_fitCS




# # CS using nlme
# 
# fitCS_1 <- lme(fixed = yield ~ 1 + year,
#                random = list(gen = pdCompSymm(~ year - 1)),
#                dat)
# 
# varcor <- VarCorr(fitCS_1)


# Diagonal - heterogenous variance and zero correlation

diagMM <- model.matrix(~ -1 + year, dat)
dat1 <- cbind(dat, diagMM)

ranef_terms <- paste0("(0 + ", colnames(diagMM), "|gen)")
formDIAG <- reformulate(c("year", ranef_terms), response = "yield")

fitDIAG <- lmer(formDIAG, dat1)

(var_df <- as.data.frame(VarCorr(fitDIAG)))

varG <- subset(var_df, grp != "Residual", vcov, drop = TRUE)


# What does the G matrix look like for environments?
G_fitDIAG <- G_blank
diag(G_fitDIAG) <- varG
G_fitDIAG


# fitDIAG_1 <- lme(fixed = yield ~ 1 + year,
#                  random = list(gen = pdDiag(~ year - 1)),
#                  dat)
# 
# VarCorr(fitDIAG_1)


# Unstructured - heterogenous variance and heterogeneous covariance

fitUS <- lmer(yield ~ 1 + year + (0 + year|gen), dat)

(var_df <- as.data.frame(VarCorr(fitUS)))

varG <- subset(var_df, grp == "gen" & is.na(var2), vcov, drop = TRUE)

covGE_mat <- subset(var_df, grp == "gen" & !is.na(var2)) %>% 
  select(var1, var2, vcov) %>% 
  spread(var2, vcov) %>% 
  column_to_rownames("var1") %>% 
  as.matrix()

G_fitUS <- G_blank

G_fitUS[upper.tri(G_fitUS)] <- covGE_mat[upper.tri(covGE_mat, diag = TRUE)]
G_fitUS[lower.tri(G_fitUS)] <- G_fitUS[upper.tri(G_fitUS)]

diag(G_fitUS) <- varG

G_fitUS









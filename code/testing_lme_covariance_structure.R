## Testing nlme
##
##

library(agridat)
library(tidyverse)
library(nlme)

# Load a dataset
dat <- minnesota.barley.yield %>%
  as_tibble() %>%
  filter(year %in% 1929:1930) %>%
  mutate(env = interaction(site, year, sep = ".", drop = TRUE),
         dummy = factor(1)) %>%
  droplevels()




# Fit a model
fit <- lme(fixed = yield ~ env, random = list(gen = pdIdent(~ env - 1)), data = dat)

# Variance-covariance matrix of genotypes in environments
G <- getVarCov(fit)

# compount symmetry (constant variance and covariance)
fit1 <- lme(fixed = yield ~ env, random = list(gen = pdCompSymm(~ env - 1)), data = dat)
G1 <- getVarCov(fit1)

# Diagonal
fit2 <- lme(fixed = yield ~ env, random = list(gen = pdDiag(~ env - 1)), data = dat)
G2 <- getVarCov(fit2)

# Diagonal plus CS
fit3 <- lme(fixed = yield ~ env, random = list(gen = pdDiag(~ env - 1), dummy = pdIdent(~ gen:env)), data = dat)


# Unstructured
fit3 <- lme(fixed = yield ~ env, random = ~ env | gen, data = dat)
fit3 <- lme(fixed = yield ~ env, random = list(gen = pdNatural(~ env - 1)), data = dat)
G3 <- getVarCov(fit2)




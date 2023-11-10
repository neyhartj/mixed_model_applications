
##################################################################
#
# R script used to fit linear mixed models in Tolhurst et al. 2022
#
# Script author: D.J. Tolhurst
#
##################################################################

# Assume that data.df is a MET data-frame with the following columns
# Yld - numeric, n = 10,368 observations
# Env - factor, p = 24 levels
# Trial - factor, t = 72 levels
# Block - factor, b = 2 levels
# Column - factor, c levels
# Row - factor, r levels
# Gkeep - factor, v = 204 levels, and corresponds to those genotypes with marker data
# Gdrop - factor, vd = 4 levels, and corresponds to those genotypes without marker data
# (Gdrop is not required if all genotypes have marker data)

# Also assume that:
# rc is an object which contains the names of environments in which random column effects are fitted
# rr is an object which contains the names of environments in which random row effects are fitted

# Lastly, assume that Gg is a v x v genomic relationship with column and row names ordered according to the levels in Gkeep

############################################
#
# Baseline models
#
############################################

# Independent model
id.asr <- asreml(Yld ~ Env,
                 random= ~ idv(Env):vm(Gkeep, Gg) +
                           diag(Env):Trial +
                           diag(Env):Trial:Block +
                           at(Env, rc):Column +
                           at(Env, rr):Row,
                 residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                 sparse= ~Env:Gdrop,
                 data.df= data.df,
                 na.action = na.method(x="include", y = "include"),
                 workspace = 6e8)
id.asr <- update(id.asr)


# Diagonal model
diag.asr <- asreml(Yld ~ Env,
                   random= ~ diag(Env):vm(Gkeep, Gg) +
                             diag(Env):Trial +
                             diag(Env):Trial:Block +
                             at(Env, rc):Column +
                             at(Env, rr):Row,
                   residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                   sparse= ~Env:Gdrop,
                   data.df= data.df,
                   na.action = na.method(x="include", y = "include"),
                   workspace = 6e8)
diag.asr <- update(diag.asr)


#  -- Models with simple main effects --
# compound symmetric model
cs.asr <- asreml(Yld ~ Env,
                 random= ~ vm(Gkeep, Gg) + idv(Env):vm(Gkeep, Gg) +
                           diag(Env):Trial +
                           diag(Env):Trial:Block +
                           at(Env, rc):Column +
                           at(Env, rr):Row,
                  residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                  sparse= ~Env:Gdrop,
                  data.df = data.df,
                  na.action = na.method(x="include", y = "include"),
                  workspace = 6e8)
cs.asr <- update(cs.asr)


# main effects plus diagonal model
mdiag.asr <- asreml(Yld ~ Env,
                    random= ~ vm(Gkeep, Gg) + diag(Env):vm(Gkeep, Gg) +
                              diag(Env):Trial +
                              diag(Env):Trial:Block +
                              at(Env, rc):Column +
                              at(Env, rr):Row,
                    residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                    sparse= ~Env:Gdrop,
                    data.df = data.df,
                    na.action = na.method(x="include", y = "include"),
                    workspace = 6e8)
mdiag.asr <- update(mdiag.asr)


############################################
#
# Regressions on latent covariates
#
############################################

#  -- Models with simple main effects --
# Factor analytic (FAMk) model with intercepts explicitly fitted, Equation 5 in the manuscript
# Note that the factor analytic model is fitted in terms of its two components, that is a reduced rank (rr) + diagonal (diag) model
FAM1.asr <- asreml(Yld ~ Env,
                   random= ~  vm(Gkeep, Gg) +
                              rr(Env, 2):vm(Gkeep, Gg) + diag(Env):vm(Gkeep, Gg) +
                              diag(Env):Trial +
                              diag(Env):Trial:Block +
                              at(Env, rc):Column +
                              at(Env, rr):Row,
                   residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                   sparse= ~Env:Gdrop,
                   data.df = data.df,
                   na.action = na.method(x="include", y = "include"),
                   workspace = 6e8)
FAM1.asr <- update(FAM1.asr)


#  -- Models with generalised main effects --
# Conventional factor analytic (FAk) model, Equation 10 in the manuscript
# Note that the factor analytic model is fitted in terms of its two components, that is a reduced rank (rr) + diagonal (diag) model
FA1.asr <- asreml(Yld ~ Env,
                  random= ~  rr(Env, 1):vm(Gkeep, Gg) + diag(Env):vm(Gkeep, Gg) +
                             diag(Env):Trial +
                             diag(Env):Trial:Block +
                             at(Env, rc):Column +
                             at(Env, rr):Row,
                  residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                  sparse= ~Env:Gdrop,
                  data.df = data.df,
                  na.action = na.method(x="include", y = "include"),
                  workspace = 6e8)
FA1.asr <- update(FA1.asr)


############################################
#
# Regressions on known covariates
#
############################################

# Assume S is a pxq matrix of known environmental covariates
# Env.names is a p-vector of environment names

#  -- Models without translational invariance --
# Random regression model in Jarquin et al. (2014), Equation 13 in the manuscript with variance matrix in Equation 15
S.1 <- data.df.frame(Env.names, S)
rreg1.asr <- asreml(Yld ~ Env,
                    random= ~ vm(Gkeep, Gg) +
                              str(~mbf(S):vm(Gkeep, Gg), vmodel = ~ idv(18):vm(Gkeep, Gg)) +
                              # str is used for general variance structures
                              diag(Env):vm(Gkeep, Gg) +
                              diag(Env):Trial +
                              diag(Env):Trial:Block +
                              at(Env, rc):Column +
                              at(Env, rr):Row,
                    residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                    sparse= ~Env:Gdrop,
                    mbf= list(S = list(key= c("Env","Env"), cov = "S.1")),
                    data.df = data.df,
                    na.action = na.method(x="include", y = "include"),
                    workspace = 6e8)
rreg1.asr <- update(rreg1.asr)


# Random regression model in Heslot et al. (2014), Equation 13 in the manuscript with variance matrix in Equation 14
S.1 <- data.df.frame(Env.names, S)
rreg2.asr <- asreml(Yld ~ Env,
                    random= ~ vm(Gkeep, Gg) +
                              str(~mbf(S):vm(Gkeep, Gg), vmodel = ~ diag(18):vm(Gkeep, Gg)) +
                              diag(Env):vm(Gkeep, Gg) +
                              diag(Env):Trial +
                              diag(Env):Trial:Block +
                              at(Env, rc):Column +
                              at(Env, rr):Row,
                    residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                    mbf= list(S = list(key= c("Env","Env"), cov = "S.1")),
                    sparse= ~Env:Gdrop,
                    data.df = data.df,
                    na.action = na.method(x="include", y = "include"),
                    workspace = 6e8)
rreg2.asr <- update(rreg2.asr)


#  -- Models with translational invariance --
# Factor analytic regression (FARk) model, Equation 18 in the manuscript
A.1 <- data.df.frame(Env.names, 1, S, 0)
# There is an issue in ASReml-R when using rr within the mbf argument.
# To resolve this issue, a column of zeros is required in A
FAR1.asr <- asreml(Yld ~ Env,
                   random= ~ str(~mbf(A):vm(Gkeep, Gg), vmodel = ~ rr(19, 1):vm(Gkeep, Gg)) +
                             diag(Env):vm(Gkeep, Gg) +
                             diag(Env):Trial +
                             diag(Env):Trial:Block +
                             at(Env, rc):Column +
                             at(Env, rr):Row,
                   residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                   mbf= list(A = list(key= c("Env","Env"), cov = "A.1")),
                   sparse= ~Env:Gdrop,
                   data.df = data.df,
                   na.action = na.method(x="include", y = "include"),
                   workspace = 6e8)
FAR1.asr <- update(FAR1.asr)


############################################
#
# Regressions on known and latent covariates
#
############################################

# Integrated factor analytic (IFAk) model

#  -- Models with generalised main effects and translational invariance--
# Dependent form of the IFAk model, Equation 23 in the manuscript
# Construct Gamma as a (p-q)xq orthogonal projection matrix, such that t(S) %*% Gamma = 0
Gamma <- diag(p) - S %*% solve(t(S) %*% S) %*% t(S)
B <- cbind(S, Gamma[,1:(p-q)])
B.1 <- data.df.frame(Env.names, B, 0)
# There is an issue in ASReml-R when using rr within the mbf argument.
# To resolve this issue, a column of zeros is required in B

dIFA1.asr <- asreml(Yld ~ Env,
                    random= ~ str(~mbf(B):vm(Gkeep, Gg), vmodel = ~ rr(24, 1):vm(Gkeep, Gg)) +
                              diag(Env):vm(Gkeep, Gg) +
                              diag(Env):Trial +
                              diag(Env):Trial:Block +
                              at(Env, rc):Column +
                              at(Env, rr):Row,
                    residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                    mbf= list(B = list(key= c("Env","Env"), cov = "B.1")),
                    sparse= ~Env:Gdrop,
                    data.df = data.df,
                    na.action = na.method(x="include", y = "include"),
                    workspace = 6e8)
dIFA1.asr <- update(dIFA1.asr)

# The factor loadings can be obtained via:
dLam.IFA1 <- matrix(dIFA1.asr$vparameters[grep("mbf.*fa", names(dIFA1.asr$vparameters))], ncol=1)

# Alternatively, the factor loadings in the IFAk model can be
# obtained directly from the fit of the conventional FAk model
Lam.FA1 <- matrix(FA1.asr$vparameters[grep("rr.*fa", names(FA1.asr$vparameters))], ncol= 1)
dLam.IFA1 <- solve(B) %*% Lam.FA1


#  -- Models with generalised main effects, but without translational invariance--
# Independent form of the IFAk model
# variant 1, IFAk model with orthogonal known and latent covariates, Equation 49 in the Supplementary Material
S.1 <- data.df.frame(Env.names, S, 0)
Gamma.1 <- data.df.frame(Env.names, Gamma[,1:(p-q)], 0)
iIFA1.v1.asr <- asreml(Yld ~ Env,
                       random= ~ str(~mbf(S):vm(Gkeep, Gg), vmodel = ~ rr(18, 1):vm(Gkeep, Gg)) +
                                 str(~mbf(Gamma):vm(Gkeep, Gg), vmodel = ~ rr(6, 1):vm(Gkeep, Gg)) +
                                 diag(Env):vm(Gkeep, Gg) +
                                 diag(Env):Trial +
                                 diag(Env):Trial:Block +
                                 at(Env, rc):Column +
                                 at(Env, rr):Row,
                       residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                       mbf= list(S = list(key= c("Env","Env"), cov = "S.1"),
                                 Gamma = list(key= c("Env","Env"), cov = "Gamma.1")),
                       sparse= ~Env:Gdrop,
                       data.df = data.df,
                       na.action = na.method(x="include", y = "include"),
                       workspace = 6e8)
iIFA1.v1.asr <- update(iIFA1.v1.asr)


# variant 2, IFAk model without orthogonal known and latent covariates, Equation 51 in The Supplementary Material
iIFA1.v2.asr <- asreml(Yld ~ Env,
                       random= ~ str(~mbf(S):vm(Gkeep, Gg), vmodel = ~ rr(18, 1):vm(Gkeep, Gg)) +
                                 rr(Env, 1):vm(Gkeep, Gg) +
                                 diag(Env):vm(Gkeep, Gg) +
                                 diag(Env):Trial +
                                 diag(Env):Trial:Block +
                                 at(Env, rc):Column +
                                 at(Env, rr):Row,
                       residual = ~  dsum(~ar1(Column):ar1(Row)| Env),
                       mbf= list(S = list(key= c("Env","Env"), cov = "S.1")),
                       sparse= ~Env:Gdrop,
                       data.df = data.df,
                       na.action = na.method(x="include", y = "include"),
                       workspace = 6e8)
iIFA1.v2.asr <- update(iIFA1.v2.asr)



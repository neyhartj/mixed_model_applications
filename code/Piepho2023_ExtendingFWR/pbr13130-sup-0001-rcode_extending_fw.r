
#################### Extending Finlay-Wilkinson regression with environmental covariates #########################################################################################################

#################### Initialization #########################################################################################################

#### Loading packages

library(asreml)
library(pls)
library(reshape2)

#### Preparing data: *Lettuce data from van Eeuwijk 1992 TAG 85:92-100*

nitrate=structure(list(ENV=structure(1:18,levels=c("1","2","3",
"4","5","6","7","8","9","10","11","12","13","14","15",
"16","17","18"),class="factor"),Pa=c(3.113,3.379,3.067,
3.202,3.921,4.153,4.851,4.547,3.721,3.581,3.312,3.439,
3.195,2.89,2.7,3.143,2.746,3.273),DM=c(2.835,3.222,
2.326,2.663,3.365,3.97,4.512,4.203,3.505,3.298,3.13,
3.329,3.047,2.297,2.43,2.71,2.47,2.384),Pi=c(2.629,
2.848,2.511,2.23,3.028,3.444,4.01,3.429,3.337,3.287,
2.959,3.254,2.948,2.295,2.172,2.429,2.226,2.555),GT=c(1.988,
2.823,2.12,1.638,2.653,2.813,3.504,2.944,2.425,2.389,
2.28,2.561,2.696,2.237,2.004,2.26,2.126,2.167),RW=c(2.199,
3.002,2.692,2.187,2.935,2.865,3.135,2.616,2.177,2.159,
1.797,2.843,2.61,1.93,2.194,2.406,2.332,2.545),Wi=c(2.414,
2.95,2.598,2.171,2.931,3.232,3.624,3.052,2.525,2.681,
2.152,3.035,2.902,2.414,2.392,2.438,2.185,2.386),Tr=c(1.248,
2.176,1.032,1.062,2.007,2.341,3.08,2.817,1.917,1.744,
1.365,1.927,1.914,1.462,1.374,1.536,1.287,1.616),Ls=c(2.38,
3.196,2.355,1.599,2.942,3.289,3.612,3.07,2.83,2.726,
2.178,3.058,3.138,2.274,2.144,2.464,2.621,2.813)),row.names=c(NA,
-18L),class="data.frame")

environment=structure(list(ENV=structure(1:18,levels=c("1","2","3",
"4","5","6","7","8","9","10","11","12","13","14","15",
"16","17","18"),class="factor"),x1=c(2.1,2.1,2.1,2.1,
2.2,2,2,1.9,2.2,2,2,1.6,1.5,1.5,1.5,2.3,1.5,1.5),
x2=c(1136L,1345L,1700L,1076L,960L,316L,145L,109L,
555L,641L,676L,1951L,1651L,2281L,1244L,1398L,2041L,
1326L),x3=c(993L,1277L,2191L,1090L,779L,482L,117L,
93L,504L,663L,666L,2427L,1789L,2359L,1456L,1852L,
1515L,1416L),x4=c(911L,1250L,2586L,1323L,539L,421L,
102L,127L,415L,596L,541L,2413L,1276L,2376L,1604L,
2719L,1350L,1779L),x5=c(881L,1815L,2556L,1065L,457L,
556L,42L,42L,383L,780L,546L,2286L,1518L,2514L,1398L,
2975L,988L,1580L),x6=c(14.75,11.78,16.14,14.81,13.61,
12.81,11.91,10.96,9.64,8.45,7.7,9.22,11.78,13.15,
14.07,14.51,14.93,15.26),x7=c(11.23,13.28,16.48,
13.61,12.46,11.23,10.29,8.71,8.5,7.98,9.04,12.05,
13.01,14.45,15.37,15.96,16.23,16.39),x8=c(13.15,
14.93,16.37,12.74,10.89,9.04,8.05,7.73,9.9,11.78,
12.6,14.39,15.21,15.62,16.23,16.45,16.48,16.41)),row.names=c(NA,
-18L),class="data.frame")

lettuce0=melt(nitrate,id.vars = "ENV")
colnames(lettuce0)[2:3]=c("GENOTYPE","nitrate")
lettuce0=merge(lettuce0,cbind(ENV=environment[,1],scale(environment[,-1])),by="ENV")
lettuce0=lettuce0[order(lettuce0$ENV,lettuce0$GENOTYPE),]


# environmental covariates matrices

C=scale(environment[,-1]) # C is identical to 'environment' data.frame but scaled

# A column of zeros is required in for each rank of rr in asreml-R
C1 <- data.frame(nitrate$ENV,C,0) # rr(,1)
C2 <- data.frame(nitrate$ENV,C,0,0) # rr(,2)

asreml.options(maxit=100,extra=10)




#################### FA approach #########################################################################################################

##### Initial models : M1 to M4

#### Model M1 :

mdl1=asreml(fixed = nitrate~-1+GENOTYPE+ENV,
            na.action = na.method(x="include", y = "include"),
            workspace = "0.2gb",
            data = lettuce0)

md1_alt <- lm(nitrate ~ -1 + GENOTYPE + ENV, data = lettuce0)

# Equals results in Table 1
logLik(md1_alt, REML = TRUE) * -2


#### Model M2 :

mdl2=asreml(fixed = nitrate~-1+GENOTYPE+ENV,
            residual = ~ar1(ENV):id(GENOTYPE),
            na.action = na.method(x="include", y = "include"),
            workspace = "0.2gb",
            data = lettuce0)



#### Model M3 :

mdl_M3=asreml(fixed = nitrate~GENOTYPE+ENV,
              random=~str(~mbf(Co):GENOTYPE, vmodel = ~ rr(8,1):id(GENOTYPE)),
              mbf= list(Co = list(key= c("ENV","nitrate.ENV"), cov = "C1")),
              residual = ~ar1(ENV):id(GENOTYPE),
              na.action = na.method(x="include", y = "include"),
              workspace = "0.2gb",
              data = lettuce0)







# extract lambdas
lambda=summary(mdl_M3)$varcomp[grep(x = rownames(summary(mdl_M3)$varcomp),pattern = "fa",fixed = T),"component"]

# Computing z1
z1=t(lambda %*% t(C))

lettuce=merge(x = lettuce0,y = data.frame(ENV=unique(lettuce0$ENV),z1=z1))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#### FA1

mdl_fa1=asreml(fixed = nitrate~-1+GENOTYPE+z1:GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",
               data = lettuce)

effects=mdl_fa1$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept"],pred)

# get the slopes
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope"
lettuce=merge(lettuce[,colnames(lettuce)!="slope"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

summary(mdl_fa1)


#### Model M4 :
rm(lettuce)
# Need to specify M3 variance components as starting values to help convergence
mdl_M4=asreml(fixed = nitrate~GENOTYPE+ENV,
              random=~str(~mbf(Co):GENOTYPE, vmodel = ~ rr(8,2):id(GENOTYPE)),

              mbf= list(Co = list(key= c("ENV","nitrate.ENV"), cov = "C2")),
              residual = ~ar1(ENV):id(GENOTYPE),
              na.action = na.method(x="include", y = "include"),
              workspace = "0.2gb",
              data = lettuce0,G.param = mdl_M3$G.param,R.param = mdl_M3$R.param)

# extract lambdas
lambda=summary(mdl_M4)$varcomp[grep(x = rownames(summary(mdl_M4)$varcomp),pattern = "fa",fixed = T),"component"]
lambda1=lambda[1:8]
lambda2=lambda[9:16]

# Computing z1 and z2
z1=t(lambda1 %*% t(C))
z2=t(lambda2 %*% t(C))

lettuce=merge(x = lettuce0,y = data.frame(ENV=unique(lettuce0$ENV),z1=z1,z2=z2))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#### FA2

mdl_fa2=asreml(fixed = nitrate~-1+GENOTYPE+z1:GENOTYPE+z2:GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",
               data = lettuce)

summary(mdl_fa2)


effects=mdl_fa2$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept"],pred)

# get the slopes for z1
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope1"
lettuce=merge(lettuce[,colnames(lettuce)!="slope1"],pred)

# get the slopes for z2
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z2:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z2:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope2"
lettuce=merge(lettuce[,colnames(lettuce)!="slope2"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]




#################### PLS approach #########################################################################################################
rm(lettuce)
# run the PLS to get z1 and z2

Y=as.matrix(nitrate[,-1])
X=as.matrix(C)

pls=mvr(Y~X,scale = F,center = F,ncomp = 3,method = "oscorespls")

z1=pls$scores[,1]
z2=pls$scores[,2]

#### PLS1

lettuce=merge(x = lettuce0,y = data.frame(ENV=unique(lettuce0$ENV),z1=z1))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

mdl_pls1=asreml(fixed = nitrate~-1+GENOTYPE+z1:GENOTYPE,
                random=~ENV,
                residual = ~ar1(ENV):id(GENOTYPE),
                na.action = na.method(x="include", y = "include"),
                workspace = "0.2gb",
                data = lettuce)

effects=mdl_pls1$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept"],pred)

# get the slopes
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope"
lettuce=merge(lettuce[,colnames(lettuce)!="slope"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#### PLS2
rm(lettuce)
lettuce=merge(x = lettuce0,y = data.frame(ENV=unique(lettuce0$ENV),z1=z1,z2=z2))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

mdl_pls2=asreml(fixed = nitrate~-1+GENOTYPE+z1:GENOTYPE+z2:GENOTYPE,
                random=~ENV,
                residual = ~ar1(ENV):id(GENOTYPE),
                na.action = na.method(x="include", y = "include"),
                workspace = "0.2gb",
                data = lettuce)

effects=mdl_pls2$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept"],pred)

# get the slopes for z1
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope1"
lettuce=merge(lettuce[,colnames(lettuce)!="slope1"],pred)

# get the slopes for z2
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z2:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z2:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope2"
lettuce=merge(lettuce[,colnames(lettuce)!="slope2"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#################### FR approach #########################################################################################################
rm(lettuce)
#### FR model
lettuce=lettuce0
mdl_FR=asreml(fixed = nitrate~-1+GENOTYPE+x1:GENOTYPE+x2:GENOTYPE+x3:GENOTYPE+x4:GENOTYPE+x5:GENOTYPE+x6:GENOTYPE+x7:GENOTYPE+x8:GENOTYPE,
              random=~ENV,
              residual = ~ar1(ENV):id(GENOTYPE),
              na.action = na.method(x="include", y = "include"),
              workspace = "0.2gb",
              data = lettuce0)

summary(mdl_FR)

#################### RA approach #########################################################################################################

# Use FR model to get fitted values for the multiplicative terms

lettuce$FV=predict(object = mdl_FR,classify = "GENOTYPE:x1:x2:x3:x4:x5:x6:x7:x8",ignore = "GENOTYPE",levels=list(GENOTYPE=lettuce$GENOTYPE,x1=lettuce$x1,x2=lettuce$x2,x3=lettuce$x3,x4=lettuce$x4,x5=lettuce$x5,x6=lettuce$x6,x7=lettuce$x7,x8=lettuce$x8),parallel = T)$pvals[,"predicted.value"]
Y=reshape::cast(data = lettuce[,c("GENOTYPE","ENV","FV")],formula = GENOTYPE~ENV,value = "FV")[,-1]
row.names(Y)=levels(lettuce$GENOTYPE)

# Then get SVD of Y matrix

z1=svd(Y)$v[,1]
z2=svd(Y)$v[,2]



#### RA1
rm(lettuce)
lettuce=merge(x = lettuce0,y = data.frame(ENV=unique(lettuce0$ENV),z1=z1))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

mdl_ra1=asreml(fixed = nitrate~-1+GENOTYPE+z1:GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",
               data = lettuce)

effects=mdl_ra1$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept"],pred)

# get the slopes
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope"
lettuce=merge(lettuce[,colnames(lettuce)!="slope"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#### RA2
rm(lettuce)
lettuce=merge(x = lettuce0,y = data.frame(ENV=unique(lettuce0$ENV),z1=z1,z2=z2))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

mdl_ra2=asreml(fixed = nitrate~-1+GENOTYPE+z1:GENOTYPE+z2:GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",
               data = lettuce)

effects=mdl_ra2$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept"],pred)

# get the slopes for z1
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope1"
lettuce=merge(lettuce[,colnames(lettuce)!="slope1"],pred)

# get the slopes for z2
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z2:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z2:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope2"
lettuce=merge(lettuce[,colnames(lettuce)!="slope2"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#################### CCR approach #########################################################################################################
rm(lettuce)
# initializing

z1_it=rowMeans(nitrate[,-1]) # Initial value for z1 : Env means from raw data

lettuce=merge(x = lettuce0[,colnames(lettuce0)!="z1_it"],y = data.frame(ENV=unique(lettuce0$ENV),z1_it=z1_it))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

mdl_CCR_init=asreml(fixed = nitrate~-1+GENOTYPE+z1_it:GENOTYPE,
                    # random=~ENV,
                    residual = ~ar1(ENV):id(GENOTYPE),
                    na.action = na.method(x="include", y = "include"),
                    workspace = "0.2gb",maxit=100,
                    data = lettuce)

v_e=0.1
v_r=summary(mdl_CCR_init)$varcomp["ENV:GENOTYPE!R",1]
cor_r=summary(mdl_CCR_init)$varcomp["ENV:GENOTYPE!ENV!cor",1]

effects=mdl_CCR_init$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept_it"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept_it"],pred)

# get the slopes for z1
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1_it:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1_it:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope_it"
pred$slope_it=pred$slope_it/mean(pred$slope_it) ### standardizing slopes to mean(slopes)=1
lettuce=merge(lettuce[,colnames(lettuce)!="slope_it"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

# initialization done


# Start iterating

maxit=50 # maximum number of iterations
conv_crit=10^-4 # convergence criterion
dev_it=10^6 # initial value to check convergence
conv=F
i=1

while (i<=maxit & !conv){

  #### criss step

  lettuce$nitrate_it=lettuce$nitrate-lettuce$intercept_it # remove intercept from observed values

  # get starting values for criss model
  criss_i=asreml(fixed = nitrate_it~-1+x1:slope_it+x2:slope_it+x3:slope_it+x4:slope_it+x5:slope_it+x6:slope_it+x7:slope_it+x8:slope_it,
                 random=~ENV,
                 residual = ~ar1(ENV):id(GENOTYPE),
                 na.action = na.method(x="include", y = "include"),
                 workspace = "0.2gb",maxit=100,
                 data = lettuce,start.values = T)

  # fix variance parameters in criss step
  criss_i$vparameters.table["ENV",2:3]=c(v_e,"F")
  criss_i$vparameters.table["ENV:GENOTYPE!R",2:3]=c(v_r,"F")
  criss_i$vparameters.table["ENV:GENOTYPE!ENV!cor",2:3]=c(cor_r,"F")
  criss_i$vparameters.table$Value=as.numeric(criss_i$vparameters.table$Value)

  # run criss model
  criss=asreml(fixed = nitrate_it~-1+x1:slope_it+x2:slope_it+x3:slope_it+x4:slope_it+x5:slope_it+x6:slope_it+x7:slope_it+x8:slope_it,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb", maxit=100,
               data = lettuce,G.param = criss_i$vparameters.table,R.param = criss_i$vparameters.table)

  # get z1
  z1_it=t(rev(criss$coefficients$fixed[1:8]) %*% t(C))
  lettuce=merge(x = lettuce[,colnames(lettuce)!="z1_it"],y = data.frame(ENV=unique(lettuce0$ENV),z1_it=z1_it))
  lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



  #### cross step

  # run cross model
  cross=asreml(fixed = nitrate~-1+GENOTYPE+z1_it:GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",
               data = lettuce)

  # variance components estimated once for each iteration, in the cross step
  v_e=summary(cross)$varcomp["ENV",1]
  v_r=summary(cross)$varcomp["ENV:GENOTYPE!R",1]
  cor_r=summary(cross)$varcomp["ENV:GENOTYPE!ENV!cor",1]

  effects=cross$coefficients


  # get the intercepts
  pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
  pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
  colnames(pred)[1]="intercept_it"
  lettuce=merge(lettuce[,colnames(lettuce)!="intercept_it"],pred)

  # get the slopes for z1
  pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1_it:GENOTYPE_",value = T,fixed = T),]))
  pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1_it:GENOTYPE_",replacement = "")
  colnames(pred)[1]="slope_it"
  pred$slope_it=pred$slope_it/mean(pred$slope_it) ### standardizing slopes to mean(slopes)=1
  lettuce=merge(lettuce[,colnames(lettuce)!="slope_it"],pred)

  lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



  #### Check convergence

  diff=abs(dev_it-(-2*cross$loglik))

  if (diff<conv_crit){
    conv=T
  }else{
    dev_it=(-2*cross$loglik)
    i=i+1
  }
}
print(i)
print(diff)

# save last iteration estimates
lettuce$intercept=lettuce$intercept_it
lettuce$z1=lettuce$z1_it

lettuce=lettuce[,!(colnames(lettuce) %in% c("nitrate_it","intercept_it","z1_it","slope_it"))]

# get the slopes for z1
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1_it:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1_it:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope"
lettuce=merge(lettuce[,colnames(lettuce)!="slope"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



#### CCR2
rm(lettuce)
# initializing

mdl_add=asreml(fixed = nitrate~-1+GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",
               data = lettuce0)

# Initial value for z1 and z2 from residuals of additive model
Y=matrix(resid(mdl_add),nrow = 8,ncol = 18,byrow = F)
row.names(Y)=levels(lettuce0$GENOTYPE)

# Then get SVD of Y matrix
z1_it=svd(Y)$v[,1]
z2_it=svd(Y)$v[,2]

lettuce=merge(x = lettuce0[,colnames(lettuce0)!="z1_it"],y = data.frame(ENV=unique(lettuce0$ENV),z1_it=z1_it))
lettuce=merge(x = lettuce[,colnames(lettuce)!="z2_it"],y = data.frame(ENV=unique(lettuce0$ENV),z2_it=z2_it))
lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

mdl_CCR_init=asreml(fixed = nitrate~-1+GENOTYPE+z1_it:GENOTYPE+z2_it:GENOTYPE,
                    random=~ENV,
                    residual = ~ar1(ENV):id(GENOTYPE),
                    na.action = na.method(x="include", y = "include"),
                    workspace = "0.2gb",maxit=100,
                    data = lettuce)

v_e=summary(mdl_CCR_init)$varcomp["ENV",1]
v_r=summary(mdl_CCR_init)$varcomp["ENV:GENOTYPE!R",1]
cor_r=summary(mdl_CCR_init)$varcomp["ENV:GENOTYPE!ENV!cor",1]

effects=mdl_CCR_init$coefficients

# get the intercepts
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
colnames(pred)[1]="intercept_it"
lettuce=merge(lettuce[,colnames(lettuce)!="intercept_it"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

# Predict fitted values for multiplicative terms
lettuce$FV=predict(object = mdl_CCR_init,classify = "GENOTYPE:z1_it:z2_it",ignore = c("GENOTYPE"),levels=list(GENOTYPE=lettuce$GENOTYPE,z1_it=lettuce$z1_it,z2_it=lettuce$z2_it),parallel = T)$pvals[,"predicted.value"]
Y=reshape::cast(data = lettuce[,c("GENOTYPE","ENV","FV")],formula = GENOTYPE~ENV,value = "FV")[,-1]
row.names(Y)=levels(lettuce$GENOTYPE)

# Then get SVD of Y matrix
slope1_it=svd(Y)$u[,1]*(svd(Y)$d[1])
slope2_it=svd(Y)$u[,2]*(svd(Y)$d[2])

# get the slopes for z1
pred=as.data.frame(cbind(slope1_it,levels(lettuce$GENOTYPE)))
colnames(pred)=c("slope1_it","GENOTYPE")
pred$slope1_it=as.numeric(pred$slope1_it)
lettuce=merge(lettuce[,colnames(lettuce)!="slope1_it"],pred)

# get the slopes for z2
pred=as.data.frame(cbind(slope2_it,levels(lettuce$GENOTYPE)))
colnames(pred)=c("slope2_it","GENOTYPE")
pred$slope2_it=as.numeric(pred$slope2_it)
lettuce=merge(lettuce[,colnames(lettuce)!="slope2_it"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]

# initialization done



# Start iterating

maxit=50 # maximum number of iterations
conv_crit=10^-4 # convergence criterion
dev_it=10^6 # initial value to check convergence
conv=F
i=1

while (i<=maxit & !conv){

  #### criss step

  lettuce$nitrate_it=lettuce$nitrate-lettuce$intercept_it # remove intercept from observed values

  # get starting values for criss model
  criss_i=asreml(fixed = nitrate_it~-1+x1:slope1_it+x2:slope1_it+x3:slope1_it+x4:slope1_it+x5:slope1_it+x6:slope1_it+x7:slope1_it+x8:slope1_it+
                   x1:slope2_it+x2:slope2_it+x3:slope2_it+x4:slope2_it+x5:slope2_it+x6:slope2_it+x7:slope2_it+x8:slope2_it,
                 random=~ENV,
                 residual = ~ar1(ENV):id(GENOTYPE),
                 na.action = na.method(x="include", y = "include"),
                 workspace = "0.2gb",maxit=100,
                 data = lettuce,start.values = T)

  # fix variance parameters in criss step
  criss_i$vparameters.table["ENV",2:3]=c(v_e,"F")
  criss_i$vparameters.table["ENV:GENOTYPE!R",2:3]=c(v_r,"F")
  criss_i$vparameters.table["ENV:GENOTYPE!ENV!cor",2:3]=c(cor_r,"F")
  criss_i$vparameters.table$Value=as.numeric(criss_i$vparameters.table$Value)

  # run criss model
  criss=asreml(fixed = nitrate_it~-1+x1:slope1_it+x2:slope1_it+x3:slope1_it+x4:slope1_it+x5:slope1_it+x6:slope1_it+x7:slope1_it+x8:slope1_it+
                 x1:slope2_it+x2:slope2_it+x3:slope2_it+x4:slope2_it+x5:slope2_it+x6:slope2_it+x7:slope2_it+x8:slope2_it,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",maxit=100,
               data = lettuce,G.param = criss_i$vparameters.table,R.param = criss_i$vparameters.table)

  # Predict fitted values for multiplicative terms
  lettuce$FV=predict(object = criss,classify = "slope1_it:slope2_it:x1:x2:x3:x4:x5:x6:x7:x8",levels=list(x1=lettuce$x1,x2=lettuce$x2,x3=lettuce$x3,x4=lettuce$x4,x5=lettuce$x5,x6=lettuce$x6,x7=lettuce$x7,x8=lettuce$x8,slope1_it=lettuce$slope1_it,slope2_it=lettuce$slope2_it),parallel = T)$pvals[,"predicted.value"]
  Y=reshape::cast(data = lettuce[,c("GENOTYPE","ENV","FV")],formula = GENOTYPE~ENV,value = "FV")[,-1]
  row.names(Y)=levels(lettuce$GENOTYPE)

  # Then get SVD of Y matrix
  z1_it=svd(Y)$v[,1]
  z2_it=svd(Y)$v[,2]

  # get z1 and z2
  lettuce=merge(x = lettuce[,colnames(lettuce)!="z1_it"],y = data.frame(ENV=unique(lettuce0$ENV),z1_it=z1_it))
  lettuce=merge(x = lettuce[,colnames(lettuce)!="z2_it"],y = data.frame(ENV=unique(lettuce0$ENV),z2_it=z2_it))
  lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



  #### cross step

  # run cross model
  cross=asreml(fixed = nitrate~-1+GENOTYPE+z1_it:GENOTYPE+z2_it:GENOTYPE,
               random=~ENV,
               residual = ~ar1(ENV):id(GENOTYPE),
               na.action = na.method(x="include", y = "include"),
               workspace = "0.2gb",maxit=100,
               data = lettuce)

  # variance components estimated once for each iteration, in the cross step
  v_e=summary(cross)$varcomp["ENV",1]
  v_r=summary(cross)$varcomp["ENV:GENOTYPE!R",1]
  cor_r=summary(cross)$varcomp["ENV:GENOTYPE!ENV!cor",1]

  effects=cross$coefficients

  # get the intercepts
  pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "^GENOTYPE_",value = T,fixed = F),]))
  pred$GENOTYPE=gsub(x = row.names(pred),pattern = "GENOTYPE_",replacement = "")
  colnames(pred)[1]="intercept_it"
  lettuce=merge(lettuce[,colnames(lettuce)!="intercept_it"],pred)

  lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]


  # Predict fitted values for multiplicative terms
  lettuce$FV=predict(object = cross,classify = "GENOTYPE:z1_it:z2_it",ignore = c("GENOTYPE"),levels=list(GENOTYPE=lettuce$GENOTYPE,z1_it=lettuce$z1_it,z2_it=lettuce$z2_it),parallel = T)$pvals[,"predicted.value"]
  Y=reshape::cast(data = lettuce[,c("GENOTYPE","ENV","FV")],formula = GENOTYPE~ENV,value = "FV")[,-1]
  row.names(Y)=levels(lettuce$GENOTYPE)

  # Then get SVD of Y matrix
  slope1_it=svd(Y)$u[,1]*(svd(Y)$d[1])
  slope2_it=svd(Y)$u[,2]*(svd(Y)$d[2])

  # get the slopes for z1
  pred=as.data.frame(cbind(slope1_it,levels(lettuce$GENOTYPE)))
  colnames(pred)=c("slope1_it","GENOTYPE")
  pred$slope1_it=as.numeric(pred$slope1_it)
  lettuce=merge(lettuce[,colnames(lettuce)!="slope1_it"],pred)

  # get the slopes for z2
  pred=as.data.frame(cbind(slope2_it,levels(lettuce$GENOTYPE)))
  colnames(pred)=c("slope2_it","GENOTYPE")
  pred$slope2_it=as.numeric(pred$slope2_it)
  lettuce=merge(lettuce[,colnames(lettuce)!="slope2_it"],pred)

  lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]



  ### Check convergence

  diff=abs(dev_it-(-2*cross$loglik))
  if (diff<conv_crit){
    conv=T
  }else{
    dev_it=(-2*cross$loglik)
    i=i+1
  }
}
print(i)
print(diff)

# save last iteration estimates
lettuce$intercept=lettuce$intercept_it
lettuce$z1=lettuce$z1_it
lettuce$z2=lettuce$z2_it

lettuce=lettuce[,!(colnames(lettuce) %in% c("nitrate_it","intercept_it","z1_it","z2_it","slope1_it","slope2_it"))]

# get the slopes for z1
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z1_it:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z1_it:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope1"
lettuce=merge(lettuce[,colnames(lettuce)!="slope1"],pred)

# get the slopes for z2
pred=as.data.frame(cbind(effects$fixed[grep(x = row.names(effects$fixed),pattern = "z2_it:GENOTYPE_",value = T,fixed = T),]))
pred$GENOTYPE=gsub(x = row.names(pred),pattern = "z2_it:GENOTYPE_",replacement = "")
colnames(pred)[1]="slope2"
lettuce=merge(lettuce[,colnames(lettuce)!="slope2"],pred)

lettuce=lettuce[order(lettuce$ENV,lettuce$GENOTYPE),]




#################### END #########################################################################################################







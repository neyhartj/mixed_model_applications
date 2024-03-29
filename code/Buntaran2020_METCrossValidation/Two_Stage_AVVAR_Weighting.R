######################################
##### Two Stage AVVAR Weigthing #####
####################################

# R setup ###########
require(data.table) # data.table instead of data.frame format for coding efficiency
require(plyr)       # plyr::ldply() to collapse list after 1st stage
require(asreml)     # Fit model with ASReml-R 4
require(stringr)    # Split the character to provide factor's name
#####################

options("scipen"=100,"digits" =4 ) # set numbering format

##### Data input ######
ww <- read.csv("Winter Wheat 2016.csv",h=T)
ww <- na.omit(ww)

##### Change to factor #####
cols <- c("Alpha","Rep","Cultivar")
ww[cols] <- lapply(ww[cols], factor)

ww <- data.table(ww)

##### Make column Zone_Loc #####
ww$Zone_Loc <- factor(paste(ww$Zone, ww$Location, ww$Year,sep = "_"))

trials     <- nlevels(ww$Zone_Loc)
Envs <- levels(ww$Zone_Loc)

##### Make data list for Stage I #####
data_list <- matrix(data=list(), nrow=length(Envs), ncol=1, 
                    dimnames=list(Envs, c("data_Set")))

##### Make a list of Trials #####
for(i in 1:trials){
  print(i)
  b <- levels(ww$Zone_Loc)
  c <- b[i]
  Zone_Loc <- as.factor(c)
  Zone_Loc <- data.table(Zone_Loc)
  f <- merge(ww,Zone_Loc,by="Zone_Loc")
  assign(paste0("data_", b[i]), f)
  data_list[[i, "data_Set" ]] <- f
  
  rm(b, c, f, Zone_Loc)
}

data_list <- data.table(ldply(data_list[, "data_Set"], data.frame, .id="Zone_Loc"))

stgI_list <- matrix(data=list(), nrow=length(Envs), ncol=1, 
                    dimnames=list(Envs, c("lsmeans")))

asreml.options(maxit=100) # Set asreml iteration

############################
##### Stage I LSMEANS #####
##### per location   #####
for (i in Envs){
  
  Edat <- droplevels(subset(data_list, Zone_Loc==i))
  
  print(i)
  
  mod.1 <- asreml(fixed     = Yield ~ Cultivar,
                  random      = ~ Rep + Rep:Alpha,
                  data        = Edat,
                  predict     = predict.asreml(classify = "Cultivar"))
  
  update.asreml(mod.1)
  print(summary.asreml(mod.1)$varcomp)
  
  blue<- predict(mod.1, classify="Cultivar", levels=levels(Edat$Cultivar), vcov=TRUE)
  blue.1 <- data.table(blue$pvals)[, c(1:3)] 
  names(blue.1) <- c("Cultivar", "Yield_LSM", "se")
  
  V.j <- blue$vcov # get the variance-covariance matrix of Stage I analysis
  
  t.j <- nlevels(blue.1$Cultivar) # get the number of Varieties
  One <- matrix(rep(1, t.j), nrow = t.j) # Matrix consists of value 1
  
  trace.V.j <- sum(diag(blue$vcov)) # calculate the trace
  
  Vdiff.j <- (t.j*trace.V.j - t(One)%*%V.j%*%One)/(t.j*(t.j - 1)/2) # calculate the VDIFF
  
  R.j = as.numeric(0.5*(Vdiff.j))*diag(rep(1, t.j), nrow = t.j) # final calculation of AVVAR
  
  blue.1[ , ':='(var=se^2, avvar=diag(1/R.j))] # put the AVVAR on the list
  
  stgI_list[[i, "lsmeans" ]] <- blue.1 # put the all results on the list
  
  #rm(Edat,mod.1, blue, blue.1, V.j, t.j, trace.V.j, sed.j, R.j)
}

#######################################################
##### Preparing dataset of Stage I for Stage II ######

##### Unlist the results of Stage I and format as data.table #####
stgII_list <- data.table(ldply(stgI_list[, "lsmeans"], data.frame, .id="Zone_Loc"))

stgII_list$Zone <- factor(str_split_fixed(stgII_list$Zone_Loc, "_", 2)[,1]) # Make Zone column by split the record in Zone_Loc column
stgII_list$Location <- factor(str_split_fixed(stgII_list$Zone_Loc, "_", 2)[,2])   # Make Location by split the record in Zone_Loc column
stgII_list$Year <- factor(str_split_fixed(stgII_list$Zone_Loc, "_", 3)[,3])   # Make Year by split the record in Zone_Loc column


############################
##### Stage II BLUP ######
##### Zone analysis #####
mod.2 <- asreml(Yield_LSM  ~ Zone,
                random    = ~Cultivar + Zone:Location + Zone:Cultivar + Cultivar:Zone:Location,
                weights   = avvar,
                family    = asr_gaussian(dispersion=1.0),
                data      = stgII_list,
                predict   = predict.asreml(classify = "Cultivar:Zone"))

update.asreml(mod.2)

print(summary.asreml(mod.2)$varcomp) # print the variance components

blup.1<- data.table((mod.2$predictions$pvals[1:4])) # set the BLUP results as data.table

blup.1.a <- blup.1[order(Zone,-predicted.value),]  # Sort the results to see the highest yield in each zone
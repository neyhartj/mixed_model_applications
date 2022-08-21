######################################
##### Two Stage Fully-efficient #####
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
ww$Zone_Loc <- factor(paste(ww$Zone, ww$Location,sep = "_"))

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

stgI_list <- matrix(data=list(), nrow=length(Envs), ncol=2,
                    dimnames=list(Envs, c("lsmeans","Rvar")))

asreml.options(maxit=100) # Set asreml iteration

############################
##### Stage I LSMEANS #####
##### per location   #####
for (i in Envs){

  Edat <- droplevels(subset(ww, Zone_Loc==i))

  print(i)

  mod.1 <- asreml(fixed      = Yield ~ Cultivar,
                  random      = ~ Rep + Rep:Alpha,
                  data        = Edat,
                  predict     = predict.asreml(classify = "Cultivar"))

  update.asreml(mod.1)
  print(summary.asreml(mod.1)$varcomp)

  blue<- predict(mod.1, classify="Cultivar", levels=levels(Edat$Cultivar), vcov=TRUE)
  blue.1 <- data.table(blue$pvals)[, c(1:3)]
  names(blue.1) <- c("Cultivar", "Yield_LSM", "se")

  vcov  <- as.matrix(blue$vcov) # get the full variance-covariance matrix
  dimnames(vcov) <- list(blue.1$Cultivar,blue.1$Cultivar) # put the name of the variety in variance-covariance matrix
  vtab  <- melt(vcov)
  Rvar <- data.table(name=paste0("Zone_Loc_",i,"!Cultivar_",vtab$Var1,":",vtab$Var2), value=vtab$value) # variance-covariance matrix for stage II

  stgI_list[[i, "lsmeans" ]] <- blue.1  # put the all results on the list
  stgI_list[[i, "Rvar"]] <- Rvar        # put variance-covariance matrix of adj. means per genotype as named vector

}

#######################################################
##### Preparing dataset of Stage I for Stage II ######

##### Unlist the results of Stage I and format as data.table #####
lsm_stageI <- data.table(ldply(stgI_list[, "lsmeans"], data.frame, .id="Zone_Loc"))
lsm_stageI <- lsm_stageI[order(lsm_stageI$Zone_Loc,lsm_stageI$Cultivar),]

lsm_stageI$Zone <- factor(str_split_fixed(lsm_stageI$Zone_Loc, "_", 2)[,1])   # Make Zone column by split the record in Zone_Loc column
lsm_stageI$Location <- factor(str_split_fixed(lsm_stageI$Zone_Loc, "_", 2)[,2])    # Make Location by split the record in Zone_Loc column


# This part does not fit the model but  obtain the list object "R.param" to insert the variance-covariance from Stage I.
initmod2 <- asreml(fixed    = Yield_LSM ~ Zone,
                   random   = ~Cultivar + Zone:Location + Zone:Cultivar + Cultivar:Location,
                   residual = ~dsum(~us(Cultivar)|Zone_Loc),
                   data     = lsm_stageI,
                   family   = asr_gaussian(dispersion=1.0), # fix residual variance to 1
                   start.values = T)                       # ask for R.param and exit before model fitting


R.param <- initmod2$R.param #  Manipulate the object by putting in the stage I variance-covariance

for (i in Envs){
  DT.Rp   <- data.table(name=names(R.param[[i]]$Cultivar$initial))
  DT.Rp   <- stgI_list[[i,"Rvar"]][DT.Rp, on="name", nomatch=0] # Insert estimates from stage I variance-covaraince
  v.Rp    <- DT.Rp$value
  names(v.Rp) <- DT.Rp$name
  R.param[[i]]$Cultivar$initial <- v.Rp                   # replace default initial values
  R.param[[i]]$Cultivar$con     <- rep("F", length(v.Rp)) # set all components to "Fixed"
  rm(DT.Rp, v.Rp)
}

############################
##### Stage II BLUP ######
##### Zone analysis #####
##### Fit stage II model with R.param that filled with variance-covariance from the first stage ####
asreml.options(pworkspace="6gb",workspace="6gb") # Since we have large variance-covariance matrix, we need more space
mod.2 <- asreml(Yield_LSM  ~ Zone,
               random    = ~Cultivar + Zone:Location + Zone:Cultivar + Cultivar:Zone:Location,
               residual = ~dsum(~us(Cultivar)|Zone_Loc),
               data     = lsm_stageI,
               R.param  = R.param,
               predict=predict.asreml(classify = "Cultivar:Zone")) # Note that the error variance fixation to 1 is included in R.param

summary.asreml(mod.2)$varcomp

blup.1<- data.table((mod.2$predictions$pvals[1:4])) # set the BLUP results as data.table

blup.1.a <- blup.1[order(Zone,-predicted.value),]  # Sort the results to see the highest yield in each zone
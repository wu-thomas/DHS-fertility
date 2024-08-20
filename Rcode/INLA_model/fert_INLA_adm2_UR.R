#################################################
###### Load Packages
#################################################
#library(DHS.rates)

rm(list = ls())

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

library(INLA)
library(rgdal)
library(sp)

library(SUMMER)


library(ggplot2)

library(mapproj)
library(spdep)
library(sp)
library(RColorBrewer)
library(data.table)
library(survival)

library(rgdal)
library(dplyr)

################################################################
#########   set parameters
################################################################
country <- 'Tanzania'

## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-4)], collapse = "/")

data_dir<-paste0(home_dir,'/Data/')



setwd(paste(data_dir,country,sep=''))


info.name <- paste0(country, "_general_info.Rdata")

load(file = paste0(info.name, sep=''))


pop.year <- beg.year:end.year
pop.abbrev <- tolower(gadm.abbrev)

rw.order = 2

###############################################################################
#########   create folder to store results
###############################################################################

setwd(paste0(res.dir))

dir.create(file.path('.', 'INLA_model'))

setwd(paste0(res.dir,'/INLA_model/'))

dir.create(file.path('.', 'admin-1'))

dir.create(file.path('.', 'admin-2'))


################################################################
#########   load aggregation functions
################################################################

setwd(paste0(code.dir,'/INLA_model/'))
source('helper_functions_fert_post.R')



################################################################
######### load polygons
################################################################


setwd(paste0(data.dir,'/shapeFiles_gadm'))


country_shp_analysis <- readRDS('country_shp_analysis.rds')
country_shp_smoothed <- readRDS('country_shp_smoothed.rds')
admin1_info <- readRDS('admin1_info.rds')
admin2_info <- readRDS('admin2_info.rds')

poly.adm0 = country_shp_analysis[['National']]
poly.adm1 = country_shp_analysis[['Admin-1']]
poly.adm2 = country_shp_analysis[['Admin-2']]

################################################################
#########   prepare neighborhood graph
################################################################

Amat <- getAmat(poly.adm2, poly.adm2$REGNAME)

## check whether there is island, if so, attach it to the nearest mainland admin region
## by changing the affinity matrix
which(rowSums(Amat)==0)

### for example attach mafia to Rufiji in Tanzania 
Amat[131,133] <-  Amat[133,131] <- 1


sp.prec.mx <- as.matrix(-Amat)


diag(sp.prec.mx) <- -rowSums(sp.prec.mx)

#quick check, sum to zero
apply(sp.prec.mx,1,sum)
apply(sp.prec.mx,2,sum)

N.area <- dim(sp.prec.mx)[1]
scaled.sp.prec.mx <- INLA::inla.scale.model(sp.prec.mx, constr = list(A = matrix(1, 
                                                                                 1, dim(sp.prec.mx)[1]), e = 0))



###############################################################################
#########   load prepared data
###############################################################################

setwd(paste0(data.dir))
load(file = paste0('prepared_IR_dat/fert_cluster_dat_',survey_year,'.rda'))

fert_pyears$agegrp.int<- as.integer(as.numeric(substr(fert_pyears$agegr,1,2))/5-2)

admin2.names <- admin2_info$data
admin2.names <- admin2.names[,c('admin2.name','admin2.char')]

fert_pyears <- merge(fert_pyears,admin2.names,by='admin2.name')

fert_pyears$region.num <- as.integer(substr(fert_pyears$admin2.char,8,9))

# model setup
natl_fert_clust <- fert_pyears
natl_fert_clust <- natl_fert_clust[natl_fert_clust$period %in% c(beg.year:end.year),]

natl_fert_clust$region.int <- as.integer(natl_fert_clust$region.num)
natl_fert_clust$time.int <-  as.integer(as.character(natl_fert_clust$period))-
     min(as.integer(as.character(natl_fert_clust$period)))+1



natl_fert_clust$cutoff_bias <- 0
natl_fert_clust[natl_fert_clust$period==(survey_year-6),]$cutoff_bias <- (-1)
natl_fert_clust[natl_fert_clust$period==(survey_year-5),]$cutoff_bias <- 1



N.time <- length(unique(natl_fert_clust$period))
N.agegrp <- length(unique(natl_fert_clust$agegrp))



###############################################################################
#########   prepare space agegroup grid
###############################################################################



x <- expand.grid(1:N.agegrp, 1:N.area)
agegrp.area <- data.frame(agegrp.int = x[, 1], region.int = x[, 2], agegrp.area = c(1:nrow(x)))


natl_fert_clust <- natl_fert_clust %>% right_join(agegrp.area, c("agegrp.int", "region.int"))

natl_fert_clust <-natl_fert_clust[order(natl_fert_clust$agegrp.int,
                                        natl_fert_clust$time.int, natl_fert_clust$region.int),]




###############################################################################
#########   prepare space time grid
###############################################################################



x <- expand.grid(1:N.time, 1:N.area)
time.area <- data.frame(time.int = x[, 1], region.int = x[, 2], time.area = c(1:nrow(x)))


natl_fert_clust <- natl_fert_clust %>% right_join(time.area, c("time.int", "region.int"))

natl_fert_clust <-natl_fert_clust[order(natl_fert_clust$agegrp.int,
                                        natl_fert_clust$time.int, natl_fert_clust$region.int),]


###############################################################################
#########   prepare agegroup time grid
###############################################################################



x <- expand.grid(1:N.agegrp, 1:N.time)
agegrp.time <- data.frame(agegrp.int = x[, 1], time.int = x[, 2], agegrp.time = c(1:nrow(x)))

natl_fert_clust <- natl_fert_clust %>% right_join(agegrp.time, c("agegrp.int", "time.int"))

natl_fert_clust <-natl_fert_clust[order(natl_fert_clust$agegrp.int,
                                        natl_fert_clust$time.int, natl_fert_clust$region.int),]

natl_fert_clust$agegrp.time <- as.integer(natl_fert_clust$agegrp.time)


natl_fert_clust_U <- natl_fert_clust[natl_fert_clust$urban=='urban',] ## check urban/rural vs U/R
natl_fert_clust_R <- natl_fert_clust[natl_fert_clust$urban=='rural',]

natl_fert_clust_U <- natl_fert_clust_U[complete.cases(natl_fert_clust_U[,c(1:10)]),]
natl_fert_clust_R <- natl_fert_clust_R[complete.cases(natl_fert_clust_U[,c(1:10)]),]


###############################################################################
#########   prepare temporal and mother's age group rw
###############################################################################



inla.rw = utils::getFromNamespace("inla.rw", 
                                  "INLA")
RW.time.prec <- inla.rw(n = N.time, order = rw.order, scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                        sparse = TRUE)
RW.agegrp.prec <- inla.rw(n = N.agegrp, order = rw.order, scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                          sparse = TRUE)

scaled.RW.time.prec <- inla.rw(n = N.time, order = rw.order, scale.model = T, # set scale.model  = F because we'll scale in the formula
                               sparse = TRUE)
scaled.RW.agegrp.prec <- inla.rw(n = N.agegrp, order = rw.order , scale.model = T, # set scale.model  = F because we'll scale in the formula
                                 sparse = TRUE)
### hyperparameters
mu0 <- log(0.1/0.9)
sig2.0 <- 10000
# Prior parameters for the gamma on the precisions.
a <- 1
b <- 0.01

BYM2Prior <- list(
  phi = list(
    prior = "pc",
    param = c(0.5, 2/3)),
  prec = list(
    prior = "pc.prec",
    param = c(a, b)))

nSamp <- 1000

mu.bias <- 0.05
sig.bias <- sqrt(0.05) # consider a smaller variance if model not converging 


##################################################################################################################
######  Different spatial field for U/R, typeIV space and agegroup * UR
##################################################################################################################





### space x agegroup 
# Kronecker product between ICAR x RW2
R.space.age <- scaled.RW.agegrp.prec %x% scaled.sp.prec.mx

int.eigen.space.age <- eigen(R.space.age)
n.constr.space.age <- N.agegrp*N.area -(N.agegrp - rw.order)*(N.area - 1)

A.space.age  <- t(matrix(int.eigen.space.age$vectors[,(dim(R.space.age)[1]-n.constr.space.age+1):dim(R.space.age)[1]],
                         ncol = n.constr.space.age))

constr.space.age  <- list(A = A.space.age, e = rep(0, dim(A.space.age)[1]))

### space x time 
# Kronecker product between ICAR x RW2
R.space.time <- scaled.RW.time.prec %x% scaled.sp.prec.mx

int.eigen.space.time <- eigen(R.space.time)
n.constr.space.time  <- N.time*N.area -(N.time - rw.order)*(N.area - 1)

A.space.time   <- t(matrix(int.eigen.space.time$vectors[,(dim(R.space.time )[1]-n.constr.space.time +1):dim(R.space.time )[1]],
                         ncol = n.constr.space.time ))

constr.space.time   <- list(A = A.space.time , e = rep(0, dim(A.space.time )[1]))

### agegroup x time 
# Kronecker product between RW2 x RW2
R.age.time <- scaled.RW.agegrp.prec %x% scaled.RW.time.prec 

int.eigen.age.time <- eigen(R.age.time)
n.constr.age.time  <- N.time*N.agegrp -(N.time - rw.order)*(N.agegrp - rw.order)

A.age.time   <- t(matrix(int.eigen.age.time$vectors[,(dim(R.age.time )[1]-n.constr.age.time +1):dim(R.age.time )[1]],
                           ncol = n.constr.age.time ))

constr.age.time   <- list(A = A.age.time , e = rep(0, dim(A.age.time )[1]))

#pc.u = 1
#pc.alpha = 0.01
#hyperpc1 <- list(prec = list(prior = "pc.prec", param = c(pc.u , pc.alpha)))
#natl_fert_clust_R$time.int.rep <- natl_fert_clust_R$time.int

formula.full.int <- event ~ 1 + cutoff_bias +
  f(time.int, model = "bym", scale.model = TRUE, constr = TRUE,
    rankdef = rw.order, # good practice to specify rank deficiency of precision matrix
    graph = RW.time.prec)+
  #  hyper= c(0.0005,3)) +\
  #f(time.int,model='rw2', constr = TRUE,  extraconstr = NULL, hyper = hyperpc1)      +
  #f(time.int.rep,model='iid', constr = TRUE,  extraconstr = NULL, hyper = hyperpc1)      +
  f(region.int, model = "bym2", 
    scale.model = TRUE,
    constr = TRUE,
    rankdef = 1, # good practice to specify rank deficiency of precision matrix
    graph = sp.prec.mx,
    hyper = BYM2Prior) +
  f(agegrp.int, model = "bym", 
    scale.model = TRUE,
    constr = TRUE,
    rankdef = rw.order, # good practice to specify rank deficiency of precision matrix
    graph = RW.agegrp.prec)+
  f(agegrp.area, 
    model = "generic0", Cmatrix = R.space.age, extraconstr = constr.space.age,
    rankdef = n.constr.space.age, 
    param = c(1, 0.01))+
  f(time.area, 
    model = "generic0", Cmatrix = R.space.time, extraconstr = constr.space.time,
    rankdef = n.constr.space.time, 
    param = c(1, 0.01))+
  f(agegrp.time, 
    model = "generic0", Cmatrix = R.age.time, extraconstr = constr.age.time,
    rankdef = n.constr.age.time,
    param = c(1, 0.01))
  
  


mod.R.full.int <- inla(formula.full.int,
                     data = natl_fert_clust_R,
                     family = "nbinomial",
                     E = pyears,
                     control.fixed = list(mean.intercept = c(default=mu0),
                                          mean=list(cutoff_bias=mu.bias, default=mu0),
                                          prec.intercept = c(default=1/sig2.0),
                                          prec=list(cutoff_bias=1/sig.bias, default=1/sig2.0)),
                           control.predictor = list(compute = TRUE),
                           control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
                           verbose=T)



mod.U.full.int <- inla(formula.full.int,
                       data = natl_fert_clust_U,
                       family = "nbinomial",
                       E = pyears,
                       control.fixed = list(mean.intercept = c(default=mu0),
                                            mean=list(cutoff_bias=mu.bias, default=mu0),
                                            prec.intercept = c(default=1/sig2.0),
                                            prec=list(cutoff_bias=1/sig.bias, default=1/sig2.0)),
                       #control.fixed = list(mean.intercept = c(mu0),
                        #                    prec.intercept = c(1/sig2.0)),
                       control.predictor = list(compute = TRUE),
                       control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
                       verbose=T)

U.res.inla <- mod.U.full.int
R.res.inla <- mod.R.full.int

summary(mod.U.full.int)
summary(mod.R.full.int)


setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))
saveRDS(mod.U.full.int,file='fitted_U_typeIV_rw2.rds')
saveRDS(mod.R.full.int,file='fitted_R_typeIV_rw2.rds')

################################################################
#########   set up urban fraction data frame
################################################################

setwd(paste0(data_dir,country,'/UR/Fractions'))

agegrp_admin2_urban<- readRDS(paste0('agegrp_admin2_urban_weights.rds'))
colnames(agegrp_admin2_urban)  <- c('region','year','urban_frac','agegrp','agegrp.int','rural_frac')

agegrp_admin2_urban <- agegrp_admin2_urban[agegrp_admin2_urban$year %in% c(beg.year:end.year),]
agegrp_admin2_urban$time.int = as.integer(as.factor(agegrp_admin2_urban$year))
agegrp_admin2_urban$region.int = as.integer(word(agegrp_admin2_urban$region,2,sep = "_"))



N.agegrp <- length(unique(agegrp_admin2_urban$agegrp.int))
N.area <- length(unique(agegrp_admin2_urban$region.int))
N.time <- length(unique(agegrp_admin2_urban$time.int))

x <- expand.grid(1:N.agegrp, 1:N.area)
agegrp.area <- data.frame(agegrp.int = x[, 1], region.int = x[, 2], agegrp.area = c(1:nrow(x)))

agegrp_admin2_urban <- agegrp_admin2_urban %>% right_join(agegrp.area, c("agegrp.int", "region.int"))



x <- expand.grid(1:N.time, 1:N.area)
time.area <- data.frame(time.int = x[, 1], region.int = x[, 2], time.area = c(1:nrow(x)))
agegrp_admin2_urban <- agegrp_admin2_urban %>% right_join(time.area, c("time.int", "region.int"))


x <- expand.grid(1:N.agegrp, 1:N.time)
agegrp.time <- data.frame(agegrp.int = x[, 1], time.int = x[, 2], agegrp.time = c(1:nrow(x)))
agegrp_admin2_urban <- agegrp_admin2_urban %>% right_join(agegrp.time, c("agegrp.int", "time.int"))



################################################################
#########   aggregation
################################################################


res.admin2.aggre <- aggUR.diffmod.typeIV(U.res.inla=mod.U.full.int,R.res.inla=mod.R.full.int, 
                                         nSamp = 1000,urb.frac=agegrp_admin2_urban)

setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

saveRDS(res.admin2.aggre,file='aggreUR_diffmod_typeIV_rw2.rds')

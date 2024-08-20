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

admin2.names <- admin2_info$data
admin2.names <- admin2.names[,c('admin2.name','admin2.char')]
admin2.names$GADM <- admin2.names$admin2.name
admin2.names$Internal <- admin2.names$admin2.char

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
admin2.names$GADM <- admin2.names$admin2.name
admin2.names$Internal <- admin2.names$admin2.char

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
#########   aggregate TFR, admin2, 3-year period 
###############################################################################

setwd(data_dir)
#adm1.pop.frame
load(file= paste(country,'/','worldpop/pop_frac/', 'admin2_tf_pop_frac', '.rda', sep = ''))  

## load U/R results
setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

UR.res <- readRDS(paste0('aggreUR_diffmod_typeIV_rw2.rds'))



#tmp.check2 <- UR.res2$Combined.est
## set up grid

admin.list<-admin2.names$admin2.char

est.grid <- UR.res$Combined.est[,c('region','year','agegrp')]
est.grid <- inner_join(est.grid, adm2.pop.frame, by = c("region" = "admin2_idx", "year" = "year"))
est.grid$wt <- est.grid$admin2_pop


res.adm2.tfr.3yrs <-data.frame()


for (last.year in c(survey_year,survey_year-3,survey_year-6)) {
  
  last.period <- paste0(last.year-2,'_',last.year,sep='')
  last.years <- c((last.year-2):last.year)
  
  res.overall.all <-data.frame()
  res.overall.draws <- matrix (0, ncol = length(admin.list) ,nrow= nSamp)
  
  ## prepare TFR
  
  for ( i in 1:length(admin.list)){
    
    
    overall.aggre.trf <- cal.aggre.TFR.wt( aggre.draws= UR.res$Combined.draws,nSamp = 1000,
                                           est.grid= est.grid,period= last.years,alpha=0.1)
    
    res.overall.all <- rbind(res.overall.all, overall.aggre.trf$Combined.est[i,])
    res.overall.draws[,i] <- overall.aggre.trf$Combined.draws[,i]
    
    
    
  }
  
  res.overall.all$period <- paste0(last.year-2,'-',last.year,sep='')
  res.adm2.tfr.3yrs <- rbind(res.adm2.tfr.3yrs,res.overall.all)
  
  
}


res.adm2.tfr.3yrs %>%
  group_by(period) %>%
  summarise_at(vars(p_Med), list(name = mean))




setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

write.csv(res.adm2.tfr.3yrs,row.names = F,
          paste0(country, '_stratified_admin2_tfr.csv'))



setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

saveRDS(res.adm2.tfr.3yrs,file=paste0('admin2_tfr_UR_diff_res.rds'))


###############################################################################
#########   visualize TFR, admin2
###############################################################################

#library(sf)
toplot.poly.adm2 <-  poly.adm2


toplot.poly.adm2@data$Internal<-admin2.names$admin2.names[c(1:dim(admin2.names)[1])]


g.adm2.tfr.3yrs <- mapPlot(data = res.adm2.tfr.3yrs, geo = toplot.poly.adm2, variables = c("period"), 
                           size=0.1,ylim= c(min(res.adm2.tfr.3yrs$p_Med),max(res.adm2.tfr.3yrs$p_Med)),
                           values = c("p_Med"), by.data = "region", by.geo = "Internal",
                           direction = -1, legend.label = "TFR",
                           is.long = TRUE, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='TFR',
                                label.position = "bottom"))

setwd(paste0(res_dir,country,'/Visualization'))

ggsave(g.adm2.tfr.3yrs, width=11, height = 10, file = paste0("admin2_TFR_3yrs.pdf"))




###############################################################################
#########   aggregate UR specific TFR, admin2
###############################################################################

setwd(data_dir)
#adm1.pop.frame
load(file= paste(country,'/','worldpop/pop_frac/', 'admin2_tf_pop_frac', '.rda', sep = ''))  

## load U/R results
setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

UR.res <- readRDS(paste0('aggreUR_diffmod_typeIV_rw2.rds'))

## set up grid

admin.list<-admin2.names$Internal

est.grid <- UR.res$Combined.est[,c('region','year','agegrp')]
est.grid <- inner_join(est.grid, adm2.pop.frame, by = c("region" = "admin2_idx", "year" = "year"))
est.grid$wt <- est.grid$admin2_pop


res.UR.diff.3yrs <-data.frame()


for (last.year in c(survey_year,survey_year-3,survey_year-6)) {
  
  last.period <- paste0(last.year-2,'_',last.year,sep='')
  last.years <- c((last.year-2):last.year)
  
  res.U.all <-data.frame()
  res.U.draws <- matrix (0, ncol = length(admin.list) ,nrow= nSamp)
  
  res.R.all <-data.frame()
  res.R.draws <- matrix (0, ncol = length(admin.list) ,nrow= nSamp)
  
  res.UR.diff.all <-data.frame()
  res.UR.diff.draws <- matrix (0, ncol = length(admin.list) ,nrow= nSamp)
  
  
  ## prepare TFR
  
  for ( i in 1:length(admin.list)){
    
    
    U.aggre.trf <- cal.aggre.TFR.wt( aggre.draws= UR.res$Urban.draws,nSamp = 1000,
                                     est.grid= est.grid,period= last.years,alpha=0.1)
    R.aggre.trf <- cal.aggre.TFR.wt( aggre.draws= UR.res$Rural.draws,nSamp = 1000,
                                     est.grid= est.grid,period= last.years,alpha=0.1)
    UR.diff.trf <- cal.aggre.TFR.wt( aggre.draws= UR.res$Rural.draws-UR.res$Urban.draws,nSamp = 1000,
                                     est.grid= est.grid,period= last.years,alpha=0.1)
    
    res.U.all <- rbind(res.U.all, U.aggre.trf$Combined.est[i,])
    res.U.draws[,i] <- U.aggre.trf$Combined.draws[,i]
    
    res.R.all <- rbind(res.R.all, R.aggre.trf$Combined.est[i,])
    res.R.draws[,i] <- R.aggre.trf$Combined.draws[,i]
    
    res.UR.diff.all <- rbind(res.UR.diff.all, UR.diff.trf$Combined.est[i,])
    res.UR.diff.draws[,i] <- UR.diff.trf$Combined.draws[,i]
    
    #cv.res.all <- rbind(cv.res.all, cv.aggre.trf$Combined.est[i,])
    #cv.res.draws[,i] <- cv.aggre.trf$Combined.draws[,i]
    
    
  }
  
  res.UR.diff.all$period <- paste0(last.year-2,'-',last.year,sep='')
  res.UR.diff.3yrs <- rbind(res.UR.diff.3yrs,res.UR.diff.all)
  
  
}


res.UR.diff.3yrs %>%
  group_by(period) %>%
  summarise_at(vars(p_Med), list(name = mean))





setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

saveRDS(res.UR.diff.3yrs,file=paste0('admin2_tfr_UR_diff_res.rds'))

###############################################################################
#########   visualize UR difference TFR, admin2
###############################################################################


toplot.poly.adm2 <- poly.adm2

toplot.poly.adm2@data$Internal<-admin2.names$Internal[c(1:dim(admin2.names)[1])]


g.tfr.ur.diff.3yrs <- mapPlot(data = res.UR.diff.3yrs, geo = toplot.poly.adm2, variables = c("period"), 
                              size=0.1,ylim= c(min(res.UR.diff.3yrs$p_Med),max(res.UR.diff.3yrs$p_Med)),
                              values = c("p_Med"), by.data = "region", by.geo = "Internal",
                              direction = -1, legend.label = "Difference in U/R TFR",
                              is.long = TRUE, ncol = 3)+
  theme (legend.position = 'bottom',legend.key.height=unit(0.5,'cm'),
         legend.text=element_text(size=14),
         legend.key.width = unit(2,'cm'),legend.title = element_text(size=18),
         strip.text.x = element_text(size = 15))+
  guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                title='Difference in U/R TFR (rural minus urban)',
                                label.position = "bottom"))

setwd(paste0(res_dir,country,'/Visualization'))

ggsave(g.tfr.ur.diff.3yrs, width=11, height = 10, file = paste0("admin2_TFR_UR_diff_3yrs.pdf"))


###############################################################################
#########   visualize 15-19 ASFR, admin2
###############################################################################

setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

adm2.asfr.overall.3yrs <- readRDS(file=paste0('asfr_overall_3yrs.rds'))

adm2.15.19.overall.3yrs <- adm2.asfr.overall.3yrs[adm2.asfr.overall.3yrs$agegrp=='15-19',]

#library(sf)
toplot.poly.adm2 <-  poly.adm2


toplot.poly.adm2@data$Internal<-admin2.names$Internal[c(1:dim(admin2.names)[1])]



###############################################################################
#########   visualize 15-19 ASFR, admin2, trend
###############################################################################

setwd(paste0(res_dir,country,'/INLA_model/admin-2/'))

adm2.asfr.yearly.res <- readRDS(file=paste0('aggreUR_diffmod_typeIV_rw2.rds'))
adm2.asfr.yearly <- adm2.asfr.yearly.res$Combined.est
adm2.asfr.15.19.yearly <- adm2.asfr.yearly[adm2.asfr.yearly$agegrp=='15-19',]
adm2.asfr.15.19.yearly <- merge(adm2.asfr.15.19.yearly,admin2.names,by.x='region',by.y='Internal',all.x=T)
adm2.asfr.15.19.yearly$region.name <- adm2.asfr.15.19.yearly$GADM

my.dodge <- position_dodge(width = 0.15)

g.adm2.asfr.15.19.yearly <- ggplot(aes(x = year, color = region.name), data = adm2.asfr.15.19.yearly)
g.adm2.asfr.15.19.yearly <- g.adm2.asfr.15.19.yearly + geom_point(aes(y = p_Med), position = my.dodge)

g.adm2.asfr.15.19.yearly <-  g.adm2.asfr.15.19.yearly+ 
  geom_line(aes(y = p_Med, color = region.name,group = region.name), position = my.dodge, alpha = 0.5)

g.adm2.asfr.15.19.yearly <- g.adm2.asfr.15.19.yearly+  theme_bw()+theme (legend.position = 'bottom',
                           legend.text=element_text(size=15),
                           legend.title = element_text(size=18),
                           strip.text.x = element_text(size = 15),
                           axis.text.x = element_text(size=15),
                           axis.title.x = element_text(size=16),
                           axis.title.y = element_text(size=16),
                           axis.text.y = element_text(size=15))+
  ylab("ASFR (per women-year)") + xlab("year") + 
  guides(color=guide_legend("Admin-2 regions"))+
  scale_x_continuous(breaks=c(2010:2020))



setwd(paste0(res_dir,country,'/Visualization'))

ggsave(g.adm2.asfr.15.19.yearly, width=12, height = 9, file = paste0("admin2_15_19_ASFR_yearly.pdf"))


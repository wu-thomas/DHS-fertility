rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.

# Load libraries and info ----------------------------------------------------------

options(gsubfn.engine = "R")
#library(rgdal)
library(raster)
#library(rgeos)
library(sqldf)
library(geosphere)
library(Matrix)
library(openxlsx)
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


svy_year <- survey_year


pop.year <- beg.year:end.year
pop.abbrev <- tolower(gadm.abbrev)


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
##############################################################################
#########   create population file folders
##############################################################################

setwd(paste0(data_dir,country))

dir.create(file.path('.', 'UR'))

dir.create(file.path('.', 'UR/Fractions'))




##############################################################################
######### Function to calculate subnational urban fractions
##############################################################################
#' Function to calculate subnational urban fractions
#' 
#' @param natl.grid A national grid
#' must contain a column indicating completeness
#' @examples
#' 

get.urb.frac <- function(natl.grid,pop.ras,urb.ras,admin.poly=NULL,admin.level=NULL,
                         admin.name=NULL){
  

  # assign population density
  natl.grid$pop_den<-raster::extract(pop.ras, natl.grid[c('x','y')])
  
  # assign urban rural predicted probability
  natl.grid$urban_prob<-raster::extract(urb.ras, natl.grid[c('x','y')])
  
  # exclude pixels without classification
  natl.grid[is.na(natl.grid$urban_prob),]$pop_den <- NA    # omit pixels without classification
  
  # calculate national fraction
  if(is.null(admin.poly)){
    natl.frac <- sum(natl.grid$urban_prob*natl.grid$pop_den,na.rm = TRUE)/
      sum(natl.grid$pop_den,na.rm=TRUE)
    return(natl.frac)  
  }
  
  
 
  natl.grid$adm.char<- natl.grid[,c(paste0(admin.level,'.char')),] 

  # calculate fractions
  adm.urb.frac <- natl.grid %>%
    group_by(adm.char) %>%
    summarise(Frac = sum(pop_den*urban_prob, na.rm = TRUE)/sum(pop_den,na.rm = TRUE))
  
  adm.urb.frac <- adm.urb.frac[!is.na(adm.urb.frac$adm.char),]
  
  # assign GADM names
  if(!is.null(admin.name)){
    adm.urb.frac <- adm.urb.frac[!is.na(adm.urb.frac$adm.char),]
    
    match_id <- match(admin.name$Internal,adm.urb.frac$adm.char)
    adm.urb.frac <- adm.urb.frac[match_id,]
    adm.urb.frac$admin.name <- admin.name$GADM
  }
  
  return(adm.urb.frac)
}




##############################################################################
#########   load classifications
##############################################################################

setwd(paste0(data_dir,country,'/UR',sep=''))

uncrc.urb.ras <- raster(paste0(country,'_UR_ind_1km.tif'))

setwd(paste(data_dir,country,sep=''))


##############################################################################
#########   load total population surface
##############################################################################

setwd(paste0(data_dir,country,'/worldpop'))

# UNadjusted population counts
worldpop <- raster(paste0(country.abbrev,
                          '_ppp_',frame_year,
                          '_1km_Aggregated_UNadj.tif',sep=''))


##############################################################################
#########   create national grid
##############################################################################

if(FALSE){
## set up grid
urb_dat<-as.data.frame(coordinates(worldpop))
colnames(urb_dat)<-c('x','y')


## admin 1 region 
points.frame <- as.data.frame(urb_dat[,c("x", "y")])
points.frame <- SpatialPoints(points.frame)
poly.over.adm1 <- SpatialPolygons(poly.adm1@polygons)
admin1.key <- over(points.frame, poly.over.adm1)
urb_dat$admin1<-admin1.key
urb_dat$admin1<-as.factor(urb_dat$admin1)
urb_dat$admin1.char<-paste0('admin1_',as.character(urb_dat$admin1))
urb_dat[is.na(urb_dat$admin1),]$admin1.char <- NA
urb_dat$admin1.name <- as.character(poly.adm1@data$NAME_1)[urb_dat$admin1]

## admin 2 region 
points.frame <- as.data.frame(urb_dat[,c("x", "y")])
points.frame <- SpatialPoints(points.frame)
poly.over.adm2 <- SpatialPolygons(poly.adm2@polygons)
admin2.key <- over(points.frame, poly.over.adm2)
urb_dat$admin2<-admin2.key
urb_dat$admin2<-as.factor(urb_dat$admin2)
urb_dat$admin2.char<-paste0('admin2_',as.character(urb_dat$admin2))
urb_dat[is.na(urb_dat$admin2),]$admin2.char <- NA
urb_dat$admin2.name <-  as.character(poly.adm2@data$NAME_2)[urb_dat$admin2]

urb_dat$pixel_idx <- c(1:dim(urb_dat)[1])


## save data
setwd(paste0(data_dir,country,sep=''))

dir.create(file.path('.', 'prepared_dat'))

# save(urb_dat,file="prepared_dat/natl_grid.rda")
}

setwd(paste0(data_dir,country,sep=''))
load(file="prepared_dat/natl_grid.rda")


##############################################################################
#########  calculate subnational urban fractions
##############################################################################

years <- c(beg.year:end.year)
adm1.weight.frame <- data.frame()
adm2.weight.frame <- data.frame()

age_group <- c(3:9)*5


for ( t in 1:length(years)){
  print(t)
  
  for(age_int in 1:7){
  print(age_int)
  year <- years[t]
  
  # load U5 population at year t
  setwd(paste0(data_dir,country,'/worldpop/pop_1km'))
  
  age <- age_group[age_int]
  f_age_pop<-raster(paste0(pop.abbrev,'_', age, '_',year,'_1km.tif'))
  
  # admin1 urban fraction for U5 population at year t
  f_age_urb_admin1<-get.urb.frac(natl.grid = urb_dat,
                              pop.ras=f_age_pop,
                              urb.ras=uncrc.urb.ras,
                              admin.poly=poly.adm1,
                              admin.level='admin1',
                              admin.name = admin1.names)
  
  f_age_urb_admin1$years <- year
  f_age_urb_admin1$agegrp <- paste0(age,'-',age+4,sep='')
  f_age_urb_admin1$agegrp.int <- age_int
  
  adm1.weight.frame <- rbind(adm1.weight.frame,f_age_urb_admin1)
  
  # admin2 urban fraction for U5 population at year t
  f_age_urb_admin2<-get.urb.frac(natl.grid = urb_dat,
                                 pop.ras=f_age_pop,
                              urb.ras=uncrc.urb.ras,
                              admin.poly=poly.adm2,
                              admin.level='admin2',
                              admin.name = admin2.names)
  
  f_age_urb_admin2$years <- year
  f_age_urb_admin2$agegrp <- paste0(age,'-',age+4,sep='')
  f_age_urb_admin2$agegrp.int <- age_int
  adm2.weight.frame <- rbind(adm2.weight.frame,f_age_urb_admin2)
  }
  #setwd(res_dir)
  
  # save calculated urban fractions
  #saveRDS(u5_urb_admin1,file=paste0('UR/U5_fraction/','admin1_',
  #                                  year, '_urban_frac_bart.rds'))
 # saveRDS(u5_urb_admin2,file=paste0('UR/U5_fraction/','admin2_',
  #                                  year, '_urban_frac_bart.rds'))
  
}


# process admin 1 urban rural weights data frame
adm1.weight.frame <- adm1.weight.frame[,c('adm.char','years','Frac','agegrp','agegrp.int')]
colnames(adm1.weight.frame) <- c('region','year','urban','agegrp','agegrp.int')
adm1.weight.frame$rural <- 1 - adm1.weight.frame$urban

# process admin 2 urban rural weights data frame
adm2.weight.frame <- adm2.weight.frame[,c('adm.char','years','Frac','agegrp','agegrp.int')]
colnames(adm2.weight.frame) <- c('region','year','urban','agegrp','agegrp.int')
adm2.weight.frame$rural <- 1 - adm2.weight.frame$urban

# save weights frames

setwd(paste0(data_dir,country,'/UR/Fractions'))

colnames(adm1.weight.frame)  <- c('region','year','urban_frac','agegrp','agegrp.int','rural_frac')
colnames(adm2.weight.frame)  <- c('region','year','urban_frac','agegrp','agegrp.int','rural_frac')

saveRDS(adm1.weight.frame,paste0('agegrp_admin1_urban_weights.rds'))
saveRDS(adm2.weight.frame,paste0('agegrp_admin2_urban_weights.rds'))


#adm1.weight.frame$year <- adm1.weight.frame$years
#adm2.weight.frame$year <- adm2.weight.frame$years
#colnames(adm1.weight.frame)  <- c('region','year','urban_frac','agegrp','agegrp.int','rural_frac')
#colnames(adm2.weight.frame)  <- c('region','year','urban_frac','agegrp','agegrp.int','rural_frac')
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

################################################################
#########   prepare function and dictionary (each row is one admin2 region)
################################################################

# link accross admin1 and admin2
#adm_link <- as.data.frame(poly.adm2)
#adm_link <- adm_link[,c('NAME_1','NAME_2')]
#colnames(adm_link)<-c('admin1_name','admin2_name')

# get unique ID for admin2
adm_link <- admin2_info$data[,c('admin1.name','admin2.name','admin2.char')]
colnames(adm_link) <- c('admin1_name','admin2_name','admin2_idx')
#adm_link$admin2_idx<- admin2_info$data$admin2.char

# get admin1 index
adm1_match<-match(adm_link$admin1_name, admin1_info$data$admin1.name)
adm_link$admin1_idx<- admin1_info$data$admin1.char[adm1_match]

# not using below because of repeated admin2 name
# adm2_match<-match(adm_link$admin2_name, admin2.names$GADM)#
# adm_link$admin2_idx<- admin2.names$Internal[adm2_match]

# function that calculates population in each admin2 area
pop_adm2<-function(adm2.shp, wp,admin_pop_dat){
  
  # make sure polygons have same crs as population raster
  adm2.shp <- spTransform(adm2.shp, crs(wp))
  
  # list of admin2 regions
  adm2.names <- adm2.shp$NAME_2
  
  # admin 2 level population
  wp.adm2.list <- lapply(1:nrow(adm2.shp), function(x) {
    list(state_id = x, state_raster = mask(crop(wp,adm2.shp[x,]),
                                           adm2.shp[x,]))
  })
  
  # store total population at admin 2
  pop.adm2<-vector()
  for ( j in 1:nrow(adm2.shp)){
    pop_j<-wp.adm2.list[[j]]
    pop.adm2[j]<-sum(values(pop_j$state_raster),na.rm=TRUE)
    
  }
  
  # add admin2 population 
  admin_pop_dat$admin2_pop<-pop.adm2
  # not using matching because repeated admin2 names
  # assume order in admin.link is the same as polygon
  
  return (admin_pop_dat)
}
##############################################################################
#########   create population file folders
##############################################################################

setwd(paste0(data_dir,country))

dir.create(file.path('.', 'worldpop'))
dir.create(file.path('worldpop', 'pop_100m'))
dir.create(file.path('worldpop', 'pop_1km'))
dir.create(file.path('worldpop', 'pop_frac'))


setwd(paste0(data_dir,country,'/worldpop'))

options(timeout = 3000) # adjust this time, should be longer than each download

file <- paste0(pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')
if(!file.exists(file)){
  
  url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                frame_year, "/", toupper(pop.abbrev),"/",      
                pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')
  
  download.file(url, file, method = "libcurl",mode="wb")
}

# UNadjusted population counts
worldpop <- raster(paste0(country.abbrev,
                          '_ppp_',frame_year,
                          '_1km_Aggregated_UNadj.tif',sep=''))


##############################################################################
#########   download age group specific population and aggregate
##############################################################################


setwd(paste0(data_dir,country,'/worldpop'))

pop.year <- beg.year:min(end.year,2020)
pop.abbrev <- tolower(gadm.abbrev)

options(timeout = 30000) # adjust this time, should be longer than each download
rigorousFileTest = F # set to TRUE after files have been downloaded to test 
# if files were downloaded correctly, i.e. if they can be loaded into R
for(year in pop.year){
  
  print(year)
  # includes ages 0-1 years and 1-5 years
  
  for(age in c(3:9)*5){
    
    ### down load 100m files
    setwd(paste0(data_dir,country,'/worldpop/pop_100m'))
    
    sex<- 'f'
    
      file <- paste0(pop.abbrev,'_', sex, '_', age, '_', year,'_100m.tif')
      
      # check if the raster file exists. If rigorousFileTest == TRUE, also check 
      # if the file can be successfully loaded
      goodFile = file.exists(file)
      if(goodFile && rigorousFileTest) {
        goodFile = tryCatch(
          {
            test = raster(file)
            TRUE
          }, 
          error = function(e) {FALSE}
        )
      }
      
      if(!goodFile){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/", 
                      year, "/", toupper(pop.abbrev), "/", pop.abbrev, "_", 
                      sex, "_", age, "_", year, ".tif")
        download.file(url, file, method = "libcurl",mode="wb")
      }
      
      
      
      pop_f_agegrp<-raster(file)
      
      proj4string(pop_f_agegrp) <- sf::st_crs(poly.adm1)
      
      pop_grid<-as.data.frame(coordinates(worldpop))
      colnames(pop_grid)<-c('x','y')
      
      
      
      pop_f_aggregate <- aggregate(pop_f_agegrp, fact=10,sum)
      
      pop_grid$f_pop <- raster::extract(pop_f_aggregate, 
                                         pop_grid[c('x','y')])
      
      f_agegrp_pop<-worldpop
      values(f_agegrp_pop)<-pop_grid$f_pop 
      
      setwd(paste0(data_dir,country,'/worldpop/pop_1km'))
      
      writeRaster(f_agegrp_pop, overwrite=TRUE,
                  paste0(pop.abbrev,'_', age, '_',year,'_1km.tif'))
  }
}


### currently cannot download pop files after 2020 directly from worldpop, use 2020 surface for 2021, 2022

if(end.year > 2020){
  for(tmp.year in c(2021:end.year)){
    
    setwd(paste0(data_dir,country,'/worldpop/pop_1km'))
    
    for(age in c(3:9)*5){
    file_2020 <- paste0(pop.abbrev,'_', age, '_',2020,'_1km.tif')
    file_later <- paste0(pop.abbrev,'_', age, '_',tmp.year,'_1km.tif')
    
    file.copy(file_2020, file_later)
    }
    
  }
}



################################################################
#########  get population fractions for total fertility rate
################################################################

pop.year <- c(beg.year:end.year)

# read worldpop rasters, summarize and save admin-2 population
for(year in pop.year){
  
  print(year)
  
  f_15_name = paste0(pop.abbrev,'_', 15, '_',year,'_1km.tif')
  f_20_name = paste0(pop.abbrev,'_', 20, '_',year,'_1km.tif')
  f_25_name = paste0(pop.abbrev,'_', 25, '_',year,'_1km.tif')
  f_30_name = paste0(pop.abbrev,'_', 30, '_',year,'_1km.tif')
  f_35_name = paste0(pop.abbrev,'_', 35, '_',year,'_1km.tif')
  f_40_name = paste0(pop.abbrev,'_', 40, '_',year,'_1km.tif')
  f_45_name = paste0(pop.abbrev,'_', 45, '_',year,'_1km.tif')
  
  setwd(paste0(data_dir,country,'/worldpop/pop_1km'))
  
  pop_f_15<-raster(f_15_name)
  pop_f_20<-raster(f_20_name)
  pop_f_25<-raster(f_25_name)
  pop_f_30<-raster(f_30_name)
  pop_f_35<-raster(f_35_name)
  pop_f_40<-raster(f_40_name)
  pop_f_45<-raster(f_45_name)
  
  
  proj4string(pop_f_15) <- proj4string(pop_f_20) <- 
    proj4string(pop_f_25) <- proj4string(pop_f_30) <- 
    proj4string(pop_f_35) <- proj4string(pop_f_40) <- 
    proj4string(pop_f_45) <- proj4string(poly.adm1)
  
  
  pop_f_15_49<- pop_f_15+pop_f_20+
    pop_f_25+pop_f_30+
    pop_f_35+pop_f_40+
    pop_f_45
  
  # admin 2 population fraction 
  adm2_pop<-pop_adm2(adm2.shp=poly.adm2,
                     wp=pop_f_15_49,
                     admin_pop_dat=adm_link)
  # admin 1 population 
  adm2_pop<-adm2_pop %>% 
    group_by(admin1_idx) %>% 
    mutate(admin1_pop = sum(admin2_pop))
  
  
  # fraction of admin2 w.r.t. admin1
  adm2_pop$admin2_frac<-adm2_pop$admin2_pop/
    adm2_pop$admin1_pop
  
  setwd(data_dir)
  
  save(adm2_pop, file = paste(country,'/','worldpop/pop_frac/', 'admin2_tf_pop_frac_', year, '.rda', sep = ''))
  
  # sanity check, fraction for admin2 in each admin1 sum up to 1
  #print(aggregate(admin2_frac~admin1_idx, data = adm2_pop, FUN = sum))    
  
  
}

# summarize
# load admin-2 population weights (with respect to each admin-1)
adm1.pop.frame <- data.frame()
adm2.pop.frame <- data.frame()

for (i in 1:length(pop.year)){

  year <- pop.year[i]
  setwd(data_dir)
  
  load(paste(country,'/','worldpop/pop_frac/', 'admin2_tf_pop_frac_', year, '.rda', sep = ''))  
  
  ### admin-1 level population 
  
  adm1.pop<-adm2_pop[!duplicated(adm2_pop[,c('admin1_idx')]),]
  # create an ordered admin1 list
  match.order = match( adm1.pop$admin1_idx,paste("admin1", 1: nrow(adm1.pop), 
                                                 sep = "_"))
  adm1.pop = adm1.pop[match.order, ]
  
  adm1.pop<-adm1.pop[,c('admin1_name','admin1_idx','admin1_pop')]
  adm1.pop$year <- year
  adm1.pop.frame <- rbind (adm1.pop.frame,adm1.pop)
  
  adm2.pop <- adm2_pop[,c('admin2_name','admin2_idx','admin2_pop')]
  adm2.pop$year <- year
  adm2.pop.frame <- rbind (adm2.pop.frame,adm2.pop)
  
}

setwd(data_dir)

save(adm1.pop.frame,file= paste(country,'/','worldpop/pop_frac/', 'admin1_tf_pop_frac', '.rda', sep = ''))  
save(adm2.pop.frame,file= paste(country,'/','worldpop/pop_frac/', 'admin2_tf_pop_frac', '.rda', sep = ''))  

################################################################
#########  get population fractions for age specific fertility rate
################################################################

pop.year <- c(beg.year:end.year)

# read worldpop rasters, summarize and save admin-2 population
for(age in c(3:9)*5){
  
  for(year in pop.year){
    
    print(year)
    
    f_name = paste0(pop.abbrev,'_', age, '_',year,'_1km.tif')
    
    setwd(paste0(data_dir,country,'/worldpop/pop_1km'))
    
    pop_f<-raster(f_name)
    
    
    proj4string(pop_f) <- sf::st_crs(poly.adm1)
    
    
    # admin 2 population fraction 
    adm2_pop<-pop_adm2(adm2.shp=poly.adm2,
                       wp=pop_f,
                       admin_pop_dat=adm_link)
    # admin 1 population 
    adm2_pop<-adm2_pop %>% 
      group_by(admin1_idx) %>% 
      mutate(admin1_pop = sum(admin2_pop))
    
    
    # fraction of admin2 w.r.t. admin1
    adm2_pop$admin2_frac<-adm2_pop$admin2_pop/
      adm2_pop$admin1_pop
    
    setwd(data_dir)
    
    save(adm2_pop, file = paste(country,'/','worldpop/pop_frac/', 'admin2_f_',age,'_pop_frac_', year, '.rda', sep = ''))
    
    # sanity check, fraction for admin2 in each admin1 sum up to 1
    #print(aggregate(admin2_frac~admin1_idx, data = adm2_pop, FUN = sum))    
    
    
  }
}



# summarize
# load admin-2 population weights (with respect to each admin-1)

adm1.pop.frame <- data.frame()
adm2.pop.frame <- data.frame()

for(age in c(3:9)*5){
  
#adm1.agegrp.pop.frame <- data.frame()

for (i in 1:length(pop.year)){
  
  year <- pop.year[i]
  setwd(data_dir)
  
  load(paste(country,'/','worldpop/pop_frac/', 'admin2_f_',age,'_pop_frac_', year, '.rda', sep = ''))  
  
  ### admin-1 level population 
  
  adm1.pop<-adm2_pop[!duplicated(adm2_pop[,c('admin1_idx')]),]
  # create an ordered admin1 list
  match.order = match( adm1.pop$admin1_idx,paste("admin1", 1: nrow(adm1.pop), 
                                                 sep = "_"))
  adm1.pop = adm1.pop[match.order, ]
  
  adm1.pop<-adm1.pop[,c('admin1_name','admin1_idx','admin1_pop')]
  adm1.pop$year <- year
  adm1.pop$agegr <- paste0(age,'-',age+4)
  
  adm1.pop.frame <- rbind (adm1.pop.frame,adm1.pop)
  
  
  
  adm2_pop$year <- year
  adm2_pop$agegr <- paste0(age,'-',age+4)
  adm2.pop.frame <- rbind (adm2.pop.frame,adm2_pop)
  
}
}


setwd(data_dir)

save(adm1.pop.frame,file= paste(country,'/','worldpop/pop_frac/', 'admin1_asfr_pop_frac', '.rda', sep = ''))  
save(adm2.pop.frame,file= paste(country,'/','worldpop/pop_frac/', 'admin2_asfr_pop_frac', '.rda', sep = ''))  


################################################################
#########   load libraries
################################################################
#### Libraries ####

rm(list = ls())

options(warn = -1)
library(SUMMER)
library(dplyr)
library(rgdal)
library(INLA)
library(ggplot2)
library(maptools)

options(gsubfn.engine = "R")
library(rgeos)
library(sqldf)
library(raster)


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
#########   load worldpop
################################################################


### automated downloading, if not working, try manually download
setwd(paste0(data_dir,'/',country,'/worldpop/'))
pop.abbrev <- tolower(gadm.abbrev)


# UNadjusted population counts
worldpop <- raster(paste0(country.abbrev,
                          '_ppp_',frame_year,
                          '_1km_Aggregated_UNadj.tif',sep=''))

################################################################
#########  Load cluster data
################################################################

# load DHS cluster information
setwd(paste(data_dir,country,sep=''))
load( file = paste0('prepared_IR_dat/fert_cluster_dat_',survey_year,'.rda'))

# clusters
cluster_list<-fert_pyears[!duplicated(fert_pyears[c('cluster','survey_year',
                                            'LONGNUM','LATNUM')]),]

# keep only the clusters from DHS surveys using the recent census
cluster_list<-cluster_list[cluster_list$survey %in% survey_year,]



################################################################
#########   Correct urban cluster
################################################################

# The function below corrects urban clusters that are misclassified to be rural due to jittering. Jittered location is not the exact location 
# but a randomly shifted location for the purpose of confidentiality. This may results in some urban clusters to be jittered to some rural areas, 
# which is unexpected to our classification algorithm and therefore correcting the misclassified clusters is of interest. 

# The codes below fulfills the process above by assigning the possibly misclassified urban clusters to the nearest most densely populated areas.
# It's generally true that the urban areas tend to have a higher population density so this process can alleviate the side effect
# of jittering.
constr_prior <- function(obs,jitter_r,prob_r,poly_admin,pop_ras){
  
  # for the cluster, find its coordinates and admin1 area
  admin1_index<-obs$admin1
  sp_xy<-SpatialPoints(as.data.frame(obs)[,c('LONGNUM','LATNUM')])
  proj4string(sp_xy) <- CRS(sf::st_crs(poly_admin)$proj4string)

  #pt<-as.data.frame(obs)[,c('LONGNUM','LATNUM')]
  #colnames(pt)<-c('x','y')
  
  # generate gitter 
  #jitter_r<-2000
  cluster_buffer<-buffer(sp_xy, width=jitter_r)
  
  # extract pixels within the buffer
  temp_pop_cir<-mask(crop(pop_ras,cluster_buffer),
                     cluster_buffer)
  
  # put admin area restriction
  admin1_poly<-poly_admin[poly_admin$NAME_1== admin1_index,]
  temp_pop_admin<-mask(crop(temp_pop_cir,admin1_poly),
                       admin1_poly)
  
  
  # check whether need to adjust for constraint 
  cir_val<-values(temp_pop_cir)
  admin_val<-values(temp_pop_admin)
  
  admin_adj<-length(which(!is.na(cir_val)))!=length(which(!is.na(admin_val)))
  
  if(admin_adj){
    #normc<-admin1_normc(pt,jitter_r,admin1_poly,ntrial=1000)
    normc<-1
  }else{  normc<-1}
  
  
  
  ## prepare sample frame
  
  temp_val<-values(temp_pop_admin)
  pop_index<-which(!is.na(temp_val))
  
  
  temp_frame<-as.data.frame(coordinates(temp_pop_admin))
  pixel_candidate<-temp_frame[pop_index,]
  pixel_candidate$pop_den<-temp_val[pop_index]
  
  pixel_candidate$center_x<-obs$LONGNUM
  pixel_candidate$center_y<-obs$LATNUM
  pixel_candidate$dist<-diag(distm(pixel_candidate[,c('x','y')], 
                                   pixel_candidate[,c('center_x','center_y')]))
  
  pixel_candidate$unn_w<-pixel_candidate$pop_den*
    1/(2*pi * 2 * pixel_candidate$dist)*normc*prob_r
  
  pixel_candidate$normc<-normc
  return(pixel_candidate[,c("x","y",'normc','unn_w')])
  #return(pixel_candidate)
  
}

# only correct urban clusters 
# check if the strata are named 'urban' and 'rural' or 'U', 'R'
urban_clus<-cluster_list[cluster_list$urban=='urban',]
rural_clus<-cluster_list[cluster_list$urban=='rural',]

urban_clus$x_adj<-NA
urban_clus$y_adj<-NA
rural_clus$x_adj<-rural_clus$LONGNUM
rural_clus$y_adj<-rural_clus$LATNUM

# points have to stay within the same admin2 region 
for( i in 1:dim(urban_clus)[1]){
  
  print(i)
  temp_frame<-constr_prior(urban_clus[i,],2000,1,poly.adm1,worldpop)
  p_mode = sqldf("SELECT * FROM temp_frame GROUP BY x,y ORDER BY SUM(unn_w) DESC LIMIT 1")
  urban_clus[i,]$x_adj<-p_mode$x
  urban_clus[i,]$y_adj<-p_mode$y
  
}


prep_dat<-rbind(urban_clus,rural_clus)
#xy <- as.matrix(prep_dat[c('x_adj','y_adj')])


# create directory to store cluster data
setwd(paste0(data_dir,'/',country))

if(!dir.exists(paths = paste0('prepared_dat/'))){
  dir.create(path = paste0('prepared_dat/'))
}


save(prep_dat,file='prepared_dat/prep_dat.rda')



################################################################
######### adding covariates 
################################################################

crc_dat<-prep_dat

# set up corrected xy for clusters
xy_crc <- as.matrix(crc_dat[c('x_adj','y_adj')])
crc_dat$x<-crc_dat$x_adj # x_adj and y_adj: corrected coordinates
crc_dat$y<-crc_dat$y_adj

# extract covariates
crc_dat$pop_den<-raster::extract(worldpop,xy_crc)

# only retain part of the columns to reduce redundancy
col_select<-c('clusterIdx','urban','admin1','admin2',
              'admin1.name','admin2.name',
              'admin1.char','admin2.char',
              'survey','pop_den','x','y')

crc_dat_final<-crc_dat[,col_select]

setwd(paste0(data_dir,'/',country))

save(crc_dat_final,file='prepared_dat/crc_dat.rda')




################################################################
######### prepare data without urban correction
################################################################

uncrc_dat<-prep_dat

# set up uncorrected xy for clusters
xy_uncrc <- as.matrix(uncrc_dat[c('LONGNUM','LATNUM')])
uncrc_dat$x<-uncrc_dat$LONGNUM # x_adj and y_adj: corrected coordinates
uncrc_dat$y<-uncrc_dat$LATNUM


# extract covariates
uncrc_dat$pop_den<-raster::extract(worldpop,xy_uncrc)


# keep columns
col_select<-c('clusterIdx','urban','admin1','admin2',
              'admin1.name','admin2.name',
              'admin1.char','admin2.char',
              'survey','pop_den','x','y')

uncrc_dat_final<-uncrc_dat[,col_select]

setwd(paste0(data_dir,'/',country))

save(uncrc_dat_final,file='prepared_dat/uncrc_dat.rda')



################################################################
#########   prepare national grid
################################################################
## set up grid
urb_dat<-as.data.frame(coordinates(worldpop))
colnames(urb_dat)<-c('x','y')

urb_dat$pop_den<-raster::extract(worldpop,urb_dat[c('x','y')])

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

save(urb_dat,file="prepared_dat/natl_grid.rda")
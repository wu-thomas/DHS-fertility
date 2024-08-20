################################################################
#########   load libraries
################################################################
rm(list = ls())
options(warn = -1)
#### Libraries ####
library(SUMMER)
#library(classInt)
#library(RColorBrewer)
library(dplyr)
#library(tidyr)
#library(rgdal)
library(scales)
library(INLA)
#library(survey)
library(ggplot2)
library(raster)
library(maptools)
library(gridExtra)
#library(mgcv)
library(caret)
#library(geosphere)
#library(rgeos)
library(haven)
library(labelled)
library(data.table)
options(gsubfn.engine = "R")
library(sqldf)

library(gstat)
library(stringdist)
library(openxlsx)

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


admin2.names <- admin2_info$data
admin2.names <- admin2.names[,c('admin2.name','admin2.char')]
admin2.names$GADM <- admin2.names$admin2.name
admin2.names$Internal <- admin2.names$admin2.char


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


worldpop.df <-
  as.data.frame(worldpop, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit() 

worldpop.df$pop <- worldpop.df$rwa_ppp_2012_1km_Aggregated_UNadj
worldpop.df <- worldpop.df[,c('x','y','pop')]

g.pop.base <- ggplot(data = worldpop.df) +geom_raster(aes(x = x, y = y, fill = pop))+
  scale_fill_gradientn(colours = rev(terrain.colors(100)),limits=c(0,20000),
                       values=rescale(seq(0,8,1)),
                       space = "Lab",na.value = "grey")+
  theme_bw()+xlab('lon')+ylab('lat')+
  guides(fill = guide_colourbar(barheight=unit(6, "cm"),
    title.position = "top",
                                title.hjust = 0.5,title.vjust=3.5,
                                title='Pop',
                                label.position = "left"))
  
overlay.adm1 <- fortify(poly.adm1, region="NAME_1")
overlay.adm2 <- fortify(poly.adm2, region="NAME_2")

g.pop.base <- g.pop.base + geom_polygon(data=overlay.adm2, aes(x=long, y=lat, group=group),linewidth=0.3, color="darkgrey", alpha=0) 
g.pop.base <- g.pop.base + geom_polygon(data=overlay.adm1, aes(x=long, y=lat, group=group),linewidth=0.5, color="black", alpha=0) 

setwd(paste0(data_dir,'/',country,'/worldpop/'))

ggsave(g.pop.base, width=10, height = 8, file = paste0("Rwanda_2012_population.pdf"))

################################################################
#########   sampling frame urban proportion at admin1 
################################################################


setwd(paste(data_dir,country,sep=''))

# read the .txt or .xlsx file containing urban population fraction at admin1 level.
if (file.exists(paste(country.abbrev, "frame_urb_prop.txt", sep = "_"))){
  frame <- read.delim(paste(country.abbrev, "frame_urb_prop.txt", sep = "_"), 
                      header = FALSE,sep=' ')
}
if (file.exists(paste(country.abbrev, "frame_urb_prop.xlsx", sep = "_"))){
  frame <- read.xlsx(paste(country.abbrev, "frame_urb_prop.xlsx", sep = "_"))
}


# # identify column for fraction (need additional processing in general)
frame[,c(2,4)] <- lapply(frame[,c(2,4)],   ## function to remove comma in numbers
                         function(x){as.numeric(gsub(",", "", x))})
frame$frac <- frame[, 2]/frame[, 4]

# greedy algorithm to match admin names 
adm.ref <- expand.grid(tolower(frame[, 1]),        
                        tolower(admin2.names$GADM)) # Distance matrix in long form
names(adm.ref) <- c("frame_name","gadm_name")
### string distance,  jw=jaro winkler distance, try 'dl' if not working
adm.ref$dist <- stringdist(adm.ref$frame_name,
                            adm.ref$gadm_name, method="jw") 

greedyAssign <- function(a,b,d){
  x <- numeric(length(a)) # assgn variable: 0 for unassigned but assignable, 
  # 1 for already assigned, -1 for unassigned and unassignable
  while(any(x==0)){
    min_d <- min(d[x==0]) # identify closest pair, arbitrarily selecting 1st if multiple pairs
    a_sel <- a[d==min_d & x==0][1] 
    b_sel <- b[d==min_d & a == a_sel & x==0][1] 
    x[a==a_sel & b == b_sel] <- 1
    x[x==0 & (a==a_sel|b==b_sel)] <- -1
  }
  cbind(a=a[x==1],b=b[x==1],d=d[x==1])
}

match_order<-data.frame(greedyAssign(adm.ref$frame_name,
                                     adm.ref$gadm_name,
                                     adm.ref$dist))

# create reference table 
ref.tab <- admin2.names
ref.tab$matched_name <- frame$V1[match_order$a] ### check!!!
ref.tab$urb_frac <- frame$frac[match_order$a] 




################################################################
#########   admin 1 threshold 
################################################################

## load grid
setwd(paste(data_dir,country,sep=''))
load(file='prepared_dat/natl_grid.rda')

# index the grid
urb_dat$index <- c(1:nrow(urb_dat))
adm2_dat <- split( urb_dat , f = urb_dat$admin2 )


# This function computes the urban population threshold for a given admin1 area.
# This is done by keep counting the urban locations until the urban population fraction in the reference table is reached.
thresh_urb<-function(adm_grid,ref_tab){
  
  # sort grid population
  vals <- adm_grid$pop_den
  vals[is.na(vals)] <- 0
  sort.obj <- sort.int(vals, decreasing = TRUE, index.return = TRUE, method = 'shell')
  svals <- sort.obj$x
  svals.int <- sort.obj$ix
  
  # extract cutoff proportion based on admin1
  adm.idx <- adm_grid$admin2.char[1]
  cutoff <- ref_tab[ref_tab$Internal==adm.idx,]$urb_frac
  
  # determine population threshold and urban rural
  csvals <- cumsum(svals)/sum(svals)
  is.urb <- csvals <= cutoff
  org.isurb <- is.urb[invPerm(svals.int)]
  threshold <- min(vals[org.isurb == 1]) #cutoff
  
  # prepare return object (grid with urban/rural)
  adm_grid$threshold <- threshold
  adm_grid$urban <- as.numeric(org.isurb)
  #adm_grid[is.na(adm_grid$pop_den),]$urban<-NA
  
  return(adm_grid)
  
}

urb_list<-lapply(adm2_dat, FUN=thresh_urb,ref_tab=ref.tab)

urb_class <- do.call("rbind", urb_list)


urb_grid <- urb_dat
urb_grid$urb_ind <-NA
urb_grid[urb_class$index,]$urb_ind <- urb_class$urban


urb_surf<-worldpop
values(urb_surf)<-urb_grid$urb_ind

setwd(paste(data_dir,country,sep=''))

if(!dir.exists(paths = paste0('UR/'))){
  dir.create(path = paste0('UR/'))
}

setwd(paste(data_dir,country,'/UR',sep=''))

writeRaster(urb_surf, overwrite=TRUE,
            paste0(country,'_UR_ind_1km.tif'))

## save reference table along with calculated threshold 
thresh_ref <- urb_class[!duplicated(urb_class[,c('admin2')]),]
ref.tab$threshold <- thresh_ref$threshold # check whether the thresholds are sensible (shouldn't be NA or all 0)

setwd(paste(data_dir,country,sep=''))

write.xlsx(ref.tab, file='prepared_dat/reference_table.xlsx',
           row.names = FALSE)



################################################################
#########   U/R plot
################################################################

setwd(paste(data_dir,country,'/UR',sep=''))

urb_surf <- raster(paste0(country,'_UR_ind_1km.tif'))

urb.ras.df <-
  as.data.frame(urb_surf, xy = TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit() 

urb.ras.df$urban <- urb.ras.df$Rwanda_UR_ind_1km
urb.ras.df <- urb.ras.df[,c('x','y','urban')]
urb.ras.df$urban <- as.factor(as.character(urb.ras.df$urban))

g.urb.base <- ggplot(data = urb.ras.df) +geom_raster(aes(x = x, y = y, fill = urban))+
  scale_fill_manual(values = rev(terrain.colors(2)))+
  theme_bw()+xlab('lon')+ylab('lat')

overlay.adm1 <- fortify(poly.adm1, region="NAME_1")
overlay.adm2 <- fortify(poly.adm2, region="NAME_2")

g.urb.base <- g.urb.base + geom_polygon(data=overlay.adm2, aes(x=long, y=lat, group=group),linewidth=0.3, color="darkgrey", alpha=0) 
g.urb.base <- g.urb.base + geom_polygon(data=overlay.adm1, aes(x=long, y=lat, group=group),linewidth=0.5, color="black", alpha=0) 

setwd(paste0(data_dir,'/',country,'/UR/'))

ggsave(g.urb.base, width=10, height = 8, file = paste0("Rwanda_UR_classification_map.pdf"))


################################################################
#########   check classification accuracy based on clusters
################################################################

setwd(paste(data_dir,country,sep=''))

load('prepared_dat/crc_dat.rda')
load('prepared_dat/uncrc_dat.rda')

### remove rows with missing covariates, could also build model with missing data
crc_dat<-crc_dat_final[complete.cases(crc_dat_final), ]
uncrc_dat<-uncrc_dat_final[complete.cases(uncrc_dat_final), ]


xy_crc <- as.matrix(crc_dat[c('x','y')])
xy_uncrc <- as.matrix(uncrc_dat[c('x','y')])

# extract the urban/rural prediction
crc_dat$urb_pred<-raster::extract(urb_surf,xy_crc)
uncrc_dat$urb_pred<-raster::extract(urb_surf,xy_uncrc)
pred_crc <- factor( ifelse(crc_dat$urb_pred ==1 ,"U","R" )) # make sure U/R or urban/rural
pred_crc  <- relevel(pred_crc, "U") # make sure levels are same 
pred_uncrc <- factor( ifelse(uncrc_dat$urb_pred ==1 ,"U","R" ))
pred_uncrc  <- relevel(pred_uncrc, "U") # make sure levels are same 


# compute the confusion to evaluate the accuracy
confmatrix_crc<-caret::confusionMatrix(
  data = pred_crc,
  reference = as.factor(crc_dat$urban)
)

confmatrix_crc

setwd(paste(data_dir,country,'/UR',sep=''))

save(confmatrix_crc,file='confmatrix_crc.rda')

confmatrix_uncrc<-caret::confusionMatrix(
  data = pred_uncrc,
  reference = as.factor(uncrc_dat$urban)
)

confmatrix_uncrc

setwd(paste(data_dir,country,'/UR',sep=''))

save(confmatrix_uncrc,file='confmatrix_uncrc.rda')


##############################################################################
#########   load packages
##############################################################################

rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.


# Load libraries and info ----------------------------------------------------------
library(dplyr)
library(labelled)
library(survival)
library(haven)
library(rdhs)
library(surveyPrev)

################################################################
#########   set parameters
################################################################
country <- 'Tanzania'

## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-4)], collapse = "/")

data.dir<-paste0(home.dir,'/Data/')

setwd(paste(data.dir,country,sep=''))


info.name <- paste0(country, "_general_info.Rdata")

load(file = paste0(info.name, sep=''))

##############################################################################
#########   load helper functions 
##############################################################################

setwd(paste(code.dir))
source('DataProcessing_helper.R')


##############################################################################
#########   create folders for prepared data
##############################################################################

setwd(paste(data.dir))

dir.create(file.path('.', 'prepared_IR_dat'))

##############################################################################
###### RDHS Configration 
##############################################################################


## login
# set API to get DHS data -- you will need to change this to your information!
rdhs::set_rdhs_config(email = "xxxxx",
                project = "xxxxxxxxxxxxxxxxxxxxxxxxxx")

rdhs::update_rdhs_config(email = "xxxxxxx", password = T,
                   project = "xxxxxxxxxxxxxxxxxxxxxxxxxxxx")

##############################################################################
#########   load survey meta data
##############################################################################


setwd(paste0(data.dir))

### run for the first time use saved object later
if(!file.exists('DHS_meta.rda')){
  DHS.country.meta <- rdhs::dhs_countries()
  DHS.survey.meta <- rdhs::dhs_surveys()
  DHS.dataset.meta <- rdhs::dhs_datasets()
  
  save(DHS.country.meta,
       DHS.survey.meta,
       DHS.dataset.meta,
       file='DHS_meta.rda')
}else{
  
  ### use saved object if not first run
  setwd(paste0(data.dir))
  load('DHS_meta.rda')
}


##############################################################################
#########   load polygon files
##############################################################################

setwd(paste0(data.dir,'/shapeFiles_gadm'))


####  !!!! 
### need to customize based on country specific information
### i.e. what GADM levels corresponds to admin-1 and admin-2 of our interest
###
  
  
### only run for the first time and use stored results later
if(!file.exists('country_shp_analysis.rds')|
   !file.exists('country_shp_smoothed.rds')|
   !file.exists('admin1_info.rds')|
   !file.exists('admin2_info.rds')){
  
  ### need to customize based on country specific information
  ### i.e. what GADM levels corresponds to admin-1 and admin-2 of our interest
  
  country_shp_analysis <- get_country_GADM(country=country,resolution=1)
  country_shp_smoothed <- get_country_GADM(country=country,resolution=2)
  country_shp_analysis <- lapply(country_shp_analysis, function(x) {
    sf::st_set_crs(x, 4326)
  })
  country_shp_smoothed <- lapply(country_shp_smoothed, function(x) {
    sf::st_set_crs(x, 4326)
  })
  
  
  admin1_info <- surveyPrev::adminInfo(poly.adm = country_shp_analysis[['Admin-1']],
                                       admin = 1,by.adm="NAME_1") ### DHS admin-1 is GADM admin-2!!
  admin1_info$data$admin1.char <- paste0("admin1_", 1:dim(admin1_info$data)[1])
  
  admin2_info <- surveyPrev::adminInfo(poly.adm = country_shp_analysis[['Admin-2']],
                                       by.adm.upper='NAME_1',
                                       admin = 2,by.adm="NAME_2") ### DHS admin-2 is GADM admin-3!!
  admin2_info$data$admin2.char <- paste0("admin2_", 1:dim(admin2_info$data)[1])
  
  saveRDS(country_shp_analysis,'country_shp_analysis.rds')
  saveRDS(country_shp_smoothed,'country_shp_smoothed.rds')
  saveRDS(admin1_info,'admin1_info.rds')
  saveRDS(admin2_info,'admin2_info.rds')

  
  #save(country_shp_analysis,
  #     country_shp_smoothed,
  #     admin1_info,admin2_info,
  #     file='shape_info.rda')
}else{
  
  ### use stored results for later runs 
  setwd(paste0(data.dir,'/shapeFiles_gadm'))
  
  country_shp_analysis <- readRDS('country_shp_analysis.rds')
  country_shp_smoothed <- readRDS('country_shp_smoothed.rds')
  admin1_info <- readRDS('admin1_info.rds')
  admin2_info <- readRDS('admin2_info.rds')
}


##############################################################################
#########   process DHS surveys
##############################################################################


### only run for the first time and use stored results later

setwd(paste0(data.dir,'/DHS_data'))

if( !file.exists(paste0(country,'_',survey_year,'_GPS.rds')) | !file.exists(paste0(country,'_',survey_year,'_IR.rds')) ){
  # Find DHS surveys ----------------------------------------------------------
  
  #get country ID
  countryId <- DHS.country.meta[DHS.country.meta$ISO3_CountryCode==toupper(gadm.abbrev),]
  
  
  potential_surveys <-  DHS.dataset.meta %>% dplyr::filter(CountryName==country &
                                                             SurveyYear == survey_year-survey_year_span &
                                                             ((FileType == 'Individual Recode' &
                                                                 FileFormat=='Stata dataset (.dta)') |
                                                                (FileType == 'Geographic Data')))
  
  #only keep surveys with both IR recode and a geographic dataset
  dhs_survey_ids <- as.numeric(unique(potential_surveys$SurveyNum)[sapply(unique(potential_surveys$SurveyNum),
                                                                          function(num){
                                                                            if(sum(c("Individual Recode","Geographic Data") %in% (potential_surveys %>% filter(SurveyNum==num))$FileType) ==2){return(T)
                                                                            }else(return(F))})])
  
  surveys <- potential_surveys %>% filter(SurveyNum %in% dhs_survey_ids) %>% 
    group_by(SurveyYear) %>% arrange(SurveyYear,DatasetType)
  
  
  # CHECK THAT SURVEYS FOR CORRECT COUNTRY HAVE BEEN CHOSEN
  #unique(surveys$CountryName)
  
  data.paths.tmp <- rdhs::get_datasets(surveys[surveys$SurveyYear==survey_year-survey_year_span,]$FileName, clear_cache = T)
  
  
  raw.geo.dat <- readRDS(paste0(data.paths.tmp[1]))
  raw.IR.dat <- readRDS(paste0(data.paths.tmp[2]))
  
  
  setwd(paste0(data.dir,'/DHS_data'))
  
  saveRDS(raw.geo.dat,file=paste0(country,'_',survey_year,'_GPS.rds'))
  saveRDS(raw.IR.dat,file=paste0(country,'_',survey_year,'_IR.rds'))
}else{
  
  ### use stored results for later run
  setwd(paste0(data.dir,'/DHS_data'))
  
  raw.geo.dat <- readRDS(file=paste0(country,'_',survey_year,'_GPS.rds'))
  raw.IR.dat <- readRDS(file=paste0(country,'_',survey_year,'_IR.rds'))
}


##############################################################################
#########   process fertility data
##############################################################################


# Process data for each DHS survey year ----------------------------------------------------------

# The codes below first loads the raw DHS data, then it assigns the GPS coordinates 
# to each sampling cluster and admin regions where
# the sampling is conducted and assigns the admin regions where the clusters are located.

### prepare person-year data
fert_cluster_dat <- prepare_fert_cluster(data = raw.IR.dat,
                                 country = country,
                                 survey_year = survey_year,
                                 by.individual=F)

### prepare cluster info
cluster.info <- surveyPrev::clusterInfo(geo = raw.geo.dat,
                           poly.adm1 = country_shp_analysis[['Admin-1']],
                           by.adm1="NAME_1", ### very important to check whether levels in GAMD is consistent with DHS admin levels
                           poly.adm2 = country_shp_analysis[['Admin-2']],
                           by.adm2="NAME_2")

### merge cluster info
fert_cluster_dat <- left_join(fert_cluster_dat,cluster.info$data,by=c("cluster"))
fert_cluster_dat<- fert_cluster_dat[!(is.na(fert_cluster_dat$LONGNUM)), ]

### Save processed fertility data  ----------------------------------------------------------

setwd(paste0(data.dir))
fert_pyears = fert_cluster_dat
save(fert_pyears, file = paste0('prepared_IR_dat/fert_cluster_dat_',survey_year,'.rda'))


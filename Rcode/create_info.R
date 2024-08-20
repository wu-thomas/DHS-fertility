##################################################################
##################################################################
# This script is used to generate the info file for a given country
##################################################################
##################################################################
rm(list = ls())

################################################################
#########   set parameters
################################################################

# Files info (For those lines with ### xxx ### above, please fill in as commented)
country <- 'Tanzania'

### please fill in the country abbreviation in all upper case of gadm files ### (e.g. fill in SEN for gadm36_SEN_3.shp)
gadm.abbrev <- "TZA"
DHS.abbrev <- 'TZ'

### please fill in the following information ####

survey_year <- 2022
frame_year <- 2012 ### year of the sampling frame creation
survey_year_span <- 0 ### for single year survey, this should be 0, for surveys like Rwanda 2019-2020, this is 1


### setting rest of parameters using info from above
country.abbrev <- tolower(gadm.abbrev)           # lower the country gadm abbreviation 
beg.year <- survey_year-9 # the first year of the interest
end.year <- survey_year # last year we would like to project to 




################################################################
#########   set parameters
################################################################

## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

data.dir<-paste0(home.dir,'/Data/',country)
res.dir<-paste0(home.dir,'/Results/',country)
#pop_dir<-paste0(main_dir,'/Population/')

code.dir <- paste0(home.dir,'/Rcode/',country)



setwd(paste(data.dir))

info.name <- paste0(country, "_general_info.Rdata")

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

save.image(file = paste0(info.name, sep=''))


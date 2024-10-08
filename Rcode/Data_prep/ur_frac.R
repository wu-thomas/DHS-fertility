################################################################
#########   load libraries
################################################################
rm(list = ls())

#### Libraries ####
library(SUMMER)
library(classInt)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(rgdal)
library(scales)
library(INLA)
library(survey)
library(ggplot2)
library(raster)
library(maptools)
library(gridExtra)
library(mgcv)
library(caret)
library(geosphere)
library(rgeos)
library(haven)
library(labelled)
library(data.table)
library(sqldf)
library(sp)
library(gstat)
library(stringdist)
library(openxlsx)

# extract file location of this script
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home_dir <- paste(code.path.splitted[1: (length(code.path.splitted)-4)], collapse = "/")
countries <- countries <- scan(paste0(home_dir, "/countries_implemented.txt"), character(), quote = "")
country <- countries[length(countries)] # retrieve the country being analyzed
info.name <- paste0(country, "_general_info.Rdata")

load(file = paste0(info.name, sep=''))


################################################################
#########   set directories
################################################################

data_dir <- paste0(home_dir,'/Data/', country) # set the directory to store the data
res_dir <- paste0(home_dir,'/Results/', country) # set the directory to store the results (e.g. fitted R objects, figures, tables in .csv etc.)
pop_dir <- paste0(data_dir,'/Population') # set the directory to store the population surface files

################################################################
#########   load files
################################################################

setwd(paste0(data.dir,'/shapeFiles_gadm'))

country_shp_analysis <- readRDS('country_shp_analysis.rds')
country_shp_smoothed <- readRDS('country_shp_smoothed.rds')
admin1_info <- readRDS('admin1_info.rds')
admin2_info <- readRDS('admin2_info.rds')

poly.adm0 = country_shp_analysis[['National']]
poly.adm1 = country_shp_analysis[['Admin-1']]
poly.adm2 = country_shp_analysis[['Admin-2']]

admin1.names = admin1_info$data
admin1.names$GADM=admin1.names$admin1.name

################################################################
#########   generate U/R population table
################################################################
table_generator = function(){
  pop_frac = as.data.frame(matrix(NA, nrow = length(admin1.names$GADM), ncol = 4))
  colnames(pop_frac) = c("admin1", "urban", "rural", "total")
  pop_frac$admin1 = admin1.names$GADM
  table_choice = menu(c("R","Excel"), 
                      title = "You would like to create the U/R summary table via? (Choose option 1 or 2 from below)")
  if (table_choice == 1){
    choice = menu(c("Exact urban and rural population of each admin-1", 
                 "Average household size and number of households of each admin-1"), 
               title = "What information do you know? (Choose option 1 or 2 from below)")
                  
    if (choice == 1){
      for (i in 1: nrow(pop_frac)) {
        urban_pop = as.numeric(readline(prompt = paste0("Please enter the urban population of ", pop_frac$admin1[i], " (no comma) :")))
        total_pop = as.numeric(readline(prompt = paste0("Please enter the total population of ", pop_frac$admin1[i], " (no comma) :")))
        pop_frac[i, -1] = c(urban_pop, NA, total_pop, NA)
      }
      while (sum(is.na(pop_frac[, c(2, 4)])) > 0) {
        missing_admin = which(is.na(pop_frac[, c(2, 4)]), arr.ind = T)
        for (i in 1: nrow(missing_admin)) {
          val = as.numeric(readline(prompt = paste0("The ", c("urban", "total")[missing_admin[i, 2]], " population of ", pop_frac$admin1[missing_admin[i, 1]], " is missing, please enter (no comma) :")))
          pop_frac[missing_admin[i, 1], c("urban", "total")[missing_admin[i, 2]]] = val
        }
      }
      pop_frac$rural = pop_frac$total - pop_frac$urban
    } 
    
    if(choice == 2){
      household_table = as.data.frame(matrix(NA, nrow = length(admin1.names$GADM), ncol = 5))
      colnames(household_table) = c("admin1", "urbanHH", "urbanHH_avg", "totalHH", "totalHH_avg")
      household_table$admin1 = admin1.names$GADM
      for (i in 1: nrow(pop_frac)) {
        urban_household = as.numeric(readline(prompt = paste0("Please enter the number of urban households of ", pop_frac$admin1[i], " (no comma) :")))
        urban_household_size = as.numeric(readline(prompt = paste0("Please enter the averaged urban household size of ", pop_frac$admin1[i], " (no comma) :")))
        total_household = as.numeric(readline(prompt = paste0("Please enter the number of total households of ", pop_frac$admin1[i], " (no comma) :")))
        total_household_size = as.numeric(readline(prompt = paste0("Please enter the averaged total household size of ", pop_frac$admin1[i], " (no comma) :")))
        
        household_table[i, -1] = c(urban_household, urban_household_size, total_household, total_household_size)
      }
      temp = c("averaged urban household size", "number of urban households", 
               "averaged total household size", "number of total households")
      while (sum(is.na(household_table)) > 0) {
        missing_admin = which(is.na(household_table), arr.ind = T)
        for (i in 1: nrow(missing_admin)) {
          val = as.numeric(readline(prompt = paste0("The ", temp[missing_admin[i, 2]], " of ", pop_frac$admin1[missing_admin[i, 1]], " is missing, please enter (no comma) :")))
          household_table[missing_admin[i, 1], c("urban", "total")[missing_admin[i, 2]]] = val
        }
      }
      pop_frac$urban = household_table$urbanHH * household_table$urbanHH_avg
      pop_frac$total = household_table$totalHH * household_table$totalHH_avg
      pop_frac$rural = pop_frac$total - pop_frac$urban
    }
    print("This is the U/R population summary table you've created")
    print(pop_frac)
    write.xlsx(pop_frac, file = paste(country.abbrev, "frame_urb_prop.xlsx", sep = "_"))
    print(paste0("U/R population summary table saved in ", data_dir))
  }
  if (table_choice == 2){
    write.xlsx(pop_frac, file = paste(country.abbrev, "frame_urb_prop.xlsx", sep = "_"))
    print(paste0("The table has been saved into ", data_dir, ". Please fill in the table before we proceed"))
  }
}

table_generator()






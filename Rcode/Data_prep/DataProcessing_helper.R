###############################################################
###  load GADM files
###############################################################

get_country_GADM <- function(country,resolution=1) {
  
  country_iso3 <- DHS.country.meta[DHS.country.meta$CountryName==country,'ISO3_CountryCode']
  
  gadm_list <- list()
  levels <- 0
  repeat {
    tmp.gadm <- geodata::gadm(country = country_iso3, resolution=resolution,
                              level = levels,
                              path = tempdir())
    if (is.null(tmp.gadm)) {
      break
    } else {
      tmp.gadm <- sf::st_as_sf(tmp.gadm)
      tmp.gadm <- sf::st_set_crs(tmp.gadm, 4326)
      
      n_region <- dim(tmp.gadm)[1]
      #message(paste0('n region: ',n_region))
      if(n_region >1000){break}
      
      
      if(levels==0){      gadm_list[['National']]  <- tmp.gadm
      }else{
        gadm_list[[paste0('Admin-',levels)]]  <- tmp.gadm}
      levels <- levels + 1
    }
  }
  
  
  return(gadm_list)
}


###############################################################
###  prepare fertility
###############################################################

#' Calculate Person-Years of Exposure for Fertility
#'
#' This function calculates the person-years of exposure (women-years contributed)
#' (adapted from demogsurv package)
#'
#' @param data.IR DHS IR recode raw data.
#' @param survey_year The year of the survey.
#' @param country The country of the survey. Special adjustment for Ethiopia.
#' @param begin_year Use data starting from which year, default is 10 year before survey
#' @param end_year Use data till which year, default is survey year
#' @param by.individual whether prepare data for individual level (for individual level covariate) or aggregate to cluster level


prepare_fert_cluster <- function(data.IR, survey_year, country, begin_year =NULL, end_year = NULL,
                                 tips= NULL,by.individual= T) {
  
  # Country-specific CMC adjustment
  cmc.adj <- ifelse(country == 'Ethiopia', 92, 0)
  
  if(is.null(begin_year)){
    begin_year <- as.numeric(survey_year) - 11
  }
  if(is.null(end_year)){
    end_year <- as.numeric(survey_year) +1
  }
  
  
  # Set up parameters for calculating person-years
  agegr = 3:10*5
  period = c(begin_year:end_year) ### !!important, the grid should include one more year after survey
  tips=tips
  clusters=~v001
  strata=~v024+v025
  id="caseid"
  dob="v011"
  intv = "v008"
  weight= "v005"
  varmethod = "lin"
  bvars = grep("^b3\\_[0-9]*", names(data.IR), value=TRUE)
  birth_displace = 1e-6
  origin=1900
  scale=12
  counts=FALSE
  clustcounts = FALSE
  
  # Prepare person-years
  data.IR$id <- data.IR[[id]]
  data.IR$dob <- data.IR[[dob]]
  data.IR$intv <- data.IR[[intv]]
  data.IR$weights <- data.IR[[weight]] / mean(data.IR[[weight]])
  
  
  # Define grouping variables and create model frame
  
  ### individual level data (for maternal education covariate)
  if(!by.individual){
    
    by <- ~1
    vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
    f <- as.formula(paste("~", paste(vars, collapse = "+")))
    
  }else{
    
    by= ~caseid
    f <- formula('~ caseid')
    
  }
  
  mf <- model.frame(formula = f, data = data.IR, na.action = na.pass,
                    id = id, weights = weights, dob = dob, intv = intv)
  
  # Ensure birth variables are provided
  if (!is.character(bvars) || length(bvars) == 0) {
    stop("`bvars` must be specified as a variable or list of variables containing child DOB. DHS default values b3_01, b3_02, ... were not found.")
  }
  
  # Prepare birth data
  births <- model.frame(paste("~", paste(bvars, collapse="+")),
                        data.IR, na.action=na.pass, id=id)
  births <- reshape(births,
                    idvar="(id)", timevar="bidx",
                    varying=bvars, v.names="bcmc", direction="long")
  
  births <- births[!is.na(births$bcmc), ]
  births$bcmc <- births$bcmc + births$bidx * birth_displace
  
  # Rename for tmerge
  names(mf)[names(mf) == "(id)"] <- "id_"
  names(births)[names(births) == "(id)"] <- "id_"
  
  # Merge event data
  epis <- tmerge(mf, mf, id=id_, tstart=`(dob)`, tstop=`(intv)`)
  epis <- tmerge(epis, births, id=id_, birth = event(bcmc))
  
  # Apply CMC adjustment
  epis$tstart <-   epis$tstart +cmc.adj
  epis$tstop <-   epis$tstop +cmc.adj 
  epis$`(dob)` <-     epis$`(dob)` +cmc.adj 
  epis$`(intv)` <-     epis$`(intv)` +cmc.adj 
  
  mod.data <- epis
  formula <- f
  
  agegr=agegr
  tips=tips
  event="birth"
  weights="(weights)"
  origin=origin
  scale=scale
  dob="(dob)"
  intv="(intv)"
  tstart="tstart"
  tstop="tstop"
  
  # Function to create labels for intervals
  .epis_labels <- function(x){
    if("labels" %in% attributes(x))
      return(labels(x)[-length(x)])
    lower <- x[-length(x)]
    upper <- x[-1]-1
    val <- ifelse(lower==upper, lower, paste0(lower, "-", upper))
    gsub("-Inf", "+", val)
  }
  
  
  if(!is.null(period)){
    mod.data$period <- tcut(mod.data[[tstart]], (period-origin)*scale, .epis_labels(period))
    formula <- update(formula, ~. + period)
  }
  
  if(!is.null(agegr)){
    mod.data$agegr <- tcut(mod.data[[tstart]] - mod.data[[dob]], agegr*scale, .epis_labels(agegr))
    formula <- update(formula, ~. + agegr)
  }
  
  if(!is.null(tips)){
    mod.data$tips <- tcut(mod.data[[tstart]] - mod.data[[intv]], -rev(tips)*scale, rev(.epis_labels(tips)))
    formula <- update(formula, ~. + tips)
  }
  
  formula <- update(formula, bquote(Surv(.(as.name(tstop)) - .(as.name(tstart)), .(as.name(event))) ~ .))
  
  # Handle weights if provided
  if (!is.null(weights)) {
    mod.data <- mod.data %>%
      mutate(weights = mod.data[[weights]])
  } else {
    mod.data <- mod.data %>%
      mutate(weights = 1)
  }
  
  # Calculate person-years
  fert_pyears <- survival::pyears(formula, mod.data, scale = scale, data.frame = TRUE)$data
  
  
  if(by.individual){
    vars.to.merge <- c('caseid','v001','v005','v022','v023','v024','v025','v133','v010')
    fert_pyears <- fert_pyears %>% left_join(data.IR[,vars.to.merge],by='caseid')
    
  }else{
    vars.to.merge <- c('v001','v005','v022','v023')
    fert_pyears <- fert_pyears %>% left_join(data.IR[!duplicated(data.IR[,c('v001')]),vars.to.merge],by='v001')
  }
  
  pre <- ""
  strat <- attr(fert_pyears[, paste0(pre, "v025")], which='labels')
  names(strat) <- tolower(names(strat))
  
  fert_pyears[, paste0(pre, "v025")] <- ifelse(fert_pyears[, paste0(pre, "v025")] == strat["urban"][[1]],'urban','rural')
  fert_pyears[, "urban"] <- factor(fert_pyears[, paste0(pre, "v025")], levels = c('urban','rural'))
  fert_pyears[, paste0(pre, "v024")] <- factor(labelled::unlabelled(fert_pyears[, paste0(pre, "v024")]))
  fert_pyears[, paste0(pre, "v023")] <- factor(labelled::unlabelled(fert_pyears[, paste0(pre, "v023")]))
  fert_pyears[, paste0(pre, "v022")] <- factor(labelled::unlabelled(fert_pyears[, paste0(pre, "v022")]))
  
  
  fert_pyears <- fert_pyears %>%
    mutate(cluster = v001)
  
  fert_pyears$weights=fert_pyears$v005/100000
  
  fert_pyears$survey_year <-survey_year

  return(fert_pyears)
}







###############################################################
###  prepare maternal education
###############################################################

#' Calculate Person-Years of Exposure for Fertility
#'
#' This function calculates the person-years of exposure (women-years contributed)
#' (adapted from demogsurv package)
#'
#' @param data.IR DHS IR recode raw data.
#' @param survey_year The year of the survey.
#' @param country The country of the survey. Special adjustment for Ethiopia.
#' @param begin_year Use data starting from which year, default is 10 year before survey
#' @param end_year Use data till which year, default is survey year
#' @param by.individual whether prepare data for individual level (for individual level covariate) or aggregate to cluster level


prepare_educ <- function(data.IR) {
  
  raw.dat.tmp <-data.IR
  pre <- ""
  
  strat <- attr(raw.dat.tmp[, paste0(pre, "v025")], which='labels')
  names(strat) <- tolower(names(strat))
  
  raw.dat.tmp[, paste0(pre, "v025")] <- ifelse(raw.dat.tmp[, paste0(pre, "v025")] == strat["urban"][[1]],'urban','rural')
  raw.dat.tmp[, paste0(pre, "v025")] <- factor(raw.dat.tmp[, paste0(pre, "v025")], levels = c('urban','rural'))
  raw.dat.tmp[, paste0(pre, "v024")] <- factor(labelled::unlabelled(raw.dat.tmp[, paste0(pre, "v024")]))
  raw.dat.tmp[, paste0(pre, "v023")] <- factor(labelled::unlabelled(raw.dat.tmp[, paste0(pre, "v023")]))
  raw.dat.tmp[, paste0(pre, "v022")] <- factor(labelled::unlabelled(raw.dat.tmp[, paste0(pre, "v022")]))
  
  dat.tmp<-  raw.dat.tmp %>%
    mutate (contra.use= as.numeric(v361==1))%>%
    dplyr::  select(c(cluster= paste0(pre, "v001"),
                      householdID= paste0(pre, "v002"),
                      v022= paste0(pre, "v022"),
                      v023= paste0(pre, "v023"),
                      v024= paste0(pre, "v024"),
                      weight= paste0(pre, "v005"),
                      urban= paste0(pre, "v025"),
                      educ.yrs = v133,
                      contra.use= contra.use
                      ))  
  
  return(dat.tmp)
  
}

  
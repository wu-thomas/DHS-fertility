###############################################################################
# functions.R                                                                 #
#    Help functions for aggregate                       #
###############################################################################

# Required libraries
library(survey)
library(INLA)
library(ggplot2)
library(raster)
library(readxl)
library(stringr)

#library(foreach)
#library(doMC)

#registerDoMC(detectCores()-1)
################################################################
#########  Aggregate draws of estimates, fixed/sampled fractions 
################################################################

#' @u.draws Urban specific draws at native scale, dim: n_samp*n_admin, ordered by admin index
#' @r.draws Rural specific draws at native scale, dim: n_samp*n_admin, ordered by admin index
#' @fixed.urb Indicator whether is fixed 
#' @nSamp number of samples 
#' @admin.ref reference table for admin names and admin index



aggre_Adm_samp  = function(u.draws,r.draws,urb.frac,nSamp=1000){
  
  n_comb <- dim(urb.frac)[1]  
  
  ## fixed urban proportions
  
    
    aggre.samp <- matrix(0, nrow = nSamp, ncol = n_comb)
    
    for (i in 1:n_comb){
      
      frac <- urb.frac$urban_frac[i]
      
      if(frac == 1){
        urban.samp <- u.draws[,i]
        aggre.samp[,i] <- urban.samp
        
      }else{
        
        if(frac == 0){
          rural.samp <- r.draws[,i]
          aggre.samp[,i]  <- rural.samp
          
        }else{
          
          urban.samp <- u.draws[,i]
          rural.samp <- r.draws[,i]
          
          aggre.samp[,i]  <- urban.samp*frac+ rural.samp*(1-frac)
        }
        
      }
      
    }
    
    all.q = matrix(0, nrow = n_comb, ncol = 3)
    
    for(i in 1:n_comb){
      
      all.q[i,] = quantile(aggre.samp[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      
    }
    
    overall.res = urb.frac
    overall.res$p_Low = all.q[,1]
    overall.res$p_Med = all.q[,2]
    overall.res$p_Upp = all.q[,3]

  
  
  
  # Calculate quantiles
  U.q = matrix(0, nrow = n_comb, ncol = 3)
  R.q = matrix(0, nrow = n_comb, ncol = 3)
  
  for(i in 1:n_comb){
    U.q[i,] = quantile(u.draws[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    R.q[i,] = quantile(r.draws[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  Urban.res = urb.frac
  Urban.res$p_Low = U.q[,1]
  Urban.res$p_Med = U.q[,2]
  Urban.res$p_Upp = U.q[,3]
  
  Rural.res = urb.frac
  Rural.res$p_Low = U.q[,1]
  Rural.res$p_Med = U.q[,2]
  Rural.res$p_Upp = U.q[,3]
  
  row.names(Urban.res) <- NULL
  row.names(Rural.res) <- NULL
  row.names(overall.res) <- NULL
  
  return(list(Combined.est = overall.res,
              Urban.res = Urban.res,
              Rural.res = Rural.res,
              Combined.draws = aggre.samp,
              Urban.draws = u.draws,
              Rural.draws = r.draws ))
  
  
  
  
  #row.names(overall.res) <- NULL
  
  #return (overall.res)
  
  
  
  
}



################################################################
#########  Direct estimates admin1
################################################################

aggreDirect = function(urb.res,urb.frac){
  
  aggre.est <- vector()
  
  for (i in 1:dim(urb.frac)[1]){
    
    adm.name <- urb.frac$admin.name[i]
    frac <- urb.frac$Frac[i]
    
    if(frac == 1){
      urban.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='U',]$p_Med
      aggre.est[i] <- urban.est
      
    }else{
      
      if(frac == 0){
        rural.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='R',]$p_Med
        aggre.est[i] <- rural.est
        
      }else{
        
        urban.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='U',]$p_Med     
        rural.est <- urb.res[urb.res$admin1==adm.name & urb.res$urb =='R',]$p_Med
        
        aggre.est[i] <- urban.est*frac+ rural.est*(1-frac)
      }
      
    }
    
  }
  
  return (data.frame(admin1 = urb.frac$admin.name,
                     aggre.est = aggre.est, 
                     urb.frac =  urb.frac$Frac))
}



################################################################
#########  INLA, main effect (BYM2), UR, same spatial field
################################################################


aggAdminStrat = function(res.inla, nSamp = 1000,urb.frac){
  ###!!! urb.frac must match the year/agegrp/region in the model
  

  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_comb = dim(urb.frac)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.R =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.comb = matrix(0, nrow = nSamp, ncol = n_comb)
  
  
  # Get number of variations
  N.area = length(unique(urb.frac$region))
  N.time = length(unique(urb.frac$year))
  N.agegrp = length(unique(urb.frac$agegrp))
  
  # Get indicies in the sample
  tags = res.inla$misc$configs$contents$tag
  
  time.start.idx = res.inla$misc$configs$contents$start[which(grepl('time', tags))]
  agegrp.start.idx = res.inla$misc$configs$contents$start[which(grepl('agegrp', tags))]
  area.start.idx = res.inla$misc$configs$contents$start[which(grepl('region', tags))]
  
  Intercept.idx = res.inla$misc$configs$contents$start[which(grepl('Intercept', tags))]
  urb.idx  = res.inla$misc$configs$contents$start[which(grepl('urban', tags))]
  
  
  # Sample urban and rural
  #nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  
  for(i in 1:nSamp){
    for(j in 1:n_comb){
      
      ### get index for specific draw
      area.idx = as.integer(urb.frac[j,'region.int'])
      time.idx = as.integer(urb.frac[j,'time.int'])
      agegrp.idx = as.integer(urb.frac[j,'agegrp.int'])
      
      
      post.sample.latent = post.sample[[i]]$latent
      
      etaRur.tmp = post.sample.latent[time.start.idx + time.idx -1] + 
        post.sample.latent[area.start.idx + area.idx -1] +
        post.sample.latent[agegrp.start.idx + agegrp.idx -1] +
        post.sample.latent[Intercept.idx] 
      
      etaUrb.tmp = etaRur.tmp +  post.sample.latent[urb.idx] 
      
      
      # Nugget is overdispersion
      #cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      #est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      #est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
      est.U[i,j] = etaUrb.tmp
      est.R[i,j] = etaRur.tmp
    }
  }
  # calculate combined results 
  est.comb <- aggre_Adm_samp(u.draws=exp(est.U),r.draws=exp(est.R),urb.frac=urb.frac,
                                 nSamp=1000)
  
  return(est.comb)
  
  
  
  
}


################################################################
#########  INLA, same UR spatial field, Type IV interaction
################################################################

aggAdmin.sameUR.typeIV = function(res.inla, nSamp = 1000,urb.frac){
  
  
  
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_comb = dim(urb.frac)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.R =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.comb = matrix(0, nrow = nSamp, ncol = n_comb)
  
  
  # Get number of variations
  N.area = length(unique(urb.frac$region))
  N.time = length(unique(urb.frac$year))
  N.agegrp = length(unique(urb.frac$agegrp))
  
  # Get indicies in the sample
  tags = res.inla$misc$configs$contents$tag
  
  time.start.idx = res.inla$misc$configs$contents$start[which(grepl('time', tags))]
  agegrp.start.idx = res.inla$misc$configs$contents$start[which(grepl('agegrp.int', tags))]
  
  area.start.idx = res.inla$misc$configs$contents$start[which(grepl('region.int', tags))]

  area.agegrp.start.idx = res.inla$misc$configs$contents$start[which(grepl('agegrp.area', tags))]

  Intercept.idx = res.inla$misc$configs$contents$start[which(grepl('Intercept', tags))]
  urb.idx  = res.inla$misc$configs$contents$start[which(grepl('urban', tags))]
  
  
  # Sample urban and rural
  #nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  
  for(i in 1:nSamp){
    for(j in 1:n_comb){
      
      ### get index for specific draw
      area.idx = as.integer(urb.frac[j,'region.int'])
      time.idx = as.integer(urb.frac[j,'time.int'])
      agegrp.idx = as.integer(urb.frac[j,'agegrp.int'])
      agegrp.area.idx = as.integer(urb.frac[j,'agegrp.area'])
      
      
      post.sample.latent = post.sample[[i]]$latent
      
      etaRur.tmp = post.sample.latent[time.start.idx + time.idx -1] + 
        post.sample.latent[area.start.idx + area.idx -1] +
        post.sample.latent[agegrp.start.idx + agegrp.idx -1] +
        post.sample.latent[area.agegrp.start.idx + agegrp.area.idx -1] +
        post.sample.latent[Intercept.idx] 
      
      etaUrb.tmp = post.sample.latent[time.start.idx + time.idx -1] + 
        post.sample.latent[area.start.idx + area.idx -1] +
        post.sample.latent[agegrp.start.idx + agegrp.idx -1] +
        post.sample.latent[area.agegrp.start.idx + agegrp.area.idx -1] +
        post.sample.latent[Intercept.idx]+
        post.sample.latent[urb.idx] 
      
      
      # Nugget is overdispersion
      #cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      #est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      #est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
      est.U[i,j] = etaUrb.tmp
      est.R[i,j] = etaRur.tmp
    }
  }
  # calculate combined results 
  est.comb <- aggre_Adm_samp(u.draws=exp(est.U),r.draws=exp(est.R),urb.frac=urb.frac,
                             nSamp=1000)
  
  return(est.comb)
  
  
  
}

################################################################
#########  INLA, different UR spatial field, Type IV interaction
################################################################

aggAdmin.diffUR.typeIV = function(res.inla, nSamp = 1000,urb.frac){
  

  
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  n_comb = dim(urb.frac)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.R =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.comb = matrix(0, nrow = nSamp, ncol = n_comb)
  
  
  # Get number of variations
  N.area = length(unique(urb.frac$region))
  N.time = length(unique(urb.frac$year))
  N.agegrp = length(unique(urb.frac$agegrp))
  
  # Get indicies in the sample
  tags = res.inla$misc$configs$contents$tag
  
  time.start.idx = res.inla$misc$configs$contents$start[which(grepl('time', tags))]
  agegrp.start.idx = res.inla$misc$configs$contents$start[which(grepl('agegrp.int', tags))]
  
  area.U.start.idx = res.inla$misc$configs$contents$start[which(grepl('region.int.U', tags))]
  area.R.start.idx = res.inla$misc$configs$contents$start[which(grepl('region.int.R', tags))]
  
  area.agegrp.U.start.idx = res.inla$misc$configs$contents$start[which(grepl('agegrp.area.U', tags))]
  area.agegrp.R.start.idx = res.inla$misc$configs$contents$start[which(grepl('agegrp.area.R', tags))]
  
  Intercept.idx = res.inla$misc$configs$contents$start[which(grepl('Intercept', tags))]
  urb.idx  = res.inla$misc$configs$contents$start[which(grepl('urban', tags))]
  
  
  # Sample urban and rural
  #nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  
  for(i in 1:nSamp){
    for(j in 1:n_comb){
      
      ### get index for specific draw
      area.idx = as.integer(urb.frac[j,'region.int'])
      time.idx = as.integer(urb.frac[j,'time.int'])
      agegrp.idx = as.integer(urb.frac[j,'agegrp.int'])
      agegrp.area.idx = as.integer(urb.frac[j,'agegrp.area'])
      
      
      post.sample.latent = post.sample[[i]]$latent
      
      etaRur.tmp = post.sample.latent[time.start.idx + time.idx -1] + 
        post.sample.latent[area.R.start.idx + area.idx -1] +
        post.sample.latent[agegrp.start.idx + agegrp.idx -1] +
        post.sample.latent[area.agegrp.R.start.idx + agegrp.area.idx -1] +
        post.sample.latent[Intercept.idx] 
      
      etaUrb.tmp = post.sample.latent[time.start.idx + time.idx -1] + 
        post.sample.latent[area.U.start.idx + area.idx -1] +
        post.sample.latent[agegrp.start.idx + agegrp.idx -1] +
        post.sample.latent[area.agegrp.U.start.idx + agegrp.area.idx -1] +
        post.sample.latent[Intercept.idx]+
        post.sample.latent[urb.idx] 
      
      
      # Nugget is overdispersion
      #cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      #est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      #est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
      est.U[i,j] = etaUrb.tmp
      est.R[i,j] = etaRur.tmp
    }
  }
  # calculate combined results 
  est.comb <- aggre_Adm_samp(u.draws=exp(est.U),r.draws=exp(est.R),urb.frac=urb.frac,
                             nSamp=1000)
  
  return(est.comb)
  
  
  
}






################################################################
#########  INLA, different UR models, all Type IV interaction
################################################################

aggUR.diffmod.typeIV = function(U.res.inla,R.res.inla, nSamp = 1000,urb.frac){
  
  
  
  # Draw posterior samples
  U.post.sample = inla.posterior.sample(n = nSamp, result = U.res.inla)
  R.post.sample = inla.posterior.sample(n = nSamp, result = R.res.inla)
  
  # Initialize storage
  n_comb = dim(urb.frac)[1]
  est.U =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.R =  matrix(0, nrow = nSamp, ncol = n_comb)
  est.comb = matrix(0, nrow = nSamp, ncol = n_comb)
  
  
  # Get number of variations
  N.area = length(unique(urb.frac$region))
  N.time = length(unique(urb.frac$year))
  N.agegrp = length(unique(urb.frac$agegrp))
  
  # Get indicies in the sample for urban and rural
  
  ## urban
  tags = U.res.inla$misc$configs$contents$tag
  
  U.time.start.idx = U.res.inla$misc$configs$contents$start[which(grepl('time.int', tags))]
  U.agegrp.start.idx = U.res.inla$misc$configs$contents$start[which(grepl('agegrp.int', tags))]
  U.area.start.idx = U.res.inla$misc$configs$contents$start[which(grepl('region.int', tags))]

  U.agegrp.area.start.idx = U.res.inla$misc$configs$contents$start[which(grepl('agegrp.area', tags))]
  U.time.area.start.idx = U.res.inla$misc$configs$contents$start[which(grepl('time.area', tags))]
  U.agegrp.time.start.idx = U.res.inla$misc$configs$contents$start[which(grepl('agegrp.time', tags))]
  
  U.Intercept.idx = U.res.inla$misc$configs$contents$start[which(grepl('Intercept', tags))]
  
  ## rural
  tags = R.res.inla$misc$configs$contents$tag
  
  R.time.start.idx = R.res.inla$misc$configs$contents$start[which(grepl('time.int', tags))]
  R.agegrp.start.idx = R.res.inla$misc$configs$contents$start[which(grepl('agegrp.int', tags))]
  R.area.start.idx = R.res.inla$misc$configs$contents$start[which(grepl('region.int', tags))]
  
  R.agegrp.area.start.idx = R.res.inla$misc$configs$contents$start[which(grepl('agegrp.area', tags))]
  R.time.area.start.idx = R.res.inla$misc$configs$contents$start[which(grepl('time.area', tags))]
  R.agegrp.time.start.idx = R.res.inla$misc$configs$contents$start[which(grepl('agegrp.time', tags))]
  
  R.Intercept.idx = R.res.inla$misc$configs$contents$start[which(grepl('Intercept', tags))]
  
  
  # Sample urban and rural
  #nugSimStd = rnorm(1e4, mean = 0, sd = 1)
  #foreach(i = 1:100) %dopar% {
  for(i in 1:nSamp){
    #print(i)
    U.post.sample.latent = U.post.sample[[i]]$latent
    R.post.sample.latent = R.post.sample[[i]]$latent
    
    for(j in 1:n_comb){
      
      ### get index for specific draw
      area.idx = as.integer(urb.frac[j,'region.int'])
      time.idx = as.integer(urb.frac[j,'time.int'])
      agegrp.idx = as.integer(urb.frac[j,'agegrp.int'])
      agegrp.area.idx = as.integer(urb.frac[j,'agegrp.area'])
      time.area.idx = as.integer(urb.frac[j,'time.area'])
      agegrp.time.idx = as.integer(urb.frac[j,'agegrp.time'])
      
      
      etaRur.tmp = R.post.sample.latent[R.time.start.idx + time.idx -1] + 
        R.post.sample.latent[R.area.start.idx + area.idx -1] +
        R.post.sample.latent[R.agegrp.start.idx + agegrp.idx -1] +
        R.post.sample.latent[R.agegrp.area.start.idx + agegrp.area.idx -1] +
        R.post.sample.latent[R.time.area.start.idx + time.area.idx -1] +
        R.post.sample.latent[R.agegrp.time.start.idx + agegrp.time.idx -1] +
        R.post.sample.latent[R.Intercept.idx] 
      
      etaUrb.tmp = U.post.sample.latent[U.time.start.idx + time.idx -1] + 
        U.post.sample.latent[U.area.start.idx + area.idx -1] +
        U.post.sample.latent[U.agegrp.start.idx + agegrp.idx -1] +
        U.post.sample.latent[U.agegrp.area.start.idx + agegrp.area.idx -1] +
        U.post.sample.latent[U.time.area.start.idx + time.area.idx -1] +
        U.post.sample.latent[U.agegrp.time.start.idx + agegrp.time.idx -1] +
        U.post.sample.latent[U.Intercept.idx] 
      
      
      # Nugget is overdispersion
      #cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      #est.U[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      #est.R[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))    }
      est.U[i,j] = etaUrb.tmp
      est.R[i,j] = etaRur.tmp
    }
  }
  # calculate combined results 
  est.comb <- aggre_Adm_samp(u.draws=exp(est.U),r.draws=exp(est.R),urb.frac=urb.frac,
                             nSamp=1000)
  
  return(est.comb)
  
  
  
}



################################################################
#########  ASFR to TFR
################################################################

cal.aggre.TFR = function(aggre.draws, nSamp = 1000,urb.frac,
                                  period=NULL) {
  
  
  if(is.null(period)){
    period <- c(min(urb.frac$year):max(urb.frac$year))
    
    if(length(unique(urb.frac$year))==1){period <- c(min(urb.frac$year))}
  }
  
  
  urb.frac$ID <- c(1:dim(urb.frac)[1])
  
  region.vec <- unique(urb.frac$region)
  N.area = length(region.vec)
  
  TFR.samp = matrix(0, nrow = nSamp, ncol = N.area)
  
  
  for (i in 1:N.area){
    
    id <- urb.frac[urb.frac$region==region.vec[i] & urb.frac$year %in% period,]$ID
    
    
    TFR.samp[,i] <- rowSums(aggre.draws[,id])*5/length(period)
    
  }
  
  TFR.samp.q = matrix(0, nrow = N.area, ncol = 3)
  
  
  for(i in 1:N.area){
    TFR.samp.q[i,] = quantile(TFR.samp[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  res.frame <- data.frame(region=region.vec,
                          p_Low = TFR.samp.q[,1],
                          p_Med = TFR.samp.q[,2],
                          p_Upp = TFR.samp.q[,3])
  
  return(res.frame)
  
  
  
}

################################################################
#########  aggregate ASFR over year
################################################################

cal.aggre.ASFR.wt = function(aggre.draws, nSamp = 1000,est.grid,
                            period=NULL,alpha=0.05) {
  
  
  if(is.null(period)){
    period <- c(min(est.grid$year):max(est.grid$year))
    
    if(length(unique(est.grid$year))==1){period <- c(min(est.grid$year))}
  }
  
  
  est.grid$ID <- c(1:dim(est.grid)[1])
  
  region.vec <- unique(est.grid$region)
  N.area = length(region.vec)
  
  agegrp.vec <- unique(est.grid$agegrp)
  N.agegrp = length(agegrp.vec)
  
  ASFR.samp.list <- list()
  ASFR.res.frame.list <- list()
    
    
  for (age.idx in 1:N.agegrp){
  
  ASFR.samp = matrix(0, nrow = nSamp, ncol = N.area)

  for (i in 1:N.area){
    
    id <- est.grid[est.grid$region==region.vec[i] & est.grid$year %in% period &
                     est.grid$agegrp==agegrp.vec[age.idx],]$ID
    
    
    #TFR.samp[,i] <- rowSums(aggre.draws[,id])*5/length(period)
    ASFR.samp[,i] <- rowSums(aggre.draws[,id]*est.grid$wt[id])/sum(est.grid$wt[id])
    
  }
  
  ASFR.samp.q = matrix(0, nrow = N.area, ncol = 3)
  
  
  for(i in 1:N.area){
    ASFR.samp.q[i,] = quantile(ASFR.samp[,i], probs = c(alpha/2, 0.50, 1-alpha/2), na.rm = TRUE)
  }
  
  ASFR.res.frame <- data.frame(region=region.vec,
                          p_Low = ASFR.samp.q[,1],
                          p_Med = ASFR.samp.q[,2],
                          p_Upp = ASFR.samp.q[,3])
  
  ASFR.samp.list[[age.idx]] <- ASFR.samp
  ASFR.res.frame.list[[age.idx]] <- ASFR.res.frame
  
  }
  
  names(ASFR.samp.list) <- agegrp.vec
  names(ASFR.res.frame.list) <- agegrp.vec
  
  #return(res.frame)
  return(list(Combined.est = ASFR.res.frame.list,
              Combined.draws = ASFR.samp.list))
  
  
}

################################################################
#########  ASFR to TFR, yearly
################################################################

cal.aggre.TFR.wt = function(aggre.draws, nSamp = 1000,est.grid,
                         period=NULL,alpha=0.05) {
  
  
  if(is.null(period)){
    period <- c(min(est.grid$year):max(est.grid$year))
    
    if(length(unique(est.grid$year))==1){period <- c(min(est.grid$year))}
  }
  
  
  est.grid$ID <- c(1:dim(est.grid)[1])
  
  region.vec <- unique(est.grid$region)
  N.area = length(region.vec)
  
  TFR.samp = matrix(0, nrow = nSamp, ncol = N.area)
  
  
  for (i in 1:N.area){
    
    id <- est.grid[est.grid$region==region.vec[i] & est.grid$year %in% period,]$ID
    
    
    #TFR.samp[,i] <- rowSums(aggre.draws[,id])*5/length(period)
    TFR.samp[,i] <- rowSums(aggre.draws[,id]*est.grid$wt[id])*5/sum(est.grid$wt[id])*7
    
  }
  
  TFR.samp.q = matrix(0, nrow = N.area, ncol = 3)
  
  
  for(i in 1:N.area){
    TFR.samp.q[i,] = quantile(TFR.samp[,i], probs = c(alpha/2, 0.50, 1-alpha/2), na.rm = TRUE)
  }
  
  res.frame <- data.frame(region=region.vec,
                          p_Low = TFR.samp.q[,1],
                          p_Med = TFR.samp.q[,2],
                          p_Upp = TFR.samp.q[,3])
  
  #return(res.frame)
  return(list(Combined.est = res.frame,
              Combined.draws = TFR.samp))
  
  
}
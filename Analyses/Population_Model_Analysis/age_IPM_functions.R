## Title: Plant-Microbe meta-population model with a bayesian framework
## Purpose: functions for building matrix population model from vital rate estimates adjusted by estimates of microbial effects from the 30-bald greenhouse experiment
## Authors: Joshua Fowler
#############################################################

invlogit<-function(x){exp(x)/(1+exp(x))}

# Parameter assembly function ---------------------------------------------

# make_mods <- function(seedling_surv, surv,
#                       seedling_flw, flw,
#                       recruitment){#flw, fert, seeds_per_stem, seed_mortality, seed_germ1, seed_germ2, seedling_surv, seedling_size){
#   models <- list()
#   models$surv <- surv
#   models$seedling_surv <- seedling_surv
#   models$flw <- flw
#   models$seedling_flw <- seedling_flw
#   models$recruitment <- recruitment
#   models$seed_mortality <- seed_mortality
#   
#   
#   
#   
#   # models$seed_mortality <- seed_mortality
#   # models$seed_germ1 <- seed_germ1
#   # models$seed_germ2 <- seed_germ2  
# 
#   return(models)
# }
make_params <- function(post_draws = NA,
                        iter = NA,
                        surv_par,  sdlg_surv_par,
                        flw_par, sdlg_flw_par,
                        flw_count_par,
                        recruit_par,
                        seed_mortality,
                        recruitment_adjustment,
                        season,
                        
                        bald.rfx=F,bald = NULL,
                        year.rfx=F,year=NULL,
                        preddata,
                        microbe_off=0,germ_microbe = 0, grow_microbe = 0, flw_microbe = 0,
                        size_bounds){
                        #surv_par,surv_sdlg_par,grow_par,grow_sdlg_par,flow_par,fert_par,spike_par,seed_par,recruit_par){
params <- c()
draw <- post_draws[iter]
params$draw <- draw 
params$newdata <- preddata
params$season <- season

  if(year.rfx==F){
    year.rfx_surv <- year.rfx_sdlg.surv <- year.rfx_flw <- year.rfx_sdlg.flw <- year.rfx_rct <-  0
    params$newdata$year.t1 = NA}

  if(year.rfx==T){
    params$newdata$year.t1 <-  year}

  
  if(bald.rfx==F){
    bald.rfx_surv <- bald.rfx_sdlg.surv <- bald.rfx_flw <- bald.rfx_sdlg.flw <- bald.rfx_rct <-  0
    }
  
  if(bald.rfx==T){
    params$newdata <- preddata[preddata$bald == bald,]

     if(bald %in% parse_number(colnames(surv_par[,grepl("r_bald", colnames(surv_par))]))){
     bald.rfx_surv <- surv_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_sdlg.surv <- sdlg_surv_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_flw <- flw_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_sdlg.flw <- sdlg_flw_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_rec <- recruit_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     

     }
     
     if(!bald %in% parse_number(colnames(surv_par[,grepl("r_bald", colnames(surv_par))]))){
       bald.rfx_surv <- rnorm(n = 1, mean = 0, sd = surv_par[draw, "sd_bald__Intercept"]);
       bald.rfx_sdlg.surv <- rnorm(n = 1, mean = 0, sd = sdlg_surv_par[draw, "sd_bald__Intercept"]);
       bald.rfx_flw <- rnorm(n = 1, mean = 0, sd = flw_par[draw, "sd_bald__Intercept"]);
       bald.rfx_sdlg.flw <- rnorm(n = 1, mean = 0, sd = sdlg_flw_par[draw, "sd_bald__Intercept"]);
       bald.rfx_rec <- rnorm(n = 1, mean = 0, sd = recruit_par[draw, "sd_bald__Intercept"]);
     }
}
  
  params$year.t1 <- year
  params$bald <- bald
  params$bald.rfx_surv <- bald.rfx_surv
  params$bald.rfx_sdlg.surv <- bald.rfx_sdlg.surv
  params$bald.rfx_flw <- bald.rfx_flw
  params$bald.rfx_sdlg.flw <- bald.rfx_sdlg.flw
  params$bald.rfx_rec <- bald.rfx_rec
  
  
  
  
  #tack on size bounds
  params$max_age <- size_bounds$max_age
  params$min_age <- size_bounds$min_age
  

  
  # params$surv_fire <-  surv_par[draw,"b_time_since_fire"]
  # params$surv_elev <-  surv_par[draw,"b_rel_elev"]
  params$surv_int <- surv_par[draw, paste0("b_season",season)] + surv_par[draw,"b_time_since_fire"] * params$newdata$time_since_fire + surv_par[draw,"b_rel_elev"] * params$newdata$rel_elev + bald.rfx_surv
  params$surv_slope <-   surv_par[draw, "b_log_age.t"]
  
  params$sdlg.surv_int <- sdlg_surv_par[draw, paste0("b_season",season)] + sdlg_surv_par[draw,"b_time_since_fire"] * params$newdata$time_since_fire + sdlg_surv_par[draw,"b_rel_elev"] * params$newdata$rel_elev + bald.rfx_sdlg.surv
  

  if(season %in% c("winter", "spring")){
    params$flw_int <- params$flw_slope <-  -Inf
    params$sdlg.flw_int <- -Inf
  }

  if(!(season %in% c("winter", "spring"))){
  params$flw_int <- flw_par[draw, paste0("b_season",season)] + flw_par[draw,"b_time_since_fire"] * params$newdata$time_since_fire + flw_par[draw,"b_rel_elev"] * params$newdata$rel_elev + bald.rfx_flw
  params$flw_slope <-   flw_par[draw, "b_log_age.t"]
  
  params$sdlg.flw_int <- sdlg_flw_par[draw, paste0("b_season",season)] + sdlg_flw_par[draw,"b_time_since_fire"] * params$newdata$time_since_fire + sdlg_flw_par[draw,"b_rel_elev"] * params$newdata$rel_elev + bald.rfx_sdlg.flw
  }
  
  params$flw_count_int <- flw_count_par[draw, "b_Intercept"] +  flw_count_par[draw, "b_live_sterilesterile"]*microbe_off

  params$recruit_int <- (recruit_par[draw, paste0("b_season",season)] + recruit_par[draw,"b_time_since_fire"] * params$newdata$time_since_fire + recruit_par[draw,"b_rel_elev"] * params$newdata$rel_elev + bald.rfx_rec)
  params$recruit_adjust <- recruitment_adjustment
  params$seed_mortality <- seed_mortality

  params$microbe_off <- microbe_off
  # params$germ_microbe <- germ_microbe
  # params$grow_microbe <- grow_microbe
  # params$flw_microbe <- flw_microbe
  if(is.data.frame(germ_microbe)){
  params$germ_microbe <- germ_microbe[germ_microbe$bald == bald,][iter,]}
  if(is.data.frame(grow_microbe)){
  params$grow_microbe <- grow_microbe[grow_microbe$bald == bald,][iter,]}
  if(is.data.frame(flw_microbe)){
  params$flw_microbe <- flw_microbe[flw_microbe$bald == bald,][iter,]}
  
  if(microbe_off==0&!is.data.frame(germ_microbe)){
    params$germ_microbe$rel_diff <- 0
  }
  if(microbe_off==0&!is.data.frame(grow_microbe)){
    params$grow_microbe$rel_diff <- 0
  }
  if(microbe_off==0&!is.data.frame(flw_microbe)){
    params$flw_microbe$rel_diff <- 0
  }


  return(params)
}




# Vital rate functions ----------------------------------------------------

sx<-function(x,params){
  xb <- (pmax(pmin(x,params$max_age), params$min_age))

  int <- params$surv_int + params$surv_int*params$microbe_off*params$grow_microbe$rel_diff
  mu <- invlogit(int + params$surv_slope*log(xb))

  return(mu)
}

sdlg.s<-function(params){
  int <- params$sdlg.surv_int + params$sdlg.surv_int*params$microbe_off*params$grow_microbe$rel_diff
  mu <- invlogit(int)
  
  return(mu)
}





fx<-function(x,params){
  xb <- (pmax(pmin(x,params$max_age), params$min_age))
  if(params$season %in% c("winter", "spring")){
    flw_prob <- 0
  }else{
    int <- (params$flw_int + params$flw_int*params$microbe_off*params$flw_microbe$rel_diff)
  
  flw_prob <- invlogit(int + params$flw_slope*log(xb))
  }
  flw_count <- exp(params$flw_count_int)
  # note that this is really flw counts, but I assume there is a high success from flw to seed
  seeds <- flw_prob*flw_count # I'm pulling flw counts from the greenhouse; currently this does not depend on size/age, but I have that data in the greenhouse
  return(seeds)
}

sdlg.fx<-function(params){
  if(params$season %in% c("winter", "spring")){
    flw_prob <- 0
  }else{
  int <- params$sdlg.flw_int + params$sdlg.flw_int*params$microbe_off*params$flw_microbe$rel_diff
  flw_prob <- invlogit(int)
  }
  flw_count <- exp(params$flw_count_int)
  seeds <- flw_prob*flw_count # I'm pulling flw counts from the greenhouse
  return(seeds)
}


recruitment <- function(params){
  int <- params$recruit_int/params$recruit_adjust
  rec <- exp(params$recruit_int + params$recruit_int*params$microbe_off*params$germ_microbe$rel_diff)
  return(rec)
}



seed_bank <- function(params){
  bank <- (1-params$seed_mortality)*(1-(recruitment(params)))             
  return(bank)
}







# Bigmatrix function ------------------------------------------------------
bigmatrix<-function(params,models, extension=1){  
  matdim <- params$max_age+extension # size range of the data, extension here adds sizes larger than max size which accounts for probability density lost at predicted sizes larger than our upper bound
  y <- 1:matdim 
    #fertility transition
  Fmat <- matrix(0,matdim+1,matdim+1) # +2 dimension for  seedbank 
  

  Fmat[1,2] <- sdlg.fx(params) - sdlg.fx(params)*recruitment(params) 
  Fmat[1,3:(matdim+1)] <- fx(y[-1],params) - fx(y[-1],params)*recruitment(params) 
  Fmat[2,2] <-   sdlg.fx(params)*recruitment(params) # recruiting directly as seedling
  Fmat[2,3:(matdim+1)] <- fx(y[-1],params)*recruitment(params) # recruiting directly as seedling
    
  Fmat[1,1] <- seed_bank(params)
  Fmat[2,1] <- 1- recruitment(params)
  

  #growth/survival transition
  Tmat <-matrix(0,matdim+1,matdim+1)
  Tmat[2,1] <- sdlg.s(params)
  
  # assign survival probability to the subdiagonal
  diag(Tmat[-1,-ncol(Tmat)])[-1] <- sx(y[-1], params)

  

  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}


# Annual transition matrix as product of seasonal (spring, summer, autumn, winter) matrices
# P is spring to summer, Q is summer to autumn, R is autumn to winter, S is winter to spring
annual.matrix <- function(P, Q, R, S, interval){
  if(interval == "spring"){
    Annual <- S %*% R %*% Q %*% P
  }
  if(interval == "summer"){
    Annual <- P %*% S %*% R %*% Q
  }
  if(interval == "fall"){
    Annual <- Q %*% P %*% S %*% R
  }
  if(interval == "winter"){
    Annual <- R %*% Q %*% P %*% S
  }
  
  return(Annual)
} 

# checking for eviction
# plot(y, sx(y,models,params), xlab="size", type="l",
#      ylab="Survival Probability", lwd = 2)
# points(y,apply(Tmat[2:matdim+1,1:matdim+1],2,sum), col="red",lwd=1,cex = .5,pch=19)

# for ERYCUN, I found that 25 size bins, and a relatively small extension (2) is sufficient
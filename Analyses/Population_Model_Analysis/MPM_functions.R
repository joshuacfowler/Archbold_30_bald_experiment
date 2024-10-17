## Title: Plant-Microbe meta-population model with a bayesian framework
## Purpose: functions for building matrix population model from vital rate estimates adjusted by estimates of microbial effects from the 30-bald greenhouse experiment
## Authors: Joshua Fowler
#############################################################

invlogit<-function(x){exp(x)/(1+exp(x))}

# Parameter assembly function ---------------------------------------------
make_mods <- function(grow, surv, flw, fert, head, recruit){
  models <- list()
  models$grow <- grow
  models$surv <- surv
  models$flw <- flw
  models$fert <- fert
  models$head <- head
  models$recruit <- recruit
  return(models)
}
make_params <- function(bald.rfx=F,year.rfx=F,year=NULL,size_bounds){
                        #surv_par,surv_sdlg_par,grow_par,grow_sdlg_par,flow_par,fert_par,spike_par,seed_par,recruit_par){

  # newdata <- data.fram(log_size.t = 
  # if(plot.rfx==F){
  #   
  #   plot.rfx_surv <- plot.rfx_surv_sdlg <- plot.rfx_grow <- plot.rfx_grow_sdlg <- plot.rfx_flow <- plot.rfx_fert <- plot.rfx_spike <- plot.rfx_rct <-  0}
  # # 
  # # if(bald.rfx==T){
  # #   ## timing and survival and growth (size_t / y_t1) is meant to line up with reproduction (size_t1 / y_t1)
  # #   plot.rfx_surv <- surv_par$tau_year[draw,species,(endo_var_U+1),(year)]; 
  # #   plot.rfx_surv_sdlg <-surv_sdlg_par$tau_year[draw,species,(endo_var_U+1),(year)];
  # #   plot.rfx_grow <- grow_par$tau_year[draw,species,(endo_var_U+1),(year)];
  # #   plot.rfx_grow_sdlg <- grow_sdlg_par$tau_year[draw,species,(endo_var_U+1),(year)];
  # #   plot.rfx_flow <- flow_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)]; # fitting 
  # #   plot.rfx_fert <- fert_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)]; 
  # #   plot.rfx_spike <- spike_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)];
  # #   plot.rfx_rct <- recruit_par$tau_year[draw,species,(endo_var_F+1),year];
  # # }
  # # 
  # 
  params <- c()
  # #survival
  # params$surv_int <- surv_par$beta0[draw,species] + 
  #   endo_mean_U * surv_par$betaendo[draw,species] + 
  #   original * surv_par$betaorigin[draw] +
  #   rfx_surv
  # params$surv_spei <- spei_surv
  # params$surv_slope <- surv_par$betasize[draw,species]
  # params$surv_slope_2 <- surv_par$betasize_2[draw,species]
  # 
  # # seedling survival
  # params$surv_sdlg_int <- surv_sdlg_par$beta0[draw,species] + 
  #   endo_mean_U * surv_sdlg_par$betaendo[draw,species] + 
  #   rfx_surv_sdlg
  # params$surv_sdlg_spei <- spei_surv_sdlg
  # 
  # #growth
  # params$grow_int <- grow_par$beta0[draw,species] + 
  #   endo_mean_U * grow_par$betaendo[draw,species] + 
  #   original * grow_par$betaorigin[draw] +
  #   rfx_grow
  # params$grow_spei <- spei_grow
  # params$grow_slope <- grow_par$betasize[draw,species] 
  # params$grow_slope_2 <- grow_par$betasize_2[draw,species]  
  # 
  # params$grow_sigma <- grow_par$sigma[draw] 
  # # seedling growth
  # params$grow_sdlg_int <- grow_sdlg_par$beta0[draw,species] + 
  #   endo_mean_U * grow_sdlg_par$betaendo[draw,species] + 
  #   rfx_grow_sdlg
  # params$grow_sdlg_spei <- spei_grow_sdlg
  # params$grow_sdlg_sigma <- grow_sdlg_par$sigma[draw] 
  # 
  # #flowering
  # params$flow_int <- flow_par$beta0[draw,species] + 
  #   endo_mean_F * flow_par$betaendo[draw,species] + 
  #   original * flow_par$betaorigin[draw] +
  #   spei_flow + rfx_flow
  # params$flow_spei <- spei_flow
  # params$flow_slope <- flow_par$betasize[draw,species]  
  # params$flow_slope_2 <- flow_par$betasize_2[draw,species]  
  # 
  # #fertility
  # params$fert_int <- fert_par$beta0[draw,species] +
  #   endo_mean_F * fert_par$betaendo[draw,species] +
  #   original * fert_par$betaorigin[draw] +
  #   rfx_fert
  # params$fert_spei <- spei_fert
  # params$fert_slope <- fert_par$betasize[draw,species]
  # params$fert_slope_2 <- fert_par$betasize_2[draw,species]
  # 
  # 
  # #spikelets
  # params$spike_int <- spike_par$beta0[draw,species]  +
  #   endo_mean_F * spike_par$betaendo[draw,species] +
  #   original * spike_par$betaorigin[draw] +
  #   rfx_spike
  # params$spike_spei <- spei_spike
  # params$spike_slope <- spike_par$betasize[draw,species]  
  # params$spike_slope_2 <- spike_par$betasize_2[draw,species]  
  # 
  # 
  # #seeds per spikelet
  # params$seeds_per_spike <- seed_par$beta0[draw,species] + 
  #   endo_mean_F * seed_par$betaendo[draw,species]
  # #recruits per seed
  # params$recruits_per_seed <- recruit_par$beta0[draw,species] + 
  #   endo_mean_F * recruit_par$betaendo[draw,species] +
  #   rfx_rct
  # params$recruits_spei <- spei_rct
  params$year.t1 <- year
  params$bald <- bald.rfx
  
  #tack on size bounds
  params$max_size <- size_bounds$max_size
  params$min_size <- size_bounds$min_size
  
  return(params)
}


# Vital rate functions ----------------------------------------------------
sx<-function(x,models,params){
  xb<-pmin(x,params$max_size) # any predicted plants larger than max size will set to be max size
  newdata <- data.frame(log_size.t = xb, bald = params$year.t1, bald = params$bald)
  posterior_epred(object = models$surv, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,]
  # invlogit(params$surv_int + params$surv_slope*log(xb) + params$surv_slope_2*(log(xb)^2)*quadratic)
}


gxy <- function(x,y,models,params){
  xb<-pmin(x,params$max_size)
  newdata <- data.frame(log_size.t = xb, year.t1 = params$year.t1, bald = params$bald)
  pred_mu <- posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1, dpar = c("mu"))[1,]
  pred_sigma <- posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1, dpar = c("sigma"))[1,]
  return(dnorm(x=y, 
                  mean=pred_mu,
                  sd = exp(pred_sigma)))
}


pxy<-function(x,y,models,params){
  sx(x,models,params) * gxy(x,y,models,params)
}



fx<-function(x,models,params){
  xb<-pmin(x,params$max_size) # any predicted plants larger than max size will set to be max size
  newdata <- data.frame(log_size.t = xb, year.t1 = params$year.t1, bald = params$bald)
  flw_prob <- posterior_epred(object = models$flw, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,]
  flw_stem <- exp(posterior_linpred(object = models$fert, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,])
  # flw_head <- posterior_linpred(object = models$head, newdata = newdata, re_formula = NULL,, allow_new_levels = TRUE, ndraws = 1)[1,]
  flw_head <- models$head
  recruit <- models$recruit
  seedlings <- flw_prob * flw_stem * flw_head * recruit
  return(seedlings)
}



# Bigmatrix function ------------------------------------------------------
bigmatrix<-function(params,models,matdim, extension=1){   
  Lb <- params$min_size - extension # size range of the data, extension here adds sizes smaller than min size which accounts for probability density lost at predicted sizes smaller than our lower bound
  Ub <- params$max_size + extension # size range of the data, extension here adds sizes larger than max size which accounts for probability density lost at predicted sizes larger than our upper bound
  h <- (Ub-Lb)/matdim 
  y <- Lb + (1:matdim)*h - h/2# implementing the midpoint rule to create a vector of sizes from min to max (+ extension)
  #fertility transition
  Fmat <- matrix(0,matdim,matdim)
  Fmat[1,1:(matdim)]<-fx(x = y, models,params)
  
  #growth/survival transition
  Tmat <-matrix(0,matdim,matdim)
  Tmat[1:(matdim),1:(matdim)] <- t(outer(y,y, FUN = pxy,models=models,params=params))
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}




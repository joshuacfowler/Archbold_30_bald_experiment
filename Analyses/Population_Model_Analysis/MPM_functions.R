## Title: Plant-Microbe meta-population model with a bayesian framework
## Purpose: functions for building matrix population model from vital rate estimates adjusted by estimates of microbial effects from the 30-bald greenhouse experiment
## Authors: Joshua Fowler
#############################################################

invlogit<-function(x){exp(x)/(1+exp(x))}

# Parameter assembly function ---------------------------------------------

make_mods <- function(grow, surv, flw, fert, seeds_per_stem, seed_mortality, seed_germ1, seed_germ2, seedling_surv, seedling_size){
  models <- list()
  models$grow <- grow
  models$surv <- surv
  models$flw <- flw
  models$fert <- fert
  models$seeds_per_stem <- seeds_per_stem
  models$seed_mortality <- seed_mortality
  models$seed_germ1 <- seed_germ1
  models$seed_germ2 <- seed_germ2  
  models$seedling_surv <- seedling_surv
  models$seedling_size <- seedling_size
  return(models)
}
make_params <- function(bald.rfx=F,year.rfx=F,year=NULL,preddata,
                        microbe_off=0,germ_microbe, grow_microbe, flw_microbe,
                        size_bounds){
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

  params$year.t1 <- year
  params$bald <- bald.rfx
  
  
  params$microbe_off <- microbe_off
  params$germ_microbe <- germ_microbe[germ_microbe$bald == bald.rfx,]
  params$grow_microbe <- grow_microbe[grow_microbe$bald == bald.rfx,]
  params$flw_microbe <- flw_microbe[flw_microbe$bald == bald.rfx,]
  
  #tack on size bounds
  params$max_size <- size_bounds$max_size
  params$min_size <- size_bounds$min_size
  
  params$newdata <- preddata[preddata$bald == bald.rfx,]
  params$newdata$year.t1 <-  year
  return(params)
}


# Vital rate functions ----------------------------------------------------
sx<-function(x,models,params){
  xb<-data.frame(log_size.t = pmin(x,params$max_size))
  newdata <- merge(params$newdata, xb)
  mu <- posterior_epred(object = models$surv, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,]
  # invlogit(params$surv_int + params$surv_slope*log(xb) + params$surv_slope_2*(log(xb)^2)*quadratic)
 return(mu)
}


gxy <- function(x,y,models,params){
  xb<-data.frame(log_size.t = pmin(x,params$max_size))
  newdata <- merge(params$newdata, xb)
  pred_mu <- posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1, dpar = c("mu"))[1,]
  pred_sigma <- posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1, dpar = c("sigma"))[1,]
  grow <- dnorm(x=y, 
                mean=pred_mu,
                sd = exp(pred_sigma))
  grow1 <- grow+grow*params$microbe_off*params$grow_microbe$rel_diff              
  return(grow1)
}


pxy<-function(x,y,models,params){
  sx(x,models,params) * gxy(x,y,models,params)
}




fx<-function(x,models,params){
  xb<-data.frame(log_size.t = pmin(x,params$max_size))
  newdata <- merge(params$newdata, xb)
  flw_prob1 <- posterior_epred(object = models$flw, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,]
  flw_prob <- flw_prob1+flw_prob1*params$microbe_off*params$flw_microbe$rel_diff
  flw_stem <- exp(posterior_linpred(object = models$fert, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,])
  # flw_head <- posterior_linpred(object = models$head, newdata = newdata, re_formula = NULL,, allow_new_levels = TRUE, ndraws = 1)[1,]
  # recruit <- models$recruit + models$recruit*params$microbe_off*params$germ_microbe$rel_diff
  germination <- models$seed_germ1 + models$seed_germ1*params$microbe_off*params$germ_microbe$rel_diff
  
  recruits1 <- flw_prob * flw_stem * models$seeds_per_stem * germination
  return(recruits1)
}


seed_production <- function(x,models,params){
  xb<-data.frame(log_size.t = pmin(x,params$max_size))
  newdata <- merge(params$newdata, xb)
  flw_prob1 <- posterior_epred(object = models$flw, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,]
  flw_prob <- flw_prob1+flw_prob1*params$microbe_off*params$flw_microbe$rel_diff
  flw_stem <- exp(posterior_linpred(object = models$fert, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,])
  # flw_head <- posterior_linpred(object = models$head, newdata = newdata, re_formula = NULL,, allow_new_levels = TRUE, ndraws = 1)[1,]
  # recruit <- models$recruit + models$recruit*params$microbe_off*params$germ_microbe$rel_diff
  germination <- models$seed_germ1 + models$seed_germ1*params$microbe_off*params$germ_microbe$rel_diff
  seeds <- flw_prob * flw_stem * models$seeds_per_stem * (1-germination)*(1-models$seed_mortality)
  return(seeds)
}



seed_bank <- function(models, params){
  bank <- (1-models$seed_mortality)*(1-(models$seed_germ2 + models$seed_germ2*params$microbe_off*params$germ_microbe$rel_diff))             
  return(bank)
}



# Filling out the recruitment from seed bank size distribution. Currently I am just filling these with the prediction for the minimum size, but I can do this better.
recruit_surv <- function(models, params){
  xb<-data.frame(log_size.t = params$min_size)
  newdata <- merge(params$newdata, xb)
  mu <- posterior_epred(object = models$seedling_surv, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1)[1,]
  # invlogit(params$surv_int + params$surv_slope*log(xb) + params$surv_slope_2*(log(xb)^2)*quadratic)
  #return(mu)
}
# 
recruit_size <- function(y,models, params){
  xb<-data.frame(log_size.t = params$min_size)
  newdata <- merge(params$newdata, xb)
  pred_mu <- posterior_linpred(object = models$seedling_size, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1, dpar = c("mu"))[1,]
  pred_sigma <- posterior_linpred(object = models$seedling_size, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 1, dpar = c("sigma"))[1,]
  grow <- dnorm(x=y, 
                mean=pred_mu,
                sd = exp(pred_sigma))
  grow1 <- grow+grow*params$microbe_off*params$grow_microbe$rel_diff              
  return(grow1)
}

p_rec_y <- function(y, models, params){
  recruit_surv(models, params) * recruit_size(y, models, params)
}

emergence <- function(y,models, params){
  p_rec_y(y,models, params) * (models$seed_germ2 + models$seed_germ2*params$microbe_off*params$germ_microbe$rel_diff)
}





# Bigmatrix function ------------------------------------------------------
bigmatrix<-function(params,models,matdim, extension=.1){   
  Lb <- params$min_size # size range of the data, extension here adds sizes smaller than min size which accounts for probability density lost at predicted sizes smaller than our lower bound
  Ub <- params$max_size*(1+extension) # size range of the data, extension here adds sizes larger than max size which accounts for probability density lost at predicted sizes larger than our upper bound
  h <- (Ub-Lb)/matdim 
  y <- Lb + (1:matdim)*h - h/2# implementing the midpoint rule to create a vector of sizes from min to max (+ extension)
  #fertility transition
  Fmat <- matrix(0,matdim+1,matdim+1) # +1 dimension for seedbank
  
  
  Fmat[2:(matdim+1), 2:matdim+1]<-fx(x = y, models, params) # seeds recruiting immediately
  Fmat[1,2:(matdim+1)] <- seed_production(x = y, models, params)
  Fmat[1,1] <- seed_bank(models, params)
  
  #growth/survival transition
  Tmat <-matrix(0,matdim+1,matdim+1)
  Tmat[2:(matdim+1),1] <- emergence(y,models=models,params=params)

  
  Tmat[2:(matdim+1),2:(matdim+1)] <- h*t(outer(y,y, FUN = pxy,models=models,params=params))
  
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}




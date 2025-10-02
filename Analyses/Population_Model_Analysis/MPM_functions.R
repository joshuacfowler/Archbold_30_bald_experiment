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
make_params <- function(draw = NA,
                        surv_par, grow_par, flw_par, fert_par,
                        bald.rfx=F,bald = NULL,
                        year.rfx=F,year=NULL,
                        preddata,
                        microbe_off=0,germ_microbe = 0, grow_microbe = 0, flw_microbe = 0,
                        size_bounds){
                        #surv_par,surv_sdlg_par,grow_par,grow_sdlg_par,flow_par,fert_par,spike_par,seed_par,recruit_par){
params <- c()
params$draw <- draw
params$newdata <- preddata

  if(year.rfx==F){
    year.rfx_surv <- year.rfx_grow <- year.rfx_flw <- year.rfx_fert <- year.rfx_spike <- year.rfx_rct <-  0
    params$newdata$year.t1 = NA}

  if(year.rfx==T){
    params$newdata$year.t1 <-  year}

  
  if(bald.rfx==F){
    bald.rfx_surv <- bald.rfx_grow <- bald.rfx_flw <- bald.rfx_fert <- bald.rfx_spike <- bald.rfx_rct <-  0}
  
  if(bald.rfx==T){
    params$newdata <- preddata[preddata$bald == bald,]

     if(bald %in% parse_number(colnames(surv_par[,grepl("r_bald", colnames(surv_par))]))){
     bald.rfx_surv <- surv_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_grow <- grow_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_flw <- flw_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     bald.rfx_fert <- fert_par[draw,paste0("r_bald[", bald, ",Intercept]")]; 
     
     # bald.rfx_grow <- grow_par$tau_year[draw,species,(endo_var_U+1),(year)];
     # bald.rfx_flow <- flow_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)]; # fitting 
     # bald.rfx_fert <- fert_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)]; 
     # bald.rfx_spike <- spike_par$tau_year[draw,species,(endo_var_F+1),(year-repro_offset)];
     # bald.rfx_rct <- recruit_par$tau_year[draw,species,(endo_var_F+1),year];
     }
     
     if(!bald %in% parse_number(colnames(surv_par[,grepl("r_bald", colnames(surv_par))]))){
       bald.rfx_surv <- rnorm(n = 1, mean = 0, sd = surv_par[draw, "sd_bald__Intercept"]);
       bald.rfx_grow <- rnorm(n = 1, mean = 0, sd = grow_par[draw, "sd_bald__Intercept"]);
       bald.rfx_flw <- rnorm(n = 1, mean = 0, sd = flw_par[draw, "sd_bald__Intercept"]);
       bald.rfx_fert <- rnorm(n = 1, mean = 0, sd = fert_par[draw, "sd_bald__Intercept"]);
       
     }
}
  
  params$year.t1 <- year
  params$bald <- bald
  params$bald.rfx_surv <- bald.rfx_surv
  params$bald.rfx_grow <- bald.rfx_grow
  params$bald.rfx_flw <- bald.rfx_flw
  params$bald.rfx_fert <- bald.rfx_fert
  
  
  
  
  #tack on size bounds
  params$max_size <- size_bounds$max_size
  params$min_size <- size_bounds$min_size
  

  


  params$surv_fire <- surv_par[draw,"b_time_since_fire"]
  params$surv_elev <-  surv_par[draw,"b_rel_elev"]
  params$surv_int <- surv_par[draw,"b_Intercept"] 
  # compiling GAM related parameters for size-relationship with help from BRMS helper functions in each function
  # surv_prep <- brms:::prepare_predictions(surv_mod, newdata = new_dat, allow_new_levels = FALSE)
  # lin_pred_prep = surv_prep$dpars$mu$fe$b[draw,"b_Intercept"] + surv_prep$dpars$mu$sm$fe$Xs*mean(surv_prep$dpars$mu$sm$fe$bs) + surv_prep$dpars$mu$sm$re$sx$Zs[[1]] %*% colMeans(surv_prep$dpars$mu$sm$re$sx$s[[1]])
  
  params$grow_fire <- grow_par[draw,"b_time_since_fire"]
  params$grow_elev <-  grow_par[draw,"b_rel_elev"]
  params$grow_int <- grow_par[draw,"b_Intercept"] 
  
  params$flw_fire <- flw_par[draw,"b_time_since_fire"]
  params$flw_elev <-  flw_par[draw,"b_rel_elev"]
  params$flw_int <- flw_par[draw,"b_Intercept"] 
  
  params$fert_fire <- fert_par[draw,"b_time_since_fire"]
  params$fert_elev <-  fert_par[draw,"b_rel_elev"]
  params$fert_int <- fert_par[draw,"b_Intercept"] 
  params$fert_size <- fert_par[draw,"b_log_size.t"] 
  params$fert_shape <- fert_par[draw,"shape"]
  
  
  params$microbe_off <- microbe_off
  # params$germ_microbe <- germ_microbe
  # params$grow_microbe <- grow_microbe
  # params$flw_microbe <- flw_microbe
  if(is.data.frame(germ_microbe)){
  params$germ_microbe <- germ_microbe[germ_microbe$bald == bald,]}
  if(is.data.frame(grow_microbe)){
  params$grow_microbe <- grow_microbe[grow_microbe$bald == bald,]}
  if(is.data.frame(flw_microbe)){
  params$flw_microbe <- flw_microbe[flw_microbe$bald == bald,]}
  

  return(params)
}


# Vital rate functions ----------------------------------------------------

sx<-function(x,models,params){
  xb <- data.frame(log_size.t = pmax(pmin(x,params$max_size), params$min_size), surv.t1 = 0) # note: the prepare_predictions function requires a response variable, in this case has to be bernoulli, not just NA
  newdata <- merge(params$newdata, xb)
  prep <- brms:::prepare_predictions(models$surv, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw)
  
  # lin_pred_prep = surv_prep$dpars$mu$fe$b[draw,"b_Intercept"] + surv_prep$dpars$mu$sm$fe$Xs*mean(surv_prep$dpars$mu$sm$fe$bs) + surv_prep$dpars$mu$sm$re$sx$Zs[[1]] %*% colMeans(surv_prep$dpars$mu$sm$re$sx$s[[1]])
  
  int <- params$surv_int + params$surv_elev*newdata$rel_elev + params$surv_fire*newdata$time_since_fire  
  linpred <- int + prep$dpars$mu$sm$fe$Xs[,1]*prep$dpars$mu$sm$fe$bs[,1] + prep$dpars$mu$sm$re$slog_size.t$Zs[[1]] %*% prep$dpars$mu$sm$re$slog_size.t$s[[1]][1,]
  mu <- c(invlogit(linpred))
  # invlogit(params$surv_int + params$surv_slope*log(xb) + params$surv_slope_2*(log(xb)^2)*quadratic)
 return(mu)
}



gxy <- function(x,y,models,params){
  xb <- data.frame(log_size.t = pmax(pmin(x,params$max_size), params$min_size), log_size.t1 = NA_real_) # note: the prepare_predictions function requires a response variable, in this case has to be bernoulli, not just NA
  newdata <- merge(params$newdata, xb)
  prep <- brms:::prepare_predictions(models$grow, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw)
  
  int <- params$grow_int + params$grow_elev*newdata$rel_elev + params$grow_fire*newdata$time_since_fire  
  # adjust the intercept by the relative effect of microbiome for a given bald
  int_1 <- int + int*params$microbe_off*params$grow_microbe$rel_diff
  linpred <- int_1 + prep$dpars$mu$sm$fe$Xs[,1]*prep$dpars$mu$sm$fe$bs[,1] + prep$dpars$mu$sm$re$slog_size.t$Zs[[1]] %*% prep$dpars$mu$sm$re$slog_size.t$s[[1]][1,]
  
  # pred_mu <- mean(posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 500, dpar = c("mu"), transform = TRUE))
  pred_sigma <- posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, draw_ids = params$draw, dpar = c("sigma"))
  grow <- dnorm(x=y, 
                mean=linpred,
                sd = exp(pred_sigma))
  return(grow)
}
# matplot(t(grow), type = "l")
# matplot(t(grow1), type = "l")


pxy<-function(x,y,models,params){
  sx(x,models,params) * gxy(x,y,models,params)
}




fx<-function(x,models,params){
  
  xb <- data.frame(log_size.t = pmax(pmin(x,params$max_size), params$min_size), flw_status.t = 0) # note: the prepare_predictions function requires a response variable, in this case has to be bernoulli, not just NA
  newdata <- merge(params$newdata, xb)
  prep <- brms:::prepare_predictions(models$flw, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw, sample_new_levels = "gaussian")
  
  # lin_pred_prep = surv_prep$dpars$mu$fe$b[draw,"b_Intercept"] + surv_prep$dpars$mu$sm$fe$Xs*mean(surv_prep$dpars$mu$sm$fe$bs) + surv_prep$dpars$mu$sm$re$sx$Zs[[1]] %*% colMeans(surv_prep$dpars$mu$sm$re$sx$s[[1]])
  
  int.flw <- params$flw_int + params$flw_elev*params$newdata$rel_elev + params$flw_fire*params$newdata$time_since_fire  
  # adjust the intercept by the relative effect of microbiome for a given bald
  int.flw_1 <- int.flw + int.flw*params$microbe_off*params$flw_microbe$rel_diff
  
  linpred.flw <- c(int.flw_1) + 
    c(prep$dpars$mu$sm$fe$Xs[,1]*prep$dpars$mu$sm$fe$bs[,1]) + c(prep$dpars$mu$sm$re$slog_size.t$Zs[[1]] %*% prep$dpars$mu$sm$re$slog_size.t$s[[1]][1,] )+ # adding in the smooth terms
    c(prep$dpars$mu$re$Z$bald[1]*prep$dpars$mu$re$r$bald) # adding in the bald random effects matrix
    
  flw_prob <- as.numeric(invlogit(linpred.flw))
  
  # implementing the truncation here by rescaling the expectation by the probability of Y=0, given by (1/1+shape*mu)^1/shape
  trunc <- newdata[1,];trunc$log_size.t <- 0
  fert_mu <- exp(posterior_linpred(models$fert, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw))
  
  flw_stem <- fert_mu/(1- ((1/(1+params$fert_shape*fert_mu))^(1/params$fert_shape)))

  # flw_head <- posterior_linpred(object = models$head, newdata = newdata, re_formula = NULL,, allow_new_levels = TRUE, ndraws = 1)[1,]
  # recruit <- models$recruit + models$recruit*params$microbe_off*params$germ_microbe$rel_diff
  germination <- models$seed_germ1 + models$seed_germ1*params$microbe_off*params$germ_microbe$rel_diff
  
  recruits1 <- c(flw_prob * flw_stem * models$seeds_per_stem * germination)
  return(recruits1)
}



seed_production <- function(x,models,params){
  xb <- data.frame(log_size.t = pmax(pmin(x,params$max_size), params$min_size), flw_status.t = 0) # note: the prepare_predictions function requires a response variable, in this case has to be bernoulli, not just NA
  newdata <- merge(params$newdata, xb)
  prep <- brms:::prepare_predictions(models$flw, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw, sample_new_levels = "gaussian")
  
  # lin_pred_prep = surv_prep$dpars$mu$fe$b[draw,"b_Intercept"] + surv_prep$dpars$mu$sm$fe$Xs*mean(surv_prep$dpars$mu$sm$fe$bs) + surv_prep$dpars$mu$sm$re$sx$Zs[[1]] %*% colMeans(surv_prep$dpars$mu$sm$re$sx$s[[1]])
  
  int.flw <- params$flw_int + params$flw_elev*params$newdata$rel_elev + params$flw_fire*params$newdata$time_since_fire  
  # adjust the intercept by the relative effect of microbiome for a given bald
  int.flw_1 <- int.flw + int.flw*params$microbe_off*params$flw_microbe$rel_diff
  
  linpred.flw <- c(int.flw_1) + 
    c(prep$dpars$mu$sm$fe$Xs[,1]*prep$dpars$mu$sm$fe$bs[,1]) + c(prep$dpars$mu$sm$re$slog_size.t$Zs[[1]] %*% prep$dpars$mu$sm$re$slog_size.t$s[[1]][1,] )+ # adding in the smooth terms
    c(prep$dpars$mu$re$Z$bald[1]*prep$dpars$mu$re$r$bald) # adding in the bald random effects matrix

  flw_prob <- as.numeric(invlogit(linpred.flw))
  
  # implementing the truncation here by rescaling the expectation by the probability of Y=0, given by (1/1+shape*mu)^1/shape
  trunc <- newdata[1,];trunc$log_size.t <- 0
  fert_mu <- exp(posterior_linpred(models$fert, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw))
  
  flw_stem <- fert_mu/(1- ((1/(1+params$fert_shape*fert_mu))^(1/params$fert_shape)))
  
  # flw_head <- posterior_linpred(object = models$head, newdata = newdata, re_formula = NULL,, allow_new_levels = TRUE, ndraws = 1)[1,]
  # recruit <- models$recruit + models$recruit*params$microbe_off*params$germ_microbe$rel_diff
  germination <- models$seed_germ1 + models$seed_germ1*params$microbe_off*params$germ_microbe$rel_diff
  seeds <- c(flw_prob * flw_stem * models$seeds_per_stem * (1-germination)*(1-models$seed_mortality))
  return(seeds)
}



seed_bank <- function(models, params){
  bank <- (1-models$seed_mortality)*(1-(models$seed_germ2 + models$seed_germ2*params$microbe_off*params$germ_microbe$rel_diff))             
  return(bank)
}



# Filling out the recruitment from seed bank size distribution. Currently I am just filling these with the prediction for the minimum size, but I can do this better.
recruit_surv <- function(models, params){
  
  xb <- data.frame(log_size.t = params$min_size, surv.t1 = 0) # note: the prepare_predictions function requires a response variable, in this case has to be bernoulli, not just NA
  newdata <- merge(params$newdata, xb)
  prep <- brms:::prepare_predictions(models$seedling_surv, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw)
  
  # lin_pred_prep = surv_prep$dpars$mu$fe$b[draw,"b_Intercept"] + surv_prep$dpars$mu$sm$fe$Xs*mean(surv_prep$dpars$mu$sm$fe$bs) + surv_prep$dpars$mu$sm$re$sx$Zs[[1]] %*% colMeans(surv_prep$dpars$mu$sm$re$sx$s[[1]])
  
  int <- params$surv_int + params$surv_elev*newdata$rel_elev + params$surv_fire*newdata$time_since_fire  
  linpred <- c(int) + prep$dpars$mu$sm$fe$Xs[,1]*prep$dpars$mu$sm$fe$bs[,1] + prep$dpars$mu$sm$re$slog_size.t$Zs[[1]] %*% prep$dpars$mu$sm$re$slog_size.t$s[[1]][1,]
  mu <- c(invlogit(linpred))
  # invlogit(params$surv_int + params$surv_slope*log(xb) + params$surv_slope_2*(log(xb)^2)*quadratic)
  return(mu)
}
# 
recruit_size <- function(y,models, params){
  
  xb <-   xb<-data.frame(log_size.t = params$min_size, log_size.t1 = NA_real_) # note: the prepare_predictions function requires a response variable, in this case has to be bernoulli, not just NA
  newdata <- merge(params$newdata, xb)
  prep <- brms:::prepare_predictions(models$seedling_size, newdata = newdata, allow_new_levels = TRUE, draw_ids = params$draw)
  
  int <- params$grow_int + params$grow_elev*newdata$rel_elev + params$grow_fire*newdata$time_since_fire  
  # adjust the intercept by the relative effect of microbiome for a given bald
  int_1 <- int + int*params$microbe_off*params$grow_microbe$rel_diff
  linpred <- c(int_1) + prep$dpars$mu$sm$fe$Xs[,1]*prep$dpars$mu$sm$fe$bs[,1] + prep$dpars$mu$sm$re$slog_size.t$Zs[[1]] %*% prep$dpars$mu$sm$re$slog_size.t$s[[1]][1,]
  
  # pred_mu <- mean(posterior_linpred(object = models$grow, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, ndraws = 500, dpar = c("mu"), transform = TRUE))
  pred_sigma <- posterior_linpred(object = models$seedling_size, newdata = newdata, re_formula = NULL, allow_new_levels = TRUE, draw_ids = params$draw, dpar = c("sigma"))
  grow <- dnorm(x=y, 
                mean=linpred,
                sd = exp(pred_sigma))
  return(grow)
}

p_rec_y <- function(y, models, params){
  recruit_surv(models, params) * recruit_size(y, models, params)
}

emergence <- function(y,models, params){
  p_rec_y(y,models, params) * (models$seed_germ2 + models$seed_germ2*params$microbe_off*params$germ_microbe$rel_diff)
}





# Bigmatrix function ------------------------------------------------------
bigmatrix<-function(params,models,matdim, extension=.1){  
  Lb <- params$min_size-extension # size range of the data, extension here adds sizes smaller than min size which accounts for probability density lost at predicted sizes smaller than our lower bound
  Ub <- params$max_size+extension # size range of the data, extension here adds sizes larger than max size which accounts for probability density lost at predicted sizes larger than our upper bound
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

  
  Tmat[2:(matdim+1),2:(matdim+1)] <- h*t(outer(y,y, pxy,models=models,params=params))
  
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}

# checking for eviction
# plot(y, sx(y,models,params), xlab="size", type="l",
#      ylab="Survival Probability", lwd = 2)
# points(y,apply(Tmat[2:matdim+1,1:matdim+1],2,sum), col="red",lwd=1,cex = .5,pch=19)

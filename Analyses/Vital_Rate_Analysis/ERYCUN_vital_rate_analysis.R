####################################################################################
# Purpose: Analyse cleaned Archbold vital rate data from Eryngium cuneifolium
# Author: Joshua Fowler
# Date: Oct 8, 2024
####################################################################################
##### Set up #####

library(renv) # track package versions
# renv::record("renv@1.0.7")
# renv::init()
# renv::snapshot()
# # ren
# renv::restore()
# # renv::update()

library(tidyverse)
library(readxl)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)

library(moments) # to calculate skew
library(patchwork) # to add ggplots together

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


Lkurtosis=function(x) log(kurtosis(x)); 


# non-parametric measures of skew and kurtosis
NPsd = function(x) {
  u = diff(quantile(x,c(0.25,0.75))/(qnorm(0.75)-qnorm(0.25)))
  as.numeric(u)
}	

NPskewness = function(x,p=0.05) {
  q = quantile(x,c(p,0.5,1-p))
  u = (q[3]+q[1]-2*q[2])/(q[3]-q[1]);
  return(as.numeric(u)); 
  
}	

NPkurtosis=function(x,p=0.05) {
  q = quantile(x,c(p,0.25,0.75,1-p))
  qN = qnorm(c(p,0.25,0.75,1-p))
  u = (q[4]-q[1])/(q[3]-q[2]);
  uN = (qN[4]-qN[1])/(qN[3]-qN[2]);
  return (as.numeric(u/uN-1)) 
}


####################################################################################
###### sourcing in the cleaned vital rate data #####################################
####################################################################################
source(paste0(getwd(),"/Analyses/data_processing.R"))
ERYCUN <- ERYCUN_covariates %>% 
  mutate(log_size.t = log(ros_diameter.t),
         log_size.t1 = log(ros_diameter.t1),
         time_since_fire = time_since_fire_actual) %>% 
  filter(Site_tag != "royce_ranch" & bald != "200")%>% 
  filter(year.t >1989)
  
####################################################################################
###### filtering out NAs for each vital rate #####################################
####################################################################################
# ERYCUN

ERYCUN_seedling_surv.df <- ERYCUN %>% 
  filter(ARCHBOLD_surv == 5) %>% 
  filter(!is.na(surv.t1) & !is.na(log_size.t)) 
ERYCUN_seedling_grow.df <- ERYCUN %>% 
  filter(!is.na(ros_diameter.t1) & !is.na(log_size.t)) 


ERYCUN_surv.df <- ERYCUN %>% 
  filter(ARCHBOLD_surv == 5) %>% 
  filter(!is.na(surv.t1) & !is.na(log_size.t)) 

ERYCUN_flw_status.df <- ERYCUN %>% 
  filter(ARCHBOLD_surv != 5) %>% 
  filter(!is.na(flw_status.t) & !is.na(log_size.t)) 
ERYCUN_flw_stem.df <- ERYCUN %>% 
  filter(ARCHBOLD_surv != 5) %>% 
  filter(flw_status.t == 1 & !is.na(flw_stem.t) & flw_stem.t !=0 & !is.na(log_size.t))
ERYCUN_flw_head.df <- ERYCUN %>% 
  filter(ARCHBOLD_surv != 5) %>% 
  filter(flw_status.t == 1 & !is.na(flw_head.t) & flw_head.t>0 & !is.na(log_size.t))
ERYCUN_growth.df <- ERYCUN %>% 
  filter(ARCHBOLD_surv != 5) %>% 
  filter(!is.na(ros_diameter.t1) & !is.na(log_size.t)) 


################################################################################
######### Setting up MCMC parameters ###########################################
################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 2500, 
  iter = 5000, 
  thin = 1, 
  chains = 3
)



####################################################################################
# Survival models ------------------------------------------------------------------
####################################################################################
### ERYCUN  ----
# starting first with a model without environmental covariates
erycun.survival <- brm(surv.t1~ 1 + log_size.t +I(log_size.t^2) + rel_elev + time_since_fire + (1|year.t1) + (1|bald), 
                       data = ERYCUN_surv.df,
                       family = "bernoulli",
                       prior = c(set_prior("normal(0,1)", class = "b")),
                       warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)

saveRDS(erycun.survival, "erycun.survival.rds")
erycun.survival <- readRDS("erycun.survival.rds")






####################################################################################
###### Growth models ###############################################################
####################################################################################
# started first with a model with gaussian distribution, but did a poor job fitting
# trying a gamma distribution to account for positive only values and more flexible mean-variance
# after fitting the gamma, the mean looks much better, but sd and skew are still not ideal. 
# trying to model size dependence in variance with distributional model formula option of brms and thin-plate splines 

# reading in the functions to access the sinhasinh family which allows separately modeling mean, variance, skew, and kurtosis.
source('Analyses/Vital_Rate_Analysis/stan_funs.R')
stanvars <- stanvar(scode = stan_funs(name = "sinhasinh"), block = "functions")



# formula <- bf(log_size.t1 ~ 1 + s(log_size.t, bs = "tp") + (1|year.t1) + (1|bald),
#               sigma ~ 1 + s(log_size.t, bs = "tp"),
#               eps ~  1 + s(log_size.t, bs = "tp"),
#               delta ~  1 + s(log_size.t, bs = "tp"))
# erycun.growth <- brm(formula,
#                      data = ERYCUN_growth.df,
#                      family = sinhasinh,
#                      stanvars = stanvars,
#                      prior = c(set_prior("normal(0,5)", class = "b"),
#                                set_prior("normal(0,5)", class = "Intercept")),
#                                # set_prior("cauchy(0,5)", class = "sd")),
#                      #set_prior("inv_gamma(0.4, 0.3)", class = "shape")),
#                      warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)


formula <- bf(log_size.t1 ~ 1 + s(log_size.t, bs = "tp") + rel_elev + time_since_fire + (1|year.t1) + (1|bald),
              sigma ~ 1 + s(log_size.t, bs = "tp"))
erycun.growth <- brm(formula,
                     data = ERYCUN_growth.df,
                     family = gaussian,
                     prior = c(set_prior("normal(0,5)", class = "b"),
                               set_prior("normal(0,5)", class = "Intercept")),
                     # set_prior("cauchy(0,5)", class = "sd")),
                     #set_prior("inv_gamma(0.4, 0.3)", class = "shape")),
                     warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)

saveRDS(erycun.growth, "erycun.growth.rds")
erycun.growth <- readRDS("erycun.growth.rds")


####################################################################################
###### flowering models ###############################################################
####################################################################################
### ERYCUN  ----
# starting first with a model for flowering probability
erycun.flw_status <- brm(flw_status.t~ 1 + s(log_size.t, bs = "tp") + rel_elev + time_since_fire+ (1|year.t1) + (1|bald), data = ERYCUN_flw_status.df,
                       family = "bernoulli",
                       prior = c(set_prior("normal(0,1)", class = "b")),
                       warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)

saveRDS(erycun.flw_status, "erycun.flw_status.rds")

source('Analyses/Vital_Rate_Analysis/stan_funs.R')

stanvars <- stanvar(scode = stan_funs(name = "PIG"), block = "functions")


erycun.flw_stem <- brm(flw_stem.t  ~ 1 + log_size.t + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = ERYCUN_flw_stem.df,
                         family = negbinomial,
                         # stanvars = stanvars,
                         prior = c(set_prior("normal(0,5)", class = "b"),
                                   set_prior( "normal(0,1)", class = "shape")),
                         warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)
saveRDS(erycun.flw_stem, "erycun.flw_stem.rds")

# trying to implement a Poisson-inverse gaussian distribution as a mixture model of inverse gaussian and poisson in BRMS
PIG_mix <- mixture(poisson, poisson)
prior <- c(
  prior(normal(0.5), Intercept, dpar = mu1),
  prior(normal(0,5), Intercept, dpar = mu2)
)
formula_mix <- bf(flw_stem.t|trunc(lb = 0) ~ 1, mu1 ~ log_size.t, mu2 ~ mu1)
erycun.flw_stem <- brm(formula_mix, data = ERYCUN_flw_stem.df,
                       family = PIG_mix,
                       prior = prior,
                       warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)




formula <- bf(flw_head.t | trunc(lb = 0)  ~ 1 + (1|year.t1),
              shape ~ 1 + (1|year.t1))

erycun.flw_head <- brm(formula, data = ERYCUN_flw_head.df,
                       family = negbinomial(),
                       prior = c(set_prior("normal(0,5)", class = "Intercept")),
                       warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)
saveRDS(erycun.flw_head, "erycun.flw_head.rds")





####################################################################################
###### Posterior predictive checks #################################################
####################################################################################

# Function for looking at binned size_t fits, particularly important for the growth kernel as this determines the transitions through the matrix model
# plots the mean, sd, skew and kertosis of the posteriors (grey) as well as the mean of the posteriors for each moment (black) and the data (red) for size bins
size_moments_ppc <- function(data,y_name,sim, n_bins, title = NA){
  require(tidyverse)
  require(patchwork)
  data$y_name <- data[[y_name]]
  bins <- data %>%
    ungroup() %>% 
    arrange(log_size.t) %>% 
    mutate(size_bin = cut_number(log_size.t, n_bins)) %>% 
    group_by(size_bin)  %>% 
    dplyr::summarize(mean_t1 = mean(y_name),
                     sd_t1 = sd(y_name),
                     skew_t1 = skewness(y_name),
                     kurt_t1 = kurtosis(y_name),
                     bin_mean = mean(log_size.t),
                     bin_n = n())
  sim_moments <- bind_cols(enframe(data$log_size.t), as_tibble(t(sim))) %>%
    rename(log_size.t = value) %>%
    arrange(log_size.t) %>%
    mutate(size_bin = cut_number(log_size.t, n_bins)) %>%
    pivot_longer(., cols = starts_with("V"), names_to = "post_draw", values_to = "sim") %>%
    group_by(size_bin, post_draw) %>%
    summarize( mean_sim = mean((sim)),
               sd_sim = sd((sim)),
               skew_sim = skewness((sim)),
               kurt_sim = kurtosis((sim)),
               bin_mean = mean(log_size.t),
               bin_n = n())
  sim_medians <- sim_moments %>%
    group_by(size_bin, bin_mean) %>%
    summarize(median_mean_sim = median(mean_sim),
              median_sd_sim = median(sd_sim),
              median_skew_sim = median(skew_sim),
              median_kurt_sim = median(kurt_sim))
  meanplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = mean_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_mean_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = mean_t1), shape = 1, color = "firebrick2") +
    theme_classic()
  sdplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = sd_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_sd_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = sd_t1), shape = 1, color = "firebrick2") + theme_classic()
  skewplot <-  ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = skew_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_skew_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = skew_t1), shape = 1, color = "firebrick2") + theme_classic()
  kurtplot <- ggplot(data = bins)+
    geom_point(data = sim_moments, aes(x = bin_mean, y = kurt_sim), color = "gray72") +
    geom_point(data = sim_medians, aes(x = bin_mean, y = median_kurt_sim),shape = 1, color = "black") +
    geom_point(aes(x = bin_mean, y = kurt_t1), shape = 1, color = "firebrick2") + theme_classic()
  size_ppc_plot <- meanplot+ sdplot+skewplot+ kurtplot+plot_annotation(title = title)
  return(size_ppc_plot)
}


#### survival ppc ####

y_sim <- posterior_predict(erycun.survival, ndraws = 500)

ppc_dens_overlay(y = ERYCUN_surv.df$surv.t1, yrep = y_sim)

mean_s_plot <-   ppc_stat(ERYCUN_surv.df$surv.t1, y_sim, stat = "mean")
sd_s_plot <- ppc_stat(ERYCUN_surv.df$surv.t1, y_sim, stat = "sd")
skew_s_plot <- ppc_stat(ERYCUN_surv.df$surv.t1, y_sim, stat = "NPskewness")
kurt_s_plot <- ppc_stat(ERYCUN_surv.df$surv.t1, y_sim, stat = "Lkurtosis")
surv_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Survival")
surv_moments



ERYCUN_surv.df <- ERYCUN_surv.df %>% 
  mutate(logsize_t = log_ros_diameter.t)


surv_size_ppc <- size_moments_ppc(data = ERYCUN_surv.df,
                                  y_name = "surv.t1",
                                  sim = y_sim, 
                                  n_bins = 10, 
                                  title = "Survival")
surv_size_ppc






#### growth ppc ####

y_sim <- posterior_predict(erycun.growth, ndraws = 500)

ppc_dens_overlay(y = ERYCUN_growth.df$log_size.t1, yrep = y_sim)

mean_s_plot <-   ppc_stat(ERYCUN_growth.df$log_size.t1, y_sim, stat = "mean")
sd_s_plot <- ppc_stat(ERYCUN_growth.df$log_size.t1, y_sim, stat = "NPsd")
skew_s_plot <- ppc_stat(ERYCUN_growth.df$log_size.t1, y_sim, stat = "NPskewness")
kurt_s_plot <- ppc_stat(ERYCUN_growth.df$log_size.t1, y_sim, stat = "NPkurtosis")
grow_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Growth")
grow_moments



grow_size_ppc <- size_moments_ppc(data = ERYCUN_growth.df,
                                  y_name = "log_size.t1",
                                  sim = y_sim, 
                                  n_bins = 10, 
                                  title = "Growth")
grow_size_ppc




#### flowering status ppc ####

y_sim <- posterior_predict(erycun.flw_status, ndraws = 500)

ppc_dens_overlay(y = ERYCUN_flw_status.df$flw_status.t, yrep = y_sim)

mean_s_plot <-   ppc_stat(ERYCUN_flw_status.df$flw_status.t, y_sim, stat = "mean")
sd_s_plot <- ppc_stat(ERYCUN_flw_status.df$flw_status.t, y_sim, stat = "sd")
skew_s_plot <- ppc_stat(ERYCUN_flw_status.df$flw_status.t, y_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(ERYCUN_flw_status.df$flw_status.t, y_sim, stat = "kurtosis")
flw_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Flowering Probability")
flw_moments



flw_status_size_ppc <- size_moments_ppc(data = ERYCUN_flw_status.df,
                                  y_name = "flw_status.t",
                                  sim = y_sim, 
                                  n_bins = 10, 
                                  title = "Growth")
flw_status_size_ppc



#### flowering stem ppc ####

y_sim <- posterior_predict(erycun.flw_stem, ndraws = 500)

ppc_dens_overlay(y = ERYCUN_flw_stem.df$flw_stem.t, yrep = y_sim) + xlim(0,50)

mean_s_plot <-   ppc_stat(ERYCUN_flw_stem.df$flw_stem.t, y_sim, stat = "mean")
sd_s_plot <- ppc_stat(ERYCUN_flw_stem.df$flw_stem.t, y_sim, stat = "sd")
skew_s_plot <- ppc_stat(ERYCUN_flw_stem.df$flw_stem.t, y_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(ERYCUN_flw_stem.df$flw_stem.t, y_sim, stat = "Lkurtosis")
flw_stem_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Flowering Stems")
flw_stem_moments



flw_stem_size_ppc <- size_moments_ppc(data = ERYCUN_flw_stem.df,
                                        y_name = "flw_stem.t",
                                        sim = y_sim, 
                                        n_bins = 10, 
                                        title = "Growth")
flw_stem_size_ppc




#### flowering head ppc ####

y_sim <- posterior_predict(erycun.flw_head, ndraws = 500)

ppc_dens_overlay(y = ERYCUN_flw_head.df$flw_head.t, yrep = y_sim) + xlim(0,500)

mean_s_plot <-   ppc_stat(ERYCUN_flw_head.df$flw_head.t, y_sim, stat = "mean")
sd_s_plot <- ppc_stat(ERYCUN_flw_head.df$flw_head.t, y_sim, stat = "sd")
skew_s_plot <- ppc_stat(ERYCUN_flw_head.df$flw_head.t, y_sim, stat = "skewness")
kurt_s_plot <- ppc_stat(ERYCUN_flw_head.df$flw_head.t, y_sim, stat = "Lkurtosis")
flw_head_moments <- mean_s_plot+sd_s_plot+skew_s_plot+kurt_s_plot + plot_annotation(title = "Flowering heads")
flw_head_moments




####################################################################################
###### Visualizing variation across the balds ######################################
####################################################################################

preddata <- expand.grid(log_size.t = mean(ERYCUN_surv.df$log_size.t), year.t1 = NA, bald = unique(ERYCUN_surv.df$bald))

surv_pred <- fitted(erycun.survival, newdata = preddata, probs = c(0.025,0.50,0.975))

surv_pred_df <- bind_cols(preddata, surv_pred)

surv_binned <- ERYCUN_surv.df %>% 
  group_by(bald) %>% 
  summarize(mean_surv = mean(surv.t1, na.rm = T),
            sample_size = n())
ggplot(surv_pred_df)+
  # geom_point(data = ERYCUN_surv.df, aes(x = bald, y = surv.t1, color = bald))+
  geom_point(data = surv_binned, aes(x = bald, y = mean_surv, size = sample_size, color = bald))+
  geom_point(aes(x = bald, y = Estimate, group = bald), color = "black")+
  geom_linerange(aes(x = bald, ymin = Q2.5, ymax = Q97.5))+
  lims(y = c(0,1)) + 
  theme_classic() + labs(size = "Sample Size", color = "Bald")

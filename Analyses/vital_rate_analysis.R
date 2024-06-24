####################################################################################
# Purpose: Analyse cleaned Archbold vital rate data from eleven (10?) species
# Author: Joshua Fowler
# Date: Jun 24, 2024
####################################################################################
##### Set up #####

library(renv) # track package versions
# renv::init()
# renv::snapshot()
renv::restore()


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


####################################################################################
###### sourcing in the cleaned vital rate data #####################################
####################################################################################
source("data.processing.R")
ERYCUN <- ERYCUN_covariates


####################################################################################
###### filtering out NAs for each vital rate #####################################
####################################################################################

ERYCUN_surv.df <- ERYCUN %>% 
  filter(!is.na(surv.t1) & !is.na(ros_diameter.t)) %>% 
  mutate(log_ros_diameter.t = log(ros_diameter.t)) %>% 
  filter(Site_tag != "royce_ranch" & bald != "200")




################################################################################
######### Setting up MCMC parameters ###########################################
################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  warmup = 1000, 
  iter = 2500, 
  thin = 1, 
  chains = 3
)



####################################################################################
###### Survival models #############################################################
####################################################################################

# starting first with a model without environmental covariates
erycun.survival <- brm(surv.t1~ 1 + log_ros_diameter.t + I(log_ros_diameter.t^2) + fire_frequency_actual + (1|bald), data = ERYCUN_surv.df,
              family = "bernoulli",
              prior = c(set_prior("normal(0,1)", class = "b")),
              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)




# Making prediction dataframe

prediction_df <- expand.grid(log_ros_diameter.t = c(min(ERYCUN_surv.df$log_ros_diameter.t, na.rm = T), median(ERYCUN_surv.df$log_ros_diameter.t, na.rm = T), max(ERYCUN_surv.df$log_ros_diameter.t, na.rm  = T)),
                             fire_frequency_actual = seq(from = min(ERYCUN_surv.df$fire_frequency_actual, na.rm = T), to = max(ERYCUN_surv.df$fire_frequency_actual, na.rm = T), length.out = 25),
                             bald = NA)

preds <- fitted(erycun.survival, newdata = prediction_df)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


# Making a binned survival df for visualization
erycun.survival.binned <- ERYCUN_surv.df %>% 
  ungroup() %>% 
         mutate(size_bin = cut_number(log_ros_diameter.t, 10)) %>% 
  group_by(size_bin) %>% 
  summarize(mean_size = mean(log_ros_diameter.t, na.rm = T),
            mean_surv = mean(surv.t1, na.rm = T),
            samplesize = n())




# now we can plot the model predictions

ggplot(data = prediction_df)+
  geom_point(data = ERYCUN_surv.df, aes( x= fire_frequency_actual, y = surv.t1), color = "grey55", alpha = .2)+
  # geom_point(data = erycun.survival.binned, aes( x= mean_size, y = mean_surv, size = samplesize), color = "skyblue3", alpha = .5)+
  geom_ribbon(aes(x = fire_frequency_actual, ymax = upr, ymin = lwr,fill = as.factor(log_ros_diameter.t),  group = log_ros_diameter.t), alpha = .4)+
  geom_line(aes(x = fire_frequency_actual, y = fit, color = as.factor(log_ros_diameter.t), group = log_ros_diameter.t), size = 1) +
  # geom_linerange(aes(x = live_sterile, ymin = lwr, ymax = upr))+ 
  labs(y = "Survival (%)") + theme_minimal()





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
    arrange(logsize_t) %>% 
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>% 
    group_by(size_bin)  %>% 
    dplyr::summarize(mean_t1 = mean(y_name),
                     sd_t1 = sd(y_name),
                     skew_t1 = skewness(y_name),
                     kurt_t1 = Lkurtosis(y_name),
                     bin_mean = mean(logsize_t),
                     bin_n = n())
  sim_moments <- bind_cols(enframe(data$logsize_t), as_tibble(t(sim))) %>%
    rename(logsize_t = value) %>%
    arrange(logsize_t) %>%
    mutate(size_bin = cut_number(logsize_t, n_bins)) %>%
    pivot_longer(., cols = starts_with("V"), names_to = "post_draw", values_to = "sim") %>%
    group_by(size_bin, post_draw) %>%
    summarize( mean_sim = mean((sim)),
               sd_sim = sd((sim)),
               skew_sim = skewness((sim)),
               kurt_sim = Lkurtosis((sim)),
               bin_mean = mean(logsize_t),
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
skew_s_plot <- ppc_stat(ERYCUN_surv.df$surv.t1, y_sim, stat = "skewness")
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

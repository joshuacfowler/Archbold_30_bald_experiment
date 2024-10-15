####################################################################################
# Purpose: Analyse cleaned Archbold vital rate data from eleven (10?) species
# Author: Joshua Fowler
# Date: Oct 8, 2024
####################################################################################
##### Set up #####

library(tidyverse)
library(readxl)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(tidybayes)

library(moments) # to calculate skew
library(patchwork) # to add ggplots together

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


####################################################################################
###### sourcing in the vital rate model objects ####################################
####################################################################################

erycun.survival <- readRDS("erycun.survival.rds")
erycun.growth <- readRDS("erycun.growth.rds")
erycun.flw_status <- readRDS("erycun.flw_status.rds")
erycun.flw_stem <- readRDS("erycun.flw_stem.rds")
erycun.flw_head <- readRDS("erycun.flw_head.rds")


####################################################################################
###### sourcing in the MPM functions script ####################################
####################################################################################
source("Analyses/Population_Model_Analysis/MPM_functions.R")

# calculating max size
source(paste0(getwd(),"/Analyses/data_processing.R"))
size_bounds_df <- ERYCUN%>% 
  group_by() %>% 
  filter(!is.na(log_size.t1)) %>%
  summarise(actual_max_size = max(log_size.t1),
            max_size = quantile(log_size.t1,probs=0.975),
            max_size_99 = quantile(log_size.t1,probs=0.99),
            actual_min_size = min(log_size.t1),
            min_size = quantile(log_size.t1,probs=0.025),
            min_size_99 = quantile(log_size.t1,probs=0.01)
            )

# setting up models and params to feed in to MPM_functions
models <- make_mods(grow = erycun.growth, surv = erycun.survival, flw = erycun.flw_status, fert = erycun.flw_stem, head = 183, recruit = .001)
params <- make_params(bald.rfx = 1, year = 2000, size_bounds = size_bounds_df)



####################################################################################
###### Calculating population summaries ####################################
####################################################################################
# we can calculate lambda, but we might also consider later looking at effects of microbes on quantities like the stable stage distribution etc.
ndraws <- 10
nbalds <- 10
lambda <- array(NA, dim = c(ndraws, nbalds))
for(i in 1:ndraws){
  for(b in 1:nbalds){
  # lambda[i] <- Re(eigen(bigmatrix(params = params, models = models, matdim = 100, extension = 0)$MPMmat)$values[1])
  lambda[i,b] <- lambda(bigmatrix(params = make_params(bald.rfx = b, year = 2000, size_bounds = size_bounds_df), 
                                models = models, matdim = 100, extension = 0)$MPMmat)
  }
}
lambda


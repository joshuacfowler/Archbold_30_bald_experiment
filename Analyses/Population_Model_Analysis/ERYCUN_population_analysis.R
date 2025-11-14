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

library(parallel) # used for parallel computing of matrix operations


library(popbio)

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
# erycun.survival <- readRDS("/Volumes/Macintosh HD - Data/Users/joshuacfowler/Documents/R_projects/Archbold_30_bald_experiment/erycun.survival.rds")
erycun.survival <- readRDS("erycun.survival.rds")
erycun.growth <- readRDS("erycun.growth.rds")
erycun.flw_status <- readRDS("erycun.flw_status.rds")
erycun.flw_stem <- readRDS("erycun.flw_stem.rds")
# erycun.flw_head <- readRDS("erycun.flw_head.rds")


surv_par <- as.data.frame(erycun.survival)
grow_par <- as.data.frame(erycun.growth)
flw_par <- as.data.frame(erycun.flw_status)
flw_stem_par <- as.data.frame(erycun.flw_stem)
# erycun.flw_head <- as.data.frame(erycun.flw_head)

# reading in the predicted microbial effects from the 30 bald experiment
germ_30bald_predictions <- read_csv("germ_30bald_predictions.csv") %>% 
  select(bald, spp_code, Iter, rel_diff) %>% 
  # summarize(rel_diff = mean(rel_diff))%>% 
  drop_na() %>% 
  filter(spp_code == "ERYCUN")


grow_30bald_predictions <- read_csv("grow_30bald_predictions.csv") %>% 
  select(bald, spp_code, Iter, rel_diff) %>% 
  # summarize(rel_diff = mean(rel_diff)) %>% 
  drop_na() %>% 
  filter(spp_code == "ERYCUN")

flw_30bald_predictions <- read_csv("flw_30bald_predictions.csv") %>% 
  select(bald, spp_code, Iter, rel_diff) %>% 
  # summarize(rel_diff = mean(rel_diff))%>% 
  drop_na() %>% 
  filter(spp_code == "ERYCUN")
####################################################################################
###### sourcing in the MPM functions script ####################################
####################################################################################
source("Analyses/Population_Model_Analysis/size_IPM_functions.R")

# calculating max size
# source(paste0(getwd(),"/Analyses/data_processing.R"))
ERYCUN_covariates <- read_csv("/Users/joshuacfowler/Dropbox/UofMiami/Demographic Data/cleaned_data/ERYCUN_covariates.csv")

ERYCUN <- ERYCUN_covariates %>% 
  mutate(log_size.t = log(ros_diameter.t),
         log_size.t1 = log(ros_diameter.t1),
         time_since_fire = time_since_fire_actual) %>% 
  filter(Site_tag != "royce_ranch" & bald != "200")%>% 
  filter(year.t >1989)



ERYCUN_surv.df <- ERYCUN %>% 
  filter(!is.na(surv.t1) & !is.na(log_size.t)) 
ERYCUN_flw_status.df <- ERYCUN %>% 
  filter(!is.na(flw_status.t) & !is.na(log_size.t)) 
ERYCUN_flw_stem.df <- ERYCUN %>% 
  filter(flw_status.t == 1 & !is.na(flw_stem.t) & flw_stem.t !=0 & !is.na(log_size.t))
ERYCUN_flw_head.df <- ERYCUN %>% 
  filter(flw_status.t == 1 & !is.na(flw_head.t) & flw_head.t>0 & !is.na(log_size.t))
ERYCUN_growth.df <- ERYCUN %>% 
  filter(!is.na(ros_diameter.t1) & !is.na(log_size.t)) 

size_bounds_df <- ERYCUN%>% 
  group_by() %>% 
  filter(!is.na(log_size.t1)) %>%
  summarise(max_size = max(log_size.t1),
            max_size_97 = quantile(log_size.t1,probs=0.975),
            max_size_99 = quantile(log_size.t1,probs=0.99),
            min_size = min(log_size.t1),
            min_size_97 = quantile(log_size.t1,probs=0.025),
            min_size_99 = quantile(log_size.t1,probs=0.01)
            )

# setting up models and params to feed in to MPM_functions


####################################################################################
###### setting up data across balds for predictions ################################
####################################################################################

# setting up the covariates for each bald for prediction

fire_to2022 <- read_xlsx(path = "~/Dropbox/UofMiami/Balds2009_FireIntensityArea_Through2022.xlsx", sheet = "Balds2009_FireIntensityArea", guess_max = 1048576)# guess_max makes the function look deeper in the columns to assign type

fire_summary <- fire_to2022 %>% 
  group_by(Bald_U, Bald_) %>% 
  filter(INTENSITY != 0) %>% 
  summarize(last_fire = max(as.numeric(Year), na.rm = TRUE),
            time_since_fire = 2023-max(as.numeric(Year), na.rm = TRUE),
            fire_list = toString(unique(as.numeric(Year))),
            fire_frequency = length(unique(as.numeric(Year)))) %>% 
  ungroup() %>% 
  add_row(Bald_U = "1S", Bald_ = 1, last_fire = NA, time_since_fire = NA, fire_list = NA, fire_frequency = 0)

# Now getting the relative elevation from the old fire history file

elev.df <- read_xlsx(path = "~/Dropbox/UofMiami/Experiment Set up/firehistory_thru2018.xlsx", sheet = "Rx_Freq", guess_max = 1048576) %>% # guess_max makes the function look deeper in the columns to assign type
  rename(rel_elev = rel.eve) %>% 
  mutate(
    bald = case_when(bald == "01S" ~ "1S",
                     bald == "01N" ~ "1N",
                     bald == "02" ~ "2",
                     bald == "05E" ~ "5E",
                     bald == "07N" ~ "7N",
                     bald == "35N" ~ "35",
                     bald == "65E" ~ "65",
                     TRUE ~ bald)) %>% 
  select(bald, name, rel_elev)


bald_covariates <- fire_summary %>% 
  left_join(elev.df, by = join_by(Bald_U == bald)) %>% 
  filter(!is.na(rel_elev), !is.na(time_since_fire)) %>% 
  rename(bald = Bald_U)

preddata_1 <- bald_covariates %>% 
  select(bald, time_since_fire, rel_elev) %>% 
  mutate(total_seeds = 1, soil_source = bald) 

treatment_df <- expand_grid(spp_code = "ERYCUN", live_sterile = c("live", "sterile"), soil_source = unique(bald_covariates$bald))

preddata <- preddata_1 %>% 
  left_join(treatment_df) %>% 
  filter(live_sterile == "live") %>% select(-live_sterile)
# params <- make_params(bald.rfx = 10, year = NA, preddata = preddata, size_bounds = size_bounds_df)



####################################################################################
###### Calculating population summaries ####################################
####################################################################################
# Taking seed dynamics from Hindle et al. "The two fertility scenarios which produced dynamics best fitting to the observed population dynamics were low first year germination (c_f=0), low germination from the seed bank (c_b=0.005), and low seed mortality (d=0.3) and low first year germination (c_f=0), high germination from the seed bank (c_b=0.04), and relatively high seed mortality (d=0.7; Fig. A5). "
models <- make_mods(grow = erycun.growth, surv = erycun.survival, flw = erycun.flw_status, fert = erycun.flw_stem, 
                    seeds_per_stem = 183, seed_mortality = .3, seed_germ1 = 0, seed_germ2 = .005,
                    seedling_surv = erycun.survival, seedling_size = erycun.growth)
# params <- make_params(post_draws = post_draws,
#                       iter = i, 
#                       bald.rfx = T, bald = 12, year.rfx = F, 
#                       surv_par = surv_par, grow_par = grow_par,
#                       flw_par = flw_par, fert_par = flw_stem_par,
#                       microbe = 0, 
#                       preddata = preddata,
#                       germ_microbe = germ_30bald_predictions,
#                       grow_microbe = grow_30bald_predictions,
#                       flw_microbe = flw_30bald_predictions,
#                       size_bounds = size_bounds_df)

# gxy(0,0,models, params)
#
# plot(fx(1:10, models, params))
# 
# 


# we can calculate lambda, but we might also consider later looking at effects of microbes on quantities like the stable stage distribution etc.
ndraws <- 3
nbalds <- length(unique(preddata$bald)[1:3])
nmicrobe <- 2

balds <- unique(preddata$bald)[1:3]
post_draws <- sample.int(7500,size=ndraws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.

microbe <- c(0,1) # 0 is alive, and 1 is sterile becuase we start with the microbes in the model, but then turn off the microbes
lambda <- array(NA, dim = c(ndraws, nbalds, nmicrobe))
params_list <- matrix_list <- matrix_balds <- matrix_iter <- list()
for(i in 1:ndraws){
  for(b in 1:nbalds){
    for(m in 1:nmicrobe){
      params_list[[paste(paste0("iter",c(1:ndraws)[i]),balds[b],c("live", "sterile")[m], sep = "_")]] <- make_params(bald.rfx = T, bald = balds[b], year.rfx = F,
                                                                 post_draws = post_draws,
                                                                 iter = i,
                                                                 preddata = preddata, 
                                                                 microbe = microbe[m],
                                                                 surv_par = surv_par,
                                                                 grow_par = grow_par, 
                                                                 flw_par = flw_par,
                                                                 fert_par = flw_stem_par,
                                                                 germ_microbe = germ_30bald_predictions,
                                                                 grow_microbe = grow_30bald_predictions,
                                                                 flw_microbe = flw_30bald_predictions,
                                                                 size_bounds = size_bounds_df)
    }
    # matrix_balds[[b]] <- matrix_list 
  }
  # matrix_iter[[i]] <- matrix_balds
  print(paste("iteration", i))
}



return_MPM <- function(params,...) {
  bigmatrix(params = params,...)$MPMmat
  }


# system.time(matrix_list <- lapply(X = params_list, FUN = return_MPM, models = models, matdim = 25, extension = 2)) # note that lapply is vectorized, and takes about the same time as using for loops, at least for a few iterations
system.time(matrix_list <- mclapply(X = params_list, FUN = return_MPM, models = models, matdim = 25, extension = 2)) # note that mclapply takes advantage of parallel computation and takes about half as much time as using lapply

lambda_list <- lapply(X = matrix_list, FUN = popbio::lambda) # for already fast operations, the parallelization isn't faster, but it's trivial

lambda_df <- enframe(lambda_list) %>% 
  unnest(cols = c("name", "value")) %>% 
  separate(name, into = c("Iter", "Bald", "Microbe")) %>% 
  rename(lambda=value)
# write_csv(lambda_df, "prelim_IPM_lambdas.csv")



# saveRDS(lambda, "lambda_microbe.Rds")
# lambda <- readRDS("lambda_microbe.Rds")
# 
# # lambda
# 
# # lambda_cube <- cubelyr::as.tbl_cube(lambda)
# lambda_df <- as_tibble(lambda, rownames = "Iter") %>% 
#   pivot_longer(cols = -Iter, names_to = "name", values_to = "lambda") %>% 
#   separate(name, into = c("Bald", "Microbe")) 
# write_csv(lambda_df, "prelim_IPM_lambdas.csv")

lambda_summary <- lambda_df %>% 
  filter(!is.na(lambda)) %>% 
  group_by(Bald, Microbe) %>% 
  summarize(lambda_mean = mean(lambda),
            lambda_97.5 = quantile(lambda, 0.975),
            lambda_02.5 = quantile(lambda, 0.025))

ggplot()+
  geom_hline(aes(yintercept = 1), linetype = "dashed")+
  geom_jitter(data = lambda_df, aes(x = Bald, y = lambda), width = .05, alpha = .4)+
  geom_point(data = lambda_summary, aes(x = Bald, y = lambda_mean, color = Microbe, group = "microbe"), size = 4)+
  geom_linerange(data = lambda_summary, aes(x = Bald, ymin = lambda_02.5, ymax = lambda_97.5, color = Microbe), lwd = 1)+
  # ylim(0,10)+
  theme_classic()


lambda_effect_percent <- lambda_summary %>% 
  select(-lambda_97.5, -lambda_02.5) %>% 
  pivot_wider(names_from = Microbe, values_from = lambda_mean) %>% 
  dplyr::summarize(percent_change = ((live-sterile)/sterile))

lambda_percent_summary <- lambda_effect_percent %>% 
  summarize(mean_change = mean(percent_change),
            max_change = max(percent_change),
            min_change = min(percent_change))

# plot(y,sx(y,models,params),xlab="Size",type="l",
#      ylab="Survival Probability",lwd=12)
# points(y,apply(bigmatrix(params = params, models = models, matdim = 100, extension = 2)$Tmat,2,sum),col="red",lwd=3,cex=.1,pch=19)

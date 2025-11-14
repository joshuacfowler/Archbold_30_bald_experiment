####################################################################################
# Purpose: Construct population model for PARCHA
# Author: Joshua Fowler
# Date: Nov 12, 2025
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
parcha.adult.survival <- readRDS("parcha.adult.survival.rds")
parcha.seedling.survival <- readRDS("parcha.seedling.survival.rds")
parcha.adult.flowering <- readRDS("parcha.adult.flowering.rds")
parcha.seedling.flowering <- readRDS("parcha.seedling.flowering.rds")
parcha.recruitment <- readRDS("parcha.recruitment.rds")


adult.surv_par <- as.data.frame(parcha.adult.survival)
seedling.surv_par <- as.data.frame(parcha.seedling.survival)
adult.flw_par <- as.data.frame(parcha.adult.flowering)
seedling.flw_par <- as.data.frame(parcha.seedling.flowering)
recruit_par <- as.data.frame(parcha.recruitment)


# reading in the predicted microbial effects from the 30 bald experiment
germ_30bald_predictions <- read_csv("germ_30bald_predictions.csv") %>% 
  select(bald, spp_code, Iter, rel_diff) %>% 
  # summarize(rel_diff = mean(rel_diff))%>% 
  drop_na() %>% 
  filter(spp_code == "PARCHA")


grow_30bald_predictions <- read_csv("grow_30bald_predictions.csv") %>% 
  select(bald, spp_code, Iter, rel_diff) %>% 
  # summarize(rel_diff = mean(rel_diff)) %>% 
  drop_na() %>% 
  filter(spp_code == "PARCHA")

flw_30bald_predictions <- read_csv("flw_30bald_predictions.csv") %>% 
  select(bald, spp_code, Iter, rel_diff) %>% 
  # summarize(rel_diff = mean(rel_diff))%>% 
  drop_na() %>% 
  filter(spp_code == "PARCHA")




####################################################################################
###### sourcing in the MPM functions script ####################################
####################################################################################
source("Analyses/Population_Model_Analysis/age_IPM_functions.R")

# calculating max age
# source(paste0(getwd(),"/Analyses/data_processing.R"))
PARCHA_covariates <- read_csv("/Users/joshuacfowler/Dropbox/UofMiami/Demographic Data/cleaned_data/PARCHA_covariates.csv")

PARCHA <- PARCHA_covariates %>% 
  filter(!is.na(surv.t1)) %>% 
  group_by(plant_id, Bald_U) %>%  
  mutate(census_end = case_when(surv.t1 == 0 ~ census_date.t1),
         census_start = min(census_date.t1),
         census_number = row_number()) %>% 
  fill(census_end, .direction = "updown") %>% 
  fill(census_start, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(age.t = case_when(is.na(age_in_months.t) ~ NA, TRUE ~ (census_number-1)),
         log_age.t = case_when(is.na(age_in_months.t) ~ NA, TRUE ~ log(census_number-1)),
         time_since_fire = time_since_fire_actual) %>% 
  mutate(season = case_when(month.t1 %in% c(3,4,5) ~ "spring",
                            month.t1 %in% c(6,7,8) ~ "summer",
                            month.t1 %in% c(9,10,11) ~ "fall",
                            month.t1 %in% c(12,1,2) ~ "winter")) %>% 
  filter(!is.na(log_age.t), !is.na(time_since_fire))


size_bounds_df <- PARCHA%>% 
  group_by() %>% 
  filter(!is.na(age.t)) %>%
  summarise(max_age = max(age.t),
            max_age_97 = quantile(age.t,probs=0.975),
            max_age_99 = quantile(age.t,probs=0.99),
            min_age = min(log_age.t),
            min_age_97 = quantile(age.t,probs=0.025),
            min_age_99 = quantile(age.t,probs=0.01)
  )


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

treatment_df <- expand_grid(spp_code = "PARCHA", live_sterile = c("live", "sterile"), soil_source = unique(bald_covariates$bald))

preddata <- preddata_1 %>% 
  left_join(treatment_df) %>% 
  filter(live_sterile == "live") %>% select(-live_sterile)
# params <- make_params(bald.rfx = 10, year = NA, preddata = preddata, size_bounds = size_bounds_df)

####################################################################################
###### Calculating population summaries ####################################
####################################################################################
# Taking seed dynamics from Hindle et al. "The two fertility scenarios which produced dynamics best fitting to the observed population dynamics were low first year germination (c_f=0), low germination from the seed bank (c_b=0.005), and low seed mortality (d=0.3) and low first year germination (c_f=0), high germination from the seed bank (c_b=0.04), and relatively high seed mortality (d=0.7; Fig. A5). "
models <- make_mods(seedling_surv = parcha.seedling.survival, 
                    surv = parcha.adult.survival,
                    seedling_flw = parcha.seedling.flowering, 
                    flw = parcha.adult.flowering,
                    recruitment = parcha.recruitment)
params <- make_params(post_draws = post_draws,
                      iter = i,
                      bald.rfx = T, bald = 12, year.rfx = F,
                      surv_par = adult.surv_par, sdlg_surv_par = seedling.surv_par, 
                      flw_par = adult.flw_par, sdlg_flw_par = seedling.flw_par, 
                      season = "winter",
                      recruit_par = recruit_par,
                      seed_mortality = 0.001,
                      microbe = 0,
                      preddata = preddata,
                      germ_microbe = germ_30bald_predictions,
                      grow_microbe = grow_30bald_predictions,
                      flw_microbe = flw_30bald_predictions,
                      size_bounds = size_bounds_df)




# we can calculate lambda, but we might also consider later looking at effects of microbes on quantities like the stable stage distribution etc.
ndraws <- 3
nbalds <- length(unique(preddata$bald)[1:3])
nmicrobe <- 2

balds <- unique(preddata$bald)[1:3]
post_draws <- sample.int(7500,size=ndraws) # The models except for seedling growth have 7500 iterations. That one has more (15000 iterations) to help it converge.

microbe <- c(0,1) # 0 is alive, and 1 is sterile becuase we start with the microbes in the model, but then turn off the microbes
lambda <- array(NA, dim = c(ndraws, nbalds, nmicrobe))
params_list <- matrix_list <- matrix_balds <- matrix_iter <- list()



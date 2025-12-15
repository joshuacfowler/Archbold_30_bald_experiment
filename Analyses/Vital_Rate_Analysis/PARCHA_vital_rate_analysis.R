####################################################################################
# Purpose: Analyse cleaned Archbold vital rate data from Paronychia chartaceae
# Author: Joshua Fowler
# Date: Oct 9, 2025
####################################################################################
##### Set up #####

# library(renv) # track package versions
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

####################################################################################
###### sourcing in the cleaned vital rate data #####################################
####################################################################################
## Filepath for stored demographic data
filepath <- c("/Users/joshuacfowler/Dropbox/UofMiami/Demographic Data")

PARCHA_covariates <- read_csv(paste0(filepath,"/cleaned_data", "/PARCHA_covariates.csv"))

# reading in greenhouse data which has some flower counts
plants <- read_xlsx(path = "~/Dropbox/UofMiami/Archbold_30baldexperiment.xlsx", sheet = "Census", guess_max = 1048576) # guess_max makes the function look deeper in the columns to assign type

flw_census <- plants %>% 
  mutate(across(everything(), as.character)) %>% 
  dplyr::select(c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds,reseeding, date_potted), 
                contains("Flw"), -contains("notes")) %>% 
  pivot_longer(cols = c(contains("Flw"))) %>% 
  separate(name, c("measurement", "census_number", "name")) %>% 
  filter(spp_code == "PARCHA") %>% 
  pivot_wider(id_cols = c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds, reseeding, date_potted, census_number), names_from = c(measurement,name), values_from = value, names_repair = "minimal") %>% 
  filter(FlwCount_measurement>1)


# getting size of these plants
size_census <- plants %>% 
  mutate(across(everything(), as.character)) %>% 
  dplyr::select(c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds,reseeding, date_potted), 
                contains("size"), -contains("notes")) %>% 
  pivot_longer(cols = c(contains("size"))) %>% 
  separate(name, c("measurement", "census_number", "name")) %>% 
  filter(spp_code == "PARCHA") %>% 
  pivot_wider(id_cols = c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds, reseeding, date_potted, census_number), names_from = c(measurement,name), values_from = value, names_repair = "minimal") %>% 
  filter(pot_id %in% flw_census$pot_id & Size_date == unique(flw_census$FlwCount_date))
  

flw_counts <- flw_census %>% 
  left_join(size_census, by = c("pot_id", "x_id", "y_id", "rep_id", "table_id", "soil_source", "live_sterile", "spp_code", "number_seeds", "reseeding", "date_potted")) %>% 
  select(pot_id, soil_source, live_sterile, spp_code, FlwCount_date, FlwCount_measurement, HeightSize_measurement, DiameterSize_measurement) %>% 
  mutate(across( c(FlwCount_measurement, HeightSize_measurement, DiameterSize_measurement),  ~ as.numeric(.))) %>% 
  mutate(area = pi*((DiameterSize_measurement/2)^2),
         flw_count = FlwCount_measurement)



####################################################################################
###### filtering out NAs for each vital rate #####################################
####################################################################################
# we can filter out some of the sites, but need to get more of the meta-data about that
# PARCHA

PARCHA_surv.df <- PARCHA_covariates %>% 
  filter(!is.na(surv.t1)) %>% 
  group_by(plant_id, Bald_U) %>%  
  mutate(census_end = case_when(surv.t1 == 0 ~ census_date.t1),
         census_start = min(census_date.t1),
         census_number = row_number()) %>% 
  fill(census_end, .direction = "updown") %>% 
  fill(census_start, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(log_age.t = case_when(is.na(age_in_months.t) ~ NA, TRUE ~ log(census_number-1)),
         time_since_fire = time_since_fire_actual) %>% 
  mutate(season = case_when(month.t1 %in% c(3,4,5) ~ "spring",
                            month.t1 %in% c(6,7,8) ~ "summer",
                            month.t1 %in% c(9,10,11) ~ "fall",
                            month.t1 %in% c(12,1,2) ~ "winter")) %>% 
  filter(!is.na(log_age.t), !is.na(time_since_fire))

PARCHA_surv_adult<- PARCHA_surv.df %>% 
  filter(log_age.t != 0)

PARCHA_surv_seedling <- PARCHA_surv.df %>% 
  filter(log_age.t == 0)



# %>% 
#   group_by(age.t) %>% 
#   summarize(length(unique(plant_id)))
 
PARCHA_flw.df <- PARCHA_covariates %>% 
  filter(!is.na(surv.t1)) %>%
  group_by(plant_id, Bald_U) %>%  
  mutate(census_end = case_when(surv.t1 == 0 ~ census_date.t1),
         census_start = min(census_date.t1),
         census_number = row_number()) %>% 
  fill(census_end, .direction = "updown") %>% 
  fill(census_start, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(log_age.t = case_when(is.na(age_in_months.t) ~ NA, TRUE ~ log(census_number-1)),
         time_since_fire = time_since_fire_actual) %>% 
  mutate(season = case_when(month.t1 %in% c(3,4,5) ~ "spring",
                            month.t1 %in% c(6,7,8) ~ "summer",
                            month.t1 %in% c(9,10,11) ~ "fall",
                            month.t1 %in% c(12,1,2) ~ "winter")) %>% 
  filter(!is.na(flw.t1), !is.na(log_age.t), !is.na(time_since_fire))

PARCHA_flw_adult<- PARCHA_flw.df %>% 
  filter(log_age.t != 0)

PARCHA_flw_seedling <- PARCHA_flw.df %>% 
  filter(log_age.t == 0)


PARCHA_recruit.df <- PARCHA_covariates %>% 
  mutate(time_since_fire = time_since_fire_actual) %>% 
  # filter(birth_month == month.t1  & birth_year == year.t1) %>% 
  group_by(time_since_fire, rel_elev, bald, plot, year.t1, month.t1, .drop = FALSE) %>%
  dplyr::summarize(recruits = sum(birth_month == month.t1  & birth_year == year.t1, na.rm = T),
                   reproductives = sum(flw.t == 1, na.rm = T),
                   reproductives_lag = sum(sum(flw.t == 1, na.rm = T), sum(lag(flw.t, 1) == 1, na.rm = T), sum(lag(flw.t, 1) == 2, na.rm = T), sum(lag(flw.t, 1) == 2, na.rm = T)),
                   ) %>% 
  filter(!is.na(recruits)) %>% 
  mutate(season = case_when(month.t1 %in% c(3,4,5) ~ "spring",
                            month.t1 %in% c(6,7,8) ~ "summer",
                            month.t1 %in% c(9,10,11) ~ "fall",
                            month.t1 %in% c(12,1,2) ~ "winter")) %>% 
  filter(!is.na(time_since_fire))


#   
# ggplot(PARCHA_recruit.df)+
#   geom_point(aes(y = recruits, x = month.t1))
#   facet_wrap(~bald, scales = "free")


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
# Adult Survival models ------------------------------------------------------------------
####################################################################################
### PARCHA  ----
# starting first with a model without environmental covariates

parcha.adult.survival <- brm(surv.t1~ 0 + season + log_age.t + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_surv_adult,
                       #surv.t1~ 1 + log_age.t, data = PARCHA_surv.df,
                       family = "bernoulli",
                       prior = c(set_prior("normal(0,1)", class = "b"),
                                 #set_prior("normal(0,1)", class = "Intercept"),
                                 set_prior("normal(0,1)", class = "sd")),
                                 #set_prior("normal(0,1)", class = "sds")),
                       warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)
# The model with spline wins model selection against one without
# parcha.adult.survival.spline <- brm(surv.t1~ 1 +  s(log_age.t, bs = "tp", k=6) + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_surv_adult,
#                              #surv.t1~ 1 + log_age.t, data = PARCHA_surv.df,
#                              family = "bernoulli",
#                              prior = c(set_prior("normal(0,1)", class = "b"),
#                                        set_prior("normal(0,1)", class = "Intercept"),
#                                        set_prior("normal(0,1)", class = "sd"),
#                                        set_prior("normal(0,1)", class = "sds")),
#                              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)
# # linear model wins here
# loo(parcha.adult.survival, parcha.adult.survival.spline)

pp_check(parcha.adult.survival, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))
# pp_check(parcha.survival.spline, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))

# pp_check(parcha.survival, ndraws = 100, type = "stat", stat = "mean")
# pp_check(parcha.survival, ndraws = 100, type = "stat", stat = "sd")

saveRDS(parcha.adult.survival, "parcha.adult.survival.rds")
parcha.adult.survival <- readRDS("parcha.adult.survival.rds")









prediction_df <- expand.grid(log_age.t = seq(from = min((PARCHA_surv_adult$log_age.t)), to = max((PARCHA_surv_adult$log_age.t)), length.out = 25),
                             time_since_fire = c(min(PARCHA_surv_adult$time_since_fire), max(PARCHA_surv_adult$time_since_fire)),
                             rel_elev = mean(PARCHA_surv_adult$rel_elev),
                             bald = NA,
                             year.t1 = NA)

preds <- fitted(parcha.survival, newdata = prediction_df, re_formula = NA)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


# Making a binned survival df for visualization
parcha.survival.binned <- PARCHA_surv_adult %>% 
  ungroup() %>% 
  mutate(age_bin = cut_interval(log_age.t, n=25),
         fire_bin = cut_interval(time_since_fire_actual, n = 2)) %>%
  group_by(age_bin, fire_bin) %>% 
  summarize(mean_age = mean(log_age.t, na.rm = T),
            mean_surv = mean(surv.t1, na.rm = T),
            samplesize = n())


ggplot(data = prediction_df)+
  geom_point(data = PARCHA_surv_adult, aes( x= (log_age.t), y = surv.t1), color = "grey55", alpha = .2)+
  geom_point(data = parcha.survival.binned, aes( x= (mean_age), y = mean_surv, size = samplesize, color = fire_bin), alpha = .5)+
  geom_ribbon(aes(x = (log_age.t), ymax = upr, ymin = lwr, fill = as.factor(time_since_fire)), alpha = .4)+
  geom_line(aes(x = (log_age.t), y = fit, color = as.factor(time_since_fire)), linewidth = 1) +
  lims(y = c(0,1))+
  labs(y = "Survival (%)") + theme_minimal()

ggplot(data = prediction_df)+
  geom_point(data = PARCHA_surv_adult, aes( x= exp(log_age.t), y = surv.t1), color = "grey55", alpha = .2)+
  geom_point(data = parcha.survival.binned, aes( x= exp(mean_age), y = mean_surv, size = samplesize, color = fire_bin), alpha = .5)+
  geom_ribbon(aes(x = exp(log_age.t), ymax = upr, ymin = lwr, fill = as.factor(time_since_fire)), alpha = .4)+
  geom_line(aes(x = exp(log_age.t), y = fit, color = as.factor(time_since_fire)), linewidth = 1) +
  lims(y = c(0,1))+
  labs(y = "Survival (%)") + theme_minimal()





####################################################################################
# Seedling Survival models ------------------------------------------------------------------
####################################################################################
### PARCHA  ----

parcha.seedling.survival <- brm(surv.t1~ 0 + season + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_surv_seedling,
                             #surv.t1~ 1 + log_age.t, data = PARCHA_surv.df,
                             family = "bernoulli",
                             prior = c(set_prior("normal(0,1)", class = "b"),
                                       set_prior("normal(0,1)", class = "sd")),
                             #set_prior("normal(0,1)", class = "sds")),
                             warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains,
                             control = list(adapt_delta = 0.99))
# The model with spline wins model selection against one without
# parcha.seedling.survival.spline <- brm(surv.t1~ 1 +  s(log_age.t, bs = "tp", k=6) + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_surv_seedling,
#                              #surv.t1~ 1 + log_age.t, data = PARCHA_surv.df,
#                              family = "bernoulli",
#                              prior = c(set_prior("normal(0,1)", class = "b"),
#                                        set_prior("normal(0,1)", class = "Intercept"),
#                                        set_prior("normal(0,1)", class = "sd"),
#                                        set_prior("normal(0,1)", class = "sds")),
#                              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)
# # linear model wins here
# loo(parcha.seedling.survival, parcha.seedling.survival.spline)

pp_check(parcha.seedling.survival, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))
# pp_check(parcha.survival.spline, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))

# pp_check(parcha.survival, ndraws = 100, type = "stat", stat = "mean")
# pp_check(parcha.survival, ndraws = 100, type = "stat", stat = "sd")

saveRDS(parcha.seedling.survival, "parcha.seedling.survival.rds")
parcha.seedling.survival <- readRDS("parcha.seedling.survival.rds")



ggplot(PARCHA_surv_seedling)+
  geom_tile(aes(x = year.t1, bald))





prediction_df <- expand.grid(#log_age.t = seq(from = min((PARCHA_surv_seedling$log_age.t)), to = max((PARCHA_surv_seedling$log_age.t)), length.out = 25),
                             time_since_fire = seq(from = min(PARCHA_surv_seedling$time_since_fire), to = max(PARCHA_surv_seedling$time_since_fire), length.out = 25),
                             rel_elev = max(PARCHA_surv_seedling$rel_elev),
                             season = unique(PARCHA_surv_seedling$season),
                             bald = NA,
                             year.t1 = NA)

preds <- fitted(parcha.seedling.survival, newdata = prediction_df, re_formula = NA)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


# Making a binned survival df for visualization
parcha.survival.binned <- PARCHA_surv_seedling %>% 
  ungroup() %>% 
  mutate(#age_bin = cut_interval(log_age.t, n=25),
         fire_bin = cut_interval(time_since_fire, n = 20)) %>%
  group_by(fire_bin, season) %>% 
  summarize(mean_fire = mean(time_since_fire, na.rm = T),
            mean_surv = mean(surv.t1, na.rm = T),
            samplesize = n())


ggplot(data = prediction_df)+
  geom_point(data = PARCHA_surv_seedling, aes( x= (time_since_fire), y = surv.t1), color = "grey55", alpha = .2)+
  geom_point(data = parcha.survival.binned, aes( x= (mean_fire), y = mean_surv, size = samplesize), alpha = .5)+
  geom_ribbon(aes(x = (time_since_fire), ymax = upr, ymin = lwr), alpha = .4)+
  geom_line(aes(x = (time_since_fire), y = fit), linewidth = 1) +
  facet_wrap(~season)+
  lims(y = c(0,1))+
  labs(y = "Survival (%)") + theme_minimal()






####################################################################################
# Adult Flowering models ------------------------------------------------------------------
####################################################################################

#######
# need to think carefully if I should use year t or year t1
parcha.adult.flowering <- brm(flw.t1~ 0 + season +  log_age.t + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_flw_adult,
                              #flw.t1~ 1 + log_age.t, data = PARCHA_flw.df,
                              family = "bernoulli",
                              prior = c(set_prior("normal(0,1)", class = "b"),
                                        set_prior("normal(0,1)", class = "sd")),
                              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains,
                              control = list(adapt_delta = 0.99))
# parcha.adult.flowering.spline <- brm(flw.t1~ 1 +  s(log_age.t, bs = "tp", k = 6) + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_flw_adult,
#                               #flw.t1~ 1 + log_age.t, data = PARCHA_flw.df,
#                               family = "bernoulli",
#                               prior = c(set_prior("normal(0,1)", class = "b"),
#                                         set_prior("normal(0,1)", class = "Intercept"),
#                                         set_prior("normal(0,1)", class = "sd"),
#                                         set_prior("normal(0,1)", class = "sds")),
#                               warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)

# checking difference with loo, spline, is slightly "better" but only just, so proceeding with linear model
# loo(parcha.adult.flowering,parcha.adult.flowering.spline)
pp_check(parcha.adult.flowering, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))
# pp_check(parcha.flowering.spline, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))

# pp_check(parcha.flowering, ndraws = 100, type = "stat", stat = "mean")
# pp_check(parcha.flowering, ndraws = 100, type = "stat", stat = "sd")

saveRDS(parcha.adult.flowering, "parcha.adult.flowering.rds")
parcha.adult.flowering <- readRDS("parcha.adult.flowering.rds")



# 
# ggplot(PARCHA_flw_adult)+
#   geom_tile(aes(time_since_fire, rel_elev)) +
#   facet_wrap(~year.t1)




prediction_df <- expand.grid(log_age.t = seq(from = min(PARCHA_flw_adult$log_age.t, na.rm = T), to = max(PARCHA_flw_adult$log_age.t, na.rm = T), length.out = 25),
                             time_since_fire = mean(PARCHA_flw_adult$time_since_fire),
                             rel_elev = max(PARCHA_flw_adult$rel_elev),
                             bald = NA,
                             year.t1 = NA)

preds <- fitted(parcha.adult.flowering, newdata = prediction_df, re_formula = NA)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


# Making a binned flowering df for visualization
parcha.flowering.binned <- PARCHA_flw_adult %>% 
  ungroup() %>% 
  mutate(age_bin = cut_interval(log_age.t, n=15)) %>%
  group_by(age_bin) %>% 
  summarize(mean_age = mean(log_age.t, na.rm = T),
            mean_flw = mean(flw.t1, na.rm = T),
            samplesize = n())


ggplot(data = prediction_df)+
  geom_point(data = PARCHA_flw_adult, aes( x= (log_age.t), y = flw.t), color = "grey55", alpha = .2)+
  geom_point(data = parcha.flowering.binned, aes( x= (mean_age), y = mean_flw, size = samplesize), color = "skyblue3", alpha = .5)+
  geom_ribbon(aes(x = (log_age.t), ymax = upr, ymin = lwr), alpha = .4)+
  geom_line(aes(x = (log_age.t), y = fit), linewidth = 1) +
  # lims(y = c(0,1))+
  labs(y = "Flowering (%)") + theme_minimal()






####################################################################################
# Seedling Flowering models ------------------------------------------------------------------
####################################################################################
### PARCHA  ----

parcha.seedling.flowering <- brm(flw.t1~ 0 + season + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_flw_seedling,
                                #flw.t1~ 1 + log_age.t, data = PARCHA_flw.df,
                                family = "bernoulli",
                                prior = c(set_prior("normal(0,1)", class = "b"),
                                          set_prior("normal(0,1)", class = "sd")),
                                control = list(adapt_delta = 0.99),
                                warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)


pp_check(parcha.seedling.flowering, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))
# pp_check(parcha.flowering.spline, ndraws = 100, type = "dens_overlay")#+lims(x = c(0,20))

# pp_check(parcha.flowering, ndraws = 100, type = "stat", stat = "mean")
# pp_check(parcha.flowering, ndraws = 100, type = "stat", stat = "sd")

saveRDS(parcha.seedling.flowering, "parcha.seedling.flowering.rds")
parcha.seedling.flowering <- readRDS("parcha.seedling.flowering.rds")



ggplot(PARCHA_flw_seedling)+
  geom_tile(aes(x = year.t1, bald))





prediction_df <- expand.grid(#log_age.t = seq(from = min((PARCHA_flw_seedling$log_age.t)), to = max((PARCHA_flw_seedling$log_age.t)), length.out = 25),
  time_since_fire = seq(from = min(PARCHA_flw_seedling$time_since_fire), to = max(PARCHA_flw_seedling$time_since_fire), length.out = 25),
  rel_elev = max(PARCHA_flw_seedling$rel_elev),
  bald = NA,
  year.t1 = NA)

prediction_df <- expand.grid(#log_age.t = seq(from = min((PARCHA_flw_seedling$log_age.t)), to = max((PARCHA_flw_seedling$log_age.t)), length.out = 25),
  time_since_fire = mean(PARCHA_flw_seedling$time_since_fire),
  rel_elev = seq(from = min(PARCHA_flw_seedling$rel_elev), to = max(PARCHA_flw_seedling$rel_elev), length.out = 25),
  bald = NA,
  year.t1 = NA)

preds <- fitted(parcha.seedling.flowering, newdata = prediction_df, re_formula = NA)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


# Making a binned flowering df for visualization
parcha.flowering.binned <- PARCHA_flw_seedling %>% 
  ungroup() %>% 
  mutate(#age_bin = cut_interval(log_age.t, n=25),
    elev_bin = cut_interval(rel_elev, n = 10)) %>%
  group_by(elev_bin) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_flw = mean(flw.t1, na.rm = T),
            samplesize = n())


ggplot(data = prediction_df)+
  geom_point(data = PARCHA_flw_seedling, aes( x= (rel_elev), y = flw.t1), color = "grey55", alpha = .2)+
  geom_point(data = parcha.flowering.binned, aes( x= (mean_elev), y = mean_flw, size = samplesize), alpha = .5)+
  geom_ribbon(aes(x = (rel_elev), ymax = upr, ymin = lwr), alpha = .4)+
  geom_line(aes(x = (rel_elev), y = fit), linewidth = 1) +
  lims(y = c(0,1))+
  labs(y = "Flowering (%)") + theme_minimal()



####################################################################################
# Recruitment model ------------------------------------------------------------------
####################################################################################
# I'm basically assuming that each of the balds have equivalent full seedbands, and this gives a reasonable approx of how germination varies across balds/years
parcha.recruitment <- brm(recruits~ 0 + season + rel_elev + time_since_fire + (1|year.t1) + (1|bald), data = PARCHA_recruit.df,
                              #flw.t1~ 1 + log_age.t, data = PARCHA_flw.df,
                              family = negbinomial(),
                              prior = c(set_prior("normal(0,1)", class = "b"),
                                        set_prior("normal(0,1)", class = "sd")),
                              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains,
                              control = list(adapt_delta = 0.99))
# The Negative binomial actually does a good job!
pp_check(parcha.recruitment, ndraws = 100, type = "dens_overlay")+lims(x = c(0,5))





saveRDS(parcha.recruitment, "parcha.recruitment.rds")
parcha.recruitment <- readRDS("parcha.recruitment.rds")





prediction_df <- expand.grid(season = unique(PARCHA_recruit.df$season),
                             time_since_fire = mean(PARCHA_recruit.df$time_since_fire, na.rm = T),
                             rel_elev = seq(from = min(PARCHA_recruit.df$rel_elev), to = max(PARCHA_recruit.df$rel_elev), length.out = 3),
                             bald = NA,
                             year.t1 = 2003)

preds <- fitted(parcha.recruitment, newdata = prediction_df, re_formula = "~1|year.t1")
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


# Making a binned flowering df for visualization
parcha.recruitment.binned <- PARCHA_recruit.df %>% 
  filter(year.t1 == 2003) %>% 
  ungroup() %>% 
  mutate(#age_bin = cut_interval(log_age.t, n=25),
    elev_bin = cut_interval(rel_elev, n = 3)) %>%
  group_by(season, elev_bin) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_recruit = mean(recruits, na.rm = T),
            samplesize = n()) %>% 
  mutate(rel_elev = case_when(elev_bin == "[0.65,1.52]" ~ 0.65,
                              elev_bin == "(1.52,2.38]" ~ 1.95,
                              elev_bin == "(2.38,3.25]" ~ 3.25,))


ggplot(data = prediction_df)+
  geom_point(data = PARCHA_recruit.df, aes( x= (rel_elev), y = recruits), color = "grey55", alpha = .2)+
  geom_point(data = parcha.recruitment.binned, aes( x= (mean_elev), y = mean_recruit, size = samplesize), alpha = .5)+
  geom_ribbon(aes(x = (rel_elev), ymax = upr, ymin = lwr), alpha = .4)+
  geom_line(aes(x = (rel_elev), y = fit), linewidth = 1) +
  facet_wrap(~season, scales = "free")+
  labs(y = "Recruits") + theme_minimal()

ggplot(data = prediction_df)+
  # geom_point(data = PARCHA_recruit.df, aes( x= season, y = recruits), color = "grey55", alpha = .2)+
  geom_point(data = parcha.recruitment.binned, aes( x= season, y = mean_recruit, size = samplesize, group = elev_bin), alpha = .5)+
  geom_linerange(aes(x = season, ymax = upr, ymin = lwr, fill = as.factor(rel_elev), group = rel_elev), alpha = .4)+
  geom_point(aes(x = season, y = fit, color = as.factor(rel_elev)), linewidth = 1) +
  facet_wrap(~as.factor(rel_elev))+
  labs(y = "Recruits") + theme_minimal()




#########################################################################################
#  average number of flowers for a flowering individual in greenhouse --------
#########################################################################################

parcha.flw_number <- brm(flw_count~ 1 + live_sterile ,  data = flw_counts,
                          #flw.t1~ 1 + log_age.t, data = PARCHA_flw.df,
                          family = negbinomial(),
                          prior = c(set_prior("normal(0,5)", class = "Intercept"),
                                    set_prior("normal(0,5)", class = "b")),
                          warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)
# negative binomial does a reasonably good job
pp_check(parcha.flw_number, ndraws = 100, type = "dens_overlay")

saveRDS(parcha.flw_number, "parcha.flw_number.rds")
parcha.flw_number <- readRDS("parcha.flw_number.rds")



parcha.flw_number

flw_counts






####################################################################################
# Messing around, trying to account for unknown age --------------------------------
####################################################################################

#notes: lognormal does a pretty good job, although have to make 0 month ages actually 1 month, maybe a reasonable thing to do
# neg. binomial doesn't do as well, but maybe could recalculate the units in terms of censuses, rather than months
# I think it makes sense to treat age as contniuous. 
# or I could consider using age in years and do some more binning
#okay, actually realizing taht I should mayve use the actual census_date.t


lm1 <- lm(log_age.t ~ 1 + as.numeric(census_start)*as.numeric(census_end)*as.numeric(census_number), data = PARCHA_surv.df)
lm2 <- lm(log_age.t ~ 1 + as.numeric(census_start)*as.numeric(census_number) + as.numeric(census_end)*as.numeric(census_date.t), data = PARCHA_surv.df)
lm3 <- lm(log_age.t ~ 1 + as.numeric(census_start) + as.numeric(census_end) + as.numeric(census_number), data = PARCHA_surv.df)

BIC(object = lm1,lm2,lm3, lm4)
summary(lm4)

prediction_df$fit <- predict(lm, newdata=prediction_df)

parcha.age <- brm(log_age.t ~ 1 + census_start*census_end*census_number, data = PARCHA_surv.df,
                  family = gaussian,
                  prior = c(set_prior("normal(0,.1)", class = "b"),
                            set_prior("normal(0,.1)", class = "Intercept"),
                            set_prior("normal(0,.1)", class = "sigma")),
                  warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)

pp_check(parcha.age, ndraws = 50, type = "dens_overlay", stat = "sd")+lims(x = c(0,20))
# pp_check(parcha.age, ndraws = 50, type = "stat", stat = "sd")

prediction_df <- expand.grid(census_number = seq(from = as.numeric(min(PARCHA_surv.df$census_number, na.rm = T)), to = as.numeric(max(PARCHA_surv.df$census_number, na.rm = T)), length.out = 25),
                             census_start = c(as.numeric(min(PARCHA_surv.df$census_start))),
                             census_end = c(as.numeric(min(PARCHA_surv.df$census_end, na.rm = T)),as.numeric(mean(PARCHA_surv.df$census_end, na.rm = T)),as.numeric(max(PARCHA_surv.df$census_end, na.rm = T))),
                             bald = NA)

preds <- fitted(parcha.age, newdata = prediction_df, re_formula = NA)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]


ggplot(data = prediction_df)+
  geom_point(data = PARCHA_surv.df, aes( x= as.numeric(census_number), y = log_age.t), color = "grey55", alpha = .2)+
  # geom_point(data = parcha.age.binned, aes( x= mean_census, y = mean_age, size = samplesize), color = "skyblue3", alpha = .5)+
  geom_ribbon(aes(x = as.numeric(census_number), ymax = upr, ymin = lwr,fill =as.factor(as.numeric(census_end))), alpha = .4)+
  geom_line(aes(x = as.numeric(census_number), y = fit, color = as.factor(as.numeric(census_end))), linewidth = 1) +
  lims( y = c(0,50))+
  labs(y = "Age (seasons)") + theme_minimal()








formula.surv <- bf(surv.t1~ 1 + mi(age))+ bernoulli()
formula.age <- bf(age|mi() ~ 1 + census_start + census_end + census_number) + gaussian()
parcha.survival <- brm(formula.surv + formula.age + set_rescor(FALSE), data = PARCHA_surv.df,
                       prior = c(set_prior("normal(0,1)", class = "b")),
                       warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)



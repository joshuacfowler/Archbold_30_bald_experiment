# Archbold 30-Bald Greenhouse Experimental analysis
# Assessment of germination success and growth for focal species and deriving estimates of microbial effects
# author: Joshua Fowler
# Date: Feb 28, 2024
library(tidyverse)
library(readxl)
library(lubridate)
library(lme4)
library(brms)
library(rstan)



invlogit<-function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x)/(1+(x))}


###############################################################################################
########## Reading in data ####################################################################
###############################################################################################
plants <- read_xlsx(path = "~/Dropbox/UofMiami/Archbold_30baldexperiment.xlsx", sheet = "Census", guess_max = 1048576) # guess_max makes the function look deeper in the columns to assign type

# Summarizing the germination for each pot
germ_census <- plants %>% 
  mutate(across(everything(), as.character)) %>% 
  dplyr::select(c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds,reseeding, date_potted), contains("germ"), -contains("notes")) %>% 
  # pivot_longer(cols = contains("germination"))
  pivot_longer(cols = c(contains("germ"))) %>% 
  separate(name, c("measurement", "census_number", "name")) %>% 
  pivot_wider(id_cols = c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds, reseeding, date_potted, census_number), names_from = c(measurement,name), values_from = value, names_repair = "minimal") 


per_pot_germ <- germ_census %>% 
  group_by(pot_id, live_sterile,soil_source, rep_id, spp_code, number_seeds, reseeding) %>% 
  dplyr::summarize(GermCount = max(as.numeric(GermCount_measurement), na.rm = T),
                   ExtraGerm_total = sum(as.numeric(AdditGermCount_measurement), na.rm = T),
                   TotalGerm = sum(GermCount, ExtraGerm_total, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(number_seeds = as.numeric(number_seeds),
         reseeding = case_when(is.na(reseeding) ~ 0, 
                               TRUE ~ as.numeric(reseeding)),
         total_seeds = number_seeds + reseeding,
         seed_prop = TotalGerm/total_seeds) %>% 
  filter(GermCount != -Inf) #%>% # These are the pots that have no measurements because we dropped them from the experiment 






# Summarizing the growth (final size) for each pot

size_census <- plants %>% 
  mutate(across(everything(), as.character)) %>% 
  dplyr::select(c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds,reseeding, date_potted), contains("size"), -contains("notes")) %>% 
  # pivot_longer(cols = contains("germination"))
  pivot_longer(cols = c(contains("size"))) %>% 
  separate(name, c("measurement", "census_number", "name")) %>% 
  pivot_wider(id_cols = c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds, reseeding, date_potted, census_number), names_from = c(measurement,name), values_from = value, names_repair = "minimal") %>% 
  mutate(HeightSize_measurement = as.numeric(HeightSize_measurement),
         DiameterSize_measurement = as.numeric(DiameterSize_measurement))


per_pot_size <- size_census %>% 
  mutate(Size = case_when(!is.na(HeightSize_measurement) & !is.na(DiameterSize_measurement) ~ pi*((DiameterSize_measurement/2)^2)*(HeightSize_measurement/3),
                          is.na(HeightSize_measurement) ~ DiameterSize_measurement,
                          is.na(DiameterSize_measurement) ~ HeightSize_measurement,
                          TRUE ~ NA)) %>% 
  filter(census_number == 10)  %>%  # for now just taking the final size census, although I think we could take size at earlier time points for ones that died
  group_by(pot_id, live_sterile,soil_source, rep_id, spp_code, number_seeds, reseeding) %>% 
  dplyr::summarize(finalHeight = HeightSize_measurement,
                   finalDiameter = DiameterSize_measurement,
                   finalSize = Size) %>% 
  ungroup() %>% 
  filter(!is.na(finalSize))




flw_census <- plants %>% 
  mutate(across(everything(), as.character)) %>% 
  dplyr::select(c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds,reseeding, date_potted), contains("Flw"), -contains("notes")) %>% 
  # pivot_longer(cols = contains("germination"))
  pivot_longer(cols = c(contains("Flw"))) %>% 
  separate(name, c("measurement", "census_number", "name")) %>% 
  pivot_wider(id_cols = c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds, reseeding, date_potted, census_number), names_from = c(measurement,name), values_from = value, names_repair = "minimal")

per_pot_flw_clean <- flw_census %>%  # for now just taking the final size census, although I think we could take size at earlier time points for ones that died
  group_by(pot_id, live_sterile,soil_source, rep_id, spp_code, number_seeds, reseeding) %>% 
  mutate(FLW_COUNT = as.numeric(FlwCount_measurement),
         FLW_STATUS = case_when(FLW_COUNT >= 1 ~ 1, FLW_COUNT == 0 ~ 0, is.na(FLW_COUNT) ~ 0)) 
  
  
  
per_pot_flw_summary <- per_pot_flw_clean %>% 
  filter(FLW_STATUS == 1) %>% 
    summarize(flw_ever = as.numeric(mean(FLW_STATUS>0, na.rm = T)>0),
            first_flw = min(FlwCount_date)) %>%  
  filter(spp_code %in% c("BALANG", "HYPCUM", "PARCHA", "POLBAS","LECCER", "LECDEC", "ERYCUN", "CHAFAS"))



per_pot_flw <- per_pot_flw_clean %>% 
  left_join(per_pot_flw_summary) %>% 
  slice(which.max(FLW_COUNT)) %>% # selecting the largest flower count
  mutate(flw_ever = case_when(is.na(flw_ever) ~ 0,
                              !is.na(flw_ever) ~ flw_ever))
  
  
per_pot_fert <- per_pot_flw %>% 
  slice(which.max(FLW_COUNT)) %>% # selecting the largest flower count
  filter(spp_code %in% c("HYPCUM", "PARCHA", "ERYCUN")) # BALANG, LECCER, LECDEC, and POLBAS had all less than 3 plants flower, so while I did record flower number, I don't think it's wise to model this. ERYCUN, I think I could potentially use, but but some of the individuals were marked as flowering when they were just bolting so they don't actually have accurate counts. There for I am only analyzing CHAFAS, HYPCUM, PARCHA

# Now I want to merge in the fire history and elevation data for each bald


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
  dplyr::select(bald, rel_elev) %>% 
  mutate(
         bald = case_when(bald == "01S" ~ "1S",
                          bald == "01N" ~ "1N",
                          bald == "02" ~ "2",
                          bald == "05E" ~ "5E",
                          bald == "07N" ~ "7N",
                          bald == "35N" ~ "35",
                          bald == "65E" ~ "65",
                          TRUE ~ bald))




germ.covariates <- per_pot_germ %>%
  left_join(fire_summary, by = join_by( soil_source == Bald_U)) %>% 
  left_join(elev.df, by = join_by(soil_source == bald))




size.covariates <- per_pot_size %>%
  left_join(fire_summary, by = join_by( soil_source == Bald_U)) %>% 
  left_join(elev.df, by = join_by(soil_source == bald))



flw.covariates <- per_pot_flw %>%
  left_join(fire_summary, by = join_by( soil_source == Bald_U)) %>% 
  left_join(elev.df, by = join_by(soil_source == bald))


fert.covariates <- per_pot_fert %>%
  left_join(fire_summary, by = join_by( soil_source == Bald_U)) %>% 
  left_join(elev.df, by = join_by(soil_source == bald))



# write_csv(germ.covariates, file = "germ.covariates.csv")

# write_csv(size.covariates, file = "size.covariates.csv")

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

my_palette <- c("#000000", "#E69F00", "#009E73")

################################################################################
######### Regressions of germination rate with environmental covariates ########
################################################################################

# starting first with a model without environmental covariates
germ.m <- brm(TotalGerm|trials(total_seeds) ~ -1 + spp_code* live_sterile, data = germ.covariates,
                            family = binomial(link = "logit"),
                            prior = c(set_prior("normal(0,1)", class = "b")),
                            warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)



# Making prediction dataframe

prediction_df <- expand.grid(spp_code = unique(germ.covariates$spp_code), live_sterile = unique(germ.covariates$live_sterile), total_seeds = 1)

preds <- fitted(germ.m, newdata = prediction_df)
prediction_df$fit <- preds[,"Estimate"]
prediction_df$lwr <- preds[,"Q2.5"]
prediction_df$upr <- preds[,"Q97.5"]

# now we can plot the model predictions

ggplot(data = prediction_df)+
  # geom_jitter(data = germ.covariates, aes( x= live_sterile, y = seed_prop), color = "blue", width = .2, alpha = .2)+
  geom_point(aes(x = live_sterile, y = fit), size = 1) +
  geom_linerange(aes(x = live_sterile, ymin = lwr, ymax = upr))+
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()




# now incorporating the environmental covariates and a random effect
germ.m <- brm(TotalGerm|trials(total_seeds) ~ -1 + spp_code*live_sterile*rel_elev + spp_code*live_sterile*fire_frequency + (1|soil_source), data = germ.covariates,
              family = binomial(link = "logit"),
              prior = c(set_prior("normal(0,1)", class = "b")),
              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)



# Making prediction dataframe
prediction_df.1 <- expand.grid(spp_code = unique(germ.covariates$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(germ.covariates$live_sterile),
                             rel_elev = c(median(germ.covariates$rel_elev)), fire_frequency = seq(min(germ.covariates$fire_frequency), max(germ.covariates$fire_frequency), by = .2))

prediction_df.2 <- expand.grid(spp_code = unique(germ.covariates$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(germ.covariates$live_sterile), 
                             rel_elev = seq(min(germ.covariates$rel_elev), max(germ.covariates$rel_elev), by = .2), fire_frequency = c(median(germ.covariates$fire_frequency)))


preds.1 <- fitted(germ.m, newdata = prediction_df.1, re_formula = NA)
prediction_df.1$fit <- preds.1[,"Estimate"]
prediction_df.1$lwr <- preds.1[,"Q2.5"]
prediction_df.1$upr <- preds.1[,"Q97.5"]


preds.2 <- fitted(germ.m, newdata = prediction_df.2, re_formula = NA)
prediction_df.2$fit <- preds.2[,"Estimate"]
prediction_df.2$lwr <- preds.2[,"Q2.5"]
prediction_df.2$upr <- preds.2[,"Q97.5"]


# now we can plot the model predictions

germ.binned.1 <- germ.covariates %>% 
  mutate(fire_bin = fire_frequency, 3) %>% 
  group_by(spp_code, fire_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_germ = mean(seed_prop, na.rm = T),
            mean_fire = mean(fire_frequency, na.rm = T))

germ.binned.2 <- germ.covariates %>% 
         mutate(elev_bin = cut_number(rel_elev, 10)) %>% 
  group_by(spp_code, elev_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_germ = mean(seed_prop, na.rm = T))




germ_plot.1 <- ggplot(data = prediction_df.1)+
  geom_ribbon(aes(x = fire_frequency, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = fire_frequency, y = fit, group = live_sterile, color = live_sterile)) +
  geom_point(data = germ.binned.1, aes( x= mean_fire, y = mean_germ, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code) + labs(y = "Proportion Germinated") + theme_light()

germ_plot.1



ggsave(germ_plot.1, filename = "germ_plot_fire.png",  width = 6, height = 6)


germ_plot.2 <- ggplot(data = prediction_df.2)+
  geom_ribbon(aes(x = rel_elev, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = rel_elev, y = fit, group = live_sterile, color = live_sterile)) +
  geom_point(data = germ.binned.2, aes( x= mean_elev, y = mean_germ, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code) + labs(y = "Proportion Germinated") + theme_light()

germ_plot.2

ggsave(germ_plot.2, filename = "germ_plot_elev.png", width = 6, height = 6)





  
########################################################################################
############ Regressions of size with environmental covariates ###########  
########################################################################################
  #Need to calculate ERYCUN based on diameter, rather than combining them because this is what is in the field, but could analyze a couple of ways
grow_data <- size.covariates  %>% 
  mutate(finalSize = case_when(spp_code == "ERYCUN" ~ finalDiameter,
                               TRUE~ finalSize),
         log_finalSize = log(finalSize))
  

grow.m <- brm(log_finalSize ~ -1 + spp_code*live_sterile*rel_elev + spp_code*live_sterile*fire_frequency + (1|soil_source), data = grow_data,
              family = "gaussian",
              prior = c(set_prior("normal(0,10)", class = "b")),
              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)



# Making prediction dataframe
prediction_df.1 <- expand.grid(spp_code = unique(grow_data$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(grow_data$live_sterile),
                               rel_elev = c(median(grow_data$rel_elev)), fire_frequency = seq(min(grow_data$fire_frequency), max(grow_data$fire_frequency), by = .2))

prediction_df.2 <- expand.grid(spp_code = unique(grow_data$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(grow_data$live_sterile), 
                               rel_elev = seq(min(grow_data$rel_elev), max(grow_data$rel_elev), by = .2), fire_frequency = c(median(grow_data$fire_frequency)))


preds.1 <- fitted(grow.m, newdata = prediction_df.1, re_formula = NA)
prediction_df.1$fit <- preds.1[,"Estimate"]
prediction_df.1$lwr <- preds.1[,"Q2.5"]
prediction_df.1$upr <- preds.1[,"Q97.5"]


preds.2 <- fitted(grow.m, newdata = prediction_df.2, re_formula = NA)
prediction_df.2$fit <- preds.2[,"Estimate"]
prediction_df.2$lwr <- preds.2[,"Q2.5"]
prediction_df.2$upr <- preds.2[,"Q97.5"]


# now we can plot the model predictions

grow.binned.1 <- grow_data %>% 
  mutate(fire_bin = fire_frequency, 3) %>% 
  group_by(spp_code, fire_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_grow = mean(log_finalSize, na.rm = T),
            mean_fire = mean(fire_frequency, na.rm = T))

grow.binned.2 <- grow_data %>% 
  mutate(elev_bin = cut_number(rel_elev, 10)) %>% 
  group_by(spp_code, elev_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_grow = mean(log_finalSize, na.rm = T))




grow_plot.1 <- ggplot(data = prediction_df.1)+
  geom_ribbon(aes(x = fire_frequency, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = fire_frequency, y = fit, group = live_sterile, color = live_sterile)) +
  geom_point(data = grow.binned.1, aes( x= mean_fire, y = mean_grow, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "Size") + theme_light()

grow_plot.1



ggsave(grow_plot.1, filename = "grow_plot_fire.png",  width = 6, height = 6)


grow_plot.2 <- ggplot(data = prediction_df.2)+
  geom_ribbon(aes(x = rel_elev, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = rel_elev, y = fit, group = live_sterile, color = live_sterile)) +
  geom_point(data = grow.binned.2, aes( x= mean_elev, y = mean_grow, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code) + labs(y = "Size") + theme_light()

grow_plot.2

ggsave(grow_plot.2, filename = "grow_plot_elev.png", width = 6, height = 6)





########################################################################################
######### Regressions of probability of flowering with environmental covariates ########
########################################################################################
  
  # starting first with a model without environmental covariates
flw.m <- brm(flw_ever ~ -1 + spp_code*live_sterile*rel_elev + spp_code*live_sterile*fire_frequency + (1|soil_source), data = flw.covariates,
                family = bernoulli,
                prior = c(set_prior("normal(0,1)", class = "b")),
                warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)



# Making prediction dataframe
prediction_df.1 <- expand.grid(spp_code = unique(flw.covariates$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(flw.covariates$live_sterile),
                               rel_elev = c(median(flw.covariates$rel_elev)), fire_frequency = seq(min(flw.covariates$fire_frequency), max(flw.covariates$fire_frequency), by = .2))

prediction_df.2 <- expand.grid(spp_code = unique(flw.covariates$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(flw.covariates$live_sterile), 
                               rel_elev = seq(min(flw.covariates$rel_elev), max(flw.covariates$rel_elev), by = .2), fire_frequency = c(median(flw.covariates$fire_frequency)))


preds.1 <- fitted(flw.m, newdata = prediction_df.1, re_formula = NA)
prediction_df.1$fit <- preds.1[,"Estimate"]
prediction_df.1$lwr <- preds.1[,"Q2.5"]
prediction_df.1$upr <- preds.1[,"Q97.5"]


preds.2 <- fitted(flw.m, newdata = prediction_df.2, re_formula = NA)
prediction_df.2$fit <- preds.2[,"Estimate"]
prediction_df.2$lwr <- preds.2[,"Q2.5"]
prediction_df.2$upr <- preds.2[,"Q97.5"]


# now we can plot the model predictions

flw.binned.1 <- flw.covariates %>% 
  mutate(fire_bin = fire_frequency, 3) %>% 
  group_by(spp_code, fire_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_flw = mean(flw_ever, na.rm = T),
            mean_fire = mean(fire_frequency, na.rm = T))

flw.binned.2 <- flw.covariates %>% 
  mutate(elev_bin = cut_number(rel_elev, 1)) %>% 
  group_by(spp_code, elev_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_flw = mean(flw_ever, na.rm = T))




flw_plot.1 <- ggplot(data = prediction_df.1)+
  geom_point(data = flw.covariates, aes(x = fire_frequency, y = flw_ever, color = live_sterile), position = position_jitter(width = 0.05, height =0.05), alpha = .3)+
  geom_ribbon(aes(x = fire_frequency, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = fire_frequency, y = fit, group = live_sterile, color = live_sterile)) +
  # geom_point(data = flw.binned.1, aes( x= mean_fire, y = mean_flw, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code) + labs(y = "Flowering Probability") + theme_light()

flw_plot.1



ggsave(flw_plot.1, filename = "flw_plot_fire.png",  width = 6, height = 6)


flw_plot.2 <- ggplot(data = prediction_df.2)+
  geom_point(data = flw.covariates, aes(x = rel_elev, y = flw_ever, color = live_sterile), position = position_jitter(width = 0.05, height =0.05), alpha = .3)+
  geom_ribbon(aes(x = rel_elev, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = rel_elev, y = fit, group = live_sterile, color = live_sterile)) +
  # geom_point(data = flw.binned.2, aes( x= mean_elev, y = mean_flw, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code) + labs(y = "Flowering Probability") + theme_light()

flw_plot.2

ggsave(flw_plot.2, filename = "flw_plot_elev.png", width = 6, height = 6)


########################################################################################
############ Regressions of size with environmental covariates ###########  
########################################################################################
#Need to calculate ERYCUN based on diameter, rather than combining them because this is what is in the field, but could analyze a couple of ways
fert_data <- fert.covariates  %>% 
  filter(FLW_COUNT >0)

fert.m <- brm(FLW_COUNT ~ -1 + spp_code*live_sterile*rel_elev + spp_code*live_sterile*fire_frequency + (1|soil_source), data = fert_data,
              family = "negbinomial",
              prior = c(set_prior("normal(0,10)", class = "b")),
              warmup = mcmc_pars$warmup, iter = mcmc_pars$iter, chains = mcmc_pars$chains)

pp_check(fert.m) + lims(x = c(0,100))

# Making prediction dataframe
prediction_df.1 <- expand.grid(spp_code = unique(fert_data$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(fert_data$live_sterile),
                               rel_elev = c(median(fert_data$rel_elev)), fire_frequency = seq(min(fert_data$fire_frequency), max(fert_data$fire_frequency), by = .2))

prediction_df.2 <- expand.grid(spp_code = unique(fert_data$spp_code), total_seeds = 1, soil_source = NA, live_sterile = unique(fert_data$live_sterile), 
                               rel_elev = seq(min(fert_data$rel_elev), max(fert_data$rel_elev), by = .2), fire_frequency = c(median(fert_data$fire_frequency)))


preds.1 <- fitted(fert.m, newdata = prediction_df.1, re_formula = NA)
prediction_df.1$fit <- preds.1[,"Estimate"]
prediction_df.1$lwr <- preds.1[,"Q2.5"]
prediction_df.1$upr <- preds.1[,"Q97.5"]


preds.2 <- fitted(fert.m, newdata = prediction_df.2, re_formula = NA)
prediction_df.2$fit <- preds.2[,"Estimate"]
prediction_df.2$lwr <- preds.2[,"Q2.5"]
prediction_df.2$upr <- preds.2[,"Q97.5"]


# now we can plot the model predictions

fert.binned.1 <- fert_data %>% 
  mutate(fire_bin = fire_frequency, 3) %>% 
  group_by(spp_code, fire_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_fert = mean(FLW_COUNT, na.rm = T),
            mean_fire = mean(fire_frequency, na.rm = T))

fert.binned.2 <- fert_data %>% 
  mutate(elev_bin = cut_number(rel_elev, 10)) %>% 
  group_by(spp_code, elev_bin,  live_sterile) %>% 
  summarize(mean_elev = mean(rel_elev, na.rm = T),
            mean_fert = mean(FLW_COUNT, na.rm = T))




fert_plot.1 <- ggplot(data = prediction_df.1)+
  geom_point(data = fert_data, aes( x= fire_frequency, y = FLW_COUNT, color = live_sterile), alpha = .5)+
  geom_ribbon(aes(x = fire_frequency, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = fire_frequency, y = fit, group = live_sterile, color = live_sterile)) +
  geom_point(data = fert.binned.1, aes( x= mean_fire, y = mean_fert, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "# of Flws or Stems") + theme_light()

fert_plot.1



ggsave(fert_plot.1, filename = "fert_plot_fire.png",  width = 6, height = 6)


fert_plot.2 <- ggplot(data = prediction_df.2)+
  geom_point(data = fert_data, aes(x = rel_elev, y = FLW_COUNT, color = live_sterile), alpha = .5)+
  geom_ribbon(aes(x = rel_elev, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = rel_elev, y = fit, group = live_sterile, color = live_sterile)) +
  # geom_point(data = fert.binned.2, aes( x= mean_elev, y = mean_fert, color = live_sterile), alpha = .5)+
  scale_color_manual(values = c(my_palette[1], my_palette[3]))+
  scale_fill_manual(values = c(my_palette[1], my_palette[3]))+
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "# of Flws or Stems") + theme_light()

fert_plot.2

ggsave(fert_plot.2, filename = "fert_plot_elev.png", width = 6, height = 6)





  

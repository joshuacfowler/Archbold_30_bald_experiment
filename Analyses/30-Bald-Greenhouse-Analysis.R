# Archbold 30-Bald Greenhouse Experimental analysis
# Assessment of germination success and growth for focal species and deriving estimates of microbial effects
# author: Joshua Fowler
# Date: Feb 28, 2024
library(tidyverse)
library(readxl)
library(lubridate)
library(lme4)




invlogit<-function(x){exp(x)/(1+exp(x))}


###############################################################################################
########## Reading in data ####################################################################
###############################################################################################
plants <- read_xlsx(path = "~/Dropbox/UofMiami/Archbold_30baldexperiment.xlsx", sheet = "Census", guess_max = 1048576) # guess_max makes the function look deeper in the columns to assign type

census <- plants %>% 
  mutate(across(everything(), as.character)) %>% 
  dplyr::select(c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds,reseeding, date_potted), contains("germ"), -contains("notes")) %>% 
  # pivot_longer(cols = contains("germination"))
  pivot_longer(cols = c(contains("germ"))) %>% 
  separate(name, c("measurement", "census_number", "name")) %>% 
  pivot_wider(id_cols = c(pot_id, x_id, y_id, rep_id, table_id, soil_source, live_sterile, spp_code, number_seeds, reseeding, date_potted, census_number), names_from = c(measurement,name), values_from = value, names_repair = "minimal") 


per_pot_germ <- census %>% 
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

# this is how to fit a generalized linear model, this is more "correct" for our data. 
# The -1 in the formula indicates that I want to use the spp name as the first intercept, also worth noting that the first "live_sterile" parameter listed in the summary is for BALANG
germ.m <- glm(cbind(TotalGerm, total_seeds - TotalGerm) ~ -1 + spp_code* live_sterile, data = germ.covariates, family = binomial)

summary(germ.m)
anova(germ.m, test = "Chisq")

# now we will make our model predictions. We need to set up a dataframe of all of our predictor variables
# then we use the predict function, and have to do a tranformation to get things in the scale of the proportion germinated
prediction_df <- expand.grid(spp_code = unique(germ.covariates$spp_code), live_sterile = unique(germ.covariates$live_sterile))

preds <- predict(germ.m, newdata = prediction_df, type = "link", se.fit = TRUE)


critval <- 1.96 ## approx 95% CI
prediction_df$upr <- germ.m$family$linkinv(preds$fit + (critval * preds$se.fit))
prediction_df$lwr <- germ.m$family$linkinv(preds$fit - (critval * preds$se.fit))
prediction_df$fit <- germ.m$family$linkinv(preds$fit)


# now we can plot the model predictions

ggplot(data = prediction_df)+
  # geom_jitter(data = germ.covariates, aes( x= live_sterile, y = seed_prop), color = "blue", alpha = .2)+
  geom_point(aes(x = live_sterile, y = fit)) +
  geom_linerange(aes(x = live_sterile, ymin = lwr, ymax = upr))+
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()



# Now if we want to get an effect of microbes for each species, we can take the difference between sterile and live, and then save the dataframe
micro_effect_df <- prediction_df %>% 
  pivot_wider(id_cols = c(spp_code), names_from = live_sterile, values_from = fit) %>% 
  mutate(micro_effect = live - sterile)


####### now fitting a model to test for variation across the soil source ######
# I will actually fit a glmm later to account for bald variation and use actual covariates, but I just want to look the raw variation
# The -1 in the formula indicates that I want to use the spp name as the first intercept, also worth noting that the first "live_sterile" parameter listed in the summary is for BALANG
germ.m <- glm(cbind(TotalGerm, total_seeds - TotalGerm) ~ -1 + spp_code* live_sterile*soil_source, data = germ.covariates, family = binomial)

summary(germ.m)
anova(germ.m, test = "Chisq")

# now we will make our model predictions. We need to set up a dataframe of all of our predictor variables
# then we use the predict function, and have to do a tranformation to get things in the scale of the proportion germinated
prediction_df <- expand.grid(spp_code = unique(germ.covariates$spp_code), live_sterile = unique(germ.covariates$live_sterile), soil_source = unique(germ.covariates$soil_source))

preds <- predict(germ.m, newdata = prediction_df, type = "link", se.fit = TRUE)


critval <- 1.96 ## approx 95% CI
prediction_df$upr <- germ.m$family$linkinv(preds$fit + (critval * preds$se.fit))
prediction_df$lwr <- germ.m$family$linkinv(preds$fit - (critval * preds$se.fit))
prediction_df$fit <- germ.m$family$linkinv(preds$fit)


# now we can plot the model predictions

ggplot(data = filter(prediction_df, spp_code == "HYPCUM"))+
  # geom_jitter(data = germ.covariates, aes( x= live_sterile, y = seed_prop), color = "blue", alpha = .2)+
  geom_point(aes(x = live_sterile, y = fit)) +
  geom_linerange(aes(x = live_sterile, ymin = lwr, ymax = upr))+
  facet_wrap(~spp_code + soil_source, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()



dataplot <- ggplot(data = filter(germ.covariates, spp_code == "HYPCUM"))+
  geom_jitter(aes( x= live_sterile, y = seed_prop), color = "blue", alpha = .2)+
  geom_stat(aes(x = live_sterile, ))
  # geom_point(aes(x = live_sterile, y = fit)) +
  # geom_linerange(aes(x = live_sterile, ymin = lwr, ymax = upr))+
  facet_wrap(~soil_source, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()

ggsave(dataplot, file = "soil_source_data_plot.png",width = 10, height = 10)


####### now fitting a random effects model to test for variation across the soil source ######

#dropping bald 97, which we mostly pulled due to contamination
data <- germ.covariates %>% filter(soil_source != 97)


#germ.m <- glm(cbind(TotalGerm, total_seeds - TotalGerm) ~ -1 + spp_code* live_sterile*rel_elev , data = data, family = binomial)

germ.m <- glmer(cbind(TotalGerm, total_seeds - TotalGerm) ~ -1 + spp_code* live_sterile*rel_elev + (1|soil_source), data = data, family = binomial,
                control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # adding this to try to help convergence issues

summary(germ.m)
anova(germ.m, test = "Chisq")



prediction_df <- expand.grid(spp_code = unique(germ.covariates$spp_code), live_sterile = unique(germ.covariates$live_sterile), rel_elev = seq(from = min(data$rel_elev), to = max(data$rel_elev), by = .2))

preds <- predict(germ.m, newdata = prediction_df, type = "link", se.fit = TRUE, re.form = NA)


critval <- 1.96 ## approx 95% CI
prediction_df$upr <- invlogit(preds$fit + (critval * preds$se.fit))
prediction_df$lwr <- invlogit(preds$fit - (critval * preds$se.fit))
prediction_df$fit <- invlogit(preds$fit)


# now we can plot the model predictions

ggplot(data = prediction_df)+
  # geom_point(data = germ.covariates, aes( x= rel_elev, y = seed_prop, color = live_sterile), alpha = .2)+
  geom_ribbon(aes(x = rel_elev, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
  geom_line(aes(x = rel_elev, y = fit, group = live_sterile, color = live_sterile)) +
  scale_fill_manual(values = c(""))
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()


  
  
# And now analyzing fire frequency
  
  
  germ.m <- glmer(cbind(TotalGerm, total_seeds - TotalGerm) ~ -1 + spp_code* live_sterile*fire_frequency + (1|soil_source), data = data, family = binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # adding this to try to help convergence issues
  
  summary(germ.m)
  anova(germ.m, test = "Chisq")
  
  
  
  prediction_df <- expand.grid(spp_code = unique(germ.covariates$spp_code), live_sterile = unique(germ.covariates$live_sterile), fire_frequency = seq(from = min(data$fire_frequency), to = max(data$fire_frequency), by = .2))
  
  preds <- predict(germ.m, newdata = prediction_df, type = "link", se.fit = TRUE, re.form = NA)
  
  
  critval <- 1.96 ## approx 95% CI
  prediction_df$upr <- invlogit(preds$fit + (critval * preds$se.fit))
  prediction_df$lwr <- invlogit(preds$fit - (critval * preds$se.fit))
  prediction_df$fit <- invlogit(preds$fit)
  
  
  # now we can plot the model predictions
  
  ggplot(data = prediction_df)+
    geom_point(data = germ.covariates, aes( x= fire_frequency, y = seed_prop, color = live_sterile), alpha = .2)+
    geom_ribbon(aes(x = fire_frequency, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
    geom_line(aes(x = fire_frequency, y = fit, group = live_sterile, color = live_sterile)) +
    # scale_fill_manual(values = c(""))
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()
  
  
  
# and now analyzing time since fire
  
data <- data %>% filter(!is.na(time_since_fire)) # dropping one bald which was never burned, so has NA values
  
  germ.m <- glmer(cbind(TotalGerm, total_seeds - TotalGerm) ~ -1 + spp_code* live_sterile*time_since_fire + (1|soil_source), data = data, family = binomial,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # adding this to try to help convergence issues
  
  summary(germ.m)
  anova(germ.m, test = "Chisq")
  
  
  
  prediction_df <- expand.grid(spp_code = unique(germ.covariates$spp_code), live_sterile = unique(germ.covariates$live_sterile), time_since_fire = seq(from = min(data$time_since_fire), to = max(data$time_since_fire), by = .2))
  
  preds <- predict(germ.m, newdata = prediction_df, type = "link", se.fit = TRUE, re.form = NA)
  
  
  critval <- 1.96 ## approx 95% CI
  prediction_df$upr <- invlogit(preds$fit + (critval * preds$se.fit))
  prediction_df$lwr <- invlogit(preds$fit - (critval * preds$se.fit))
  prediction_df$fit <- invlogit(preds$fit)
  
  
  # now we can plot the model predictions
  
  ggplot(data = prediction_df)+
    geom_point(data = germ.covariates, aes( x= time_since_fire, y = seed_prop, color = live_sterile), alpha = .2)+
    geom_ribbon(aes(x = time_since_fire, ymin = lwr, ymax = upr, group = live_sterile, fill = live_sterile), alpha = .3)+
    geom_line(aes(x = time_since_fire, y = fit, group = live_sterile, color = live_sterile)) +
    scale_fill_manual(values = c(""))
  facet_wrap(~spp_code, scales = "free_y") + labs(y = "Proportion Germinated") + theme_minimal()
  

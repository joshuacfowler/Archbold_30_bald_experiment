####################################################################################
# Purpose: Visualize vital rate and population results
# Author: Joshua Fowler
# Date: Dec 31, 2025
####################################################################################
##### Set up #####

library(tidyverse)
library(readxl)
library(ggmap)
library(sf)



####################################################################################
####### reading in greenhouse vital rate and population level microbial effects data #######
####################################################################################
germ_30bald_summary <- read_csv("germ_30bald_predictions.csv") %>% 
  group_by(spp_code, bald, rel_elev, time_since_fire) %>% 
  summarize(rel_diff_mean = mean(rel_diff))

grow_30bald_summary <- read_csv("grow_30bald_predictions.csv") %>% 
  group_by(spp_code, bald, rel_elev, time_since_fire) %>% 
  summarize(rel_diff_mean = mean(rel_diff))

flw_30bald_summary <- read_csv("flw_30bald_predictions.csv") %>% 
  group_by(spp_code, bald, rel_elev, time_since_fire) %>% 
  summarize(rel_diff_mean = mean(rel_diff))



ERYCUN_lambdas <- read_csv("ERYCUN_bald_lambdas.csv") %>% 
  mutate(spp_code = "ERYCUN")
PARCHA_lambdas <- read_csv("PARCHA_bald_lambdas.csv") %>% 
  mutate(spp_code = "PARCHA")




####################################################################################
####### covariate and map info #######
####################################################################################
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



bald_coords <- read_csv("~/Dropbox/UofMiami/Experiment Set up/ABS-bald coordinates_UTM + latlong.csv")  %>% 
  mutate(bald = case_when(bald_code == "02"~ "2",
                          bald_code == "01N"~"1N",
                          bald_code == "07N"~"7N",
                          bald_code == "05E"~"5E",
                          bald_code == "35N"~"35",
                          bald_code == "65E"~"65",
                          TRUE~bald_code))


center_point <-  c(mean(bald_coords$longitude), mean(bald_coords$latitude))

# register_google()
archbold_satellite <- get_googlemap(center = center_point, zoom = 12, source = 'google',
                                    maptype = "satellite")


# map of balds
Archbold_bald_map <- ggmap(archbold_satellite)+
  geom_point(data = bald_coords, aes(x = longitude, y = latitude),fill = "grey",  size = 2, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs( x = "Longitude", y = "Latitude", fill = "Microbial Effect")


Archbold_bald_map


####################################################################################
####### putting coordinates on predictions #######
####################################################################################

germ_bald_pred_coords <- germ_30bald_summary %>% 
  left_join(bald_coords, by = join_by(bald == bald  )) %>% 
  filter(!is.na(rel_elev))

grow_bald_pred_coords <- grow_30bald_summary %>% 
  left_join(bald_coords, by = join_by(bald == bald  )) %>% 
  filter(!is.na(rel_elev))

flw_bald_pred_coords <- flw_30bald_summary %>% 
  left_join(bald_coords, by = join_by(bald == bald  )) %>% 
  filter(!is.na(rel_elev))





####################################################################################
####### plotting predictions across the landscape #######
####################################################################################

##### Germination #####

germ_map.ERYCUN <- ggmap(archbold_satellite)+
  # geom_point(data = germ_bald_pred_coords, aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  geom_point(data = filter(germ_bald_pred_coords, spp_code == "ERYCUN"), aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 2.5, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs( x = "Longitude", y = "Latitude", fill = "Microbial Effect")

# germ_map.ERYCUN
ggsave(germ_map.ERYCUN, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/germ_map_ERYCUN.png", width = 5, height = 5)


germ_map.PARCHA <- ggmap(archbold_satellite)+
  # geom_point(data = germ_bald_pred_coords, aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  geom_point(data = filter(germ_bald_pred_coords, spp_code == "PARCHA"), aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 2.5, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs( x = "Longitude", y = "Latitude", fill = "Microbial Effect")

# germ_map.PARCHA
ggsave(germ_map.PARCHA, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/germ_map_PARCHA.png", width = 5, height = 5)


##### Growth #####

grow_map.ERYCUN <- ggmap(archbold_satellite)+
  geom_point(data =  filter(grow_bald_pred_coords, spp_code == "ERYCUN"), aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs(title = "Growth", x = "Longitude", y = "Latitude", fill = "Microbial Effect")

grow_map.ERYCUN
ggsave(grow_map.ERYCUN, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/grow_map_ERYCUN.png", width = 5, height = 6)


grow_map.PARCHA <- ggmap(archbold_satellite)+
  geom_point(data =  filter(grow_bald_pred_coords, spp_code == "PARCHA"), aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs(title = "Growth", x = "Longitude", y = "Latitude", fill = "Microbial Effect")

grow_map.PARCHA
ggsave(grow_map.PARCHA, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/grow_map_PARCHA.png", width = 5, height = 6)

##### Flowering #####


flw_map.ERYCUN <- ggmap(archbold_satellite)+
  geom_point(data = filter(flw_bald_pred_coords, spp_code == "ERYCUN"), aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs(title = "Flowering Probability", x = "Longitude", y = "Latitude", fill = "Microbial Effect")

flw_map.ERYCUN
ggsave(flw_map.ERYCUN, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/flw_map_ERYCUN.png", width = 5, height = 6)

flw_map.PARCHA <- ggmap(archbold_satellite)+
  geom_point(data = filter(flw_bald_pred_coords, spp_code == "PARCHA"), aes(x = longitude, y = latitude, fill = -rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs(title = "Flowering Probability", x = "Longitude", y = "Latitude", fill = "Microbial Effect")

flw_map.PARCHA
ggsave(flw_map.PARCHA, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/flw_map_PARCHA.png", width = 5, height = 6)



####################################################################################
####### now plotting population growth rates across the landscape   #######
####################################################################################

ERYCUN_bald_lambda <- read_csv("ERYCUN_bald_lambdas.csv")  %>% 
  filter(!is.na(Bald)) %>% 
  mutate(Bald = as.character(Bald)) %>% 
  mutate(bald = case_when(Bald == "45" ~ "45N", # need to double check that these are the right assignments
                          Bald == "1" ~ "1N",
                          Bald == "95" ~ "95S",
                          Bald == "62" ~ "62N",
                          Bald == "46" ~ "46E",
                          Bald == "41" ~ "41S",
                          Bald == "46" ~ "46E",
                          Bald == "21" ~ "21S",
                          Bald == "24" ~ "24N",
                          Bald == "29" ~ "29W",
                          Bald == "5" ~ "5E", TRUE ~ Bald)) %>% 
  left_join(preddata) %>% 
  left_join(bald_coords) %>% 
  mutate(spp_code = "ERYCUN")


PARCHA_bald_lambda <- read_csv("PARCHA_bald_lambdas.csv") %>% 
  filter(!is.na(Bald)) %>% 
  mutate(Bald = as.character(Bald)) %>% 
  mutate(bald = case_when(Bald == "45" ~ "45N", # need to double check that these are the right assignments
                          Bald == "1" ~ "1N",
                          Bald == "95" ~ "95S",
                          Bald == "62" ~ "62N",
                          Bald == "46" ~ "46E",
                          Bald == "41" ~ "41S",
                          Bald == "46" ~ "46E",
                          Bald == "21" ~ "21S",
                          Bald == "24" ~ "24N",
                          Bald == "29" ~ "29W",
                          Bald == "5" ~ "5E", TRUE ~ Bald)) %>% 
  left_join(preddata) %>% 
  left_join(bald_coords) %>% 
  mutate(spp_code = "PARCHA")

lambda_df <- bind_rows(ERYCUN_bald_lambda, PARCHA_bald_lambda)

lambda_summary <- lambda_df %>% 
  filter(!is.na(lambda)) %>% 
  group_by(spp_code, Bald, bald, Microbe, rel_elev, time_since_fire, longitude, latitude) %>% 
  summarize(lambda_mean = mean(lambda),
            lambda_97.5 = quantile(lambda, 0.975),
            lambda_02.5 = quantile(lambda, 0.025)) %>% 
  filter(!(spp_code == "PARCHA" & lambda_mean>1))

my_palette <- c("#000000", "#009E73")

ERYCUN_lambda_elev <- ggplot(filter(lambda_summary, spp_code == "ERYCUN"))+
  geom_hline(aes(yintercept = 1))+
  geom_smooth(aes(x = rel_elev, y = lambda_mean-1.1, color = Microbe), method = "lm", se = F)+
  geom_linerange(aes(x = rel_elev, ymin = lambda_02.5-1.1, ymax = lambda_97.5-1.1, color = Microbe), alpha = .4)+
  # geom_point(data = lambda_lambda_covariates, aes( x = rel_elev , y = lambda, color = Microbe), alpha = .03)+
  geom_point(aes(x = rel_elev, y = lambda_mean-1.1, color = Microbe), size = 2, alpha = .4)+
  scale_color_manual(values = c(my_palette[2], my_palette[1]))+
  scale_fill_manual(values = c(my_palette[2], my_palette[1]))+
  labs(y = expression(lambda),  x = "Rel. Elev. (m)", color = "Microbiome") + theme_light()
ERYCUN_lambda_elev
ggsave(ERYCUN_lambda_elev, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/ERYCUN_lambda_elev.png", width = 4, height = 4)

ERYCUN_lambda_fire <- ggplot(filter(lambda_summary, spp_code == "ERYCUN"))+
  geom_hline(aes(yintercept = 1))+
  geom_smooth(aes(x = time_since_fire, y = lambda_mean-1.1, color = Microbe), method = "lm", se = F)+
  geom_linerange(aes(x = time_since_fire, ymin = lambda_02.5-1.1, ymax = lambda_97.5-1.1, color = Microbe), alpha = .4)+
  # geom_point(data = lambda_lambda_covariates, aes( x = time_since_fire , y = lambda, color = Microbe), alpha = .03)+
  geom_point(aes(x = time_since_fire, y = lambda_mean-1.1, color = Microbe), size = 2, alpha = .4)+
  scale_color_manual(values = c(my_palette[2], my_palette[1]))+
  scale_fill_manual(values = c(my_palette[2], my_palette[1]))+
  labs(y = expression(lambda),  x = "Time Since Fire (years)", color = "Microbiome") + theme_light()
ERYCUN_lambda_fire
ggsave(ERYCUN_lambda_fire, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/ERYCUN_lambda_fire.png", width = 4, height = 4)





PARCHA_lambda_elev <- ggplot(filter(lambda_summary, spp_code == "PARCHA"))+
  geom_hline(aes(yintercept = 1))+
  geom_smooth(aes(x = rel_elev, y = lambda_mean, color = Microbe), method = "loess", se = F)+
  geom_linerange(aes(x = rel_elev, ymin = lambda_02.5, ymax = lambda_97.5, color = Microbe), alpha = .4)+
  # geom_point(data = lambda_lambda_covariates, aes( x = rel_elev , y = lambda, color = Microbe), alpha = .03)+
  geom_point(aes(x = rel_elev, y = lambda_mean, color = Microbe), size = 2, alpha = .4)+
  scale_color_manual(values = c(my_palette[2], my_palette[1]))+
  scale_fill_manual(values = c(my_palette[2], my_palette[1]))+
  labs(y = expression(lambda),  x = "Rel. Elev. (m)", color = "Microbiome") + theme_light()
PARCHA_lambda_elev
ggsave(PARCHA_lambda_elev, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/PARCHA_lambda_elev.png", width = 4, height = 4)

PARCHA_lambda_fire <- ggplot(filter(lambda_summary, spp_code == "PARCHA"))+
  geom_hline(aes(yintercept = 1))+
  geom_smooth(aes(x = time_since_fire, y = lambda_mean, color = Microbe), method = "loess", se = F)+
  geom_linerange(aes(x = time_since_fire, ymin = lambda_02.5, ymax = lambda_97.5, color = Microbe), alpha = .4)+
  # geom_point(data = lambda_lambda_covariates, aes( x = time_since_fire , y = lambda, color = Microbe), alpha = .03)+
  geom_point(aes(x = time_since_fire, y = lambda_mean, color = Microbe), size = 2, alpha = .4)+
  scale_color_manual(values = c(my_palette[2], my_palette[1]))+
  scale_fill_manual(values = c(my_palette[2], my_palette[1]))+
  labs(y = expression(lambda),  x = "Time Since Fire (years)", color = "Microbiome") + theme_light()
PARCHA_lambda_fire
ggsave(PARCHA_lambda_fire, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/PARCHA_lambda_fire.png", width = 4, height = 4)





####### Making Maps of lambda #####

lambda_diff <- lambda_df %>% 
  pivot_wider( names_from  = c(Microbe), values_from = c(lambda)) %>% 
  mutate(rel_diff = ((live-sterile)/sterile))

lambda_diff_summary <- lambda_diff %>% 
  group_by(spp_code, Bald, longitude, latitude) %>% 
  summarize(rel_diff_mean = mean(rel_diff),
            diff_02.5 = quantile(rel_diff, probs = .025),
            diff_97.5 = quantile(rel_diff, probs = .975)) 
  # rename(bald = Bald) %>% 
  # filter(bald != "59")


lambda_map.ERYCUN <- ggmap(archbold_satellite)+
  geom_point(data = filter(lambda_diff_summary, spp_code == "ERYCUN"), aes(x = longitude, y = latitude, fill = rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs(title = "Population Growth", x = "Longitude", y = "Latitude", fill = "Microbial Effect")

lambda_map.ERYCUN
ggsave(lambda_map.ERYCUN, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/lambda_map_ERYCUN.png", width = 5, height = 6)




lambda_map.PARCHA <- ggmap(archbold_satellite)+
  geom_point(data = filter(lambda_diff_summary, spp_code == "PARCHA"), aes(x = longitude, y = latitude, fill = rel_diff_mean),  size = 3, shape = 21, alpha = .8)+
  coord_sf(xlim = c(min(bald_coords$longitude)-.01,max(bald_coords$longitude) + .01), ylim = c(min(bald_coords$latitude)-.01,max(bald_coords$latitude) + .01))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  theme_light()+
  guides(size = "none")+
  labs(title = "Population Growth", x = "Longitude", y = "Latitude", fill = "Microbial Effect")

lambda_map.PARCHA
ggsave(lambda_map.PARCHA, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/lambda_map_PARCHA.png", width = 5, height = 6)




####################################################################################
####### trying to calculate metapopulation capacity   #######
####################################################################################

# trying to calculate metapopulation persistence:
bald_areas<- read_xlsx(path = "~/Dropbox/UofMiami/Balds2009_FireIntensityArea_Through2022.xlsx", sheet = "abs_RosemaryBalds2009", guess_max = 1048576) %>% # guess_max makes the function look deeper in the columns to assign type
  select(BALD_, `Bald_U...5`, Area) %>% 
  rename(Bald = Bald_U...5)



lambda_meta <- lambda_summary %>% 
  left_join(bald_areas, by = join_by(bald == Bald)) %>% 
  # mutate(bald_code = as.character(BALD_)) %>% 
  ungroup() %>% 
  mutate(lambda_overall_mean = mean(lambda_mean, na.rm = T),
         patch_quality = lambda_mean/lambda_overall_mean) %>% 
  filter(!is.na(longitude))

# ERYCUN
lambda_meta_live.ERYCUN <- data.frame(latitude = filter(lambda_meta, Microbe == "live" & spp_code == "ERYCUN")$latitude,
                               longitude = filter(lambda_meta, Microbe == "live"& spp_code == "ERYCUN")$longitude)
rownames(lambda_meta_live.ERYCUN) <- filter(lambda_meta, Microbe == "live"& spp_code == "ERYCUN")$Bald


lambda_meta_sterile.ERYCUN <- data.frame(latitude = filter(lambda_meta, Microbe == "sterile" & spp_code == "ERYCUN")$latitude,
                                  longitude = filter(lambda_meta, Microbe == "sterile" & spp_code == "ERYCUN")$longitude)
rownames(lambda_meta_sterile.ERYCUN) <- filter(lambda_meta, Microbe == "sterile" & spp_code == "ERYCUN")$Bald

patch_quality_live.ERYCUN <- filter(lambda_meta, Microbe == "live"& spp_code == "ERYCUN")$patch_quality
patch_quality_sterile.ERYCUN <-  filter(lambda_meta, Microbe == "sterile"& spp_code == "ERYCUN")$patch_quality

patch_area_live.ERYCUN  <- filter(lambda_meta, Microbe == "live"& spp_code == "ERYCUN")$Area
patch_area_sterile.ERYCUN  <- filter(lambda_meta, Microbe == "sterile"& spp_code == "ERYCUN")$Area

distance_live.ERYCUN  <- dist(lambda_meta_live.ERYCUN, diag=T, upper=T)
distance_sterile.ERYCUN  <- dist(lambda_meta_sterile.ERYCUN, diag=T, upper=T)

# PARCHA
lambda_meta_live.PARCHA <- data.frame(latitude = filter(lambda_meta, Microbe == "live" & spp_code == "PARCHA")$latitude,
                                      longitude = filter(lambda_meta, Microbe == "live"& spp_code == "PARCHA")$longitude)
rownames(lambda_meta_live.PARCHA) <- filter(lambda_meta, Microbe == "live"& spp_code == "PARCHA")$Bald


lambda_meta_sterile.PARCHA <- data.frame(latitude = filter(lambda_meta, Microbe == "sterile" & spp_code == "PARCHA")$latitude,
                                         longitude = filter(lambda_meta, Microbe == "sterile" & spp_code == "PARCHA")$longitude)
rownames(lambda_meta_sterile.PARCHA) <- filter(lambda_meta, Microbe == "sterile" & spp_code == "PARCHA")$Bald

patch_quality_live.PARCHA <- filter(lambda_meta, Microbe == "live"& spp_code == "PARCHA")$patch_quality
patch_quality_sterile.PARCHA <-  filter(lambda_meta, Microbe == "sterile"& spp_code == "PARCHA")$patch_quality

patch_area_live.PARCHA  <- filter(lambda_meta, Microbe == "live"& spp_code == "PARCHA")$Area
patch_area_sterile.PARCHA  <- filter(lambda_meta, Microbe == "sterile"& spp_code == "PARCHA")$Area

distance_live.PARCHA  <- dist(lambda_meta_live.PARCHA, diag=T, upper=T)
distance_sterile.PARCHA  <- dist(lambda_meta_sterile.PARCHA, diag=T, upper=T)




# calculating landscape matrices
M_live.ERYCUN <-  as.matrix(exp(-2*distance_live.ERYCUN)*(patch_quality_live.ERYCUN*patch_area_live.ERYCUN))
M_sterile.ERYCUN <-  exp(-2*distance_sterile.ERYCUN)*patch_quality_sterile.ERYCUN*patch_area_sterile.ERYCUN

capacity_live.ERYCUN <- Re(eigen(M_live.ERYCUN)$values[1])
capacity_sterile.ERYCUN <- Re(eigen(M_sterile.ERYCUN)$values[1])

percent_change_capacity.ERYCUN <- (capacity_live.ERYCUN-capacity_sterile.ERYCUN)/capacity_sterile.ERYCUN*100
ratio_capacity.ERYCUN <- (capacity_live.ERYCUN/capacity_sterile.ERYCUN)

# parcha
M_live.PARCHA <-  as.matrix(exp(-2*distance_live.PARCHA)*(patch_quality_live.PARCHA*patch_area_live.PARCHA))
M_sterile.PARCHA <-  exp(-2*distance_sterile.PARCHA)*patch_quality_sterile.PARCHA*patch_area_sterile.PARCHA

capacity_live.PARCHA <- Re(eigen(M_live.PARCHA)$values[1])
capacity_sterile.PARCHA <- Re(eigen(M_sterile.PARCHA)$values[1])

percent_change_capacity.PARCHA <- (capacity_live.PARCHA-capacity_sterile.PARCHA)/capacity_sterile.PARCHA*100
ratio_capacity.PARCHA <- (capacity_live.PARCHA/capacity_sterile.PARCHA)





library(ggraph)
library(igraph)
nearest_adj <- distance_live.ERYCUN
distance_df <- reshape2::melt(as.matrix(distance_live.ERYCUN), varnames = c("from", "to")) %>% 
  group_by(from) %>% 
  mutate(first_value = sort(value)[2],
         second_value = sort(value)[3],
         third_value = sort(value)[4],
         fourth_value = sort(value)[5],
         fifth_value = sort(value)[6]) %>% 
  filter(value %in% c(first_value, second_value, third_value, fourth_value, fifth_value)) %>% select(-first_value, -second_value, -third_value, -fourth_value, -fifth_value)


live_metapop_graph <- graph_from_data_frame(distance_df)
sterile_metapop_graph <- graph_from_data_frame(distance_df)

lo.live <- norm_coords(as.matrix(lambda_meta_live.ERYCUN[, c("longitude","latitude")]),
                  xmin = min(lambda_meta_live.ERYCUN$longitude)-.01,
                  xmax = max(lambda_meta_live.ERYCUN$longitude)+.01,
                  ymin = min(lambda_meta_live.ERYCUN$latitude)-.01,
                  ymax = max(lambda_meta_live.ERYCUN$latitude)+.01)

lo.sterile <- norm_coords(as.matrix(lambda_meta_sterile.ERYCUN[, c("longitude","latitude")]),
                       xmin = min(lambda_meta_sterile.ERYCUN$longitude)-.01,
                       xmax = max(lambda_meta_sterile.ERYCUN$longitude)+.01,
                       ymin = min(lambda_meta_sterile.ERYCUN$latitude)-.01,
                       ymax = max(lambda_meta_sterile.ERYCUN$latitude)+.01)


V(live_metapop_graph)$x<-lo.live[,1]
V(live_metapop_graph)$y<-lo.live[,2]
V(live_metapop_graph)$quality   = patch_quality_live.ERYCUN #filter(lambda_diff_summary, spp_code == "ERYCUN",!is.na(longitude))$rel_diff_mean


V(sterile_metapop_graph)$x<-lo.sterile[,1]
V(sterile_metapop_graph)$y<-lo.sterile[,2]
V(sterile_metapop_graph)$quality   = patch_quality_sterile.ERYCUN 




live_metapop_plot <- ggraph(live_metapop_graph, layout = lo.live, rescale = FALSE)+
  geom_edge_link(color = "grey", alpha = .3) + 
  geom_node_point(aes(fill = quality), size = 3, shape = 21, alpha = .6, show.legend = FALSE)+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  coord_sf()+
  # coord_sf()+
  theme_void()
live_metapop_plot

ggsave(live_metapop_plot, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/ERYCUN_live_metapop_graph.png",  width = 5, height = 6)

sterile_metapop_plot <- ggraph(sterile_metapop_graph, layout = lo.sterile, rescale = FALSE)+
  geom_edge_link(color = "grey", alpha = .3) + 
  geom_node_point(aes(fill = quality), size = 3, shape = 21, alpha = .6, show.legend = FALSE)+
  scale_fill_gradient(low = "white", high = "grey25")+
  coord_sf()+
  # coord_sf()+
  theme_void()
sterile_metapop_plot
ggsave(sterile_metapop_plot, filename = "~/Documents/R_projects/Archbold_30_bald_experiment/ERYCUN_sterile_metapop_graph.png",  width = 5, height = 6)





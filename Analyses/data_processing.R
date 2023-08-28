####################################################################################
# Purpose: Read in Archbold demographic data from eleven species, clean and process data to have consistent format for vital rate analyses
# Author: Joshua Fowler
# Date: Aug 16, 2023
####################################################################################

##### Set up #####

library(renv) # track package versions
# renv::snapshot()
# renv::restore()


library(tidyverse)
library(readxl)

####################################################################################
###### Loading in each data for each species #######################################
####################################################################################

## Filepath for stored demographic data
filepath <- c("/Users/joshuacfowler/Dropbox/UofMiami/Demographic Data")

# Reading in data for Eryngium cuneifolium (ERYCUN)
# for ERYCUN, there is demographic data from rosemary balds at Archbold as well as at Royce Ranch (one dataset from a firelane transect and another from .
ERYCUN_abs_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecabs") %>% 
  mutate(Site_tag = "archbold")

ERYCUN_apfi_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecapfi") %>% 
  mutate(Site_tag = "royce_ranch")

ERYCUN_apsc_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecapsc") %>% 
  mutate(Site_tag = "royce_ranch")






####################################################################################
###### Cleaning and merging together ERYCUN #######################################
####################################################################################

#### Loading in the data from the main Archbold sites
# removing columns related to XY coordinates only applicable in some plots as well as outdated info on burn history in the sites. We will merge in the most updated burn history later
# removing redundant demographic information columns (several contain the change in flowering stem number or rosette diameter between years)
# removing columns containing ground cover information about other species in the plots, cool information, but seems to mostly have been collected before 1994
# removing columns containing record of sand accretion on the metal id tag %>% This might be worth including as a covariate, but I don't know if it really tells us much and some of it is after the plant is already dead.
# Making a unique id for each plant in each patch in each bald
# Then pivoting the data to long format so that the demographic measures are in combined columns for each row
# We will also remove the assigned stage information from the dataset since we will work with the raw size measurements
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# There are a few instances where tags were not found during a census but relocated in the following year or the second year after. Recoding these to be alive during those censuses, but we will obviously be missing size information during those years.



ERYCUN_abs <- ERYCUN_abs_raw %>%
  dplyr::select(-gps18,-ltreb,-pull,-X,-Y,-x.cor,-y.cor,-xy.cor.notes, 
                -contains("byr2"), -tsf, -breg, -sdno95, -stat, pull98, -cohort, -oldX, -oldtag, -contains("burn2"), -hobo0710, -rx0710, -starts_with("rx"), -tsf2018) %>% 
  dplyr::select(-contains("cs9"), -contains("cr9"),-contains("cr0"), 
                -contains("agr0"), -contains("agr9"), -contains("annsur9"), -contains("annsur0"), -contains("annsur1"),
                -contains("hstg9"), -contains("age9")) %>% 
  dplyr::select(-ag, -ca, -cev, -cp, -cs, -ec, -hc, -ld, -lc, -pc, -pr, -pb, -sab, -sel, -qi, -qu, -licania, -ceratiol, -groundco, -perlitte, -sppshrub, -distshru, -oakht02, -oakdis02, -quad, -otherspp) %>% 
  dplyr::select(-starts_with("sa"), -contains("pull"), -master) %>% 
  mutate(row_id = row_number()) %>% 
  mutate(plant_id = paste(bald,patch,plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename(first_year = yr1) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, bald, patch, plant, TP, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                 "s" ~ "ARCHBOLD_surv",
                                 "stg" ~ "assigned_stage",
                                 "rdm" ~ "ros_diameter",
                                 "stm" ~ "flw_stem", "fl" ~ "flw_stem",
                                 "ht" ~ "max_stem_height",
                                 "hstm" ~ "herb_count",
                                 "h" ~ "flw_head",
                                 "sa" ~ "sand_accretion",
                                 "comm" ~ "comment",
                                 .default = as.character(measurement)),
         census_year = case_match(as.numeric(census_year),
                                  22 ~ 2022,
                                  21 ~ 2021,
                                  20 ~ 2020,
                                  19 ~ 2019,
                                  18 ~ 2018,
                                  17 ~ 2017,
                                  16 ~ 2016,
                                  15 ~ 2015,
                                  14 ~ 2014,
                                  13 ~ 2013,
                                  12 ~ 2012, 2012 ~ 2012,
                                  11 ~ 2011, 2011 ~ 2011,
                                  10 ~ 2010, 2010 ~ 2010,
                                  09 ~ 2009,
                                  08 ~ 2008,
                                  07 ~ 2007,
                                  06 ~ 2006,
                                  05 ~ 2005,
                                  04 ~ 2004,
                                  03 ~ 2003,
                                  02 ~ 2002,
                                  01 ~ 2001,
                                  00 ~ 2000,
                                  99 ~ 1999,
                                  98 ~ 1998,
                                  97 ~ 1997,
                                  96 ~ 1996,
                                  95 ~ 1995,
                                  94 ~ 1994,
                                  93 ~ 1993,
                                  92 ~ 1992,
                                  91 ~ 1991,
                                  90 ~ 1990,
                                  89 ~ 1989,
                                  88 ~ 1988,
                                  .default = as.numeric(census_year))) %>% 
  pivot_wider(id_cols = c(plant_id, Site_tag, bald, patch, plant, TP, row_id, first_year, census_year), names_from = measurement, values_from = value) %>% 
  mutate(ARCHBOLD_surv = case_when( ARCHBOLD_surv == 20 ~ 2, TRUE ~ as.numeric(ARCHBOLD_surv)),
         surv = case_match(as.numeric(ARCHBOLD_surv),
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ 1, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  group_by(plant_id) %>%
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          TRUE ~ surv)) %>% 
  mutate(surv = case_when(dplyr::lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 2 ~ 1,
                          dplyr::lag(ARCHBOLD_surv, n = 2) == 1 & dplyr::lag(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ 1,
                          TRUE ~ surv))

  
  
#### Loading in the fire lane data 
# Removing some extraneous location info
# removing redundant demographic information columns (columns with change in size and in number of stems)
# creating unique id for each plant
# Then pivoting to combine the columns for each vital rate measurement
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# For this dataset, there are several instances, particulary recruits in the 90s, where new seedlings have a size when they appear but NA survival, so I am correcting this
# And correcting cases where tag was not found  in multiple years. Assuming plant is dead, unless it is re-found in the following year or the next

ERYCUN_apfi <- ERYCUN_apfi_raw %>%
  dplyr::select(-pull, -quad) %>% 
  dplyr::select(-contains("rgr"), -contains("chst")) %>% 
  mutate(row_id = row_number(), patch = "Firelane") %>%
  mutate(plant_id = paste(patch, plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename(first_year = yr1) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, patch, plant, TP, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                "s" ~ "ARCHBOLD_surv", "a" ~ "ARCHBOLD_surv",
                                "stg" ~ "assigned_stage",
                                "rdm" ~ "ros_diameter", "rd" ~ "ros_diameter", "bdm" ~ "ros_diameter",
                                "stm" ~ "flw_stem", "fl" ~ "flw_stem", "st" ~ "flw_stem",
                                "ht" ~ "max_stem_height",
                                "hstm" ~ "herb_count",
                                "h" ~ "flw_head", "hd" ~ "flw_head", "he" ~ "flw_head",
                                "sa" ~ "sand_accretion",
                                "comm" ~ "comment",
                                .default = as.character(measurement)),
       census_year = case_match(as.numeric(census_year),
                                22 ~ 2022,
                                21 ~ 2021,
                                20 ~ 2020,
                                19 ~ 2019,
                                18 ~ 2018,
                                17 ~ 2017,
                                16 ~ 2016,
                                15 ~ 2015,
                                14 ~ 2014,
                                13 ~ 2013, 2013 ~ 2013,
                                12 ~ 2012, 2012 ~ 2012,
                                11 ~ 2011, 2011 ~ 2011,
                                10 ~ 2010, 2010 ~ 2010,
                                09 ~ 2009,
                                08 ~ 2008,
                                07 ~ 2007,
                                06 ~ 2006,
                                05 ~ 2005,
                                04 ~ 2004,
                                03 ~ 2003,
                                02 ~ 2002,
                                01 ~ 2001,
                                00 ~ 2000,
                                99 ~ 1999,
                                98 ~ 1998,
                                97 ~ 1997,
                                96 ~ 1996,
                                95 ~ 1995,
                                94 ~ 1994,
                                93 ~ 1993,
                                92 ~ 1992,
                                91 ~ 1991,
                                90 ~ 1990,
                                89 ~ 1989,
                                88 ~ 1988,
                                .default = as.numeric(census_year))) %>% 
  pivot_wider(id_cols = c(plant_id, Site_tag, patch, plant, TP, row_id, first_year, census_year), names_from = measurement, values_from = value) %>% 
  mutate(ARCHBOLD_surv = case_when( ARCHBOLD_surv == 20 ~ 2, TRUE ~ as.numeric(ARCHBOLD_surv)),
         surv = case_match(as.numeric(ARCHBOLD_surv),
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ 1, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  mutate(surv = case_when(is.na(surv) & !is.na(ros_diameter) & first_year == census_year ~ 1,
                          TRUE ~ surv)) %>% 
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 9 ~ 0,
                          TRUE ~ surv)) %>%
  # filter(ARCHBOLD_surv == 1 & lag(ARCHBOLD_surv) == 9) %>%
  select(plant_id, ARCHBOLD_surv, first_year, census_year, surv, ros_diameter)
  
  
  


  



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

# removing columns related to XY coordinates only applicable in some plots as well as outdated info on burn history in the sites. We will merge in the most updated burn history later
# removing redundant demographic information columns (several contain the change in flowering stem number or rosette diameter between years)
# removing columns containing ground cover information about other species in the plots, cool information, but seems to mostly have been collected before 1994
# removing columns containing record of sand accretion on the metal id tag %>% 
# Making a unique id for each plant in each patch in each bald
# Then pivoting the data to long format so that the demographic measures are in combined columns for each row
# We will also remove the assigned stage information from the dataset since we will work with the raw size measurements



ERYCUN_abs <- ERYCUN_abs_raw %>%
  dplyr::select(-gps18,-ltreb,-pull,-X,-Y,-x.cor,-y.cor,-xy.cor.notes, 
                -contains("byr2"), -tsf, -breg, -sdno95, -stat, pull98, -cohort, -oldX, -oldtag, -contains("burn2"), -hobo0710, -rx0710, -tsf2018) %>% 
  dplyr::select(-contains("cs9"), -contains("cr9"),-contains("cr0"), 
                -contains("agr0"), -contains("agr9"), -contains("annsur9"), -contains("annsur0"), -contains("annsur1"),
                -contains("hstg9"), -contains("age9")) %>% 
  dplyr::select(-ag, -ca, -cev, -cp, -cs, -ec, -hc, -ld, -lc, -pc, -pr, -pb, -sab, -sel, -qi, -qu, -licania, -ceratiol, -groundco, -perlitte, -sppshrub, -distshru, -oakht02, -oakdis02, -quad, -otherspp) %>% 
  dplyr::select(-contains("sa9"), -contains("sa0"), -contains("sa1"), -contains("pull"), -master) %>% 
  mutate(row_id = row_number()) %>% 
  mutate(plant_id = paste(bald,patch,plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename(first_year = yr1) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, bald, patch, plant, TP, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  dplyr::group_by(plant_id, Site_tag, bald, patch, plant, TP, row_id, census_year, measurement) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
  
  
  
  mutate(measurement = case_when(measurement == "s" ~ "surv",
                                 measurement == "stg" ~ "assigned_stage",
                                 measurement == "rdm" ~ "ros_diameter",
                                 measurement == "stm" ~ "stem_count",
                                 measurement == "ht" ~ "height",
                                 measurement == "hstm" ~ "herb_count",
                                 measurement == "h" ~ "flw_head",
                                 measurement == "sa" ~ "sand_accretion",
                                 measurement == "comm" ~ "comment",),
         census_year = case_when(census_year == 22 ~ 2022,
                                 census_year == 21 ~ 2021,
                                 census_year == 20 ~ 2020,
                                 census_year == 19 ~ 2019,
                                 census_year == 18 ~ 2018,
                                 census_year == 17 ~ 2017,
                                 census_year == 16 ~ 2016,
                                 census_year == 15 ~ 2015,
                                 census_year == 14 ~ 2014,
                                 census_year == 13 ~ 2013,
                                 census_year == 12 ~ 2012,census_year == 2012 ~ 2012,
                                 census_year == 11 ~ 2011,census_year == 2011 ~ 2011,
                                 census_year == 10 ~ 2010,census_year == 2010 ~ 2010,
                                 census_year == 09 ~ 2009,
                                 census_year == 08 ~ 2008,
                                 census_year == 07 ~ 2007,
                                 census_year == 06 ~ 2006,
                                 census_year == 05 ~ 2005,
                                 census_year == 04 ~ 2004,
                                 census_year == 03 ~ 2003,
                                 census_year == 02 ~ 2002,
                                 census_year == 01 ~ 2001,
                                 census_year == 00 ~ 2000,
                                 census_year == 99 ~ 1999,
                                 census_year == 98 ~ 1998,
                                 census_year == 97 ~ 1997,
                                 census_year == 96 ~ 1996,
                                 census_year == 95 ~ 1995,
                                 census_year == 94 ~ 1994,
                                 census_year == 93 ~ 1993,
                                 census_year == 92 ~ 1992,
                                 census_year == 91 ~ 1991,
                                 census_year == 90 ~ 1990,
                                 census_year == 89 ~ 1989,
                                 census_year == 88 ~ 1988)) %>% 
  pivot_wider(id_cols = c(plant_id, Site_tag, bald, patch, plant, TP, row_id, census_year), names_from = measurement, values_from = value) %>% 











  # pivot_longer(cols = starts_with(c("s2","s0", "s9", "s8")), names_to = "surv_year", values_to = "surv")
  pivot_longer(cols = !c(plant_id, Site_tag, bald, patch, plant, TP, row_id)) %>% 
  # pivot_wider(id_cols = c(plant_id, Site_tag, bald, patch, plant, TP)) 
  # dplyr::group_by(plant_id, Site_tag, bald, patch, plant, TP, name) %>%
  # dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  # dplyr::filter(n > 1L)
  dplyr::select(plant_id, Site_tag, bald, patch, plant, TP, s22)
  
  



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


quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
####################################################################################
###### Loading in each data for each species #######################################
####################################################################################

## Filepath for stored demographic data
filepath <- c("/Users/joshuacfowler/Dropbox/UofMiami/Demographic Data")

# Reading in data for Eryngium cuneifolium (ERYCUN)
# for ERYCUN, there is demographic data from rosemary balds at Archbold as well as at Royce Ranch (one dataset from a firelane transect and another from .
ERYCUN_abs_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecabs", col_types = "text") %>% 
  mutate(Site_tag = "archbold")

ERYCUN_apfi_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecapfi", col_types = "text") %>% 
  mutate(Site_tag = "royce_ranch")

ERYCUN_apsc_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecapsc", col_types = "text") %>% 
  mutate(Site_tag = "royce_ranch")


# Reading in data for Liatris ohlingerae (LIAOHL)
LIAOHL_raw <- read_csv(file = paste0(filepath, "/UM demographic models/LoDem_2017.csv"))




# Reading in data for Balduina angustifolia (BALANG)
BALANG_raw <- read_excel(path = paste0(filepath, "/BethStephens_BALANG_CHAFAS/BaCf_data6-2-12.xls"), sheet = "Badata_may2012")

# Reading in data for Chamaecrista fasciculate (CHAFAS)
CHAFAS_raw <- read_excel(path = paste0(filepath, "/BethStephens_BALANG_CHAFAS/BaCf_data6-2-12.xls"), sheet = "Cfdata_may2012")



####################################################################################
###### Cleaning and merging together ERYCUN ########################################
####################################################################################
# Eryngium cuneifolium data was collected starting in 1994, and is some of the longest-term monitoring at Archbold, initiated by Eric Menges and Pedro Quintana-Ascencio

#### Loading in the data from the main Archbold sites
# removing columns related to XY coordinates only applicable in some plots as well as outdated info on burn history in the sites. We will merge in the most updated burn history later
# removing redundant demographic information columns (several contain the change in flowering stem number or rosette diameter between years)
# removing columns containing ground cover information about other species in the plots, cool information, but seems to mostly have been collected before 1994
# removing columns containing record of sand accretion on the metal id tag %>% This might be worth including as a covariate, but I don't know if it really tells us much and some of it is after the plant is already dead.
# Making a unique id for each plant in each patch in each bald
# Then pivoting the data to long format so that the demographic measures are in combined columns for each row
# We will also remove the assigned stage information from the dataset since we will work with the raw size measurements
# found one plant that had a typo in the spreadsheet, where the 2022 census was missing survival info and columns for size and reproduction where shifted one column to the left. so recoding this
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# There are a few instances where tags were not found during a census but relocated in the following year or the second year after. Recoding these to be alive during those censuses, but we will obviously be missing size information during those years.
# There is one plant which has outlier size. I think this is likely just a decimal point typo based on size of the plant in preceding years, so moving that data point
# most of the time when there are no reproductive structures, the cell are left blank. particularly for new seedlings,  with rosette measurements, I think this is safe to assume is 0 stems. There are often times cases where this data is actually missing however

ERYCUN_abs <- ERYCUN_abs_raw %>%
  dplyr::select(-gps18,-ltreb,-pull,-X,-Y,-x.cor,-y.cor,-xy.cor.notes, 
                -contains("byr2"), -tsf, -breg, -sdno95, -stat, pull98, -cohort, -oldX, -oldtag, -contains("burn2"), -hobo0710, -rx0710, -starts_with("rx"), -tsf2018) %>% 
  dplyr::select(-contains("cs9"), -contains("cr9"),-contains("cr0"), 
                -contains("agr0"), -contains("agr9"), -contains("annsur9"), -contains("annsur0"), -contains("annsur1"),
                -contains("hstg9"), -contains("age9")) %>% 
  dplyr::select(-ag, -ca, -cev, -cp, -cs, -ec, -hc, -ld, -lc, -pc, -pr, -pb, -sab, -sel, -qi, -qu, -licania, -ceratiol, -groundco, -perlitte, -sppshrub, -distshru, -oakht02, -oakdis02, -quad, -otherspp) %>% 
  dplyr::select(-starts_with("sa"), -contains("pull"), -master) %>% 
  dplyr::select(-starts_with("stg")) %>% 
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
  mutate(ARCHBOLD_surv = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022 ~ "1", TRUE ~ ARCHBOLD_surv),
         ros_diameter = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022  ~ "4.7", TRUE ~ ros_diameter),
         flw_stem = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022  ~ "1", TRUE ~ flw_stem),
         flw_head = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022  ~ "7", TRUE ~ flw_head)) %>% 
  mutate(across(c(plant_id, Site_tag, bald, patch, plant, TP, row_id, comment), as.character)) %>% 
  mutate(across(c(first_year, census_year, ARCHBOLD_surv, flw_stem, flw_head, herb_count), function(x) suppressWarnings(as.integer(x)))) %>% 
  mutate(across(c(ros_diameter, max_stem_height), function(x) suppressWarnings(as.numeric(x)))) %>% 
  mutate(ARCHBOLD_surv = case_when(ARCHBOLD_surv == 20 ~ 2, TRUE ~ ARCHBOLD_surv),
         surv = case_match(as.numeric(ARCHBOLD_surv),
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ NA, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census, True NA's not necessarily dead))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  group_by(plant_id) %>%
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 0 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 2 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 6 & ARCHBOLD_surv == 6 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 6 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 6 ~ NA,
                          
                          dplyr::lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 2 ~ 1,
                          dplyr::lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 0 ~ 1,
                          dplyr::lag(ARCHBOLD_surv, n = 2) == 1 & dplyr::lag(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ 1,
                          TRUE ~ surv)) %>% 
  mutate(ros_diameter = case_when(plant_id == "54_9_155_NA_4531" & census_year == 2018 ~ ros_diameter*.1, TRUE ~ ros_diameter)) %>% 
  mutate(flw_stem = case_when(first_year == census_year & !is.na(ros_diameter) & is.na(flw_stem) ~ 0, TRUE ~ flw_stem)) %>% 
  mutate(flw_stem = case_when(!is.na(ros_diameter) & is.na(flw_stem) & is.na(flw_head) ~ 0,
                              !is.na(ros_diameter) & is.na(max_stem_height) & is.na(flw_head) ~ 0, TRUE ~ flw_stem)) %>% 
  dplyr::select(plant_id,Site_tag,bald, patch, plant, TP, row_id, first_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count, comment)


#### Loading in the Royce Ranch fire lane data 
# Removing some extraneous location info
# removing redundant demographic information columns (columns with change in size and in number of stems)
# removing measurement of sand accretion
# creating unique id for each plant
# Then pivoting to combine the columns for each vital rate measurement
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# For this dataset, there are several instances, particulary recruits in the 90s, where new seedlings have a size when they appear but NA survival, so I am correcting this
# And correcting cases where tag was not found  in multiple years. Assuming plant is dead, unless it is re-found in the following year or the next

ERYCUN_apfi <- ERYCUN_apfi_raw %>%
  dplyr::select(-pull) %>% 
  dplyr::select(-contains("rgr"), -contains("chst")) %>% 
  dplyr::select(-starts_with("sa")) %>% 
  mutate(row_id = row_number(), patch = "Firelane") %>% 
  mutate(plant_id = paste(patch, quad, plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename(first_year = yr1) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, patch, quad, plant, TP, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
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
  pivot_wider(id_cols = c(plant_id, Site_tag, patch, quad, plant, TP, row_id, first_year, census_year), names_from = measurement, values_from = value) %>% 
  mutate(across(c(plant_id, Site_tag, patch, plant, TP, row_id, comment), as.character)) %>% 
  mutate(across(c(first_year, census_year, ARCHBOLD_surv, flw_stem, flw_head, herb_count), as.integer)) %>% 
  mutate(across(c(ros_diameter, max_stem_height), as.numeric)) %>% 
  mutate(surv = case_match(ARCHBOLD_surv,
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
  dplyr::select(plant_id,Site_tag, patch, quad, plant, TP, row_id, first_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count, comment)

  
#### Loading in the Royce Ranch Scrub data
# removing some extraneous meta columns
# removing redundant demographic information columns (several contain the change in flowering stem number or rosette diameter between years)
# removing measurement of sand accretion
# removing comments and record of burn/mowing
# Removing the assigned stage information from the dataset since we will work with the raw size measurements
# creating a unique id for each plant
# removing the column for the first census year for this species because in this dataset, the column not complete and I'll just calculate it later
# creating unique id for each plant
# Then pivoting to combine the columns for each vital rate measurement
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# recoding a two size measurements, which are marked as 999.0 to be NA.
# Then there are a few cases where plants where missing but found alive in later years. Correcting survival to reflect this
# And calculating the first observation year for each individual for this dataset
ERYCUN_apsc <- ERYCUN_apsc_raw %>%
  dplyr::select(-microhabitat, -pull, -pull_temp, -pull98, -oldtag, -bigq, -s18, -subquad) %>% 
  dplyr::select(-contains("annsur"), -contains("rgr"), -contains("chst")) %>% 
  dplyr::select(-starts_with("sa")) %>% 
  dplyr::select(-starts_with("comm"), -contains("mow"), -contains("burn")) %>% 
  dplyr::select(-starts_with("stg")) %>% 
  mutate(row_id = row_number(), patch = "Scrub") %>% 
  mutate(plant_id = paste(patch, quad, plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, patch, quad, plant, TP, row_id), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                  "s" ~ "ARCHBOLD_surv", "a" ~ "ARCHBOLD_surv",
                                  "stg" ~ "assigned_stage",
                                  "rdm" ~ "ros_diameter", "rd" ~ "ros_diameter", "bdm" ~ "ros_diameter",
                                  "stm" ~ "flw_stem", "fl" ~ "flw_stem", "st" ~ "flw_stem",
                                  "ht" ~ "max_stem_height",
                                  "hstm" ~ "herb_count",
                                  "h" ~ "flw_head", "hd" ~ "flw_head", "he" ~ "flw_head", "hea" ~ "flw_head",
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
  pivot_wider(id_cols = c(plant_id, Site_tag, patch, quad, plant, TP, row_id, census_year), names_from = measurement, values_from = value) %>% 
  mutate(across(c(plant_id, Site_tag, patch, plant, TP, row_id), as.character)) %>% 
  mutate(across(c( census_year, ARCHBOLD_surv, flw_stem, flw_head, herb_count), function(x) suppressWarnings(as.integer(x)))) %>%  
  mutate(across(c(ros_diameter, max_stem_height), function(x) suppressWarnings(as.numeric(x)))) %>% 
  mutate(surv = case_match(ARCHBOLD_surv,
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ 1, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0, 88 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  mutate(ros_diameter = case_when(ros_diameter == 999.0 ~ NA, TRUE ~ ros_diameter)) %>% 
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 2~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 9 ~ 0,
                          lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 8 ~ 1,
                          TRUE ~ surv))  %>% 
  group_by(plant_id) %>% 
  mutate(census_temp = case_when(!is.na(surv) ~ census_year),first_year = min(census_temp, na.rm = T)) %>% 
  ungroup()  %>% 
  dplyr::select(plant_id,Site_tag, patch,  quad, plant, TP, row_id, first_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count)


# Still need to decide on the best way/if to combine the spatial hierarchy across sites. 

# combining data from the three sites
# Then creating columns for lagged census year size and reproduction 

ERYCUN <- ERYCUN_abs %>% 
  full_join(ERYCUN_apfi) %>% 
  full_join(ERYCUN_apsc) %>% 
  # filter(!is.na(surv)) %>% 
  group_by(plant_id) %>% 
  arrange(plant_id, census_year) %>% 
  mutate(year.t1 = census_year,
         surv.t1 = surv,
         ros_diameter.t1 = ros_diameter,
         max_stem_height.t1 = max_stem_height,
         flw_stem.t1 = flw_stem,
         flw_head.t1 = flw_head,
         herb_count.t1 = herb_count) %>% 
  mutate(year.t = dplyr::lag(year.t1, n = 1, default = NA),
         ros_diameter.t = dplyr::lag(ros_diameter.t1, n = 1, default = NA),
         max_stem_height.t = dplyr::lag(max_stem_height.t1, n = 1, default = NA),
         flw_stem.t = dplyr::lag(flw_stem.t1, n = 1, default = NA),
         flw_head.t = dplyr::lag(flw_head.t1, n = 1, default = NA),
         herb_count.t = dplyr::lag(herb_count.t1, n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% 
  select(plant_id,Site_tag, bald, patch, quad, plant, TP, row_id, first_year, 
         year.t1, ARCHBOLD_surv, surv.t1, ros_diameter.t1, max_stem_height.t1, flw_stem.t1, flw_head.t1, herb_count.t1,
         year.t, ros_diameter.t, max_stem_height.t, flw_stem.t, flw_head.t, herb_count.t)




####################################################################################
###### Cleaning and merging together LIAOHL #######################################
####################################################################################

# Liatris ohlingerae data was collected starting in 1997 up to 2017 and is some of the longest term monitoring data at Archbold. Collected primarily by Eric Menges and Pedro Quintana-Ascencio

### Loading in the LIAOHL data, cleaning out redundant columns and pivoting to long format
# removing columns that seem like indices related to moving data formats, as well as some columns that have meta info about pulling tags
# removing information about site burn history, as well as what I think is hurricane history?
# removing redundant demographic information or derived data columns (e.g the log of height, or survival across transition years)
# removing old id column which is empty
# removing columns that seem related to notes about location or different cohort groups and are only recorded for some of the plants
# removing columns, which I think are about the community of plants at the site (oaks, palms etc.)
# renaming columns that contain '#' 
# Then pivoting to long format
# recoding the column names. Here height is recorded as the sum of all heights of all measured stems, so a bit different from others that have height of tallest stem
# for individ without only rosettes, they recorded number of rosettes and total number of leaves
# herbivore stem counts have some funky outlier value so recoding that along with the survival code (even though I'm not sure we'll use this herbivore info anyways)
# recoding the archbold survival code and checking on plants that were not found and then found alive in later censuses
# Then recoding our survival column to get rid of multiple years recorded in a row as not found or as dead, and also recoding cases where plant was marked as dead/not found in one year and then alive in later


LIAOHL_temp <- LIAOHL_raw %>% 
  dplyr::select(-`\\outl0\\strokewidth0 \\strokec2 filter_$`, -299, -pull11, -prevpull, -contains("QAQC")) %>% 
  dplyr::select(-contains("burn"), -contains("TSF"), -pc_damage, -sumhurr, -hurrica) %>% 
  dplyr::select(-starts_with("lnht"), -starts_with("s9"), -starts_with("s0"), -starts_with("agr"), -dormant0910, -contains("stg"), -contains("stage"), -stclass)  %>% 
  dplyr::select(-s10_11, -s11_12, -s12_13, -s13_14, -s14_15, -s15_16) %>% 
  dplyr::select(-firstflw_age, -firstflwyear,-contains("flw"), -flower_ever, -seedyear, -`flw#`, -`top#`, -firstYr_pop, -origin) %>% 
  dplyr::select(-old_id, -problem, -nn, -location, -quadrant, -reappr, -cohort, -LC01group, -ltreb, -demproj, -Rx17, -Rx2015, -site, -hab2) %>% 
  dplyr::select(-n_oak, -oak_ht, -n_palm, -palm_ht, -n_rose, -rose_ht, -cons) %>% 
  mutate(row_id = row_number()) %>% 
  mutate(plant_id = paste(pop,plt_no,id,tp, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename_at(vars(contains("#")), ~str_replace(., "#", "")) %>% 
  pivot_longer(cols = !c(plant_id, pop, plt_no, id, tp, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                  "s" ~ "ARCHBOLD_surv", "surv" ~ "ARCHBOLD_surv",
                                  "stm" ~ "flw_stem", "st" ~ "flw_stem", 
                                  "top" ~ "herb_count_stem",
                                  "hgt" ~ "total_height", "totht" ~ "total_height",
                                  "hds" ~ "flw_head", "heads" ~ "flw_head", "tothds" ~ "flw_head",
                                  "hdsdamg" ~ "herb_count_head", "totdam" ~ "herb_count_head",
                                  "ros" ~ "num_rosettes",
                                  "lvs" ~ "num_leaves",
                                  "sa" ~ "sand_accretion",
                                  "comm" ~ "comment", "com" ~ "comment",
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
                                  12 ~ 2012,
                                  11 ~ 2011,
                                  10 ~ 2010,
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
                                  .default = as.numeric(census_year))) %>% 
  pivot_wider(id_cols = c(plant_id, pop, plt_no, id, tp, row_id, first_year, census_year), names_from = measurement, values_from = value) %>% 
  mutate(across(c(plant_id, pop, plt_no, id, tp, row_id, comment), as.character)) %>% 
  mutate(across(c(first_year, census_year, flw_stem, flw_head, herb_count_head), as.integer)) %>% 
  mutate(across(c(total_height), as.numeric)) %>% 
  mutate(surv = case_match(ARCHBOLD_surv,
                           "0" ~ 0,
                           "1" ~ 1,
                           "FALSE" ~ 0,
                           "TRUE" ~ 1,
                           "3" ~ 1, # 3 is new adult
                           "5" ~ 1, # 5 is new seedling
                           "7" ~ 1, # 7 is putative seedling, there are about 100 of these mostly in 2000 with no size measurements
                           "2" ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           "18" ~ 1, # census notes say that 18 is coded for plants which are tag not found and found later, so this was back corrected in the datasheet, but only about 31 of these individuals
                           "8" ~ 1, # census notes say that 8 is coded for plants which are dormant and found later, so this was back corrected in the datasheet
                           "19" ~ NA, # these are I assume a typo for 9, which is marked for previously dead plants. I think this is likely because each of these have no size data and spot checked to confirm they weren't refound.
                           "9" ~ NA, # 9 is for previously dead
                           "20" ~ NA, # I think this is a typo for 2, spot checking shows plants that have been marked tag not found for several years in a row
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  group_by(plant_id) %>% 
  mutate(surv = case_when(lead(surv) == 0 & surv == 0 ~ NA,
                          is.na(lead(surv)) & surv == 0 ~ NA,
                          lag(surv) == 1 & ARCHBOLD_surv == 2 & surv == 0 ~ 1,
                          lead(surv) == 1& ARCHBOLD_surv == 19 & is.na(surv) ~ 0,
                          lead(surv) == 1 & ARCHBOLD_surv == 9 & is.na(surv) ~ 0,
                          TRUE ~ surv)) %>% 
  ungroup() %>% 
  dplyr::select(plant_id, pop, plt_no, id, tp, row_id, first_year, census_year, ARCHBOLD_surv, surv, total_height, flw_stem, flw_head, herb_count_stem, herb_count_head, num_rosettes, num_leaves, comment)
                

# View(LIAOHL_temp)

# Now creating lagged columns and removing NA rows
  
LIAOHL <- LIAOHL_temp %>% 
  group_by(plant_id) %>% 
  arrange(plant_id, census_year) %>% 
  mutate(year.t1 = census_year,
         surv.t1 = surv,
         total_height.t1 = total_height,
         flw_stem.t1 = flw_stem,
         flw_head.t1 = flw_head,
         herb_count_stem.t1 = herb_count_stem,
         herb_count_head.t1 = herb_count_head,
         num_rosettes.t1 = num_rosettes,
         num_leaves.t1 = num_leaves) %>% 
  mutate(year.t = dplyr::lag(year.t1, n = 1, default = NA),
         total_height.t = dplyr::lag(total_height.t1, n = 1, default = NA),
         flw_stem.t = dplyr::lag(flw_stem.t1,n = 1, default = NA),
         flw_head.t = dplyr::lag(flw_head.t1,n = 1, default = NA),
         herb_count_stem.t = dplyr::lag(herb_count_stem.t1,n = 1, default = NA),
         herb_count_head.t = dplyr::lag(herb_count_head.t1,n = 1, default = NA),
         num_rosettes.t = dplyr::lag(num_rosettes.t1,n = 1, default = NA),
         num_leaves.t = dplyr::lag(num_leaves.t1, n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% 
  dplyr::select(plant_id, pop, plt_no, id, tp, row_id, first_year, census_year,
                year.t1, ARCHBOLD_surv, surv.t1, total_height.t1, flw_stem.t1, flw_head.t1, herb_count_stem.t1, herb_count_head.t1, num_rosettes.t1, num_leaves.t1,
                year.t, total_height.t, flw_stem.t, flw_head.t, herb_count_stem.t, herb_count_head.t, num_rosettes.t, num_leaves.t)







####################################################################################
###### Cleaning and merging together BALANG #######################################
####################################################################################
# Balduina angustifolia data was collected by Beth Stephens, starting in 2009. The data includes seed addition experiments into different microhabitats as well as monitoring of establishment and survival growth and reproduction during two years. she collected data monthly.
# The data file is a bit tricky because the raw demographic counts (size/repro) are in messy columns
# Many of the columns are unnamed/underneath a set of merged cells which have the date of the census. 
# There are also several columns in the BALANG sheet that are named CF which is the code for CHAFAS, but these are typos and should be bf. I'm correcting that in the list.
# The spreadsheet also includes separate sites, listed out vertically with their own column names. And so I am dropping those columns. The sites have slightly different dates for the germination census, but within the same months. I'm reformatting the months to make them more consistent


BALANG_colnames <- c("habitat", "site", 	"microsite",	"point",	"primary_shrub",	"secondary_shrub", "tag", 
                     "new Ba seedlings;5/26/09",
                     "estab Ba seedlings;6/5/09",	"new Ba seedlings;6/5/09", 
                     "estab Ba seedlings;6/9/09",	"new Ba seedlings;6/9/09", 
                     "estab Ba seedlings;6/18/09", "new Ba seedlings;6/18/09", 
                     "estab Ba seedlings;6/24/09",	"new Ba seedlings;6/24/09",
                     "estab Ba seedlings;7/23/09","new Ba seedlings;7/23/09", 
                     "estab Ba seedlings;8/18/09",	"new Ba seedlings;8/18/09", 
                     "estab Ba seedlings;9/21/09","new Ba seedlings;9/21/09", 
                     "estab Ba seedlings;10/29/09","new Ba seedlings;10/29/09", 
                     "estab Ba seedlings;11/23/09",	"new Ba seedlings;11/23/09", 
                     "estab Ba seedlings;12/16/09",	"new Ba seedlings;12/16/09",
                     "estab Ba seedlings;1/20/10",	"new Ba seedlings;1/20/10", 
                     "estab Ba seedlings;2/17/10","new Ba seedlings;2/17/10", 
                     "estab Ba seedlings;3/24/10","new Ba seedlings;3/24/10", 
                     "estab Ba seedlings;4/27/10",	"new Ba seedlings;4/27/10", 
                     "estab Ba seedlings;5/27/10",	"new Ba seedlings;5/27/10", 
                     "estab Ba seedlings;6/30/10",	"new Ba seedlings;6/30/10", 
                     "estab Ba seedlings;7/30/10",	"new Ba seedlings;7/30/10", 
                     "estab Ba seedlings;8/30/10",	"new Ba seedlings;8/30/10", 
                     "estab Ba seedlings;9/29/10",	"new Ba seedlings;9/29/10", 
                     "estab Ba seedlings;10/27/10",	"new Ba seedlings;10/27/10", 
                     "estab Ba seedlings;11/24/10",	"new Ba seedlings;11/24/10", 
                     "estab Ba seedlings;12/?/2010","new Ba seedlings;12/?/2010", 
                     "estab Ba seedlings;1/15/11",	"new Ba seedlings;1/15/11", 
                     "estab Ba sdlings;2/?/2011","new Ba sdlings;2/?/2011",   
                     "estab Ba sdlings;3/?/2011","new Ba sdlings;3/?/2011",   
                     "estab Ba sdlings;4/19/11","new Ba sdlings;4/19/11",   
                     "estab Ba sdlings;5/18/11","new Ba sdlings;5/18/11",	 
                     "estab Ba sdlings;6/16/11","new Ba sdlings;6/16/11",   
                     "estab Ba seedling;7/27/11",	"new Ba seedlings;7/27/11", 
                     "estab Ba seedling;8/31/11","new Ba seedlings;8/31/11", 
                     "estab Ba seedling;9/23/11",	"new Ba seedlings;9/23/11", 
                     "estab Ba seedling;10/26/11",	"new Ba seedlings;10/26/11", 
                     "estab Ba seedling;11/30/11",	"new Ba seedlings;11/30/11", 
                     "estab Ba seedlings;12/?/2011","new Ba seedlings;12/?/2011", 
                     "estab Ba seedlings;1/?/2012","new Ba seedlings;1/?/2012", 
                     "estab Ba seedlings;2/?/2012","new Ba seedlings;2/?/2012", 
                     "estab Ba seedlings;03/-/2012","new Ba seedlings;03/-/2012", 
                     "estab Ba seedlings;04/-/2012","new Ba seedlings;04/-/2012", 
                     "estab Ba seedlings;05/-/2012","new Ba seedlings;05/-/2012", 
                     "sum", 
                     paste0("height_", 1:4,";", "5/27/10"),
                     paste0("height_", 1:5,";", "6/30/10"),
                     paste0("height_", 1:4,";", "7/30/10"),
                     paste0("height_", 1:5,";", "8/30/10"),
                     paste0("height_", 1:4,";", "9/29/10"),
                     paste0("height_", 1:4,"_dupe",";", "9/29/10"),
                     paste0("height_", 1:4,";", "10/27/10"),
                     "height_1;11/24/10",
                     paste0("height_", 1:7,";", "12/?/2010"),
                     paste0("height_", 1:7,";", "1/26/11"),
                     paste0("height_", 1:4,";", "3/-/2011"),
                     paste0("height_", 1:8,";", "4/-/2011"),
                     paste0("height_", 1:7,";", "5/-/2011"),
                     paste0("height_", 1:7,";", "6/16/11"),
                     paste0("height_", 1:7,";", "7/27/11"),
                     paste0("height_", 1:6,";", "8/31/11"),
                     paste0("height_", 1:6,";", "09/-/2011"),
                     paste0("height_", 1:5,";", "10/-/2011"),
                     "height_1;11/-/2011",
                     paste0("height_", 1:5,";", "12/-/2011"),
                     paste0("height_", 1:5,";", "01/-/2012"),
                     paste0("height_", 1:6,";", "02/-/2012"),
                     paste0("height_", 1:5,";", "03/-/2012"),
                     paste0("height_", 1:4,";", "04/-/2012"),
                     paste0("height_", 1:4,";", "05/-/2012"))

BALANG_renamed <- BALANG_raw

colnames(BALANG_renamed) <- BALANG_colnames

# Dropping the secondary column names within the spreadsheet
# then pivoting the germination counts to have them in a single column. The germination counts start in May 2009 and the first size measurements occur in May 2010
# Each subplot has a unique tag, with multiple plants inside of it.
# The columns for height measurements for 9/29/10 seem to be duplicated. The only difference for the set of columns is a place where the tag number was copied over, so I am dropping these columns. The census for 10/2010 seems a bit funky with one measurement that might be too small, no census for 11/2010, and then 12/2010 seems like there was a lot of turnover in the plots and I'm worried that they column order may not always correspond to the plant id
# Currently have the germinants as individual rows, but realizing that it is probably best just to keep the germination data separate because I don't trust that the row order actually tracks the individual plants. It's impossible to tell which of the germinants survived/died before the height measurements start to be recorded.

BALANG_seedlings <- BALANG_renamed %>% 
  filter(habitat != "habitat", habitat != "2: Res") %>% mutate(habitat = case_when(habitat == "1: ABS" ~ "1", TRUE ~ habitat)) %>% 
  select(-sum,-contains("dupe"), -primary_shrub, -secondary_shrub, -contains("height")) %>%
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = -c(quote_bare(habitat, site, microsite, point, tag)), names_to = c("name","date"), names_pattern = "(.*)(?:[;])(.*)", values_to = "measurement") %>% 
  mutate(name = case_when(name == "new Ba seedlings" | name == "new Ba sdlings" ~ "new_sdlg_count",
                          name == "estab Ba seedlings" | name == "estab Ba sdlings" | name == "estab Ba seedling" ~ "estab_sdlg_count", TRUE ~ name)) %>% 
  mutate(month = lubridate::month(lubridate::mdy(date)),
         year = lubridate::year(lubridate::mdy(date))) %>% 
  pivot_wider(id_cols = c(quote_bare(habitat, site, microsite, point, tag, date, month, year)), names_from = name, values_from = measurement)  %>% 
  mutate(estab_sdlg_count = as.numeric(str_remove(estab_sdlg_count, "[?*]"))) %>% 
  group_by(tag, date) %>% 
  filter(!is.na(new_sdlg_count) & !is.na(estab_sdlg_count)) %>% 
  mutate(indiv_sdlg_list = case_when(!is.na(estab_sdlg_count) & estab_sdlg_count > 0 ~ paste(1:estab_sdlg_count, collapse = ";"),
                            TRUE ~ NA),
         indiv_sdlg_birth = case_when(!is.na(new_sdlg_count) & new_sdlg_count > 0 ~ paste(1:new_sdlg_count, collapse = ";"))) %>% 
  ungroup() %>% 
  mutate(lead_indiv_sdlg_list = lead(indiv_sdlg_list)) %>% 
  separate(lead_indiv_sdlg_list, ";", into = paste0("seedling",(1: max(BALANG_seedlings$estab_sdlg_count))), fill = "right") %>% 
  pivot_longer(cols = starts_with("seedling"), names_to = "seedling_id", values_to = "seedling_value") %>% 
  mutate(seed_surv = case_when(seedling_value >= 1 ~ 1, TRUE ~ 0))


  pivot_longer(cols = contains("height"), names_to = c("plant_no"), names_prefix = "height_", values_to = "height") %>% 
  mutate(plant_id = paste(tag, plant_no, sep = "_"))
  
  
# the growth data for individual plants is messy, but it is at least tracking clear individuals, so we will keep the two datasets separate. 
  #we want each transition year for each plant to be an individual row, so I am pivoting longer to get individual plants then separating to split up the height measurements . 
  # Then I will split out the string in the height column to keep track of flowering, etc.
  # Using stringr to pull out the height number which is always first, then extracting the number before the strings "br" and "bu" for "branches" and "buds" respectively.
  
  
BALANG_growth <- BALANG_renamed %>% 
  filter(habitat != "habitat", habitat != "2: Res") %>% mutate(habitat = case_when(habitat == "1: ABS" ~ "1", TRUE ~ habitat)) %>% 
  select(-sum,-contains("dupe"), -primary_shrub, -secondary_shrub, -contains("BA sdlings"), -contains("BA seedling")) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = starts_with("height"), names_sep = ";", names_to = c("id", "date")) %>% 
  mutate(height = parse_number(value)) %>% 
  mutate(branches = as.numeric(str_extract(value, "\\d+(?=\\sbr)|\\d+(?=\\?\\sBr)")),
         top_branches = as.numeric(str_extract(value, "\\d+(?=\\st\\sbr)|\\d+(?=\\stop\\sbr)")),
         mid_branches = as.numeric(str_extract(value, "\\d+(?=\\sm\\sbr)|\\d+(?=\\smid\\sbr)")),
         bottom_branches = as.numeric(str_extract(value, "\\d+(?=\\sb\\sbr)|\\d+(?=\\sbot\\sbr)|\\d+(?=\\sbottom\\sbr)"))) %>% 
  mutate_at(vars(branches, top_branches, mid_branches, bottom_branches), ~replace_na(.,  0)) %>% 
  mutate(branch_count = branches+top_branches+mid_branches+bottom_branches) %>% 
  mutate(buds = as.numeric(str_extract(value, "\\d+(?=\\sbu)")))
                                
                                
                                
                                
                          grepl("t br", value) ~ sum(as.numeric(str_extract(value, "\\d+(?=\\st\\sbr)")), as.numeric(str_extract(value, "\\d+(?=\\smid\\sbr)")), na.rm = TRUE)))
         
         
         
         
         bud = str_extract(value, "\\d+(?=\\sbu)"))


  
  

####################################################################################
# Purpose: Read in Archbold demographic data from eleven species, clean and process data to have consistent format for vital rate analyses
# Author: Joshua Fowler
# Date: Aug 16, 2023
####################################################################################

##### Set up #####

# library(renv) # track package versions
# # renv::init()
# # renv::snapshot()
# renv::restore()
 

library(tidyverse)
library(readxl)
library(lubridate)


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
# data from seed addition plots
BALANG_raw <- read_excel(path = paste0(filepath, "/BethStephens_BALANG_CHAFAS/BaCf_data6-2-12.xls"), sheet = "Badata_may2012")

# data from 2008 height and flowering census of naturally occurring plants
BALANG_raw_2008 <- read_excel(path = paste0(filepath, "/BethStephens_BALANG_CHAFAS/Baseedprod08_2-24-09summarystats.xls"), sheet = "cleaned")

# data from 2009 census of seedheads of naturally occurring plants
BALANG_raw_2009 <- read_excel(path = paste0(filepath, "/BethStephens_BALANG_CHAFAS/Ba_prod_foranalysis5-4-10.xlsx"), sheet = "data")


# Reading in data for Chamaecrista fasciculate (CHAFAS)
CHAFAS_raw <- read_excel(path = paste0(filepath, "/BethStephens_BALANG_CHAFAS/BaCf_data6-2-12.xls"), sheet = "Cfdata_may2012")


# Reading in the data for Chapmannia floridana (CHAFLO), data from Jenny Schafer
# There are census plots and seedling plots
CHAFLO_raw <- read_excel(path = paste0(filepath, "/JennySchafer_CHAFLO/Chaflo Demography for Josh.xlsx"), sheet = "Demo2015-2021")

CHAFLO_seedling_raw <- read_excel(path = paste0(filepath, "/JennySchafer_CHAFLO/Chaflo Seedling Survival for Josh(1).xlsx"), sheet = "Data")




# Reading in the data for Paronychia chartacea from Jenny Schafer

PARCHA_raw <- read_excel(path = paste0(filepath, "/JennySchafer_PARCHA_ProbablySameAsArchboldsVersion/Paronychia Demography Data.xlsx"))



# Reading in the data for Polygonella basiramia. This is likely Satya Witt's data, but was provided by Aaron David

POLBAS_raw <- read_csv(file= paste0(filepath, "/UM demographic models/pb1103.csv"))



# Reading in the data for Hypericum cumulicola. This is data provided by Aaron David. Majority of data collection led by Pedro Quintana-Ascenscio
HYPCUM_raw <- read_excel(path = paste0(filepath, "/UM demographic models/hcdem_ABS_2022.xlsx"), sheet = "hcdem_ABS_2022", col_types = "text")


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
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & ros_diameter<3 & flw_stem == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  mutate(ros_diameter = case_when(plant_id == "54_9_155_NA_4531" & census_year == 2018 ~ ros_diameter*.1, TRUE ~ ros_diameter)) %>% 
  mutate(flw_stem = case_when(first_year == census_year & !is.na(ros_diameter) & is.na(flw_stem) ~ 0, TRUE ~ flw_stem)) %>% 
  mutate(flw_stem = case_when(!is.na(ros_diameter) & is.na(flw_stem) & is.na(flw_head) ~ 0,
                              !is.na(ros_diameter) & is.na(max_stem_height) & is.na(flw_head) ~ 0, TRUE ~ flw_stem)) %>% 
  dplyr::select(plant_id,Site_tag,bald, patch, plant, TP, row_id, first_year, birth_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count, comment)


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
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & ros_diameter<3 & flw_stem == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  dplyr::select(plant_id,Site_tag, patch, quad, plant, TP, row_id, first_year,  birth_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count, comment)

  
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
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & ros_diameter<3 & flw_stem == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  ungroup()  %>% 
  dplyr::select(plant_id,Site_tag, patch,  quad, plant, TP, row_id, first_year, birth_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count)


# Still need to decide on the best way/if to combine the spatial hierarchy across sites. 

# combining data from the three sites
# Then creating columns for lagged census year size and reproduction 

ERYCUN <- ERYCUN_abs %>% 
  full_join(ERYCUN_apfi) %>% 
  full_join(ERYCUN_apsc) %>% 
  # filter(!is.na(surv)) %>% 
  mutate(flw_status = case_when(flw_stem>0 ~ 1, flw_stem==0 ~ 0,
                                surv == 0 ~ NA)) %>% 
  group_by(plant_id) %>% 
  arrange(plant_id, census_year) %>% 
  mutate(year.t1 = census_year,
         surv.t1 = surv,
         ros_diameter.t1 = ros_diameter,
         max_stem_height.t1 = max_stem_height,
         flw_status.t1 = flw_status,
         flw_stem.t1 = flw_stem,
         flw_head.t1 = flw_head,
         herb_count.t1 = herb_count,
         age.t1 = census_year-birth_year) %>%  
  mutate(year.t = dplyr::lag(year.t1, n = 1, default = NA),
         ros_diameter.t = dplyr::lag(ros_diameter.t1, n = 1, default = NA),
         max_stem_height.t = dplyr::lag(max_stem_height.t1, n = 1, default = NA),
         flw_status.t = dplyr::lag(flw_status.t1, n = 1, default = NA),
         flw_stem.t = dplyr::lag(flw_stem.t1, n = 1, default = NA),
         flw_head.t = dplyr::lag(flw_head.t1, n = 1, default = NA),
         herb_count.t = dplyr::lag(herb_count.t1, n = 1, default = NA),
         age.t = dplyr::lag(age.t1, n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% 
  select(plant_id,Site_tag, bald, patch, quad, plant, TP, row_id, first_year, birth_year,
         year.t1, ARCHBOLD_surv, surv.t1, ros_diameter.t1, max_stem_height.t1, flw_status.t1, flw_stem.t1, flw_head.t1, herb_count.t1, age.t1,
         year.t, ros_diameter.t, max_stem_height.t, flw_status.t, flw_stem.t, flw_head.t, herb_count.t, age.t)




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
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  ungroup() %>% 
  dplyr::select(plant_id, pop, plt_no, id, tp, row_id, first_year, birth_year, census_year, ARCHBOLD_surv, surv, total_height, flw_stem, flw_head, herb_count_stem, herb_count_head, num_rosettes, num_leaves, comment)
                

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
         num_leaves.t1 = num_leaves,
         age.t1 = census_year - birth_year) %>% 
  mutate(year.t = dplyr::lag(year.t1, n = 1, default = NA),
         total_height.t = dplyr::lag(total_height.t1, n = 1, default = NA),
         flw_stem.t = dplyr::lag(flw_stem.t1,n = 1, default = NA),
         flw_head.t = dplyr::lag(flw_head.t1,n = 1, default = NA),
         herb_count_stem.t = dplyr::lag(herb_count_stem.t1,n = 1, default = NA),
         herb_count_head.t = dplyr::lag(herb_count_head.t1,n = 1, default = NA),
         num_rosettes.t = dplyr::lag(num_rosettes.t1,n = 1, default = NA),
         num_leaves.t = dplyr::lag(num_leaves.t1, n = 1, default = NA),
         age.t = dplyr::lag(age.t1, n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% 
  dplyr::select(plant_id, pop, plt_no, id, tp, row_id, first_year, birth_year, census_year,
                year.t1, ARCHBOLD_surv, surv.t1, total_height.t1, flw_stem.t1, flw_head.t1, herb_count_stem.t1, herb_count_head.t1, num_rosettes.t1, num_leaves.t1, age.t1,
                year.t, total_height.t, flw_stem.t, flw_head.t, herb_count_stem.t, herb_count_head.t, num_rosettes.t, num_leaves.t, age.t)







####################################################################################
###### Cleaning and merging together BALANG #######################################
####################################################################################
# Balduina angustifolia data was collected by Beth Stephens, starting in 2009. The data includes seed addition experiments into different microhabitats as well as monitoring of establishment and survival growth and reproduction during two years. she collected data monthly.
# The data file is a bit tricky because the raw demographic counts (size/repro) are in messy columns
# Many of the columns are unnamed/underneath a set of merged cells which have the date of the census. 
# There are also several columns in the BALANG sheet that are named CF which is the code for CHAFAS, but these are typos and should be bf. I'm correcting that in the list.
# The spreadsheet also includes separate sites, listed out vertically with their own column names. And so I am dropping those columns. The sites have slightly different dates for the germination census, but within the same months. I'm reformatting the months to make them more consistent
# Note that there are 7 sites -> 
# Garbaria = Bald 36 at Archbold
# Grapevine = Bald 34 at Archbold
# Triangle = Bald 29 at Archbold
# Deer Bones = Bald 42 at Archbold
# Hicoria = Bald 46 at Archbold
# First_burn = Bald 25 at Archbold
# site 14 = ? at Archbold reserve (property to the west of main Archbold)
# site 23 = ? at Archbold reserve (property to the west of main Archbold)
# site 15 = ? at Archbold reserve (property to the west of main Archbold)
#
# This data also contains microsite, which describes bare sand, litter, and litter under shrub. Unclear what the numbers mean.


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
  filter(habitat != "habitat", habitat != "2: Res") %>% mutate(habitat = case_when(habitat == "1: ABS" ~ "1", TRUE ~ habitat)) %>% # here habitat == 1 is for Archbold proper with known balds, habitat == 2 is for archbold reserve, less clear exact locations, more disturbed sites
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
  mutate(bald = case_when(site == "1st burn" ~ "25",
                          site == "Grapevine" ~ "34",
                          site == "Hicoria" ~ "46",
                          TRUE ~ NA)) %>% 
  ungroup()
  
# the growth data for individual plants is messy, but it is at least tracking clear individuals, so I will keep the two datasets separate. 
  #I want each transition year for each plant to be an individual row, so I am pivoting longer to get individual plants then separating to split up the height measurements . 
  # Then I will split out the string in the height column to keep track of flowering, etc.
  # Using stringr to pull out the height number which is always first, then extracting the number before the strings "br" and "bu" and "fl" for "branches" and "buds" and "flowers" respectively.
  # Sent and Email to Beth Stephens to figure out what NB means. She says it is likely "No bolt". Still not 100% clear if that means the plant was too small to measure (like it was just a rosette) or if it wasn't found. There are several places throughout the data where a plant is recorded in one month, then maybe goes missing in the next month, and then re-appears. It's difficult to tell when these are the same plant or just a seedling that popped up then died, then a new plant appeared. Also related to this. Plants are mostly consistently recorded in the same column each census, but there are definitely a few instances where this isn't true.
# Overall, I think I will just need to filter this down until we have only the plants we are confident in.
toothpick_match <- c('WS', "BD", "BH", "RH", "WBC", "BWC", "GS", 
                     "gc", "GC", "OC", "OS", "RN", "YC", "YS",  
                     "YN", "YD", "WC", "NT", "no TP", 
                     "YH", "Ynub", "Wnub","tp?", "White")
toothpick_pattern <- paste0("\\b", paste(toothpick_match , collapse="\\b|\\b"), "\\b")


BALANG_growth <- BALANG_renamed %>% 
  filter(habitat != "habitat", habitat != "2: Res") %>% mutate(habitat = case_when(habitat == "1: ABS" ~ "1", TRUE ~ habitat)) %>% 
  select(-sum,-contains("dupe"), -primary_shrub, -secondary_shrub, -contains("BA sdlings"), -contains("BA seedling")) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = starts_with("height"), names_sep = ";", names_to = c("height_column", "date")) %>% 
  mutate(id = paste0(tag, "_", parse_number(height_column))) %>% 
  mutate(origin = "seed_addition") %>% 
  mutate(height = case_when(grepl("no ht recorded", value)~ NA,
                            TRUE ~ parse_number(value))) %>% 
  mutate(branches = as.numeric(str_extract(value, "\\d+(?=\\sbr)|\\d+(?=\\?\\sBr)|\\d+(?=br)")),
         top_branches = as.numeric(str_extract(value, "\\d+(?=\\st\\sbr)|\\d+(?=\\stop\\sbr)")),
         mid_branches = as.numeric(str_extract(value, "\\d+(?=\\sm\\sbr)|\\d+(?=\\smid\\sbr)")),
         bottom_branches = as.numeric(str_extract(value, "\\d+(?=\\sb\\sbr)|\\d+(?=\\sbot\\sbr)|\\d+(?=\\sbottom\\sbr)|\\d+(?=\\smain\\sbr)"))) %>% 
  mutate_at(vars(branches, top_branches, mid_branches, bottom_branches), ~replace_na(.,  0)) %>% 
  mutate(branch_count = branches+top_branches+mid_branches+bottom_branches) %>% 
  mutate(buds = as.numeric(str_extract(value, "\\d+(?=\\sbu)")),
         flowers = as.numeric(str_extract(value, "\\d+(?=\\sfl)"))) %>% 
  dplyr::filter(height>0 | is.na(height)) %>%
  mutate(date_temp = case_when(grepl( "/-/", date) ~ sub("\\-", "1", date),
                               grepl("/?/", date) ~ sub("\\?", "1", date),
                               TRUE ~ date),
         date_fixed = case_when(!grepl("/20", date_temp) ~ as_date(date_temp, format = "%m/%d/%y"),
                                grepl("/20", date_temp) ~ as_date(date_temp, format = "%m/%d/%Y")),
         census_month = lubridate::month(date_fixed),
         census_year = lubridate::year(date_fixed)) %>% 
  arrange(id, census_year, census_month) %>% 
  mutate(surv_inferred = case_when(!is.na(height) ~ "alive",
                          !is.na(dplyr::lag(height, n = 1)) & is.na(dplyr::lead(height, n = 1)) & is.na(dplyr::lead(height, n = 2)) & is.na(height) ~ "inferred dead",
                          TRUE ~ "not recorded"),
         monthly_surv = case_when(surv_inferred == "alive" ~ 1, surv_inferred == "inferred dead" ~ 0, TRUE ~ NA)) %>% 
  mutate(toothpick_id = str_extract(value, pattern = tolower(regex(toothpick_pattern, ignore_case = TRUE))),
         TP = case_when(is.na(toothpick_id) ~ id, TRUE ~ toothpick_id),
         plant_id = paste(habitat, site, microsite, point, tag, TP, sep = "_")) %>% 
  dplyr::select(habitat, site, microsite, point, tag,  id, plant_id, census_month, census_year,monthly_surv, height,branch_count, buds, flowers, value) %>% 
  # for this, I am taking only the yearly fall census, so plants that die in May are now added to the list of those that died in November
  # filter(!is.na(monthly_surv)) %>%
  group_by(plant_id) %>%
  arrange(plant_id,census_year,census_month) %>% 
  mutate(birth_month = case_when(is.na(dplyr::lag(monthly_surv)) & monthly_surv == 1 ~ census_month),
         birth_year = case_when(is.na(dplyr::lag(monthly_surv)) & monthly_surv == 1 ~ census_year)) %>% 
  fill(birth_month, .direction = "updown") %>% 
  fill(birth_year, .direction = "updown") %>% 
  mutate(birth_date = as.Date(paste0("01","/",birth_month,"/",birth_year), format = c("%d/%m/%Y")),
         census_date = as.Date(paste0("01","/",census_month,"/",census_year), format = c("%d/%m/%Y")),
         age = as.period(interval(birth_date, census_date)),
         age_in_months = age %/% months(1)) %>% 
  filter(age_in_months >=0) %>%
  filter(!is.na(monthly_surv)) %>% 
  mutate(bald = case_when(site == "1st burn" ~ "25",
                          site == "Grapevine" ~ "34",
                          site == "Hicoria" ~ "46",
                          TRUE ~ NA_character_)) %>% 
  ungroup()
# The most common birth month is April (4), so I am using that as the annual census date for now.
  # filter(birth_month == 4) %>% 



  
  BALANG_surv_prob <- BALANG_growth %>% 
    group_by(bald) %>% 
    summarize(surv_prob = mean(age=="1y 0m 0d 0H 0M 0S"))

  

  
  # Beth also has data from some plants that were naturally occurring. Need to incorporate this data, but for now I'm gonna ignore that. It is primarily useful for reproduction.

# BALANG_raw_2008_colnames <- c("habitat", "site", "tag", 
#                               "date;1", "flw_count;1", "fully_developed_flw_count;1", "height;1", "branch_length;1", "crown_diam;1", 
#                               "date;2", "flw_count;2", "fully_developed_flw_count;2", "height;2", "branch_length;2", "crown_diam;2", 
#                               "date;3", "flw_count;3", "fully_developed_flw_count;3", "height;3", "branch_length;3", "crown_diam;3",
#                               "date;4", "flw_count;4", "fully_developed_flw_count;4", "height;4", "branch_length;4", "crown_diam;4")
#                               
# BALANG_established2008_renamed <- BALANG_raw_2008
# 
# colnames(BALANG_established2008_renamed) <- BALANG_raw_2008_colnames
# 
# 
# BALANG_established2008 <- BALANG_established2008_renamed
#   
  
  ####################################################################################
  ###### Cleaning and merging together CHAFAS #######################################
  ####################################################################################
  # Chamaechrista fasculata data was collected by Beth Stephens, starting in 2009. The data includes seed addition experiments into different microhabitats as well as monitoring of establishment and survival growth and reproduction during two years. she collected data monthly.
  # There is a set of columns for 9/29 which aren't exactly copied of the following 9/29, so I believe there is something like this is supposed to be a different month, but it's also recorded in a different way, and it's unclear if this is actually good data, so I'm just going to drop these and use the following columns.
  # Some of these are like plants that were recorded, but never seem to show up again in later censuses. I might just drop the other set of 9/29 columns too
  
  CHAFAS_colnames <- c("habitat", "site", 	"microsite",	"point","density",	"primary_shrub",	"secondary_shrub", "tag", 
                       "new Cf seedlings;5/26/09",
                       "estab Cf seedlings;6/5/09",	"new Cf seedlings;6/5/09", 
                       "estab Cf seedlings;6/10/09",	"new Cf seedlings;6/10/09", 
                       "estab Cf seedlings;6/18/09", "new Cf seedlings;6/18/09", 
                       "estab Cf seedlings;6/24/09",	"new Cf seedlings;6/24/09",
                       "estab Cf seedlings;7/23/09","new Cf seedlings;7/23/09", 
                       "estab Cf seedlings;8/18/09",	"new Cf seedlings;8/18/09", 
                       "estab Cf seedlings;9/21/09","new Cf seedlings;9/21/09", 
                       "estab Cf seedlings;10/29/09","new Cf seedlings;10/29/09", 
                       "estab Cf seedlings;11/23/09",	"new Cf seedlings;11/23/09", 
                       "estab Cf seedlings;12/16/09",	"new Cf seedlings;12/16/09",
                       "estab Cf seedlings;1/20/10",	"new Cf seedlings;1/20/10", 
                       "estab Cf seedlings;2/17/10","new Cf seedlings;2/17/10", 
                       "estab Cf seedlings;3/24/10","new Cf seedlings;3/24/10", 
                       "estab Cf seedlings;4/27/10",	"new Cf seedlings;4/27/10", 
                       "estab Cf seedlings;5/31/10",	"new Cf seedlings;5/31/10", 
                       "estab Cf seedlings;6/30/10",	"new Cf seedlings;6/30/10", 
                       "estab Cf seedlings;7/30/10",	"new Cf seedlings;7/30/10", 
                       "estab Cf seedlings;8/30/10",	"new Cf seedlings;8/30/10", 
                       "estab Cf seedlings;9/29/10",	"new Cf seedlings;9/29/10", 
                       "estab Cf seedlings;10/27/10",	"new Cf seedlings;10/27/10", 
                       "estab Cf seedlings;11/24/10",	"new Cf seedlings;11/24/10", 
                       "estab Cf seedlings;12/?/2010","new Cf seedlings;12/?/2010", 
                       "estab Cf seedlings;1/26/11",	"new Cf seedlings;1/26/11", 
                       "estab Cf sdlings;2/16/2011","new Cf sdlings;2/16/2011",   
                       "estab Cf sdlings;3/23/2011","new Cf sdlings;3/23/2011",   
                       "estab Cf sdlings;4/19/11","new Cf sdlings;4/19/11",   
                       "estab Cf sdlings;5/18/11","new Cf sdlings;5/18/11",	 
                       "estab Cf sdlings;6/16/11","new Cf sdlings;6/16/11",   
                       "estab Cf seedling;7/27/11",	"new Cf seedlings;7/27/11", 
                       "estab Cf seedling;8/23/11","new Cf seedlings;8/23/11", 
                       "estab Cf seedling;9/23/11",	"new Cf seedlings;9/23/11", 
                       "estab Cf seedling;10/26/11",	"new Cf seedlings;10/26/11", 
                       "estab Cf seedling;11/30/11",	"new Cf seedlings;11/30/11", 
                       "estab Cf seedlings;12/?/2011","new Cf seedlings;12/?/2011", 
                       "estab Cf seedlings;1/?/2012","new Cf seedlings;1/?/2012", 
                       "estab Cf seedlings;2/?/2012","new Cf seedlings;2/?/2012", 
                       "estab Cf seedlings;03/-/2012","new Cf seedlings;03/-/2012", 
                       "estab Cf seedlings;04/-/2012","new Cf seedlings;04/-/2012", 
                       "estab Cf seedlings;05/-/2012","new Cf seedlings;05/-/2012", 
                       "sum", 
                       paste0("height_", 1:7,"_dupe",";", "9/29/10"), # There is a set of columns for 9/29 which aren't exactly copied of the following 9/29, so I believe there is something like this is supposed to be a different month, but it's also recorded in a different way, and it's unclear if this is actually good data, so I'm just going to drop these and use the following columns.
                       paste0("height_", 1:3,";", "9/29/10"),
                       paste0("height_", 1:4,";", "10/27/10"),
                       paste0("height_", 1:4,";", "11/24/10"),
                       paste0("height_", 1,";", "12/?/2010"),# seems like not all the plants were measured in this census
                       paste0("height_", 1:4,";", "1/26/11"), 
                       paste0("height_", 1:2,";", "2/16/11"), 
                       paste0("height_", 1:5,";", "3/-/2011"),
                       paste0("height_", 1:6,";", "4/-/2011"),
                       paste0("height_", 1:3,";", "5/-/2011"),
                       paste0("height_", 1:4,";", "6/16/11"),
                       paste0("height_", 1:3,";", "7/01/11"),
                       paste0("height_", 1:3,";", "8/01/11"),
                       paste0("height_", 1:3,";", "09/01/2011"),
                       paste0("height_", 1,";", "10/-/2011"),
                       paste0("height_", 1:5,";", "11/-/2011"),
                       paste0("height_", 1:3,";", "12/-/2011"),
                       paste0("height_", 1:3,";", "01/-/2012"),
                       paste0("height_", 1:3,";", "02/-/2012"),
                       paste0("height_", 1:3,";", "03/-/2012"),
                       paste0("height_", 1:3,";", "04/-/2012"),
                       paste0("height_", 1:3,";", "05/-/2012"))
  
  CHAFAS_renamed <- CHAFAS_raw
  
  colnames(CHAFAS_renamed) <- CHAFAS_colnames
  
  

  
  # Dropping the secondary column names within the spreadsheet
  # then pivoting the germination counts to have them in a single column. The germination counts start in May 2009 and the first size measurements occur in September 2010
  # Each subplot has a unique tag, with multiple plants inside of it.
  # The columns for height measurements for 9/29/10 seem to be copy-pasted error. There are also some censuses where there is no data, even though there are multiple measured plants on the censuses on either side.
  # Currently have the germinants as individual rows, but realizing that it is probably best just to keep the germination data separate because I don't trust that the row order actually tracks the individual plants. It's impossible to tell which of the germinants survived/died before the height measurements start to be recorded.
  
  CHAFAS_seedlings <- CHAFAS_renamed %>% 
    filter(habitat != "habitat", habitat != "2: Res") %>% mutate(habitat = case_when(habitat == "1: ABS" ~ "1", TRUE ~ habitat)) %>% 
    select(-sum,-contains("dupe"), -primary_shrub, -secondary_shrub, -density, -contains("height")) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols = -c(quote_bare(habitat, site, microsite, point, tag)), names_to = c("name","date"), names_pattern = "(.*)(?:[;])(.*)", values_to = "measurement") %>% 
    mutate(name = case_when(name == "new Cf seedlings" | name == "new Cf sdlings" ~ "new_sdlg_count",
                            name == "estab Cf seedlings" | name == "estab Cf sdlings" | name == "estab Cf seedling" ~ "estab_sdlg_count", TRUE ~ name)) %>% 
    mutate(month = lubridate::month(lubridate::mdy(date)),
           year = lubridate::year(lubridate::mdy(date))) %>% 
    pivot_wider(id_cols = c(quote_bare(habitat, site, microsite, point, tag, date, month, year)), names_from = name, values_from = measurement)  %>% 
    mutate(estab_sdlg_count = as.numeric(str_remove(estab_sdlg_count, "[?*]")),
           new_sdlg_count = as.numeric(str_remove(new_sdlg_count, "[?*]"))) %>% 
    group_by(tag, date) %>% 
    filter(!is.na(new_sdlg_count) & !is.na(estab_sdlg_count)) %>% 
    mutate(indiv_sdlg_list = case_when(!is.na(estab_sdlg_count) & estab_sdlg_count > 0 ~ paste(1:estab_sdlg_count, collapse = ";"),
                                       TRUE ~ NA),
           indiv_sdlg_birth = case_when(!is.na(new_sdlg_count) & new_sdlg_count > 0 ~ paste(1:new_sdlg_count, collapse = ";"))) %>% 
    mutate(bald = case_when(site == "1st burn" ~ "25",
                            site == "Grapevine" ~ "34",
                            site == "Hicoria" ~ "46",
                            TRUE ~ NA)) %>% 
    ungroup()

  
  
  # the growth data for individual plants is messy, but it is at least tracking clear individuals, so we will keep the two datasets separate. 
  #we want each transition year for each plant to be an individual row, so I am pivoting longer to get individual plants then separating to split up the height measurements . 
  # Then I will split out the string in the height column to keep track of flowering, etc.
  # Using stringr to pull out the height number which is always first, then extracting the number before the strings "br" and "bu" and "fl" for "branches" and "buds" and "flowers" respectively.
  # Sent and Email to Beth Stephens to figure out what NB means. She says it is likely "No bolt". Still not 100% clear if that means the plant was too small to measure (like it was just a rosette) or if it wasn't found. There are several places throughout the data where a plant is recorded in one month, then maybe goes missing in the next month, and then re-appears. It's difficult to tell when these are the same plant or just a seedling that popped up then died, then a new plant appeared. Also related to this. Plants are mostly consistently recorded in the same column each census, but there are definitely a few instances where this isn't true.
  # Overall, I think I will just need to filter this down until we have only the plants we are confident in.
  toothpick_match <- c('WS', "BD", "BH", "RH", "WBC", "BWC", "GS", 
                       "gc", "GC", "OC", "OS", "RN", "YC", "YS",  
                       "YN", "YD", "WC", "NT", "no TP", 
                       "YH", "Ynub", "Wnub","tp?", "White")
  toothpick_pattern <- paste0("\\b", paste(toothpick_match , collapse="\\b|\\b"), "\\b")
  
  
CHAFAS_growth <- CHAFAS_renamed %>% 
    filter(habitat != "habitat", habitat != "2: Res") %>% mutate(habitat = case_when(habitat == "1: ABS" ~ "1", TRUE ~ habitat)) %>% 
    select(-sum,-contains("dupe"), -primary_shrub, -secondary_shrub, -contains("Cf sdlings"), -contains("Cf seedling")) %>% 
    mutate(across(everything(), as.character)) %>% 
    pivot_longer(cols = starts_with("height"), names_sep = ";", names_to = c("height_column", "date")) %>% 
    mutate(id = paste0(tag, "_", parse_number(height_column))) %>% 
    mutate(origin = "seed_addition") %>% 
    mutate(height = case_when(grepl("no ht recorded", value)~ NA,
                              grepl("check to see which YD died previously (aug 2010, jan 2011)", value) ~ NA,
                              TRUE ~ parse_number(value))) %>% 
    mutate(branches = as.numeric(str_extract(value, "\\d+(?=\\sbr)|\\d+(?=\\?\\sBr)|\\d+(?=br)")),
           top_branches = as.numeric(str_extract(value, "\\d+(?=\\st\\sbr)|\\d+(?=\\stop\\sbr)")),
           mid_branches = as.numeric(str_extract(value, "\\d+(?=\\sm\\sbr)|\\d+(?=\\smid\\sbr)")),
           bottom_branches = as.numeric(str_extract(value, "\\d+(?=\\sb\\sbr)|\\d+(?=\\sbot\\sbr)|\\d+(?=\\sbottom\\sbr)|\\d+(?=\\smain\\sbr)"))) %>% 
    mutate_at(vars(branches, top_branches, mid_branches, bottom_branches), ~replace_na(.,  0)) %>% 
    mutate(branch_count = branches+top_branches+mid_branches+bottom_branches) %>% 
    mutate(buds = as.numeric(str_extract(value, "\\d+(?=\\sbu)")),
           flowers = as.numeric(str_extract(value, "\\d+(?=\\sfl)"))) %>% 
    mutate(date_temp = case_when(grepl( "/-/", date) ~ sub("\\-", "1", date),
                               grepl("/?/", date) ~ sub("\\?", "1", date),
                               TRUE ~ date),
         date_fixed = case_when(!grepl("/20", date_temp) ~ as_date(date_temp, format = "%m/%d/%y"),
                                grepl("/20", date_temp) ~ as_date(date_temp, format = "%m/%d/%Y")),
         census_month = lubridate::month(date_fixed),
         census_year = lubridate::year(date_fixed)) %>% 
    arrange(id, census_year, census_month) %>% 
    mutate(surv_inferred = case_when(!is.na(height) ~ "alive",
                                     !is.na(dplyr::lag(height, n = 1)) & is.na(dplyr::lead(height, n = 1)) & is.na(dplyr::lead(height, n = 2)) & is.na(height) ~ "inferred dead",
                                     TRUE ~ "not recorded"),
           monthly_surv = case_when(surv_inferred == "alive" ~ 1, surv_inferred == "inferred dead" ~ 0, TRUE ~ NA)) %>% 
    mutate(toothpick_id = str_extract(value, pattern = tolower(regex(toothpick_pattern, ignore_case = TRUE))),
           TP = case_when(is.na(toothpick_id) ~ id, TRUE ~ toothpick_id),
           plant_id = paste(habitat, site, microsite, point, tag, TP, sep = "_")) %>% 
    dplyr::select(habitat, site, microsite, point, tag,  id, plant_id, census_month, census_year,monthly_surv, height,branch_count, buds, flowers, value) %>% 
    # for this, I am taking only the yearly fall census, so plants that die in May are now added to the list of those that died in November
    # filter(!is.na(monthly_surv)) %>%
    group_by(plant_id) %>%
    arrange(plant_id,census_year,census_month) %>% 
    mutate(birth_month = case_when(is.na(dplyr::lag(monthly_surv)) & monthly_surv == 1 ~ census_month),
           birth_year = case_when(is.na(dplyr::lag(monthly_surv)) & monthly_surv == 1 ~ census_year)) %>% 
    fill(birth_month, .direction = "updown") %>% 
    fill(birth_year, .direction = "updown") %>% 
    mutate(birth_date = as.Date(paste0("01","/",birth_month,"/",birth_year), format = c("%d/%m/%Y")),
           census_date = as.Date(paste0("01","/",census_month,"/",census_year), format = c("%d/%m/%Y")),
           age = as.period(interval(birth_date, census_date)),
           age_in_months = age %/% months(1)) %>% 
    filter(age_in_months >=0) %>%
    filter(!is.na(monthly_surv)) %>% 
    mutate(bald = case_when(site == "1st burn" ~ "25",
                          site == "Grapevine" ~ "34",
                          site == "Hicoria" ~ "46",
                          TRUE ~ NA)) %>% 
    ungroup()
  





CHAFAS_surv_prob <- CHAFAS_growth %>% 
  group_by(bald) %>%
  summarize(surv_prob = mean(age=="1y 0m 0d 0H 0M 0S"))




# 


  
  ####################################################################################
  ###### Cleaning and merging together CHAFLO #######################################
  ####################################################################################

  # This data was collected by Jenny Schafer. Overall, the dataset/column names are pretty clear. One thing to note is that some plots were censused for different durations so the "survival" column reflects NA's in some cases where plants were not included in ensuing censuses
  # First I'm going to standardize the column names, then pivot longer and do a bit of clean up to clarify the survival measurements and to add NA values for the missing 2020 census
  # the first census is in June 2015, and then there is a followup census in Dec 2015. A few new plants were added during this census, But they are mostly accounted for in the 2016 census.
  # no data collected in 2020, so can only get some survival data from that year.
  
  column_string <- c("ARCHBOLD_surv", "stem_count", "max_height", "stem_heights", "dead_stems", "flw_status", "herb", "notes")
  CHAFLO_colnames <- c("Pop_id", "Plot_id", 	"Burn_Unit",	"Habitat","YSF",	"Plant_id",

                       paste0(column_string , "-June2021"),
                       paste0(column_string , "-July2019"),
                       paste0(column_string , "-May/June2018"),
                       paste0(column_string , "-June2017"),
                       paste0(column_string , "-June2016"),
                       paste0(column_string , "-Dec2015"),
                       paste0(c("ARCHBOLD_surv", "stem_count", "max_height", "stem_heights", "dead_stems", "flw_status", "flw_before", "notes"), "-June2015"))
  
  CHAFLO_renamed <- CHAFLO_raw
  
  colnames(CHAFLO_renamed) <- CHAFLO_colnames

  
# Max number of stems measured
  max_stem <- 8
  
CHAFLO_census <- CHAFLO_renamed %>% 
  dplyr::filter(Pop_id != "Pop. #") %>% 
  dplyr::select(!"YSF") %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = -c(Pop_id, Plot_id, Burn_Unit, Habitat, Plant_id), names_sep = "-", names_to = c("column_name", "verbatim_date")) %>% 
  mutate(tag_id = paste(Pop_id, Plot_id, Burn_Unit, Habitat, Plant_id, sep = "-"),
         census_month = case_when(grepl( "June", verbatim_date) ~ 6,
                           grepl( "July", verbatim_date) ~ 7,
                           grepl( "Dec", verbatim_date) ~ 12,
                           TRUE ~ NA),
         census_year = parse_number(verbatim_date)) %>% 
  pivot_wider(names_from = column_name, values_from = value) %>% 
  separate(stem_heights, into = paste("stem_height", 1:max_stem), fill = "right") %>% 
  mutate(across(ARCHBOLD_surv:dead_stems, as.numeric)) %>% 
  # group_by(tag_id) %>%
  mutate(surv = case_when(ARCHBOLD_surv == 1 ~ 1, 
                          ARCHBOLD_surv == 0 ~ 0,
                          ARCHBOLD_surv == 3 ~ 1,
                          ARCHBOLD_surv == 5 ~ 1,
                          ARCHBOLD_surv == 9 ~ NA)) %>%
  complete(nesting(Pop_id, Plot_id, Burn_Unit, Habitat, Plant_id, tag_id), census_year = full_seq(census_year, period = 1), fill= list(NA)) %>% 
  group_by(tag_id) %>% 
  mutate(surv = case_when(census_year == 2020 & lead(surv)==1~1, TRUE ~ surv)) %>% 
  arrange(tag_id,census_year,census_month)  %>% 
  mutate(birth_month = case_when(is.na(dplyr::lag(surv)) & surv == 1 ~ census_month),
         birth_year = case_when(is.na(dplyr::lag(surv)) & surv == 1 ~ census_year)) %>% 
  fill(birth_month, .direction = "updown") %>% 
  fill(birth_year, .direction = "updown") %>% 
  mutate(birth_date = as.Date(paste0("01","/",birth_month,"/",birth_year), format = c("%d/%m/%Y")),
         census_date = case_when(census_year == 2020 ~ as.Date(paste0("01","/","06","/",census_year), format = c("%d/%m/%Y")),
                                 census_year != 2020 ~ as.Date(paste0("01","/",census_month,"/",census_year), format = c("%d/%m/%Y"))),
         age = as.period(interval(birth_date, census_date)),
         age_in_months = age %/% months(1)) %>% 
  filter(census_date!= "2015-12-01") %>% # dropping this date, because there are no new births, and at least one of the deaths is later marked alive
  ungroup() 
  
  

CHAFLO <- CHAFLO_census %>% 
  group_by(tag_id) %>% 
  arrange(tag_id, census_year) %>% 
  mutate(year.t1 = census_year,
         surv.t1 = surv,
         age.t1 = age,
         age_in_months.t1 = age_in_months,
         stem_count.t1 = stem_count,
         max_height.t1 = max_height,
         stem_height1.t1 = `stem_height 1`,
         stem_height2.t1 = `stem_height 2`,
         stem_height3.t1 = `stem_height 3`,
         stem_height4.t1 = `stem_height 4`,
         stem_height5.t1 = `stem_height 5`,
         stem_height6.t1 = `stem_height 6`,
         stem_height7.t1 = `stem_height 7`,
         stem_height8.t1 = `stem_height 8`,
         dead_stems.t1 = `dead_stems`,
         flw_status.t1 = flw_status,
         herb.t1 = herb,
         notes.t1 = notes,
         ) %>% 
  mutate(year.t = dplyr::lag(year.t1, n = 1, default = NA),
         age.t = dplyr::lag(age.t1, n = 1, default = NA),
         age_in_months.t = dplyr::lag(age_in_months.t1, n = 1, default = NA),
         stem_count.t = dplyr::lag(stem_count.t1, n = 1, default = NA),
         max_height.t = dplyr::lag(max_height.t1, n = 1, default = NA),
         stem_height1.t = dplyr::lag(stem_height1.t1, n = 1, default = NA),
         stem_height2.t = dplyr::lag(stem_height2.t1, n = 1, default = NA),
         stem_height3.t = dplyr::lag(stem_height3.t1, n = 1, default = NA),
         stem_height4.t = dplyr::lag(stem_height4.t1, n = 1, default = NA),
         stem_height5.t = dplyr::lag(stem_height5.t1, n = 1, default = NA),
         stem_height6.t = dplyr::lag(stem_height6.t1, n = 1, default = NA),
         stem_height7.t = dplyr::lag(stem_height7.t1, n = 1, default = NA),
         stem_height8.t = dplyr::lag(stem_height8.t1, n = 1, default = NA),
         dead_stems.t = dplyr::lag(dead_stems.t1, n = 1, default = NA),
         flw_status.t = dplyr::lag(flw_status.t1, n = 1, default = NA),
         herb.t = dplyr::lag(herb.t1, n = 1, default = NA),
         notes.t = dplyr::lag(notes.t1, n = 1, default = NA)
         ) %>% 
  filter(!is.na(surv.t1)) %>% 
  dplyr::select(tag_id, Pop_id, Plot_id, Burn_Unit, Habitat, Plant_id, census_year, verbatim_date,
                year.t1, ARCHBOLD_surv, surv.t1, age.t1, age_in_months.t1, stem_count.t1, max_height.t1, stem_height1.t1,stem_height2.t1,stem_height3.t1,stem_height4.t1,stem_height5.t1,stem_height6.t1,stem_height7.t1,stem_height8.t1,dead_stems.t1, flw_status.t1, herb.t1,notes.t1,
                year.t, age.t, age_in_months.t, stem_count.t, max_height.t, stem_height1.t,stem_height2.t,stem_height3.t,stem_height4.t,stem_height5.t,stem_height6.t,stem_height7.t,stem_height8.t,dead_stems.t, flw_status.t, herb.t,notes.t)









# # Looking at dec and Jun start dates
# chaflo_2015 <- CHAFLO_census %>% 
#   filter( year == 2015)


# ggplot(data = chaflo_2015) +
#   geom_histogram(aes(x = surv, fill = month, group = month), position = "dodge")
#            
#   
# 
# chaflo_2015 %>% 
#   dplyr::select(tag_id, max_height, month) %>% 
#   pivot_wider(names_from = month, names_prefix = "max_height", values_from = max_height) %>% 
# ggplot() +
#   geom_jitter(aes(x = max_height6, y = max_height12), position = "dodge") + lims(x = c(0,25), y = c(0,25))



# Now adding in the data from seedling addition experiment plots. These were set up and followed between Aug 2016 and  may 2023
# Each row records the germination and survival of a single seed. The columns are called Germ/Surv. They are what I typically call survival: "is the individual present and alive during the census?" and you can figure out the first germination the first time they show up as alive.
# There is a set of columns labeled with duplicate date (9/8/16). I assume this was actually a second set of measurements in the same month but a couple weeks later, which is the case for the rest of the data, so I am labeling 9/18/16. The 2023 census also only recorded month, so making this the middle of the month
# We also have to fix the survival a little bit because it records stuff as dead for subsequent censuses

CHAFLO_seedling_colnames <- c("Habitat", "Burn_Unit", 	"Treatment",	"PVC_number",
                              "ARCHBOLD_surv;8/5/16",
                              "ARCHBOLD_surv;8/15/16",
                              "ARCHBOLD_surv;8/25/16","height;8/25/16",
                              "ARCHBOLD_surv;9/8/16","height;9/8/16",
                              "ARCHBOLD_surv;9/18/16","height;9/18/16",
                              "ARCHBOLD_surv;10/5/16","height;10/5/16","notes;10/5/16",
                              "ARCHBOLD_surv;10/20/16","height;10/20/16","notes;10/20/16",
                              "ARCHBOLD_surv;11/8/16","height;11/8/16","notes;11/8/16",
                              "ARCHBOLD_surv;11/30/16","height;11/30/16","notes;11/30/16",
                              "ARCHBOLD_surv;12/19/16","height;12/19/16","notes;12/19/16",
                              "ARCHBOLD_surv;1/21/17","height;1/21/17","notes;1/21/17",
                              "ARCHBOLD_surv;2/27/17","height;2/27/17",
                              "ARCHBOLD_surv;3/31/17","height;3/31/17",
                              "ARCHBOLD_surv;5/1/17","height;5/1/17","notes;5/1/17",
                              "ARCHBOLD_surv;6/3/17","height;6/3/17","notes;6/3/17",
                              "ARCHBOLD_surv;7/26/17","height;7/26/17","notes;7/26/17",
                              "ARCHBOLD_surv;9/30/17","height;9/30/17","notes;9/30/17",
                              "ARCHBOLD_surv;11/10/17","height;11/10/17","notes;11/10/17",
                              "ARCHBOLD_surv;12/14/17","height;12/14/17","notes;12/14/17",
                              "ARCHBOLD_surv;3/30/18","height;3/30/18","notes;3/30/18",
                              "ARCHBOLD_surv;6/12/18","height;6/12/18","notes;6/12/18",
                              "ARCHBOLD_surv;7/17/19","height;7/17/19","notes;7/17/19",
                              "ARCHBOLD_surv;6/14/21","height;6/14/21","notes;6/14/21",
                              "ARCHBOLD_surv;6/29/22","stem_count;6/29/22","height;6/29/22","notes;6/29/22",
                              "ARCHBOLD_surv;5/15/23","stem_count;5/15/23","height;5/15/23","notes;5/15/23")

CHAFLO_seedling_renamed <- CHAFLO_seedling_raw

colnames(CHAFLO_seedling_renamed) <- CHAFLO_seedling_colnames


CHAFLO_seedling_census <- CHAFLO_seedling_renamed %>% 
  dplyr::filter(Habitat != "Habitat") %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = -c(Habitat, Burn_Unit, Treatment, PVC_number), names_sep = ";", names_to = c("column_name", "date")) %>% 
  mutate(year = year(mdy(date)),
         month = month(mdy(date)),
         day = day(mdy(date))) %>% 
  mutate(tag_id = paste(Habitat, Burn_Unit, Treatment, PVC_number, sep = "-")) %>% 
  pivot_wider(names_from = column_name, values_from = value) %>% 
  mutate(across(c(ARCHBOLD_surv, height, stem_count), as.numeric)) %>% 
  arrange(PVC_number,year,month,day) %>% 
  mutate(surv = case_when(ARCHBOLD_surv == 0 ~ 0, ARCHBOLD_surv == 1 ~ 1,
                          lag(ARCHBOLD_surv, n = 1) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv, n = 1) == 1 & ARCHBOLD_surv == 0 ~ NA,
                          ARCHBOLD_surv == 9 ~ NA,
                          ARCHBOLD_surv == 2 ~ NA,
                          TRUE ~ ARCHBOLD_surv)) %>% 
  complete(nesting(Habitat, Burn_Unit, Treatment, PVC_number, tag_id), year = full_seq(year, period = 1), fill= list(NA)) %>% 
  group_by(tag_id) %>% 
  mutate(surv = case_when(year == 2020 & lead(surv)==1~1, TRUE ~ surv)) %>% ungroup()


CHAFLO_test <- CHAFLO_seedling_census %>% 
  filter(lead(ARCHBOLD_surv, n = 1) &  ARCHBOLD_surv == 0)

  
  
  
# ARCHBOLD_surv == 1 & lag(ARCHBOLD_surv = 0) ~ NA,
  

####################################################################################
###### Cleaning and merging together PARCHA #######################################
####################################################################################
# PARCHA data comes from plots where all individuals were marked in 2003-ish, and new recruits added, but only for the first few years until 2005/2006. these tagged plants were followed until they died out by 2010. 
# as a conseuence we don't have complete records of all plants in each plot. but we do have lots of individuals. Some of these individuals we know their age accurately, and others we don't know. 

# column_string <- c("ARCHBOLD_surv","resprout", "survival", "stem_count", "max_height", "stem_heights", "dead_stems", "flw_status", "herb", "notes")
PARCHA_colnames <- c("bald", "fire", 	"gap",	"site","plot",	"tag",
                     "ARCHBOLD_surv;Mar2010", "resprout;Mar2010",
                     "ARCHBOLD_surv;Dec2009",
                     "ARCHBOLD_surv;July2009", "flowering;July2009",
                     "ARCHBOLD_surv;Oct2009", 
                     "ARCHBOLD_surv;Sep2008", "flowering;Sep2008",
                     "ARCHBOLD_surv;May2008", 
                     "ARCHBOLD_surv;Feb2008", 
                     "ARCHBOLD_surv;Nov2007", "flowering;Nov2007",
                     "ARCHBOLD_surv;Aug2007", "flowering;Aug2007",
                     "ARCHBOLD_surv;May2007",
                     "ARCHBOLD_surv;Feb2007",
                     "ARCHBOLD_surv;Nov2006", "flowering;Nov2006", "resprout;Nov2006",
                     "ARCHBOLD_surv;Aug2006", "flowering;Aug2006",
                     "ARCHBOLD_surv;May2006", "resprout;May2006",
                     "ARCHBOLD_surv;Feb2006", 
                     "ARCHBOLD_surv;Nov2005", "flowering;Nov2005", 
                     "ARCHBOLD_surv;Aug2005", "flowering;Aug2005", 
                     "ARCHBOLD_surv;May2005", "resprout;May2005", 
                     "ARCHBOLD_surv;Feb2005", "resprout;Feb2005",
                     "ARCHBOLD_surv;Nov2004", "flowering;Nov2004", 
                     "ARCHBOLD_surv;Aug2004", "flowering;Aug2004",
                     "ARCHBOLD_surv;May2004", 
                     "ARCHBOLD_surv;Feb2004", 
                     "ARCHBOLD_surv;Nov2003", "flowering;Nov2003",
                     "length;Nov2003", "width;Nov2003", "height;Nov2003", 
                     "ARCHBOLD_surv;Aug2003", "flowering;Aug2003",
                     "ARCHBOLD_surv;May2003",
                     "ARCHBOLD_surv;Feb2003",
                     "resprout",	"nseed503",	"nadlt503",	"nseed803",	"nadlt803",	"nsed1103",	"nadt1103",	"nseed204",	"nadlt204",	"nseed504",	"nadlt504",	"nseed804",	"nadlt804",	"nsed1104",	"nadt1104", "nseed205",	"nadlt205",	"nseed505",	"nadlt505",	"nseed805",	"nadlt805",	"nsed1105"	,"nadt1105",	"nseed206",	"nadlt206")

PARCHA_renamed <- PARCHA_raw

colnames(PARCHA_renamed) <- PARCHA_colnames

PARCHA_dates <- PARCHA_renamed %>% 
  dplyr::select(!starts_with("n"), -resprout) %>% 
  mutate(row_id = row_number()) %>% 
  rename(fire_unit = fire) %>% 
  mutate(plant_id = paste(bald,site, fire_unit, gap, tag, plot, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character) ) %>% 
  pivot_longer(cols = !c(plant_id, bald,site, fire_unit, gap, tag, plot, row_id), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])")  %>% 
  mutate(census_month = case_when( grepl( "Feb", measurement) ~ 2,
                           grepl( "Mar", measurement) ~ 3,
                           grepl( "May", measurement) ~ 5,
                           grepl( "July", measurement) ~ 7,
                           grepl( "Aug", measurement) ~ 8,
                           grepl( "Sep", measurement) ~ 9,
                           grepl( "Oct", measurement) ~ 10,
                           grepl( "Nov", measurement) ~ 11,
                           grepl( "Dec", measurement) ~ 12),
         measurement = str_split(measurement, ";", simplify = T)[, 1]) %>% 
  pivot_wider(names_from = measurement, values_from = value) %>% 
  mutate(surv = case_when(ARCHBOLD_surv == 1 ~ 1, 
                          ARCHBOLD_surv == 0 ~ 0,
                          ARCHBOLD_surv == 3 ~ 1,
                          ARCHBOLD_surv == 5 ~ 1,
                          ARCHBOLD_surv == 7 ~ 1,
                          ARCHBOLD_surv == 8 ~ NA,
                          ARCHBOLD_surv == 9 ~ NA,
                          ARCHBOLD_surv == 2 ~ NA)) %>% 
  group_by(plant_id) %>% 
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 | ARCHBOLD_surv == 7 ~ as.numeric(census_year)),
         birth_month = case_when(ARCHBOLD_surv == 5 | ARCHBOLD_surv == 7 ~ as.numeric(census_month))) %>% 
  fill(birth_year, .direction = "updown") %>% 
  fill(birth_month, .direction = "updown") %>% 
  mutate(birth_date = as.Date(paste0("01","/",birth_month,"/",birth_year), format = c("%d/%m/%Y")),
         census_date = as.Date(paste0("01","/",census_month,"/",census_year), format = c("%d/%m/%Y")),
         age = as.period(interval(birth_date, census_date)),
         age_in_months = age %/% months(1)) %>% 
  filter(age>=0|is.na(age)) %>% 
  filter(ARCHBOLD_surv != 9 & ARCHBOLD_surv !=2 & ARCHBOLD_surv !=8)


# ggplot(PARCHA_dates)+
#   geom_histogram(aes(x = census_year, fill = factor(as.numeric(age))), position = "stack", stat = "count")
# 
# ggplot(PARCHA_dates)+
#   geom_histogram(aes(x = year(age), position = , stat = "count"))


PARCHA <- PARCHA_dates %>% 
  group_by(plant_id) %>% 
  arrange(plant_id, census_year, census_month) %>% 
  mutate(year.t1 = census_year,
         month.t1 = census_month,
         census_date.t1 = census_date,
         surv.t1 = surv,
         resprout.t1 = resprout,
         age.t1 = age,
         age_in_months.t1 = age_in_months,
         flw.t1 = flowering,
         length.t1 = length,
         width.t1 = width,
         height.t1 = height) %>% 
  mutate(year.t = dplyr::lag(year.t1,  n = 1, default = NA),
         month.t = dplyr::lag(month.t1,  n = 1, default = NA),
         census_date.t = dplyr::lag(census_date.t1,  n = 1, default = NA),
         surv.t = dplyr::lag(surv.t1,  n = 1, default = NA),
         resprout.t = dplyr::lag(resprout.t1,  n = 1, default = NA),
         age.t = dplyr::lag(age.t1,  n = 1, default = NA),
         age_in_months.t = dplyr::lag(age_in_months.t1,  n = 1, default = NA),
         flw.t = dplyr::lag(flw.t1,  n = 1, default = NA),
         length.t = dplyr::lag(length.t1,  n = 1, default = NA),
         width.t = dplyr::lag(width.t1,  n = 1, default = NA),
         height.t = dplyr::lag(height.t1,  n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% 
  select(bald, fire_unit, gap, site, plot, tag, row_id, plant_id, 
         birth_year, birth_month, birth_date,
         ARCHBOLD_surv, year.t1, month.t1, census_date.t1, surv.t1, resprout.t1, age.t1, age_in_months.t1,
         flw.t1, length.t1, width.t1 , height.t1,
         year.t, month.t, census_date.t,  surv.t,        
         resprout.t, age.t, age_in_months.t, flw.t, 
         length.t, width.t, height.t)




####################################################################################
###### Cleaning and merging together POLBAS #######################################
####################################################################################
# POLBAS data comes from nine study sites collected between 1996 and 2003. Data was collected two times a year (May and Nov.).
# The data is actually pretty clear in its current form. 
POLBAS_colnames <- c("Habitat","site", "timesincefire_2000",
                     "transect", "quadrat", "birth_year", "first_record",
                     "qsize",
                     "tag",
                     "Survival_notes;Nov2003","Diameter;Nov2003", "Length;Nov2003","Sex;Nov2003", "Flw_Stem;Nov2003",
                     "Survival_notes;May2003",
                     "Survival_notes;Nov2002","Diameter;Nov2002", "Length;Nov2002","Flw_Stem;Nov2002","Sex;Nov2002", 
                     "Survival_notes;May2002",
                     "Survival_notes;Nov2001","Survival_notes_copy;Nov2001","Diameter;Nov2001", "Length;Nov2001","Flw_Stem;Nov2001","Sex;Nov2001", 
                     "Survival_notes;May2001","Survival_notes_copy;May2001",
                     "Survival_notes;Nov2000","Survival_notes_copy;Nov2000","Diameter;Nov2000", "Length;Nov2000","Flw_Stem;Nov2000","Sex;Nov2000", 
                     "Survival_notes;May2000","Survival_notes_copy;May2000",
                     "Survival_notes;Nov1999","Survival_notes_copy;Nov1999","Diameter;Nov1999", "Standing_height;Nov1999","Flw_Stem;Nov1999","Sex;Nov1999",
                     "Survival_notes;May1999","Survival_notes_copy;May1999",
                     "Survival_notes;Nov1998","Survival_notes_copy;Nov1998","Diameter;Nov1998", "Standing_height;Nov1998","Flw_Stem;Nov1998","Sex;Nov1998",
                     "Survival_notes;Nov1997","Survival_notes_copy;Nov1997","Diameter;Nov1997", "Standing_height;Nov1997","Flw_Stem;Nov1997","Sex;Nov1997",
                     "Survival_notes;May1997", "Survival_notes_copy;May1997",
                     "Survival_notes;Nov1996","Survival_notes_copy;Nov1996","Diameter;Nov1996", "Standing_height;Nov1996","Flw_Stem;Nov1996","Sex;Nov1996",
                     "stage00","stage99","stage98","stage97","stage96","lnfs00","lnfs99","lnfs98","lnfs97", 
                     "lnfs96","chht9900","chht9899","chht9798","chht9697","size00","size99","size98","size97","size96","filter_$")

POLBAS_renamed <- POLBAS_raw

colnames(POLBAS_renamed) <- POLBAS_colnames

POLBAS_dates <- POLBAS_renamed %>% 
  dplyr::select(!starts_with("stage") & !starts_with("lnfs") & !starts_with("chht") &  !starts_with("size")) %>%  
  dplyr::select(-qsize, -`filter_$`, -timesincefire_2000) %>% 
  mutate(row_id = row_number()) %>% 
  mutate(bald = site) %>% 
  mutate(plant_id = paste(Habitat,site, transect, quadrat, tag, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character) ) %>% 
  pivot_longer(cols = !c(plant_id, bald,site, Habitat, transect, quadrat, tag, row_id, birth_year, first_record), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(census_month = case_when( grepl( "Nov", measurement) ~ 11,
                                   grepl( "May", measurement) ~ 5),
         measurement = str_split(measurement, ";", simplify = T)[, 1])  %>% 
  pivot_wider(names_from = measurement, values_from = value) %>% 
  mutate(surv_biannual = case_when(Survival_notes == "new unknown" ~ 1,
                          Survival_notes == "new seedling" ~ 1,
                          Survival_notes == "new mature" ~ 1,
                          Survival_notes == "new plant" ~ 1,
                          Survival_notes == "seedling" ~ 1,
                          Survival_notes == "survived" ~ 1, 
                          Survival_notes == "Survived" ~ 1, 
                          Survival_notes == "alive" ~ 1, 
                          Survival_notes == "4" ~ 1,
                          Survival_notes == "4.0" ~ 1,
                          Survival_notes == "dead" ~ 0, 
                          Survival_notes == "died" ~ 0, 
                          Survival_notes == "Died" ~ 0, 
                          Survival_notes == "previously dead" ~ NA,
                          Survival_notes == "prev dead" ~ NA,
                          Survival_notes == "6.0" ~ NA,
                          Survival_notes == "tag not found" ~ NA)) %>% 
# for this, I am taking only the yearly fall census, so plants that die in May are now added to the list of those that died in November
  group_by(plant_id) %>%
  arrange(plant_id,census_year,census_month) %>% 
  mutate(surv = case_when(dplyr::lag(census_month, n = 1) == 5 &  dplyr::lag(surv_biannual, n = 1) == 0 & is.na(surv_biannual) ~ 0, 
                          TRUE ~ surv_biannual)) %>% 
  mutate(birth_year = case_when(birth_year == 0 ~ NA,
                                Survival_notes == "new seedling" ~ census_year,
                                Survival_notes == "seedling" ~ census_year,
                                TRUE ~ birth_year)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  filter(census_month == 11) %>% ungroup() 

    


# POLBAS_age_summary <- POLBAS_dates %>%
#   group_by(plant_id) %>%
#   summarize(age = sum(surv))
# 
# 
# 
# POLBAS_sex_summary <- POLBAS_dates %>%
#   group_by(plant_id) %>%
#   summarize(sex_list = (unique(Sex)),
#             sex_count = length(unique(Sex)))



POLBAS <- POLBAS_dates %>% 
  mutate(across(c(census_year, birth_year, surv, Length, Standing_height, Diameter), as.numeric)) %>% 
  group_by(plant_id) %>% 
  arrange(plant_id, census_year) %>% 
  mutate(year.t1 = census_year,
         surv.t1 = surv,
         flw.t1 = Flw_Stem,
         length.t1 = Length,
         height.t1 = Standing_height,
         width.t1 = Diameter,
         age.t1 = census_year - birth_year) %>% 
  mutate(year.t = dplyr::lag(year.t1,  n = 1, default = NA),
         surv.t = dplyr::lag(surv.t1,  n = 1, default = NA),
         flw.t = dplyr::lag(flw.t1,  n = 1, default = NA),
         length.t = dplyr::lag(length.t1,  n = 1, default = NA),
         height.t = dplyr::lag(height.t1,  n = 1, default = NA),
         width.t = dplyr::lag(width.t1,  n = 1, default = NA),
         age.t = dplyr::lag(age.t1, n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% ungroup() %>% 
  select(plant_id, bald,site, Habitat, transect, quadrat, tag, row_id, birth_year, first_record,
         year.t1, surv.t1, age.t1, length.t1, height.t1, width.t1, 
         year.t, surv.t, age.t1, length.t, height.t, width.t)




####################################################################################
###### Cleaning and merging together HYPCUM  #######################################
####################################################################################
# HYPCUM data comes from 14 study sites, reduced to 9 sites after 2017. Data was collected annually in August, although some earlier censuses include recruit at different time points. collected between 1995 and 2022 
# reproductive data is messy because in some years they recorded number of reproductive structures, and in others, they only record flowering status (0/1). In some years they only did full counts for a subset of plants. At least in some cases, these are recorded within the notes for that year.
# Metadata also says that site "103" in older data is equivalent to bald 45.
HYPCUM_colnames <- c("site", "gap", "quadrat", "tag", "TP", "pulled",
                     "ARCHBOLD_surv;Aug2022","Height;Aug2022","Stem_count;Aug2022","repro_count;Aug2022","Notes;Aug2022",
                     "ARCHBOLD_surv;Jul2021","Height;Jul2021","Stem_count;Jul2021","repro_count;Jul2021","distance_to_Rosemary;July2021", "Notes;Jul2021",
                     "ARCHBOLD_surv;Jul2020","Height;Jul2020","Stem_count;Jul2020","repro_count;Jul2020","Notes;Jul2020",
                     "ARCHBOLD_surv;Jul2019", "distance_to_Rosemary;Jul2020", "cohabitants_in_20cm;Jul2020", "fire_check;Jul2020", "Height;Jul2019","Stem_count;Jul2019","Flw_count;Jul2019","Notes;Jul2019",
                     "ARCHBOLD_surv;Jul2018","Height;Jul2018","Stem_count;Jul2018","repro_count;Jul2018",
                     "ARCHBOLD_surv;Jul2017", "Notes;Jul2018",
                     "filter_$", "year","Height;Jul2017", "ltreb",
                     "Stem_count;Jul2017","Flw_count;Jul2017","Notes;Jul2017",
                     "ARCHBOLD_surv;Jul2016","Height;Jul2016","Stem_count;Jul2016","repro_count;Jul2016","Notes;Jul2016",
                     "ARCHBOLD_surv;Jul2015","Height;Jul2015","Stem_count;Jul2015","repro_count;Jul2015","Notes;Jul2015",
                     "ARCHBOLD_surv;Jul2014","Height;Jul2014","Stem_count;Jul2014","repro_count;Jul2014","Notes;Jul2014",
                     "ARCHBOLD_surv;Aug2013","Height;Aug2013","Stem_count;Aug2013","repro_count;Aug2013","Notes;Aug2013",
                     "ARCHBOLD_surv;Aug2012","Height;Aug2012","Stem_count;Aug2012","repro_count;Aug2012","Notes;Aug2012",
                     "ARCHBOLD_surv;Aug2011","Height;Aug2011","Stem_count;Aug2011","repro_count;Aug2011","Notes;Aug2011",
                     "ARCHBOLD_surv;Aug2010","Height;Aug2010","Stem_count;Aug2010","repro_count;Aug2010","Notes;Aug2010",
                     "ARCHBOLD_surv;Jul2009","Height;Jul2009","Stem_count;Jul2009","repro_count;Jul2009",
                     "ARCHBOLD_surv;Jul2008","Height;Jul2008","Stem_count;Jul2008","repro_count;Jul2008",
                     "ARCHBOLD_surv;Feb2008",
                     "ARCHBOLD_surv;Aug2007","Height;Aug2007","Stem_count;Aug2007","Flw_count;Aug2007",
                     "oldtag","burn_year", "burn97", "burn99", "burn04", "burn10", "burn15", "burn17",
                     "Notes;Jul2009",
                     "ARCHBOLD_surv;Feb2007",
                     "ARCHBOLD_surv;Aug2006", "ARCHBOLD_surv;Feb2006",
                     "ARCHBOLD_surv;Aug2005", "ARCHBOLD_surv;Feb2005",
                     "ARCHBOLD_surv;Aug2004", "ARCHBOLD_surv;Feb2004",
                     "ARCHBOLD_surv;Aug2003", "ARCHBOLD_surv;Feb2003",
                     "ARCHBOLD_surv;Aug2002", "ARCHBOLD_surv;Feb2002",
                     "ARCHBOLD_surv;Aug2001", "ARCHBOLD_surv;Feb2001",
                     "ARCHBOLD_surv;Aug2000", "ARCHBOLD_surv;Feb2000",
                     "ARCHBOLD_surv;Aug1999", "ARCHBOLD_surv;Feb1999",
                     "ARCHBOLD_surv;Aug1998", "ARCHBOLD_surv;Feb1998",
                     "ARCHBOLD_surv;Aug1997", "ARCHBOLD_surv;Feb1997",
                     "ARCHBOLD_surv;Aug1996", "ARCHBOLD_surv;Feb1996",
                     "ARCHBOLD_surv;Aug1995", "ARCHBOLD_surv;Feb1995",
                     "ARCHBOLD_surv;Aug1994",
                     "hurricane",
                     "Height;Aug1994",
                     "Height;Aug1995",
                     "Height;Aug1996",
                     "Height;Aug1997",
                     "Height;Aug1998",
                     "Height;Aug1999",
                     "Height;Aug2000",
                     "Height;Aug2001",
                     "Height;Aug2002",
                     "Height;Aug2003",
                     "Height;Aug2004",
                     "Height;Aug2005",
                     "Height;Aug2006",
                     "repro_count;Aug1994",
                     "repro_count;Aug1995",
                     "repro_count;Aug1996",
                     "repro_count;Aug1997",
                     "repro_count;Aug1998",
                     "repro_count;Aug1999",
                     "est_repro_count;Aug2000",
                     "repro_count;Aug2000",
                     "repro_count;Aug2001",
                     "est_repro_count;Aug2002",
                     "repro_count;Aug2002",
                     "repro_count;Aug2003",
                     "repro_count;Aug2004",
                     "repro_count;Aug2005",
                     "repro_count;Aug2006",
                     "Stem_count;Aug1994",
                     "Stem_count;Aug1995",
                     "Stem_count;Aug1996",
                     "Stem_count;Aug1997",
                     "Stem_count;Aug1998",
                     "Stem_count;Aug2000",
                     "Stem_count;Aug2003",
                     "Stem_count;Aug2004",
                     "Stem_count;Aug2005",
                     "Stem_count;Aug2006",
                     "Distance;Aug1994",
                     "Distance;Aug1995",
                     "Distance;Aug1996",
                     "Distance;Aug1997",
                     "sum_stem_length;Aug1995",
                     "sum_stem_length;Aug1996",
                     "classes1994","classes1995","classes1996","classes1997","classes1998","classes1999","classes2000","classes2001",
                     "errors",
                     "neig0201", "neig0802", "neig0801","roseclas", "lichen03","pb03","pc03","rosdis","rosdiscl","rosesize",
                     "oakdis", "oaksize", "herb0801",
                     colnames(HYPCUM_raw)[194:233])
                     
                    
                     
                     

HYPCUM_renamed <- HYPCUM_raw

colnames(HYPCUM_renamed) <- HYPCUM_colnames


HYPCUM_dates <- HYPCUM_renamed %>% 
  dplyr::select(!starts_with("neig") & !starts_with("classes") & !starts_with("Distance") &
                !starts_with("annsur") & !starts_with("agr") & !starts_with("sum_stem_length") &
                !starts_with("burn") & !starts_with("fire_check") & !starts_with("cohabitants_in_20cm")) %>%  
  dplyr::select(-errors, -`filter_$`, -pull, -pulled, -sc0720, 
                -roseclas, -rosdis, -rosdiscl, -oakdis, -rosesize, -oaksize,
                -lichen03, -pc03, -pb03, -herb0801, -hurricane, -oldtag, -ltreb, -year)  %>% 
  mutate(`repro_count;Aug2000` = case_when(!is.na(`repro_count;Aug2000`) ~ `repro_count;Aug2000`, 
                                                is.na(`repro_count;Aug2000`) ~ `est_repro_count;Aug2000`,
                                           TRUE ~ NA),
         `repro_count;Aug2000` = case_when(as.numeric(`est_repro_count;Aug2000`) >= as.numeric(`repro_count;Aug2000`) ~ `est_repro_count;Aug2000`,
                                                TRUE ~ `repro_count;Aug2000`)) %>%   # there is something weird in the 2000 and 2002 data with two columns for reproduction. The "estimated reproduction" is more completed, and the other includes some 0/1 values, so I am keeping only the estimated values for this year.
  dplyr::select(!starts_with("est_repro_count")) %>% 
  mutate(row_id = row_number()) %>% 
  mutate(site = case_when(site == "103" ~ "45", TRUE ~ site),
         bald = site) %>% 
  mutate(plant_id = paste(bald, site, gap, quadrat, tag, TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character) ) %>% 
  pivot_longer(cols = !c(plant_id, bald,site, gap, quadrat, tag, TP, row_id), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(census_month = case_when( grepl( "Aug", measurement) ~ 8,
                                   grepl( "Jul", measurement) ~ 7,
                                   grepl( "Feb", measurement) ~ 2),
         measurement = str_split(measurement, ";", simplify = T)[, 1])  %>% 
  pivot_wider(names_from = measurement, values_from = value) %>% 
  mutate(surv = case_when(ARCHBOLD_surv == 1 ~ 1, 
                          ARCHBOLD_surv == 0 ~ 0,
                          ARCHBOLD_surv == 3 ~ 1,
                          ARCHBOLD_surv == 5 ~ 1,
                          ARCHBOLD_surv == 9 ~ NA,
                          ARCHBOLD_surv == 2 ~ NA,
                          ARCHBOLD_surv == 6 ~ NA,
                          ARCHBOLD_surv == 12 ~ NA,
                          TRUE ~ as.numeric(ARCHBOLD_surv))) %>% 
  group_by(plant_id) %>% 
  # for this, I am keeping only the yearly Summer census, so plants that die in Feb are now added to the list of those that died in Aug/Jul
  arrange(plant_id,census_year,census_month) %>% 
  mutate(surv = case_when(dplyr::lag(census_month, n = 1) == 2 &  dplyr::lag(ARCHBOLD_surv, n = 1) == 0 & is.na(surv) ~ 0, 
                          TRUE ~ surv)) %>% 
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & Height <15 & repro_count == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  mutate(Notes = sub("/", " ", Notes),
         repro_count_true = case_when(repro_count >= 1 & !grepl("rep=", Notes) ~ as.numeric(repro_count), 
                                      grepl("rep=", Notes) ~ parse_number(Notes),
                                      repro_count == 0 ~ NA,
                                      TRUE ~ NA),
         repro_count_true = case_when(Notes == "rep=117+32+23=172 (1 2 subsample)" ~ 172,
                                      Notes == "pretty sure misentered, tag in feild is 600 (not 1000) rep=257" ~ 257, 
                                      TRUE ~ repro_count_true),
         repro_status = case_when(repro_count == 0 & is.na(repro_count_true) ~ 0,
                                  repro_count_true>=1 ~ 1, TRUE ~ NA),
         repro_count_true = case_when(census_year %in% c(2004:2009, 2011, 2012, 2014) ~ NA, 
                                      TRUE ~ repro_count_true)) %>% 
  filter(census_month != 2) %>% ungroup() 



HYPCUM <- HYPCUM_dates %>% 
    mutate(across(c(census_year, birth_year, surv, Height, Stem_count, repro_count_true, repro_status), as.numeric)) %>% 
    group_by(plant_id) %>% 
    arrange(plant_id, census_year) %>% 
    mutate(year.t1 = census_year,
           surv.t1 = surv,
           flw_status.t1 = repro_status,
           flw_count.t1 = repro_count_true,
           height.t1 = Height,
           stem_count.t1 = Stem_count,
           age.t1 = census_year - birth_year) %>% 
    mutate(year.t = dplyr::lag(year.t1,  n = 1, default = NA),
           surv.t = dplyr::lag(surv.t1,  n = 1, default = NA),
           flw_status.t = dplyr::lag(flw_status.t1,  n = 1, default = NA),
           flw_count.t = dplyr::lag(flw_count.t1,  n = 1, default = NA),
           height.t = dplyr::lag(height.t1,  n = 1, default = NA),
           stem_count.t = dplyr::lag(stem_count.t1,  n = 1, default = NA),
           age.t = dplyr::lag(age.t1, n = 1, default = NA)) %>%
    filter(!is.na(surv.t1)) %>% ungroup() 
  



####################################################################################
###### Processing the soil nutrient analysis data ##################################
####################################################################################
nutrients_raw <- read_csv(file = "/Users/joshuacfowler/Dropbox/UofMiami/Citrus Sampling Logistics/IFAS soil nutrient data/R8090.csv", skip = 14, skip_empty_rows = TRUE) %>% 
    rename(ID = `ID#`, Lab_Number = `Lab Number`)
  
nutrients_id_key <- read_xlsx(path = "/Users/joshuacfowler/Dropbox/UofMiami/Citrus Sampling Logistics/soil_analysis_id.xlsx", sheet = "soil nutrient analysis id") %>% 
  mutate(bald = case_when(project == "2024 rosemary bald sampling" ~ str_remove(Bald_id, ".+ ")),
         maehr_site = case_when(project == "citrus sampling" ~ str_remove(Bald_id, ".+ ")))


nutrients <- nutrients_raw %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  filter(ID %in% 1:154) %>% 
  select(-...14, -...15) %>% 
  mutate(across(ID:pH, ~as.numeric(.))) %>% 
  left_join(nutrients_id_key, by = c("ID" = "Soil_analysis_number"))

# ggplot(data = nutrients)+
#   geom_histogram(aes(x = pH, fill = project))
# 
# ggplot(data = nutrients)+
#   geom_histogram(aes(x = OrgMat, fill = project))
# 
# ggplot(data = nutrients)+
#   geom_histogram(aes(x = Fe, fill = project))

####################################################################################
###### Merging the datasets with the environmental covariates ######################
####################################################################################
ARCHBOLD_fire_to2022 <- read_xlsx(path = "~/Dropbox/UofMiami/Balds2009_FireIntensityArea_Through2022.xlsx", sheet = "Balds2009_FireIntensityArea", guess_max = 1048576)# guess_max makes the function look deeper in the columns to assign type

fire_summary <- ARCHBOLD_fire_to2022 %>% 
  group_by(Bald_U, Bald_) %>% 
  filter(INTENSITY != 0) %>% 
  summarize(last_fire = max(as.numeric(Year), na.rm = TRUE),
            time_since_fire = 2023-max(as.numeric(Year), na.rm = TRUE),
            fire_list = toString(unique(as.numeric(Year))),
            fire_frequency = length(unique(as.numeric(Year)))) %>% 
  ungroup() %>% 
  add_row(Bald_U = "1S", Bald_ = 1, last_fire = NA, time_since_fire = NA, fire_list = NA, fire_frequency = 0) %>% 
  mutate(bald = case_when(as.character(Bald_) %in% c("1", "45","95") ~ as.character(Bald_),
                          TRUE ~ as.character(Bald_U)),
         bald = case_when(Bald_U == "1S" ~ "1S", 
                          Bald_U == "95W" ~ "95W", Bald_U == "95N" ~ "95N", TRUE ~ bald),
         bald_simple = as.character(parse_number(bald)))
  

# Now getting the relative elevation from the old fire history file

elev.df <- read_xlsx(path = "~/Dropbox/UofMiami/Experiment Set up/firehistory_thru2018.xlsx", sheet = "Rx_Freq", guess_max = 1048576) %>% # guess_max makes the function look deeper in the columns to assign type
  rename(rel_elev = rel.eve) %>% 
  dplyr::select(bald, rel_elev) %>% 
  mutate(
    bald = case_when(bald == "01S" ~ "1S",
                     bald == "01N" ~ "1",
                     bald == "02" ~ "2",
                     bald == "05E" ~ "5E",
                     bald == "07N" ~ "7N",
                     bald == "35N" ~ "35",
                     bald == "65E" ~ "65",
                     bald == "85N" ~ "85",
                     bald == "45N" ~ "45",
                     bald == "70N"~ "70",
                     bald == "72N" ~ "72",
                     bald == "95" ~ "95S",
                     TRUE ~ bald),
    bald_simple = as.character(parse_number(bald)))



# Pulling out fire history relative to the observation year for each datapoint
# needing to change a few of the names of balds to make them line up with names in elev and fire datasets. this also should be double checked because some you have to choose N/S
ERYCUN_covariates <- ERYCUN %>% 
  left_join(elev.df) %>% 
  left_join(fire_summary, by = join_by(bald)) %>% 
    mutate(rel_elev = case_when(bald == 95 ~ 0.68, TRUE ~ rel_elev)) %>% 
  filter(Site_tag == "archbold")

# for now, we don't know if some of these balds were bald 1N or 1S, we just know Bald 1, so I'm gonna pick omang the 4 that we don't know.

PARCHA_covariates <- PARCHA %>% 
  mutate(year.t1 = as.numeric(year.t1)) %>% 
  left_join(elev.df, by = join_by(bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("1S", "24S", "46E", "49E"))) %>% select(-bald.y) %>%
  left_join(fire_summary, by = join_by(bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("24S", "46E", "49E"))) %>% select(-bald.y)



#
HYPCUM_covariates <- HYPCUM %>% 
  left_join(elev.df, by = c( "bald" = "bald_simple"), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("62N", "1S", "24S", "46E", "49E"))) %>% select(-bald.y) %>% 
  left_join(fire_summary, by = c( "bald" = "bald_simple"), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("62N", "1S", "24S", "46E", "49E"))) %>% select(-bald.y) %>% 
  mutate(rel_elev = case_when(bald == 95 ~ 0.68, TRUE ~ rel_elev))



# CHAFLO is not easy to connect yet because sites are listed as being at burn units. there is actually good information, but I need to follow up on this
# CHAFLO_covariates <- CHAFLO %>% 
#




#
CHAFAS_covariates <- CHAFAS_growth %>% 
  mutate(year.t1 = census_year) %>% 
  left_join(elev.df, by = join_by( bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("46E"))) %>% select(-bald.y) %>% 
  left_join(fire_summary, by = join_by(bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c( "46E"))) %>% select(-bald.y)
  



#
BALANG_covariates <- BALANG_growth %>% 
  mutate(year.t1 = census_year) %>% 
  left_join(elev.df, by = join_by( bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("46E"))) %>% select(-bald.y) %>% 
  left_join(fire_summary, by = join_by(bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c( "46E"))) %>% select(-bald.y)




# LIAOHL will also be hard to connect to specific balds and need to follow up
# LIAOHL_covariates <- LIAOHL %>% 


#
POLBAS_covariates <- POLBAS %>% 
  mutate(year.t1 = as.numeric(year.t1)) %>% 
  left_join(elev.df, by = join_by(bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("29E", "46E", "95W"))) %>% select(-bald.y) %>%
  left_join(fire_summary, by = join_by(bald == bald_simple), relationship = "many-to-many") %>% 
  filter(!(bald.y %in% c("29E", "46E", "95W"))) %>% select(-bald.y)





#
last_fire_fx <- function(x,y){ max(as.numeric(x)[as.numeric(x)<=y])}
fire_frequency_fx <- function(x,y){length(as.numeric(x)[as.numeric(x)<=y])/(y - 1950)}



#
ERYCUN_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(ERYCUN_covariates$fire_list), ", "), y= ERYCUN_covariates$year.t1))
ERYCUN_covariates$time_since_fire_actual <- ERYCUN_covariates$year.t1 - ERYCUN_covariates$last_fire_actual
ERYCUN_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(ERYCUN_covariates$fire_list), ", "), y= ERYCUN_covariates$year.t1))

ERYCUN_covariates <- ERYCUN_covariates %>% 
  mutate(fire_frequency_actual = case_when(bald == "1S" ~ 0/(year.t1-1950),
                                           TRUE ~ fire_frequency_actual)) %>% 
  dplyr::select(-last_fire, -fire_frequency, -time_since_fire)

#
PARCHA_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(PARCHA_covariates$fire_list), ", "), y= PARCHA_covariates$year.t1))
PARCHA_covariates$time_since_fire_actual <- PARCHA_covariates$year.t1 - PARCHA_covariates$last_fire_actual
PARCHA_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(PARCHA_covariates$fire_list), ", "), y= PARCHA_covariates$year.t1))





#
HYPCUM_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(paste0(HYPCUM_covariates$fire_list, ", 1950")), ", "), y= HYPCUM_covariates$year.t1))
HYPCUM_covariates$time_since_fire_actual <- HYPCUM_covariates$year.t1 - HYPCUM_covariates$last_fire_actual
HYPCUM_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(HYPCUM_covariates$fire_list), ", "), y= HYPCUM_covariates$year.t1))


#
CHAFAS_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(paste0(CHAFAS_covariates$fire_list, ", 1950")), ", "), y= CHAFAS_covariates$year.t1))
CHAFAS_covariates$time_since_fire_actual <- CHAFAS_covariates$year.t1 - CHAFAS_covariates$last_fire_actual
CHAFAS_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(CHAFAS_covariates$fire_list), ", "), y= CHAFAS_covariates$year.t1))



#
BALANG_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(paste0(BALANG_covariates$fire_list, ", 1950")), ", "), y= BALANG_covariates$year.t1))
BALANG_covariates$time_since_fire_actual <- BALANG_covariates$year.t1 - BALANG_covariates$last_fire_actual
BALANG_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(BALANG_covariates$fire_list), ", "), y= BALANG_covariates$year.t1))




#
POLBAS_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(paste0(POLBAS_covariates$fire_list, ", 1950")), ", "), y= POLBAS_covariates$year.t1))
POLBAS_covariates$time_since_fire_actual <- POLBAS_covariates$year.t1 - POLBAS_covariates$last_fire_actual
POLBAS_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(POLBAS_covariates$fire_list), ", "), y= POLBAS_covariates$year.t1))


#########

write_csv(ERYCUN_covariates, paste0(filepath,"/cleaned_data", "/ERYCUN_covariates.csv"))
write_csv(HYPCUM_covariates, paste0(filepath,"/cleaned_data", "/HYPCUM_covariates.csv"))
write_csv(POLBAS_covariates, paste0(filepath,"/cleaned_data", "/POLBAS_covariates.csv"))
write_csv(BALANG_covariates, paste0(filepath,"/cleaned_data", "/BALANG_covariates.csv"))
write_csv(CHAFAS_covariates, paste0(filepath,"/cleaned_data", "/CHAFAS_covariates.csv"))
write_csv(PARCHA_covariates, paste0(filepath,"/cleaned_data", "/PARCHA_covariates.csv"))
write_csv(CHAFLO, paste0(filepath,"/cleaned_data", "/CHAFLO.csv"))
write_csv(LIAOHL, paste0(filepath,"/cleaned_data", "/LIAOHL.csv"))



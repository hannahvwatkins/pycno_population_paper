#cleaning code for the small datasets within the original IUCN dataset and 
#additional small survey from SFU-Cote
library(tidyverse)
library(janitor)

#isabelle's survey data, which we'll just tack on to the IUCN small datasets
cote <- read_csv("../data/cote_survey_data.csv") %>% 
  clean_names()

#all iucn data
iucn <- read_csv("../data/IUCN_data.csv",
                      #deal with parsing issues from so many NAs
                      guess_max = 15000) %>% 
  clean_names() %>% 
  filter(region_deep == "Canada + Puget Sound") %>% 
  #remove OW since we have the more up to date version
  filter(source !="OceanWise_CitSciDive_BC") %>% 
  #remove reef data since we have the more up to date version
  filter(source != "REEF_Dive_Namerica") %>%   
  bind_rows(cote) %>% 
  mutate(date = case_when(is.na(date) ~ as.character(event_year),
                          TRUE ~ date)) %>% 
  #since there aren't consistent site names used, we'll make a location code
  #based on the lat lons
  unite("location_code_numeric", dec_lat:dec_long, remove = FALSE) %>% 
  #remove washington which is lumped in here
  filter(!grepl("WDFW",source)) %>% 
  #indicate which years the crash occurred in two different formats
  mutate(wasting = case_when(event_year < 2014 ~ 0,
                             TRUE ~ 1), #formatted for pre/post wasting
         wasting_presence = case_when(event_year == 2014 | 
                                        event_year == 2015 ~ 1,
                                      TRUE ~ 0), #formatted for years of decline
         #and create a variable to determine if multiple dives were conducted
         #at the same site on the same day
         site_ymd = paste(location_code_numeric, date, sep = "_")) %>% 
  #only use on observation per site per day to avoid pseudoreplication
  group_by(site_ymd) %>% 
  mutate(n = n(),
         multiple_surveys = case_when(n>1 ~ 1,
                                      TRUE ~ 0)) %>% 
  slice_max(pres_abs, with_ties = FALSE) %>% 
  ungroup() 

iucn_spatial <- iucn %>% 
  filter(source != "MARINe_CitSciObs_Global" & 
           source != "iNaturalist_CitSciObs_Global")
#write_csv(iucn_spatial, "../data/cleaned_iucn_data_for_spatial.csv")

iucn_2000 <- iucn_spatial %>% 
  #and because there are so few surveys pre 2000, we have to trim down the data
  #this ends up removing just the early Watson dataset and two dives from 2021
  #in the cote dataset
  filter(event_year > 1999 & event_year < 2021)
#write_csv(iucn_2000, "../data/cleaned_iucn_data.csv")

#in the above code we removed both the iNat and the MARINe global datasets 
#because there is no measure of effort for either, and in the case of iNat there
#are no absences recorded. Both datasets still provide helpful info for the 
#spatial analyses though. We have the more up to date iNat in a separate csv but
#we'll also extract the MARINe data here for use later
marine <- iucn %>% 
  filter(source == "MARINe_CitSciObs_Global")
#write_csv(marine, "../data/cleaned_marine_global.csv")

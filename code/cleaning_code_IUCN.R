#cleaning code for the small datasets within the original IUCN dataset and 
#additional small survey from SFU-Cote
library(tidyverse)
library(janitor)
library(sf)

#load in region shapefile
bioregions <- read_sf("data/pycno_bioregions.gpkg") 
#additional dive survey data used in COSEWIC assessment for presence/absence
#but not in paper for abundance estimates, as the dataset is too small to be
#modelled on its own and uses time instead of area as a metric of survey effort
#so can't be combined with the other small datasets
cote <- read_csv("data/cote_survey_data.csv") %>% 
  clean_names() %>% 
  #remove the misnamed column
  select(-totalpycnos) %>% 
  #rename the columns to match the iucn dataframe
  rename(totalpycnos = abundance) %>% 
  filter(!(is.na(totalpycnos)))

#all data originally used in the IUCN assessment, as stored in a Dryad 
#repository for Hamilton et al., 2021: https://doi.org/10.5061/dryad.9kd51c5hg
iucn <- read_csv("data/IUCN_data.csv",
                      #deal with parsing issues from so many NAs
                      guess_max = 15000) %>% 
  clean_names() %>% 
  #include only Canadian data
  filter(region_deep == "Canada + Puget Sound") %>% 
  #remove OW since we have the more up to date version
  filter(source !="OceanWise_CitSciDive_BC") %>% 
  #remove reef data since we have the more up to date version
  filter(source != "REEF_Dive_Namerica") %>%   
  bind_rows(cote) %>% 
  mutate(date = case_when(is.na(date) ~ as.character(event_year),
                          TRUE ~ date)) %>% 
  #since there aren't consistent site names used, we'll make a location code
  #based on the lat/lons
  unite("location_code_numeric", dec_lat:dec_long, remove = FALSE) %>% 
  #remove washington which is lumped in here with the Salish Sea data
  filter(!grepl("WDFW",source)) %>% 
  #and remove the presence-only data
  filter(source != "MARINe_CitSciObs_Global" & 
           source != "iNaturalist_CitSciObs_Global") %>% 
  #indicate which years the crash occurred in two different formats
  mutate(wasting = case_when(event_year < 2014 ~ 0,
                             TRUE ~ 1), #formatted for pre/post wasting
         wasting_presence = case_when(event_year == 2014 | 
                                        event_year == 2015 ~ 1,
                                      TRUE ~ 0), #formatted for years of decline
         #and create a variable to determine if multiple dives were conducted
         #at the same site on the same day
         site_ymd = paste(location_code_numeric, date, sep = "_")) %>% 
  #only use on observation per site per day to avoid pseudoreplication - would
  #prefer to model these as random effects, but there are too few to do so
  #reliably, so we'll just take the max value for a given day/site combo, since
  #it's more likely that a diver missed a pycno than made one up
  #we'll make a new variable indicating whether a day had multiple dives though
  #so we can model the possible negative bias of a site that only had one diver
  group_by(site_ymd) %>% 
  mutate(n = n(),
         multiple_surveys = case_when(n>1 ~ 1,
                                      TRUE ~ 0)) %>% 
  slice_max(totalpycnos, with_ties = FALSE) %>% 
  ungroup() %>% 
  #convert to spatial object to sort lat/lons into regions
  st_as_sf(coords = c('dec_long','dec_lat'),
         crs = st_crs(4326)) %>% 
  st_join(bioregions) %>% 
  mutate(region_numeric = case_when(region == "Haida Gwaii" ~ 1,
                                    region == "North and Central Coast" ~ 2,
                                    region == "Johnstone Strait/Northern Salish Sea" ~ 3,
                                    region == "West Coast Vancouver Island" ~ 4,
                                    region == "Salish Sea" ~ 5))  

#write_csv(iucn, "data/cleaned_iucn_data_for_spatial.csv")

iucn_2000 <- iucn %>% 
  #and because there are so few surveys pre 2000, we have to trim down the data
  #this ends up removing just the early Watson dataset and two dives from 2021
  #in the cote dataset
  filter(event_year > 1999 & event_year < 2021) %>% 
  #get the original lat/lon cols back
  mutate(lon = sf::st_coordinates(.)[,1],
         lat = sf::st_coordinates(.)[,2])
#write_csv(iucn_2000, "data/cleaned_iucn_data.csv")

#in the above code we removed both the iNat and the MARINe global datasets 
#because there is no measure of effort for either, and in the case of iNat there
#are no absences recorded. Both datasets still provide helpful info for the 
#spatial analyses though. We have the more up to date iNat in a separate csv but
#we'll also extract the MARINe data here for use later
marine <- iucn %>% 
  filter(source == "MARINe_CitSciObs_Global")
#write_csv(marine, "data/cleaned_marine_global.csv")

#cleaning code for REEF raw data
library(tidyverse)
library(janitor)
library(sf)

#load in region shapefile
bioregions <- read_sf("data/pycno_bioregions.gpkg") 

#load in REEF site locations - note this is not publicly available and will need
#to be requested from REEF if you need access
locations <- read_csv("data/REEF_locations.csv") %>% 
  filter(lat != "NULL") %>% 
  #this function is being annoying so this is going to be hella janky
  separate(lat, into = c("lat_degree","lat_minute","lat_decimal")) %>% 
  separate(lon, into = c(NA,"lon_degree","lon_minute","lon_decimal")) %>% 
  unite("lat_minute", lat_minute:lat_decimal, sep = ".") %>% 
  unite("lon_minute", lon_minute:lon_decimal, sep = ".") %>% 
  mutate(lat_min = as.numeric(lat_minute)/60,
         lat = as.numeric(lat_degree) + lat_min,
         lon_min = as.numeric(lon_minute)/60,
         lon = (as.numeric(lon_degree) + lon_min)*-1,
         geogr = geogid) %>% 
  select(lat, lon, geogr)

#dates are coded in two different ways because people who do data entry hate 
#people who have to do stats, so we'll clean them separately
pycno_presence_slash <- read_csv("data/REEF_raw_presence.csv") %>% 
  clean_names() %>% 
  filter(grepl("/", date)) %>% 
  separate(date, into = c("month", "day", "year"), sep = "/") %>% 
  #all surveys are post-1999 so reformat year to full number
  mutate(year = as.numeric(paste("20",year, sep = "")),
         month = as.numeric(month),
         day = as.numeric(day))
pycno_presence_dash <- read_csv("data/REEF_raw_presence.csv") %>% 
  clean_names() %>% 
  filter(grepl("-", date)) %>% 
  separate(date, into = c("year", "month", "day"), sep = "-") %>% 
  mutate(year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day))

#and combine both  
pycno_presence <- bind_rows(pycno_presence_slash, pycno_presence_dash) %>% 
  mutate(site_ymd = paste(geogr, day, month, year, sep = "_")) %>% 
  #select only obs from BC which start with a 1
  filter(grepl("^1", geogr))

#same issue with the dates in the abundance file
pycno_abund_slash <- read_csv("data/REEF_raw_abundance.csv") %>% 
  clean_names() %>% 
  filter(grepl("/", date)) %>% 
  separate(date, into = c("month", "day", "year"), sep = "/") %>% 
  #all surveys are post-1999 so reformat year to full number
  mutate(year = as.numeric(paste("20",year, sep = "")),
         month = as.numeric(month),
         day = as.numeric(day))
pycno_abund_dash <- read_csv("data/REEF_raw_abundance.csv") %>% 
  clean_names() %>% 
  filter(grepl("-", date)) %>% 
  separate(date, into = c("year", "month", "day"), sep = "-") %>% 
  #all surveys are post-1999 so reformat year to full number
  mutate(year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day))

#and join
pycno_abund <- bind_rows(pycno_abund_slash, pycno_abund_dash) %>% 
  dplyr::rename(geogr = geozone) %>% 
  mutate(site_ymd = paste(geogr, day, month, year, sep = "_")) %>% 
  #select only obs from BC
  filter(grepl("^1", geogr)) %>% 
  #select only the variables that aren't in the survey dataset, since there are
  #a lot of repeats and we only need the form ID to join
  dplyr::select(form, abundance)

#make final version
pycno <- full_join(pycno_abund, pycno_presence, 
                   by = c("form" = "formid"))  %>% 
  #remove true duplicates (i.e., only if the entire row is duplicated)
  distinct(.) %>% 
  #fill in NAs for abundance with zeros, since they're created by a lack of 
  #matching observation in the pycno_abund df
  mutate(abundance = case_when(is.na(abundance) ~ 0,
                               TRUE ~ abundance),
         #reformat abundance levels for Stan, which indexes starting with 1
         abundance_recode = abundance + 1,
         #code presence/absence
         presence = case_when(abundance == 0 ~ 0,
                              TRUE ~ 1),
         #reformat expertise level for Stan
         exp_binary = case_when(exp == "N" ~ 0,
                                exp == "E" ~ 1)) %>% 
  #since not many sites were surveyed more than once on the same day, including
  #a random effect of site_day is really challenging and the model can't
  #converge properly
  #so we're just going to take the highest score for a given site_day
  group_by(site_ymd) %>% 
  mutate(num_surveys = n(),
            multiple_surveys = case_when(num_surveys == 1 ~ 0,
                                         TRUE ~ 1)) %>% 
  slice_max(abundance, with_ties = FALSE) %>% 
  ungroup() %>% 
  #finally, create a variable for pre/post epidemic
  #I think part of the reason the model struggles to converge is that the
  #temporal variance is really small EXCEPT for when the crash happens and it
  #stresses Stan out - if we include a variable for pre/post wasting, this will
  #enable the annual estimates to work better
  mutate(wasting = case_when(year < 2014 | year > 2015 ~ 0,
                           TRUE ~ 1),
         geogr_region = geogr) %>% 
  separate(geogr_region, into = c('geogr_region', NA), sep = 2) %>% 
  left_join(locations, by = "geogr") 

#there are some NAs for lat/lons but we still know the regions they're in so we 
#have to split the df into two and merge after

pycno_latlon <- pycno %>% 
  filter(!is.na(lat)) %>% 
  #convert to spatial object to sort lat/lons into regions
  st_as_sf(coords = c('lon','lat'),
           crs = st_crs(4326)) %>% 
  st_join(bioregions) %>% 
  mutate(region_numeric = case_when(region == "Haida Gwaii" ~ 1,
                                    region == "North and Central Coast" ~ 2,
                                    region == "Johnstone Strait/Northern Salish Sea" ~ 3,
                                    region == "West Coast Vancouver Island" ~ 4,
                                    region == "Salish Sea" ~ 5))  
pycno_no <- pycno %>% 
  filter(is.na(lat)) %>% 
  mutate(region_numeric = case_when(geogr_region == 11 | geogr_region == 13 ~ 3,
                                    geogr_region == 19 ~ 2,
                                    geogr_region == 12 ~ 4))

pycno_regions <- bind_rows(pycno_latlon, pycno_no) %>% 
  select(-(lat:lon)) %>% 
  st_drop_geometry()
#write_csv(pycno_regions, "data/cleaned_reef_data.csv")


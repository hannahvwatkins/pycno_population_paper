library(tidyverse)
library(janitor)
library(sf)

#load in region shapefile
bioregions <- read_sf("data/pycno_bioregions.gpkg") 
#benthic habitat modeling data spanning from 2013-2021, only used for 
#presence/absence analysis in the COSEWIC report and not the paper
bhm <- read_csv("data/dfo_benthic_habitat_mapping.csv") %>% 
  clean_names() %>% 
  rename(presence = pyc) %>% 
  #and let's use the deep lat/lon as our main option (these are small 
  #transects so it doesn't really matter)
  mutate(lat = lat_deep,
         lon = lon_deep,
         #and we'll use the mean corrected depth as our main depth
         depth = mean_cor_depth_m,
         source = "bhm") %>% 
  #prevent pseudoreplication by only using one value for each date/site combo
  group_by(day, month, year, lat, lon) %>% 
  mutate(num_surveys = n(),
         multiple_surveys = case_when(num_surveys > 1 ~ 1,
                                      TRUE ~ 0)) %>% 
  slice_max(presence, with_ties = FALSE) %>% 
  ungroup()
#write_csv(bhm, "data/cleaned_dfo_benthic_habitat_mapping.csv")

#multispecies dive surveys, has counts but only from 2016 onward
multi <- read_csv("data/dfo_multispecies.csv") %>% 
  clean_names() %>% 
  rename(num_pycno = sum_of4xe) %>% 
  filter(!is.na(num_pycno))%>% 
  mutate(presence = case_when(num_pycno > 0 ~ 1,
                              TRUE ~ 0),
         source = "multi",
         lat = lat_start,
         lon = lon_start,
         depth = max_of_cor_depth_m) %>% 
  rename(num_quadrats = quadrat_count) %>% 
  group_by(day, month, year, lat, lon) %>% 
  mutate(num_surveys = n(),
         multiple_surveys = case_when(num_surveys > 1 ~ 1,
                                      TRUE ~ 0)) %>% 
  slice_max(presence, with_ties = FALSE) %>% 
  ungroup() %>%
  #convert to spatial object to sort lat/lons into regions
  st_as_sf(coords = c('lon_start','lat_start'),
           crs = st_crs(4326)) %>% 
  st_join(bioregions) %>% 
  mutate(region_numeric = case_when(region == "Haida Gwaii" ~ 1,
                                    region == "North and Central Coast" ~ 2,
                                    region == "Johnstone Strait/Northern Salish Sea" ~ 3,
                                    region == "West Coast Vancouver Island" ~ 4,
                                    region == "Salish Sea" ~ 5))  
  
#write_csv(multi, "data/cleaned_dfo_multispecies.csv")

#abalone data from 2006-2021 with some gaps
#note these data are not included in any maps as the locations are considered
#too sensitive to be released by DFO
aba1 <- read_csv("data/dfo_abalone_2006-2018.csv") %>% 
  clean_names() %>% 
  mutate(sheet_id = "A") %>% 
  unite(dummy_id, c(dummy_site_id, sheet_id))

aba2 <- read_csv("data/dfo_abalone_2019.csv") %>% 
  clean_names()%>% 
  mutate(sheet_id = "B") %>% 
  unite(dummy_id, c(dummy_id, sheet_id))

aba3 <- read_csv("data/dfo_abalone_2021.csv") %>% 
  clean_names()%>% 
  mutate(sheet_id = "C") %>% 
  unite(dummy_id, c(dummy_id, sheet_id))

aba <- bind_rows(aba1,aba2,aba3) %>% 
  filter(!is.na(num_pycno)) %>% 
  mutate(presence = case_when(num_pycno > 0 ~ 1,
                              TRUE ~ 0),
         source = "abalone") %>% 
  rename(num_quadrats = num_quads) %>% 
  #add in region based on stat areas since we don't have exact locations
  mutate(region_numeric = case_when(stat_area == 2 ~ 1,
                            stat_area > 2 & stat_area <= 11 ~ 2,
                            stat_area > 11 & stat_area <= 17 ~ 3,
                            stat_area > 22 & stat_area <= 27 ~ 4,
                            stat_area > 17 & stat_area <=20 ~ 5))
#write_csv(aba, "data/cleaned_dfo_abalone.csv")


#make one master dataset for all DFO dives
dfo <- bind_rows(aba, bhm, multi) %>% 
  dplyr::select(year, month, day, stat_area, sub_area, num_quadrats,
                num_pycno, presence, bottom_time,
                lat, lon, depth, source)
#write_csv(dfo, "data/cleaned_dfo_dives.csv")

#cleaning code for the OW Pacific Marine Life Surveys
library(tidyverse)
library(janitor)
library(lubridate)
library(sf)

#load in region shapefile
bioregions <- read_sf("data/pycno_bioregions.gpkg") 

gibbs <- read_csv("data/gibbs_survey_data.csv",
                  #max depth parsing failures can be ignored as they 
                  #automatically changed the "?"s to NAs
                  col_types = list(Comments = col_character(),
                                   `Max. depth` = col_double())) %>% 
  clean_names() %>% 
  #remove all the empty rows
  filter(!is.na(scientific_name)) %>% 
  #extract the year from the observation date
  mutate(year = year(date_observed))


#this dataset includes a lot of non-pycno info, so we're going to have to do 
#some work to get this in a format that is usable
#we essentially need to have each row correspond to a single dive, and ensure
#that dives without pycnos are counted as 0s, not removed entirely
#we'll start by extracting the info from the dives where pycnos were in fact 
#present
gibbs_pycno_only <- gibbs %>% 
  filter(scientific_name == "Pycnopodia helianthoides") %>% 
  #the data are inconsistent in that sometimes actual counts are recorded, while
  #other times an abundance "score" was recorded. We want to turn everything 
  #into this score to make it consistent across surveys
  #we want to be able to filter the actual counts into the abundance codes
  #since it's inconsistent, so we'll make a numeric col to simplify that
  #the NAs are expected, so we can ignore the warning
  mutate(tag_num = as.numeric(tags)) %>% 
  #now put everything into the appropriate abundance code
  mutate(abundance = case_when(tags == "none" ~ 0,
                               tags == "few" ~ 1,
                               tags == "present" ~ 1,
                               tags == "some" ~ 2,
                               tags == "many" ~ 3,
                               tags == "very many" ~ 4,
                               tags == "abundant" ~ 5,
                               tags == "very abundant" ~ 6,
                               tags == "5+" ~ 1,
                               tag_num < 11 ~ 1,
                               tag_num > 10 & tag_num < 26 ~ 2,
                               tag_num > 25 & tag_num < 51 ~ 3,
                               tag_num > 50 & tag_num < 101 ~ 4,
                               tag_num > 100 & tag_num < 1001 ~ 5,
                               tag_num > 1000 ~ 6),
         #and we'll make a new column here that just gives a 1 to every row
         #so that when we join it back up with the full dataset we'll be able 
         #to distinguish the dives with no pycnos from the dives with them
         presence = 1) %>% 
  #only take the cols we need to bind it back with the original dataset
  dplyr::select(date_observed, place_name, abundance, 
                presence, observer_name) %>% 
  distinct(.) %>% 
  #remove one case of a single diver recording 2 values for a single dive
  filter(!(grepl("2011-08-16", date_observed) & 
             place_name == "Sea-to-Sky, helipad mansion" & abundance == "2"))

gibbs_all <- gibbs %>% 
  #by left joining with dat, place name, and observer name, we essentially apply
  #the pycno specific info above to every row that corresponds to a single dive
  #including the rows pertaining to non-pycno species
  left_join(gibbs_pycno_only, by = c("date_observed", "place_name", 
                                     "observer_name")) %>% 
  #this means that if a row has NAs for presence, there were no pycnos on that
  #dive
  mutate(pres_abs = case_when(is.na(presence) ~ 0,
                              TRUE ~ presence),
         #and we can apply a score of 1 for dives without pycnos, NA for if 
         #there was no score/count for pycnos but we know they were present, 
         #or the original abundance score plus one if there was a count/score
         #note that we add the one because Stan indexing starts at 1, so we 
         #can't use a score of 0 for the actual zero dives and need to shift
         #everything up
         abundance_recode = case_when(pres_abs == 0 ~ 1,
                                      pres_abs == 1 & is.na(abundance) ~ 
                                        as.numeric(NA),
                                      TRUE ~ abundance + 1))

#now we'll check and see when this dataset became more standardized and when
#counts were consistently recorded
consistency <- gibbs_all %>% 
  group_by(year) %>% 
  summarize(n = n(),
            nas = sum(is.na(abundance_recode)),
            prop_na = nas/n)
#it looks like things got a lot more standardized in 2001, so
#for this analysis, we'll limit the data to 2001 for the time series analysis
#but keep the observations from pre-2001 since we can trust that all the 
#presences are real presences and this will help in our understanding of the 
#extent of occurrence 

ow_pycno <- gibbs_all %>% 
  #and we only want 1 observation per dive per diver, so we'll remove the rest
  #it doesn't matter which row we select from a given dive, since we've created
  #the same abundance code and pres_abs value for all observations within a 
  #given dive
  distinct(date_observed, place_name, observer_name, .keep_all = TRUE) %>% 
  mutate(wasting = case_when(year < 2014 ~ 0,
                             TRUE ~ 1),
         location_code = place_name) %>% 
  filter(!is.na(abundance_recode)) %>% 
  filter(!is.na(max_depth)) %>% 
  mutate(pres_abs = case_when(is.na(abundance) ~ 0,
                              abundance == 0 ~ 0,
                              TRUE ~ 1)) %>%
  #and finally, if there is more than one observation for a given dive, take 
  #the max of the abundance scores recorded and create a column to specify 
  #whether the observation is based on multiple surveys or a single survey
  #so we can include that as a covariate in the model (since it'll likely bias
  #the counts higher if there are more divers looking on a dive)
  group_by(date_observed, place_name) %>% 
  mutate(num_surveys = n(),
         multiple_surveys = case_when(num_surveys > 1 ~ 1,
                                      TRUE ~ 0)) %>% 
  slice_max(abundance_recode, with_ties = FALSE) %>% 
  ungroup() %>% 
  #get rid of the "tags" column since it doesn't necessarily refer to the pycnos
  dplyr::select(-(kingdom:common_name)) %>% 
  dplyr::select(-tags)
#write_csv(ow_pycno, "data/cleaned_oceanwise_data_for_spatial.csv")

#and filter down to the years we want for the time series analysis
ow_2001 <- ow_pycno %>% 
#remove 2022 as well since that has to be a mistake as we received this data 
#in 2021
  filter(year > 2000 & year < 2022) %>%  
  #convert to spatial object to sort lat/lons into regions
  st_as_sf(coords = c('longitude','latitude'),
         crs = st_crs(4326)) %>% 
  st_join(bioregions) %>% 
  mutate(region_numeric = case_when(region == "Haida Gwaii" ~ 1,
                                    region == "North and Central Coast" ~ 2,
                                    region == "Johnstone Strait/Northern Salish Sea" ~ 3,
                                    region == "West Coast Vancouver Island" ~ 4,
                                    region == "Salish Sea" ~ 5,
                                    #the remaining NAs fall just outside of our hand drawn boundary for region 1
                                    TRUE ~ 1)) %>% 
  #get the original lat/lon cols back
  mutate(lon = sf::st_coordinates(.)[,1],
                lat = sf::st_coordinates(.)[,2])
#write_csv(ow_2001, "data/cleaned_oceanwise_data.csv")


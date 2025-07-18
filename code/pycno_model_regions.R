
library(tidyverse) #for data wrangling
library(cmdstanr) #for running model
library(rstan) #for extracting and processing model output
library(DHARMa) #for examining residuals
library(loo) #for model comparison
source("code/custom_functions.R")

#load in data------------------------------------------------------------------
reef <- read_csv("data/cleaned_reef_data.csv")
ow <- read_csv("data/cleaned_oceanwise_data.csv",
               col_types = list(comments = col_character()))
iucn_pres_com <- read_csv("data/cleaned_iucn_data.csv") %>% 
  #remove the surveys without an area measurement since we need to be able to 
  #control for survey effort (this just removes some of the Hakai surveys and 
  #the Cote dataset)
  filter(!is.na(area)) %>% 
  filter(!is.na(depth)) %>% 
  mutate(year = event_year) %>%
  #removing pre-2009 data since there are only 7 surveys per year before that
  filter(year > 2008)
aba <- read_csv("data/cleaned_dfo_abalone.csv")
multi <- read_csv("data/cleaned_dfo_multispecies.csv") %>% 
  unite("site", c(stat_area, sub_area), remove = FALSE)#create unique site names

#prep data for Stan-------------------------------------------------------------
#these are the fixed effects matrices for the survey-level effects
#note that we need to log time and area to allow for more flexibility in the 
#relationships between survey effort and counts (e.g., a coefficient of 1 would
#represent a true offset, while a lower slope would allow for a saturating
#relationship over time/area)
X_reef <- matrix(data=c(scale(log(as.numeric(reef$btime))),
                        scale(as.numeric(reef$maxdepth)), 
                        scale(as.numeric(reef$visibility)),
                        scale(as.numeric(reef$current)),
                        reef$exp_binary, #experience level of surveyor
                        reef$multiple_surveys), #whether multiple surveys 
                 #occurred at the same time/site
                 ncol=6,nrow=nrow(reef))
X_ow <- matrix(data=c(scale(log(as.numeric(ow$bottom_time))),
                      scale(as.numeric(ow$max_depth)),
                      ow$multiple_surveys),
               ncol=3,nrow=nrow(ow))
X_iucn <- matrix(data=c(scale(as.numeric(iucn_pres_com$depth)),
                        scale(log(as.numeric(iucn_pres_com$area)))),
                 ncol=2,
                 nrow=nrow(iucn_pres_com))
X_aba <- matrix(data=c(scale(log(as.numeric(aba$num_quadrats)))),
                ncol=1,nrow=nrow(aba))
X_multi <- matrix(data=c(scale(as.numeric(multi$max_of_cor_depth_m)),
                         scale(log(as.numeric(multi$num_quadrats)))),
                  ncol=2, nrow=nrow(multi))

#create a year-level vector for whether or not a given year was pre, during, or
#post wasting
wasting_index <- reef %>% 
  select(year) %>% 
  distinct() %>% 
  arrange(year) %>% 
  mutate(wasting = case_when(year == 2014 | year == 2015 ~ 1,
                             year > 2015 ~ 2,
                             TRUE ~ 0))

#make an index linking the years in each dfo dataset to the full timespan
yr_index_reef <- reef %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

yr_index_ow <- ow %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

yr_index_iucn <- iucn_pres_com %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

yr_index_aba <- aba %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

yr_index_multi <- multi %>% 
  distinct(year) %>% 
  arrange(year) %>% 
  transmute(tt = year - 1999)

data_pycno_abun_region = list(TT=22,
                              y_reef = reef$abundance_recode,
                              N_reef = nrow(reef),
                              site_reef=as.numeric(factor(reef$geogr)),
                              N_site_reef=length(unique(reef$geogr)),
                              K_reef=length(unique(reef$abundance)),
                              X_reef=X_reef,
                              Z_reef=ncol(X_reef),
                              N_yr_reef=length(unique(reef$year)),
                              year_id_reef=as.numeric(factor(reef$year)),
                              y_ow = ow$abundance_recode,
                              N_ow = nrow(ow),
                              site_ow=as.numeric(factor(ow$place_name)),
                              N_site_ow=length(unique(ow$place_name)),
                              diver_ow=as.numeric(factor(ow$observer_name)),
                              N_dv_ow=length(unique(ow$observer_name)),
                              K_ow=length(unique(ow$abundance_recode)),
                              X_ow=X_ow,
                              Z_ow=ncol(X_ow),
                              N_yr_ow=length(unique(ow$year)),
                              year_id_ow=as.numeric(factor(ow$year)),
                              y_iucn = iucn_pres_com$totalpycnos,
                              N_iucn = nrow(iucn_pres_com),
                              site_iucn=as.numeric(factor(iucn_pres_com$location_code_numeric)),
                              N_site_iucn=length(unique(iucn_pres_com$location_code_numeric)),
                              X_iucn=X_iucn,
                              Z_iucn=ncol(X_iucn),
                              N_yr_iucn=length(unique(iucn_pres_com$event_year)),
                              year_id_iucn=as.numeric(factor(iucn_pres_com$event_year)),
                              source_iucn=as.numeric(factor(iucn_pres_com$source)),
                              N_source_iucn=length(unique(iucn_pres_com$source)),
                              wasting_index = wasting_index$wasting,
                              yr_index_reef=yr_index_reef$tt,
                              yr_index_ow = yr_index_ow$tt,
                              yr_index_iucn = yr_index_iucn$tt,
                              yr_index_aba = yr_index_aba$tt,
                              yr_index_multi = yr_index_multi$tt,
                              X_aba = X_aba,
                              X_multi = X_multi,
                              Z_aba = ncol(X_aba),
                              Z_multi = ncol(X_multi),
                              year_id_aba = as.numeric(factor(aba$year)),
                              year_id_multi = as.numeric(factor(multi$year)),
                              N_yr_aba = length(unique(aba$year)),
                              N_yr_multi = length(unique(multi$year)),
                              N_aba = nrow(aba),
                              N_multi = nrow(multi),
                              y_aba = aba$num_pycno,
                              y_multi = multi$num_pycno,
                              site_aba=as.numeric(factor(aba$stat_area)),
                              N_site_aba=length(unique(aba$stat_area)),
                              site_multi=as.numeric(factor(multi$site)),
                              N_site_multi=length(unique(multi$site)),
                              region_reef = as.numeric(factor(reef$region_numeric)),
                              region_ow = as.numeric(factor(ow$region_numeric)),
                              region_iucn = as.numeric(factor(iucn_pres_com$region_numeric)),
                              region_aba = as.numeric(factor(aba$region_numeric)),
                              region_multi = as.numeric(factor(multi$region_numeric)),
                              N_region_reef = length(unique(reef$region_numeric)),
                              N_region_ow = length(unique(ow$region_numeric)),
                              N_region_iucn = length(unique(iucn_pres_com$region_numeric)),
                              N_region_aba = length(unique(aba$region_numeric)),
                              N_region_multi = length(unique(multi$region_numeric)))

#run model with region effects to test-----------------
abun_mod_code_region_time <- cmdstan_model("code/pycno_model_with_time_varying_regions.stan")
abun_mod_region_time <- abun_mod_code_region_time$sample(data = data_pycno_abun_region,
                                                         #set seed for reproducibility - originally just
                                                         #used a random seed from cmdstanr and then 
                                                         #extracted with the $metadata() function so the
                                                         #exact results can be extracted again
                                                         #seed = 7104208, 
                                                         refresh = 50,#how often to print model progress
                                                         chains = 4,
                                                         parallel_chains = 4, 
                                                         iter_warmup = 1000,
                                                         iter_sampling = 2000,
                                                         adapt_delta = 0.95
)
#need to use save_object instead of saveRDS for cmdstanr objects
#and we also need to actually run any of the things we care about if we want 
#them to be saved which is wild but whatever
md_region_time <- abun_mod_region_time$metadata()
#check for convergence issues
di_region_time <- abun_mod_region_time$cmdstan_diagnose()
#abun_mod_region_time$save_object("model_output/pycno_model_cosewic_region_time.rds")

#and we'll convert and save as an rstan object too, just to be safe
abun_mod_region_time_rstan<- rstan::read_stan_csv(abun_mod_region_time$output_files())
#saveRDS(abun_mod_region_time_rstan, "model_output/pycno_model_cosewic_region_time_rstan.rds")

#compare looics for models with and without the region effects--------
#load in original model
abun_mod_rstan <- readRDS("model_output/pycno_model_rstan.rds")

#extract log liks from region model
log_lik_time <- extract_log_lik(abun_mod_region_time_rstan, merge_chains = FALSE)
r_eff_time <- relative_eff(exp(log_lik_time), cores = 4)
loo_time <- loo(log_lik_time, r_eff = r_eff_time, cores = 4)
print(loo_time)

#extract log liks from original model
log_lik_base <- extract_log_lik(abun_mod_rstan, merge_chains = FALSE)
r_eff_base <- relative_eff(exp(log_lik_base), cores = 4)
loo_base <- loo(log_lik_base, r_eff = r_eff_base, cores = 4)
print(loo_base)

#compare
comp <- loo_compare(loo_time, loo_base)
print(comp, simplify = FALSE)

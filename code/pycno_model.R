library(tidyverse) #for data wrangling
library(cmdstanr) #for running model
library(rstan) #for extracting and processing model output
library(DHARMa) #for examining residuals
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

#define all the data for the Stan model
data_pycno_abun = list(TT=22,
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
                       N_site_multi=length(unique(multi$site)))

#run model------------------------------------------------------------
#note that we're using the cmdstanr package for this model because it uses the
#latest version of Stan and there is a bug in the older version that rstan uses.
#the cutpoints on the ordered logistic regressions are super highly constrained, 
#and sometimes at the beginning of warmup Stan tries out values of the cutpoints
#that don't match those constraints - this is usually no problem but for some
#reason the older version of Stan causes a chain to abort if it tests out these
#inappropriate values (but only when there is an rng function in the generated
#quantities section for some reason). If you try running the model in cmdstanr
#you'll still get the same warnings, but Stan quickly figures out the 
#appropriate values to test out and runs fine - note that these warnings are 
#only cause for concern if they happen after the warmup!
abun_mod_code <- cmdstan_model("code/pycno_model.stan")
abun_mod <- abun_mod_code$sample(data = data_pycno_abun,
                                 #set seed for reproducibility - originally just
                                 #used a random seed from cmdstanr and then 
                                 #extracted with the $metadata() function so the
                                 #exact results can be extracted again
                                 seed = 9104208, 
                                 refresh = 100,#how often to print model progress
                                 chains = 4,
                                 parallel_chains = 4, 
                                 iter_warmup = 1000,
                                 iter_sampling = 3000
)
#need to use save_object instead of saveRDS for cmdstanr objects
#and we also need to actually run any of the things we care about if we want 
#them to be saved which is wild but whatever
md <- abun_mod$metadata()
#check for convergence issues
di <- abun_mod$cmdstan_diagnose()
#abun_mod$save_object("model_output/pycno_model.rds")

#and we'll convert and save as an rstan object too, just to be safe
abun_mod_rstan <- rstan::read_stan_csv(abun_mod$output_files())
#saveRDS(abun_mod_rstan, "model_output/pycno_model_rstan.rds")

#run model checks--------------------------------------------------------------
#read in model object if not already run
abun_mod <- readRDS("model_output/pycno_model.rds")
abun_mod_rstan <- readRDS("model_output/pycno_model.rds")

#check model metadata
metadata <- abun_mod$metadata()
#check for convergence issues
di <- abun_mod$cmdstan_diagnose()
#check r-hat values (should all be 1)
abundance_summary <- summary(abun_mod_rstan)$summary
#look good from a quick check but we can double check to make sure all are 
#exactly 1 when rounded to two decimal places
abundance_summary %>% 
  as_tibble() %>% 
  mutate(rhat = round(Rhat,2)) %>% 
  filter(rhat != 1)
#and there isn't a single one that's not 1 - yay! 

#posterior predictive check to compare model predictions with observed data - 
#need to use custom function to work with rstan object - this function works
#the same way as the ppc and pp_check functions available in other packages
#see custom_functions.R for details
#combine reef and ow for plotting with y_predict
all_pycno <- reef %>% 
  transmute(abundance_recode=abundance_recode, source="REEF") %>% 
  bind_rows(transmute(ow, abundance_recode=abundance_recode, source = "OW")) %>% 
  bind_rows(transmute(iucn_pres_com, abundance_recode = totalpycnos, 
                      source = "IUCN")) %>% 
  bind_rows(transmute(aba, abundance_recode = num_pycno, source="ABA")) %>% 
  bind_rows(transmute(multi, abundance_recode = num_pycno, source="MULTI"))

#to get it to work on all five datasets without needing to copy and paste the 
#code a bunch, we'll make some lists to differentiate the datasets and then 
#run a for loop to store the plots
n_source <- 5
sources <- c("REEF","OW","IUCN","ABA","MULTI")
y_predicts <- c("y_predict_reef", "y_predict_ow", "y_predict_iucn",
                "y_predict_aba", "y_predict_multi")
response <- list(reef$abundance_recode, ow$abundance_recode, 
                 iucn_pres_com$totalpycnos, aba$num_pycno, multi$num_pycno)
ppc_plots <- list()
iter <- 12000

for(i in 1:n_source){
  ppc_plots[[i]] <- ppc_stanfit(model = abun_mod_rstan,
                                iter = iter,
                                response = response[[i]],
                                plot_type = "density",
                                y_predict = y_predicts[[i]])
}
#the posterior predictive checks actually look really good!

#check residuals
for(i in 1:n_source){
  dharma_resids(model = abun_mod_rstan,
                iter = iter,
                response = response[[i]],
                y_predict = y_predicts[[i]])
}
#and the residuals aren't perfect for all five datasets but given the 
#complexity of the model and the massive variation across datasets, they look
#not too shabby
#note also that all of the DHARMa tests are based on significance thresholds, so 
#when the sample size is large, there are nearly always going to be significant 
#p-values! 

#examine pre/post crash changes-------------------------------------------------
abun_mod_rstan <- readRDS("model_output/pycno_model_rstan.rds")
params1 <- rstan::extract(abun_mod_rstan)
iter <- 12000

#use ord_to_n to compare popn estimates from different years to look at percent
#declines
pre_count <- data.frame(matrix(NA, nrow = iter, ncol =14))
for(i in 1:14){
  pre_count[,i] <- ord_to_n(params1$x[,i],params1$c_reef)
}
mean_pre <- rowMeans(pre_count)

#we won't include 2014 in either pre or post, since we don't
#know the exact date that SSWD hit
post_count <- data.frame(matrix(NA, nrow = iter, ncol = 7))
for(i in 16:22){
  post_count[,i-15] <- ord_to_n(params1$x[,i],params1$c_reef)
}
mean_post <- rowMeans(post_count)

diff <- as_tibble(cbind(mean_pre, mean_post)) %>% 
  mutate(diff = mean_post/mean_pre,
         decline = (diff - 1) * 100)
median(diff$diff - 1)
quantile(diff$diff,0.025)-1
quantile(diff$diff,0.975)-1

#examine effort----------------------------------------------------------------
#each dataset has a different metric of "survey effort" (either area or bottom
#time), and since each dataset is scaled independently, what is a "mean effort" 
#in one dataset may be very different than another. The scalar term in the model
#should account for those mean differences (e.g., if all the surveys in one 
#dataset are very small relative to another, the scalar term 'a' for that 
#dataset will likely be large). Here, we'll determine what these mean efforts 
#are.

reef_effort <- scale(log(reef$btime)) %>% 
  attributes()
ow_effort <- scale(log(ow$bottom_time)) %>% 
  attributes()
iucn_effort <- scale(log(iucn_pres_com$area)) %>% 
  attributes()
aba_effort <- scale(log(aba$num_quadrats)) %>%
  attributes()
multi_effort <- scale(log(multi$num_quadrats)) %>% 
  attributes()

effort_info <- tibble(df = c("reef","ow","iucn","aba","multi"),
                      metric = c("time","time","area","area","area"),
                      mean = NULL,
                      sd = NULL)  
efforts <- list(reef_effort,ow_effort,iucn_effort,aba_effort,multi_effort)

for(i in 1:5){
  effort_info[i,3] <- efforts[[i]][2]
  effort_info[i,4] <- efforts[[i]][3]
}

#convert back to real units
effort_info <- effort_info %>%
  mutate(unlogged = exp(`scaled:center`))


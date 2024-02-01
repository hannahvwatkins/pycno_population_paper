library(tidyverse)
library(rstan)
library(patchwork) #for plotting multiple panels
source("code/custom_functions.R")

#LOAD IN DATA-------------------------------------------------------------------
abun_mod_rstan <- readRDS("model_output/pycno_model.rds")

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

#FIGURE 1: TIME SERIES----------------------------------------------------------
#the way we convert model estimates (i.e., in link space) to real space will be
#different for each dataset (e.g., if we want to plot them on the same scale, we
#need to use the same transformation but with the dataset specific scalars, or 
#if we want to plot them on their own scales, we need to use different 
#transformations for the ordered logistic regressions vs. the negative binomial
#regressions), so we're going to extract the relevant info for each and then 
#make a for loop to help reduce repetition in the code
df_list <- list(reef, ow, iucn_pres_com, aba, multi)
n_dfs <- length(df_list)
iter <- 12000
TT <- 22
params1 <- rstan::extract(abun_mod_rstan)
#need to specify the name of the year by year estimates associated with each df
a_yr <- list(params1$a_yr_reef, params1$a_yr_ow, params1$a_yr_iucn,
             params1$a_yr_aba, params1$a_yr_multi)
a <- list(rep(0, times = iter), params1$a_ow, params1$a_iucn, 
          params1$a_aba, params1$a_multi)

#create arrays to populate
est_x <- array(NA, c(iter, TT))
x_mat <- array(NA, c(TT, 3))

#extract the estimates for the underlying state for each year in the full 
#timespan
for(i in 1:TT){
  x_coef <- data.frame(p_0.x=NA,
                       p_1.x=NA,p_2.x=NA,p_11.x=NA,
                       p_101.x=NA,
                       lambda.x=NA,
                       iter = seq(1,iter))
  x_coef[,1]<- plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,2]<- plogis(params1$c_reef[,2]-params1$x[,i])-plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,3]<- plogis(params1$c_reef[,3]-params1$x[,i])-plogis(params1$c_reef[,2]-params1$x[,i])
  x_coef[,4]<- plogis(params1$c_reef[,4]-params1$x[,i])-plogis(params1$c_reef[,3]-params1$x[,i])
  x_coef[,5]<- 1-plogis(params1$c_reef[,4]-params1$x[,i])
  x_coef[,6]<- apply(x_coef[,1:5],1,abund_tranfs)
  
  est_x[,i] <- x_coef[,6]
  
  x_mat[i,1]=median(est_x[,i])
  x_mat[i,2]=quantile(est_x[,i],0.025)
  x_mat[i,3]=quantile(est_x[,i],0.975)
}

#now make a dataframe with the actual years that we can connect the temporally-
#patchy observation estimates to
dat <- as.data.frame(x_mat) %>% 
  transmute(state_prob = V1,
            l.95.x = V2,
            u.95.x = V3,
            year = seq(from = 2000, to = 2021, by = 1))

#for each dataset, we're going to make some empty arrays to populate for each
#of the years with actual survey data
for(j in 1:n_dfs){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  #and then for each of those years, we'll backtransform based on the REEF scale
  #(i.e., after applying the relevant "a" term)
  for(i in 1:n_yr){
    a_coef <- data.frame(p_0=NA,
                         p_1=NA,p_2=NA,p_11=NA,
                         p_101=NA,
                         lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- plogis(params1$c_reef[,1]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,2]<- plogis(params1$c_reef[,2]-(a_yr[[j]][,i]-a[[j]]))-plogis(params1$c_reef[,1]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,3]<- plogis(params1$c_reef[,3]-(a_yr[[j]][,i]-a[[j]]))-plogis(params1$c_reef[,2]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,4]<- plogis(params1$c_reef[,4]-(a_yr[[j]][,i]-a[[j]]))-plogis(params1$c_reef[,3]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,5]<- 1-plogis(params1$c_reef[,4]-(a_yr[[j]][,i]-a[[j]]))
    a_coef[,6]<- apply(a_coef[,1:5],1,abund_tranfs)
    
    est_a[,i] <- a_coef[,6]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  #then connect the correct sampling year with the observation estimate
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  #then add it to the main plotting dataframe
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob = `1`,
              l.95.o = `2`,
              u.95.o = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#and we can extract the info on the number of surveys conducted each year across
#all datasets to add to the plot
num_surveys <- reef %>% 
  bind_rows(ow) %>% 
  bind_rows(transmute(iucn_pres_com, year = event_year)) %>% 
  bind_rows(aba) %>% 
  bind_rows(multi) %>% 
  group_by(year) %>% 
  summarize(n = n())

#now we can plot all the datasets together on top of the estimated underlying
#state
ggplot(data = dat) +
  geom_ribbon(aes(x = year, y = state_prob, ymin = l.95.x, ymax = u.95.x),
              fill='darkcyan', alpha = 0.2) +
  geom_line(aes(year, state_prob, colour = "Estimated\nstate"), lty=5,lwd=0.8) +
  geom_line(aes(year, obs_prob.x, colour = "REEF"), lwd=0.8) +
  geom_point(aes(year, obs_prob.x, fill = "REEF"), col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob.y, colour = "OW"),lwd=0.8) +
  geom_point(aes(year, obs_prob.y, fill = "OW"), col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob.x.x, colour = "Small datasets"),lwd=0.8) +
  geom_point(aes(year, obs_prob.x.x, fill = "Small datasets"), col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob.y.y, colour = "DFO-Abalone"), lwd=0.8) +
  geom_point(aes(year, obs_prob.y.y, fill = "DFO-Abalone"),col='white',pch=21,cex=3) +
  geom_line(aes(year, obs_prob, colour = "DFO-Multispecies"),lwd=0.8) +
  geom_point(aes(year, obs_prob, fill = "DFO-Multispecies"),col='white',pch=21,cex=3) +
  geom_text(data = num_surveys, aes(y=rep(-0.2,22), x=year, label = n), 
            size = 2.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  scale_colour_manual(name = "", breaks = c("Estimated\nstate",
                                            "REEF",
                                            "OW","DFO-Abalone",
                                            "DFO-Multispecies",
                                            "Small datasets"),
                      values = c("Estimated\nstate" = 'darkcyan', 
                                 "REEF" = "#197DC2", 
                                 "OW" ='#FFB14E',
                                 "DFO-Abalone" = "#FA8775",
                                 "DFO-Multispecies" = "#EA5F94",
                                 "Small datasets" = "#9D02D7"))+
  scale_fill_manual(name = "", breaks = c("Estimated\nstate",
                                          "REEF",
                                          "OW","DFO-Abalone",
                                          "DFO-Multispecies",
                                          "Small datasets"),
                    values = c("Estimated\nstate" = 'darkcyan', 
                               "REEF" = "#197DC2", 
                               "OW" ='#FFB14E',
                               "DFO-Abalone" = "#FA8775",
                               "DFO-Multispecies" = "#EA5F94",
                               "Small datasets" = "#9D02D7"),
                    guide = "none")+
  labs(x = "Year", y = "Abundance index") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12))
#ggsave("figures/Figure1.png", width = 8, height = 6)  

#FIGURE 2: POSTERIOR DENSITY OF REDUCTION---------------------------------------
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

ggplot(diff, aes(decline)) +
  geom_rect(data=NULL,aes(xmin=-100,xmax=-70,ymin=0,ymax=Inf),
            fill="#c96635") +
  geom_rect(data=NULL,aes(xmin=-70,xmax=-50,ymin=0,ymax=Inf),
            fill="#ca9c12") +
  geom_rect(data=NULL,aes(xmin=-50,xmax=0,ymin=0,ymax=Inf),
            fill="#386868") +
  geom_density(col = "black", fill = "grey") +
  theme_classic() +
  labs(x = "Change in abundance after SSWD epidemic (%)",
       y = "Posterior density") +
  scale_x_continuous(limits = c(-100, 0), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
#ggsave("figures/Figure2.png", height = 5, width = 7)

#FIGURE A1: ORDINAL CHECK------------------------------------------------------
abundance_summary <- summary(abun_mod_rstan)$summary
#create matrix of cutpoints from model
c_abun <- matrix(data = c(-1.55,-1.11,0.67,3.71), nrow = 1)
#create a sequence of x values (i.e., model estimates of the latent variable)
#across the observed sequence of x values
x_seq <- seq(from = -3, to = -0.4, by = 0.1)
x_seq2 <- seq(from=-5, to = 7, by = 0.1)
#backtransform the sequence
x_back <- ord_to_n(x_seq,c_abun) %>% 
  bind_cols(x_seq) %>% 
  transmute(back = `...1`,
            x = `...2`,
            #create the log-link model equivalent
            x_exp = exp(x+2.2)) 

x_back2 <- ord_to_n(x_seq2,c_abun) %>% 
  bind_cols(x_seq2) %>% 
  transmute(back = `...1`,
            x = `...2`,
            #create the log-link model equivalent
            x_exp = exp(x+2.2)) 

ggplot(x_back, aes(x, back)) +
  geom_line(col = "dodgerblue", linewidth = 1) +
  geom_line(aes(y = x_exp), col = "black", linewidth = 1) + 
  labs(x = "Untransformed true state from process model",
       y = "Estimated abundance index") +
  theme_classic()
#ggsave("figures/FigureA1.png", height = 5, width = 7)

#FIGURE A2: individual trends with error---------------------------------------
#need to specify the name of the year by year estimates associated with each df
a_yr <- list(params1$a_yr_reef, params1$a_yr_ow, params1$a_yr_iucn,
             params1$a_yr_aba, params1$a_yr_multi)
theta <- list(NA,NA,NA,NA, params1$theta)
a <- list(rep(0, times = iter), params1$a_ow, params1$a_iucn, 
          params1$a_aba, params1$a_multi)

est_x <- array(NA, c(iter, TT))
x_mat <- array(NA, c(TT, 3))

num_surveys_plotting <- list()
for(i in 1:n_dfs){
  num_surveys_plotting[[i]] <- df_list[[i]] %>% 
    group_by(year) %>% 
    summarize(n = n())
}

#for each year in the full timespan
for(i in 1:TT){
  x_coef <- data.frame(p_0.x=NA,
                       p_1.x=NA,p_2.x=NA,p_11.x=NA,
                       p_101.x=NA,
                       lambda.x=NA,
                       iter = seq(1,iter))
  x_coef[,1]<- plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,2]<- plogis(params1$c_reef[,2]-params1$x[,i])-plogis(params1$c_reef[,1]-params1$x[,i])
  x_coef[,3]<- plogis(params1$c_reef[,3]-params1$x[,i])-plogis(params1$c_reef[,2]-params1$x[,i])
  x_coef[,4]<- plogis(params1$c_reef[,4]-params1$x[,i])-plogis(params1$c_reef[,3]-params1$x[,i])
  x_coef[,5]<- 1-plogis(params1$c_reef[,4]-params1$x[,i])
  x_coef[,6]<- apply(x_coef[,1:5],1,abund_tranfs)
  
  est_x[,i] <- x_coef[,6]
  
  x_mat[i,1]=median(est_x[,i])
  x_mat[i,2]=quantile(est_x[,i],0.025)
  x_mat[i,3]=quantile(est_x[,i],0.975)
}

dat <- as.data.frame(x_mat) %>% 
  transmute(state_prob = V1,
            l.95.x = V2,
            u.95.x = V3,
            year = seq(from = 2000, to = 2021, by = 1))

#these don't need to be for loops but I got lazy and it was easier to cut and
#paste this way
#REEF only
for(j in 1:1){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(p_0=NA,
                         p_1=NA,p_2=NA,p_11=NA,
                         p_101=NA,
                         lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- plogis(params1$c_reef[,1]-a_yr[[j]][,i])
    a_coef[,2]<- plogis(params1$c_reef[,2]-a_yr[[j]][,i])-plogis(params1$c_reef[,1]-a_yr[[j]][,i])
    a_coef[,3]<- plogis(params1$c_reef[,3]-a_yr[[j]][,i])-plogis(params1$c_reef[,2]-a_yr[[j]][,i])
    a_coef[,4]<- plogis(params1$c_reef[,4]-a_yr[[j]][,i])-plogis(params1$c_reef[,3]-a_yr[[j]][,i])
    a_coef[,5]<- 1-plogis(params1$c_reef[,4]-a_yr[[j]][,i])
    a_coef[,6]<- apply(a_coef[,1:5],1,abund_tranfs)
    
    est_a[,i] <- a_coef[,6]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_r = `1`,
              l.95.r = `2`,
              u.95.r = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#OW only
for(j in 2:2){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(p_0=NA,
                         p_1=NA,p_11=NA,p_26=NA,
                         p_51=NA,p_101=NA,p_1001=NA,
                         lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- plogis(params1$c_ow[,1]-a_yr[[j]][,i])
    a_coef[,2]<- plogis(params1$c_ow[,2]-a_yr[[j]][,i])-plogis(params1$c_ow[,1]-a_yr[[j]][,i])
    a_coef[,3]<- plogis(params1$c_ow[,3]-a_yr[[j]][,i])-plogis(params1$c_ow[,2]-a_yr[[j]][,i])
    a_coef[,4]<- plogis(params1$c_ow[,4]-a_yr[[j]][,i])-plogis(params1$c_ow[,3]-a_yr[[j]][,i])
    a_coef[,5]<- plogis(params1$c_ow[,5]-a_yr[[j]][,i])-plogis(params1$c_ow[,4]-a_yr[[j]][,i])
    a_coef[,6]<- plogis(params1$c_ow[,6]-a_yr[[j]][,i])-plogis(params1$c_ow[,5]-a_yr[[j]][,i])
    a_coef[,7]<- 1-plogis(params1$c_ow[,6]-a_yr[[j]][,i])
    a_coef[,8]<- apply(a_coef[,1:7],1,abund_tranfs_ow)
    
    est_a[,i] <- a_coef[,8]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_o = `1`,
              l.95.o = `2`,
              u.95.o = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#IUCN only
for(j in 3:3){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- exp(a_yr[[j]][,i])
    
    est_a[,i] <- a_coef[,1]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_i = `1`,
              l.95.i = `2`,
              u.95.i = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}


#abalone only
for(j in 4:4){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- exp(a_yr[[j]][,i])
    
    est_a[,i] <- a_coef[,1]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_a = `1`,
              l.95.a = `2`,
              u.95.a = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}

#multi only
for(j in 5:5){
  n_yr <- length(unique(df_list[[j]]$year))
  est_a <- array(NA, c(iter, n_yr))
  a_mat <- array(NA, c(n_yr, 3)) 
  
  for(i in 1:n_yr){
    a_coef <- data.frame(lambda.o=NA,
                         iter = seq(1,iter))
    a_coef[,1]<- (1-plogis(theta[[j]])) * exp(a_yr[[j]][,i])
    
    est_a[,i] <- a_coef[,1]
    a_mat[i,1]=median(est_a[,i])
    a_mat[i,2]=quantile(est_a[,i],0.025)
    a_mat[i,3]=quantile(est_a[,i],0.975)
  }
  
  obs_years <- df_list[[j]] %>% 
    distinct(year) %>% 
    arrange(year)
  
  a_mat_yr <- cbind(a_mat, obs_years) %>% 
    transmute(obs_prob_m = `1`,
              l.95.m = `2`,
              u.95.m = `3`,
              year = year)
  
  dat <- left_join(dat, a_mat_yr, by = "year")
}


#REEF ONLY
reef_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_r, colour = "REEF"), lwd=0.8) +
  geom_point(aes(year, obs_prob_r, fill = "REEF"), col='white',pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_r, ymin = l.95.r, ymax = u.95.r,
                  fill = 'REEF'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[1]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[1]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("REEF" = "#197DC2"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("REEF" = "#197DC2"),
                    guide = "none")+
  labs(x = "", y = "",
       title = "REEF dive surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

#OW ONLY
ow_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_o, colour = "OW"), lwd=0.8) +
  geom_point(aes(year, obs_prob_o, fill = "OW"), col='white',pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_o, ymin = l.95.o, ymax = u.95.o,
                  fill = 'OW'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[2]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[2]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("OW" ='#FFB14E'), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("OW" ='#FFB14E'),
                    guide = "none")+
  labs(x = "", y = "", title = "OW Pacific Marine\nLife surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

#IUCN ONLY
iucn_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_i, colour = "Small datasets"), lwd=0.8) +
  geom_point(aes(year, obs_prob_i, fill = "Small datasets"), col='white',
             pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_i, ymin = l.95.i, ymax = u.95.i,
                  fill = 'Small datasets'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[3]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[3]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("Small datasets" = "#9D02D7"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("Small datasets" = "#9D02D7"),
                    guide = "none")+
  labs(x = "Year", y = "", title="Combined small datasets") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)


#abalone ONLY
#the geom_ribbon only works if there are consecuctive years in the model - since
#there were no surveys in 2020, the 2021 point is on it's own and we need to 
#add in a CI manually
aba_subset <- dat %>% 
  select(year, obs_prob_a:u.95.a) %>% 
  filter(year == 2021)

aba_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_a, colour = "DFO-Abalone"), lwd=0.8) +
  geom_point(aes(year, obs_prob_a, fill = "DFO-Abalone"), 
             col='white',pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_a, ymin = l.95.a, ymax = u.95.a,
                  fill = 'DFO-Abalone'),
              alpha = 0.2) +
  geom_linerange(data=aba_subset, aes(x=year,  
                                      ymin = l.95.a, ymax = u.95.a,  
                                      colour = 'DFO-Abalone'), 
                 alpha = 0.2, size = 3) +
  geom_text(data = num_surveys_plotting[[4]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[4]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("DFO-Abalone" = "#FA8775"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("DFO-Abalone" = "#FA8775"),
                    guide = "none")+
  labs(x = "Year", y = "",
       title="DFO-Abalone dive surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

#multi ONLY
multi_plot <- ggplot(data = dat) +
  geom_line(aes(year, obs_prob_m, colour = "DFO-Multispecies"), lwd=0.8) +
  geom_point(aes(year, obs_prob_m, fill = "DFO-Multispecies"), col='white',
             pch=21,cex=2) +
  geom_ribbon(aes(x = year, y = obs_prob_m, ymin = l.95.m, ymax = u.95.m,
                  fill = 'DFO-Multispecies'),
              alpha = 0.2) +
  geom_text(data = num_surveys_plotting[[5]], 
            aes(y=rep(0,nrow(num_surveys_plotting[[5]])), x=year, 
                label = n), size = 1.5) +
  theme_classic() +
  theme(plot.title = ggtext::element_markdown(face = "bold", hjust = 0.5, size=12),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  scale_colour_manual(name = "", values = c("DFO-Multispecies" = "#EA5F94"), 
                      guide = "none")+
  scale_fill_manual(name = "", values = c("DFO-Multispecies" = "#EA5F94"),
                    guide = "none")+
  labs(x = "", y = "", 
       title = "DFO-Multispecies\ndive surveys") +
  theme(legend.justification = c(1, 1), 
        legend.position = c(1, 1),
        legend.text = element_text(size = 12)) +
  ylim(0,NA)

common_y_axis <- ggplot(data.frame(l = "Annual abundance index (mean count per survey)", 
                                   x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90) + 
  theme_void() +
  #theme() +
  coord_cartesian(clip = "off")

common_y_axis+(reef_plot+ow_plot)/
  (multi_plot+aba_plot)/
  (iucn_plot + plot_spacer()) + 
  #now make the label v small relative to the plots
  plot_layout(widths = c(1,25))
#ggsave("figures/FigureA2.png", width = 7, height = 8)  

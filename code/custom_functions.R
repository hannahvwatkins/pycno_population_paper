#this is a script with a few custom functions Dan Greenberg wrote for dealing
#with the ordinal model output plus a couple of my own (with help from Helen
#Yan!)
library(tidyverse)
library(rstan)
library(DHARMa)

###Functions####
#this function converts probabilities to abundances for the reef data, for use 
#in ord_to_n below
#by changing the values used to multiply the probabilities, we can test the 
#robustness of our choice
abund_tranfs<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*2+probs[4]*11+probs[5]*101
  return(sum)
}

#we also want to know if randomly generating a number from the possible range
#yields similar results as above
abund_tranfs_sim <- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*round(runif(1,2,10))+
    probs[4]*round(runif(1,11,100))+probs[5]*round(runif(1,101,200))
  return(sum)
}

#same but for the ow categories
abund_tranfs_ow<- function(probs){
  sum<- probs[1]*0+probs[2]*1+probs[3]*11+probs[4]*26+probs[5]*51+probs[6]*101+
    probs[7]*1001
  return(sum)
}

#estimate abundances using the distribution function (which basically tells you
#how likely each score is for a given obs), and the known relationships between 
#scores and abundances from abund_tranfs
#note that it's written out slightly differently than on the Stan help file, but
#it's ultimately the same equation since 1-plogis(c-x) is the same as 
#plogis(x-c) 
ord_to_n<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
    p[,1]=plogis(c[,1]-x)
    p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
    p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
    p[,4]=plogis(c[,4]-x)-plogis(c[,3]-x)
    p[,5]=1-plogis(c[,4]-x)
    for(i in 1:length(x)){
      abund_x[i]=abund_tranfs(p[i,])  
    }
  return(abund_x)
}

#same for the randomly generated option
ord_to_n_sim<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=5,nrow=length(x))
  p[,1]=plogis(c[,1]-x)
  p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
  p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
  p[,4]=plogis(c[,4]-x)-plogis(c[,3]-x)
  p[,5]=1-plogis(c[,4]-x)
  for(i in 1:length(x)){
    abund_x[i]=abund_tranfs_sim(p[i,])  
  }
  return(abund_x)
}

#this is the equivalent function for the OW abundance categories, but since
#the latent variable (i.e. a_yr2) underlying the categories can be scaled to 
#the REEF dataset with the scalar term a, we can convert it with the REEF 
#ord_to_n above (as long as we include the a)
ord_to_n_ow<- function(x,c){
  abund_x<- numeric(length(x))
  p= matrix(ncol=7,nrow=length(x))
  p[,1]=plogis(c[,1]-x)
  p[,2]=plogis(c[,2]-x)-plogis(c[,1]-x)
  p[,3]=plogis(c[,3]-x)-plogis(c[,2]-x)
  p[,4]=plogis(c[,4]-x)-plogis(c[,3]-x)
  p[,5]=plogis(c[,5]-x)-plogis(c[,4]-x)
  p[,6]=plogis(c[,6]-x)-plogis(c[,5]-x)
  p[,7]=1-plogis(c[,6]-x)
  for(i in 1:length(x)){
    abund_x[i]=abund_tranfs(p[i,])  
  }
  return(abund_x)
}
  
#this function does a basic posterior predictive check - created it manually due
#to some issues with cmdstanr outputs, but in hindsite, it is a deeply janky 
#approach. I am repenting for my inefficient coding sins.
#to use this function on any dataset, make sure you name the output from the rng
#in the generated quantities section "y_predict" (or else change that term in
#the function to match what you used)
#model is the name of the model output, iter refers to how many 
#iterations are stored in the output (should be #chains * half the iterations in
#the model, unless you customized the burn-in or thinning), and response is the 
#location of the response variable in the original data (i.e. df$response)
#plot_type is one of "histogram" or "density"
ppc_stanfit <- function(model, iter, response, 
                        plot_type = "density", y_predict = "y_predict"){
  predictions <- as_tibble(rstan::extract(model,
                                          pars = y_predict)[[1]]) %>% 
    rownames_to_column() %>% 
    pivot_longer(!rowname, 
                 names_to = "obs", 
                 values_to = "y_predict") %>% 
    mutate(rowname = as.numeric(rowname))
  
  #select 50 random iterations to extract
  if(plot_type == "histogram"){
   predict_df<- as_tibble(sample(1:iter, 10, replace = FALSE)) %>% 
     left_join(predictions, by = c("value" = "rowname")) %>% 
     mutate(value = as.character(value))
   ppc_df <- tibble(value = "y", y_predict = response) %>% 
     rownames_to_column() %>% 
     mutate(obs = paste("V", rowname, sep = "")) %>% 
     bind_rows(predict_df)
  } else{
    ppc_df<- as_tibble(sample(1:iter, 50, replace = FALSE)) %>% 
      left_join(predictions, by = c("value" = "rowname"))
  }
 
  if(plot_type == "histogram"){
    
    fig <- ggplot() +
      geom_histogram(data = ppc_df, aes(x = y_predict, group = value)) +
      facet_wrap(~value) +
      theme_classic()
  } else { 
    fig <- ggplot() +
      geom_density(data = ppc_df, aes(x = y_predict, group = value), 
                   alpha = 0.1) +
      geom_density(aes(x = response), 
                   col = "dodgerblue", size = 1.5) +
      theme_classic()
}
  print(fig)
}

#extract residuals from a Bayesian model for use with DHARMa
dharma_resids <- function(model, iter, response, y_predict = "y_predict"){
  predictions <-
    as_tibble(rstan::extract(model,
                             pars = y_predict)[[1]]) %>%
    rownames_to_column() %>%
    pivot_longer(!rowname,
                 names_to = "obs",
                 values_to = "y_predict") %>%
    mutate(rowname = as.numeric(rowname))
  
  # randomly select 250 draws
  resid_matrix <-
    as_tibble(sample(1:iter, 250, replace = FALSE)) %>%
    left_join(predictions, by = c("value" = "rowname")) %>%
    dplyr::select(-obs) %>%
    # so we can pivot wider without pissing it off
    group_by(value) %>%
    mutate(row = row_number()) %>%
    pivot_wider(everything(), names_from = 'value', values_from = 'y_predict') %>%
    dplyr::select(-row) %>%
    as.matrix()
  
  resids <-
    createDHARMa(simulatedResponse = resid_matrix,
                 observedResponse = response,
                 fittedPredictedResponse = apply(t(resid_matrix), 2, median),
                 #note this is the one line you might need to change for other
                 #models
                 integerResponse = FALSE)
  plot(resids)
}


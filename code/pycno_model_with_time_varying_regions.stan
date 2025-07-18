//this is a custom prior for the cutpoints in the ordered logistic regression
//derived from Michael Betancourt's blog:
//https://betanalpha.github.io/assets/case_studies/ordinal_regression.html
//The code in this case study is copyrighted by Michael Betancourt and licensed 
//under the new BSD (3-clause) license: 
//https://opensource.org/licenses/BSD-3-Clause
functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
  
}
//upload the data as defined in the list in pycno_model.R
data{
  int<lower=1> N_reef;//number of observations (REEF surveys)
  int<lower=1> N_ow;
  int<lower=1> N_iucn;
  int<lower=1> N_aba; 
  int<lower=1> N_multi;
  array[N_reef] int y_reef; //abundance category for each survey (1-5)
  array[N_ow] int y_ow; //abundance value for each survey (1-7)
  array[N_iucn] int y_iucn; //counts 
  array[N_aba] int y_aba; 
  array[N_multi] int y_multi; 
  int<lower=0> N_site_reef; //number of sites
  array[N_reef] int<lower=1,upper=N_site_reef> site_reef; // array of site identities
  int<lower=0> N_site_ow; 
  array[N_ow] int<lower=1,upper=N_site_ow> site_ow; 
  int<lower=0> N_site_iucn; 
  array[N_iucn] int<lower=1,upper=N_site_iucn> site_iucn; 
  int<lower=0> N_site_aba; 
  array[N_aba] int<lower=1,upper=N_site_aba> site_aba; 
  int<lower=0> N_site_multi; 
  array[N_multi] int<lower=1,upper=N_site_multi> site_multi; 
  int<lower=0> N_source_iucn; //number of data sources
  array[N_iucn] int<lower=1,upper=N_source_iucn> source_iucn; //array of source id
  int<lower=0> N_dv_ow; //number of divers
  array[N_ow] int<lower=1,upper=N_dv_ow> diver_ow; // array of diver identities
  int Z_reef; // columns in the covariate matrix
  int Z_ow; 
  int Z_iucn;
  int Z_aba; 
  int Z_multi;
  matrix[N_reef,Z_reef] X_reef; // design matrix X
  matrix[N_ow,Z_ow] X_ow; 
  matrix[N_iucn,Z_iucn] X_iucn;
  matrix[N_aba,Z_aba] X_aba; 
  matrix[N_multi,Z_multi] X_multi;
  int K_reef; //ordinal levels
  int K_ow;
  int TT; // total timespan across all datasets
  int<lower=0> N_yr_reef; //number of years total for REEF 
  int<lower=0> N_yr_ow; //number of years with surveys
  int<lower=0> N_yr_iucn;
  int<lower=0> N_yr_aba; 
  int<lower=0> N_yr_multi;
  array[N_reef] int<lower=1,upper=N_yr_reef> year_id_reef; //array of years
  array[N_ow] int<lower=1,upper=N_yr_ow> year_id_ow; 
  array[N_iucn] int<lower=1,upper=N_yr_iucn> year_id_iucn;
  array[N_aba] int<lower=1,upper=N_yr_aba> year_id_aba; 
  array[N_multi] int<lower=1,upper=N_yr_multi> year_id_multi;
  array[TT] int wasting_index;
  array[N_yr_reef] int yr_index_reef; //index of years
  array[N_yr_ow] int yr_index_ow;
  array[N_yr_iucn] int yr_index_iucn;
  array[N_yr_aba] int yr_index_aba;
  array[N_yr_multi] int yr_index_multi;  
  
  int<lower=0> N_region_reef;
  int<lower=0> N_region_ow;
  int<lower=0> N_region_iucn;
  int<lower=0> N_region_aba;
  int<lower=0> N_region_multi;
  array[N_reef] int<lower=1,upper=N_region_reef> region_reef; // array of region identities
  array[N_ow] int<lower=1,upper=N_region_ow> region_ow; 
  array[N_iucn] int<lower=1,upper=N_region_iucn> region_iucn; 
  array[N_aba] int<lower=1,upper=N_region_aba> region_aba; 
  array[N_multi] int<lower=1,upper=N_region_multi> region_multi;
}
parameters {
  ordered[K_reef-1] c_reef; //cutpoints
  ordered[K_ow-1] c_ow; //cutpoints
  real x0; //initial popn size
  real a_ow; //scalars for time-series
  real a_iucn;
  real a_aba;
  real a_multi; 
  real u_pre; //trend pre-wasting
  real u_wasting; //trend during the crash
  real u_post; //trend post-wasting
  
  //variance on the deviance components
  real<lower = 0> sd_site_reef;
  real<lower = 0> sd_site_ow;
  real<lower = 0> sd_site_iucn;
  real<lower = 0> sd_site_aba;
  real<lower = 0> sd_site_multi;
  real<lower = 0> sd_dv_ow;
  real<lower = 0> sd_source_iucn;
  real<lower = 0> sd_r_reef; //variance on observation deviations 
  real<lower = 0> sd_r_ow; 
  real<lower = 0> sd_r_iucn; 
  real<lower = 0> sd_r_aba; 
  real<lower = 0> sd_r_multi;
  real<lower = 0> sd_q; //variance on process deviance for underlying state
  
  vector[TT] pro_dev; //process deviations
  vector[N_yr_reef] obs_dev_reef; //observation deviations 
  vector[N_yr_ow] obs_dev_ow; 
  vector[N_yr_iucn] obs_dev_iucn;
  vector[N_yr_aba] obs_dev_aba; 
  vector[N_yr_multi] obs_dev_multi;

  real<lower=0> phi_iucn; //dispersion parameters for negative binomials
  real<lower=0> phi_aba;
  real<lower=0> phi_multi;
  real theta_multi; //theta term for zero inflation

  vector[Z_reef] beta_reef; //survey-specific effort coefficients
  vector[Z_ow] beta_ow; 
  vector[Z_iucn] beta_iucn;
  vector[Z_aba] beta_aba; 
  vector[Z_multi] beta_multi; 
  vector[N_site_reef] a_site_reef; //deviation between sites
  vector[N_site_ow] a_site_ow; 
  vector[N_site_iucn] a_site_iucn; 
  vector[N_site_aba] a_site_aba; 
  vector[N_site_multi] a_site_multi; 
  vector[N_source_iucn] a_source_iucn;//deviation between sources within dataset
  vector[N_dv_ow] a_dv_ow; //deviation between divers OW
  
  //region model terms
  matrix[TT,N_region_reef] pro_dev_region_reef; //process deviations
  matrix[TT,N_region_ow] pro_dev_region_ow; //process deviations
  matrix[TT,N_region_iucn] pro_dev_region_iucn; //process deviations
  matrix[TT,N_region_aba] pro_dev_region_aba; //process deviations
  matrix[TT,N_region_multi] pro_dev_region_multi; //process deviations
  matrix[N_yr_reef,N_region_reef] obs_dev_region_reef; //observation deviations 
  matrix[N_yr_ow,N_region_ow] obs_dev_region_ow; 
  matrix[N_yr_iucn,N_region_iucn] obs_dev_region_iucn;
  matrix[N_yr_aba,N_region_aba] obs_dev_region_aba; 
  matrix[N_yr_multi,N_region_multi] obs_dev_region_multi;
  real<lower = 0> sd_r_region_reef; //variance on observation deviations 
  real<lower = 0> sd_r_region_ow; 
  real<lower = 0> sd_r_region_iucn; 
  real<lower = 0> sd_r_region_aba; 
  real<lower = 0> sd_r_region_multi;
  real<lower = 0> sd_q_region_reef;
  real<lower = 0> sd_q_region_ow;
  real<lower = 0> sd_q_region_iucn;
  real<lower = 0> sd_q_region_aba;
  real<lower = 0> sd_q_region_multi;
  vector[N_region_reef] x0_reef; 
  vector[N_region_ow] x0_ow; 
  vector[N_region_iucn] x0_iucn; 
  vector[N_region_aba] x0_aba; 
  vector[N_region_multi] x0_multi; 
  
}

transformed parameters{
  vector[TT] x; //underlying state (i.e., estimated true population)
  vector[N_yr_reef] a_yr_reef; //year by year estimates
  vector[N_yr_ow] a_yr_ow; 
  vector[N_yr_iucn] a_yr_iucn;
  vector[N_yr_aba] a_yr_aba; 
  vector[N_yr_multi] a_yr_multi;
  matrix[TT, N_region_reef] x_region_reef;
  matrix[TT, N_region_ow] x_region_ow;
  matrix[TT, N_region_iucn] x_region_iucn;
  matrix[TT, N_region_aba] x_region_aba;
  matrix[TT, N_region_multi] x_region_multi;
  matrix[N_yr_reef,N_region_reef] a_region_reef;
  matrix[N_yr_ow,N_region_ow] a_region_ow;
  matrix[N_yr_iucn,N_region_iucn] a_region_iucn;
  matrix[N_yr_aba,N_region_aba] a_region_aba;
  matrix[N_yr_multi,N_region_multi] a_region_multi;
  
  for(r in 1:N_region_reef){
    x_region_reef[1,r] = x0_reef[r] + pro_dev_region_reef[1,r]*sd_q_region_reef;
    for(t in 2:TT){
      x_region_reef[t,r] = x_region_reef[t-1,r] + pro_dev_region_reef[t,r]*sd_q_region_reef;
    }
    for(t in 1:N_yr_reef){
      a_region_reef[t,r] = x_region_reef[yr_index_reef[t],r] + obs_dev_region_reef[t,r]*sd_r_region_reef; 
    }
  }
  for(r in 1:N_region_ow){
    x_region_ow[1,r] = x0_ow[r] + pro_dev_region_ow[1,r]*sd_q_region_ow;
    for(t in 2:TT){
      x_region_ow[t,r] = x_region_ow[t-1,r] + pro_dev_region_ow[t,r]*sd_q_region_ow;
    }
    for(t in 1:N_yr_ow){
      a_region_ow[t,r] = x_region_ow[yr_index_ow[t],r] + obs_dev_region_ow[t,r]*sd_r_region_ow; 
    }
  }
  for(r in 1:N_region_iucn){
    x_region_iucn[1,r] = x0_iucn[r] + pro_dev_region_iucn[1,r]*sd_q_region_iucn;
    for(t in 2:TT){
      x_region_iucn[t,r] = x_region_iucn[t-1,r] + pro_dev_region_iucn[t,r]*sd_q_region_iucn;
    }
    for(t in 1:N_yr_iucn){
      a_region_iucn[t,r] = x_region_iucn[yr_index_iucn[t],r] + obs_dev_region_iucn[t,r]*sd_r_region_iucn; 
    }
  }
  for(r in 1:N_region_aba){
    x_region_aba[1,r] = x0_aba[r] + pro_dev_region_aba[1,r]*sd_q_region_aba;
    for(t in 2:TT){
      x_region_aba[t,r] = x_region_aba[t-1,r] + pro_dev_region_aba[t,r]*sd_q_region_aba;
    }
    for(t in 1:N_yr_aba){
      a_region_aba[t,r] = x_region_aba[yr_index_aba[t],r] + obs_dev_region_aba[t,r]*sd_r_region_aba; 
    }
  }
  for(r in 1:N_region_multi){
    x_region_multi[1,r] = x0_multi[r] + pro_dev_region_multi[1,r]*sd_q_region_multi;
    for(t in 2:TT){
      x_region_multi[t,r] = x_region_multi[t-1,r] + pro_dev_region_multi[t,r]*sd_q_region_multi;
    }
    for(t in 1:N_yr_multi){
      a_region_multi[t,r] = x_region_multi[yr_index_multi[t],r] + obs_dev_region_multi[t,r]*sd_r_region_multi; 
    }
  }
  
  
  x[1] = x0 + u_pre + pro_dev[1]*sd_q;
  for(t in 2:TT){
    if(wasting_index[t] == 0)
      x[t] = x[t-1] + u_pre + pro_dev[t]*sd_q;
    else
      if(wasting_index[t] == 1)
        x[t] = x[t-1] + u_wasting + pro_dev[t]*sd_q;
    else
      if(wasting_index[t] == 2)
        x[t] = x[t-1] + u_post + pro_dev[t]*sd_q;
  }

  for(i in 1:N_yr_reef){
    a_yr_reef[i] = x[yr_index_reef[i]] + obs_dev_reef[i]*sd_r_reef; 
  }
  for(i in 1:N_yr_ow){
    a_yr_ow[i] = a_ow + x[yr_index_ow[i]] + obs_dev_ow[i]*sd_r_ow; 
  }
  for(i in 1:N_yr_iucn){
    a_yr_iucn[i] = a_iucn + x[yr_index_iucn[i]] + obs_dev_iucn[i]*sd_r_iucn; 
  }
  for(i in 1:N_yr_aba){
    a_yr_aba[i] = a_aba + x[yr_index_aba[i]] + obs_dev_aba[i]*sd_r_aba; 
  }
  for(i in 1:N_yr_multi){
    a_yr_multi[i] = a_multi + x[yr_index_multi[i]] + obs_dev_multi[i]*sd_r_multi; 
  }
}  

model{
  //priors on sds for process and observation error
  //gamma priors on variance terms provide reasonable constraints, and the model 
  //does not appear to be strongly influenced by the choice of prior here, as I 
  //have tried a range of options with the same result. Can check the shape of 
  //these priors in R with bayesAB::plotGamma
  //note though that if the sd_q and sd_r priors do not encompass small/large 
  //enough values, with the non-centered parameterization the model will likely
  //shrink/increase these values and compensate by changing the pro_dev/obs_dev
  //values accordingly. Not a big deal in terms of actually modelling the 
  //underlying state, but this will make it impossible to directly compare sd_q
  //and sd_r. To ensure this issue isn't happening, check to make sure each set
  //of pro_devs and obs_devs is actually distributed normally around 0 with an
  //sd of 1!
  sd_q ~ gamma(2,4);    
  sd_r_reef ~ gamma(2,4);  
  sd_r_ow ~ gamma(2,4);    
  sd_r_iucn ~ gamma(2,4);
  sd_r_aba ~ gamma(2,4);  
  sd_r_multi ~ gamma(2,4);
  x0 ~ normal(0,2); //initial state
  u_pre ~ normal(0,2);
  u_wasting ~ normal(0,2);
  u_post ~ normal(0,2);
  a_ow ~ normal(0,2); //scalars
  a_iucn ~ normal(0,2);
  a_aba ~ normal(0,2);
  a_multi ~ normal(0,2);
  //varying intercepts of observation and process errors
  pro_dev ~ std_normal(); 
  obs_dev_reef ~ std_normal();
  obs_dev_ow ~ std_normal();
  obs_dev_iucn ~ std_normal();
  obs_dev_aba ~ std_normal();
  obs_dev_multi ~ std_normal();
  
  //priors on fixed effects and cutpoints
  c_reef ~ induced_dirichlet(rep_vector(1, K_reef), 0); //custom prior function
  c_ow ~ induced_dirichlet(rep_vector(1, K_ow), 0); 
  beta_reef ~ normal(0,2); 
  beta_ow ~ normal(0,2); 
  beta_iucn ~normal(0,2);
  beta_aba ~ normal(0,2); 
  beta_multi ~ normal(0,2);
  //dispersion parameters for negative binomial models - inv_gamma does a good
  //job of constraining on both the lower bound and to realistic upper values
  //and changing the choice of shape and scale does not change the model fit
  phi_iucn ~ inv_gamma(5,5);
  phi_aba ~ inv_gamma(5,5);
  phi_multi ~ inv_gamma(5,5);
  //zero inflation component of multi dataset
  theta_multi ~ normal(0,2);
  
  //sds on random effects
  sd_site_reef ~ gamma(2,4); 
  sd_site_ow ~ gamma(2,4); 
  sd_site_iucn ~ gamma(2,4);
  sd_site_aba ~ gamma(2,4); 
  sd_site_multi ~ gamma(2,4); 
  sd_dv_ow ~ gamma(2,4);   
  sd_source_iucn ~ gamma(2,4);
  
  //varying intercepts of random effects
  a_site_reef ~ std_normal();
  a_site_ow ~ std_normal();
  a_site_iucn ~ std_normal();
  a_site_aba ~ std_normal();
  a_site_multi ~ std_normal();
  a_dv_ow ~ std_normal();
  a_source_iucn ~ std_normal();
  //note that there's no a_yr prior, since it's in the transformed parameter 
  //sections
  
  //region-speficic terms
  sd_r_region_reef ~ gamma(2,4); 
  sd_r_region_ow ~ gamma(2,4); 
  sd_r_region_iucn ~ gamma(2,4); 
  sd_r_region_aba ~ gamma(2,4); 
  sd_r_region_multi ~ gamma(2,4); 
  sd_q_region_reef ~ gamma(2,4); 
  sd_q_region_ow ~ gamma(2,4); 
  sd_q_region_iucn ~ gamma(2,4); 
  sd_q_region_aba ~ gamma(2,4); 
  sd_q_region_multi ~ gamma(2,4); 
  x0_reef ~ normal(0,2); 
  x0_ow ~ normal(0,2); 
  x0_iucn ~ normal(0,2); 
  x0_aba ~ normal(0,2); 
  x0_multi ~ normal(0,2); 
  //for some reason stan doesn't like applying the priors to the full matrix
  for(t in 1:TT){
    pro_dev_region_reef[t,] ~ std_normal();
    pro_dev_region_ow[t,] ~ std_normal();
    pro_dev_region_iucn[t,] ~ std_normal();
    pro_dev_region_aba[t,] ~ std_normal();
    pro_dev_region_multi[t,] ~ std_normal();
  }
  for(t in 1:N_yr_reef){
    obs_dev_region_reef[t,] ~ std_normal();
  }
  for(t in 1:N_yr_ow){
    obs_dev_region_ow[t,] ~ std_normal();
  }
  for(t in 1:N_yr_iucn){
    obs_dev_region_iucn[t,] ~ std_normal();
  }
  for(t in 1:N_yr_aba){
    obs_dev_region_aba[t,] ~ std_normal();
  }
  for(t in 1:N_yr_multi){
    obs_dev_region_multi[t,] ~ std_normal();
  }

  //main ordered logistic/negative binomial regression models for each dataset
  //since there is no intercept in either model and all continuous variables are
  //centered around 0 (as well as the random effects of site and diver), the 
  //a_yr term tells us the expected abundance for each year at the mean value
  //of all the other variables - it basically creates a separate intercept for 
  //each year!
  for(i in 1:N_reef){
    y_reef[i] ~ ordered_logistic(a_yr_reef[year_id_reef[i]]+a_region_reef[year_id_reef[i],region_reef[i]]+a_site_reef[site_reef[i]]*sd_site_reef+X_reef[i,]*beta_reef,c_reef);
  }
  for(i in 1:N_ow){
    y_ow[i] ~ ordered_logistic(a_yr_ow[year_id_ow[i]]+a_region_ow[year_id_ow[i],region_ow[i]]+a_site_ow[site_ow[i]]*sd_site_ow+a_dv_ow[diver_ow[i]]*sd_dv_ow+X_ow[i,]*beta_ow,c_ow);
  }
  for(i in 1:N_iucn){
    y_iucn[i] ~ neg_binomial_2_log(a_yr_iucn[year_id_iucn[i]]+a_region_iucn[year_id_iucn[i],region_iucn[i]]+a_site_iucn[site_iucn[i]]*sd_site_iucn+a_source_iucn[source_iucn[i]]*sd_source_iucn+X_iucn[i,]*beta_iucn, phi_iucn);
  }
  for(i in 1:N_aba){
    y_aba[i] ~ neg_binomial_2_log(a_yr_aba[year_id_aba[i]]+a_region_aba[year_id_aba[i],region_aba[i]]+a_site_aba[site_aba[i]]*sd_site_aba+X_aba[i,]*beta_aba, phi_aba);
  }

  
  //based on posterior predictive checks of a negative binomial model, the multi
  //dataset appears to be quite strongly zero inflated (i.e., all iterations 
  //underpredict the number of zeros), so we'll add a simple zero inflation 
  //component to this part of the model to help with that
  for (i in 1:N_multi){
    if(y_multi[i] == 0) {
      target += log_sum_exp(bernoulli_logit_lpmf(1 | theta_multi),
                            bernoulli_logit_lpmf(0 | theta_multi)
                              + neg_binomial_2_log_lpmf(y_multi[i] | a_yr_multi[year_id_multi[i]]+a_region_multi[year_id_multi[i],region_multi[i]]+a_site_multi[site_multi[i]]*sd_site_multi+X_multi[i,]*beta_multi, phi_multi));
    }
    else {
       target += bernoulli_logit_lpmf(0 | theta_multi)
                  + neg_binomial_2_log_lpmf(y_multi[i] | a_yr_multi[year_id_multi[i]]+a_region_multi[year_id_multi[i],region_multi[i]]+a_site_multi[site_multi[i]]*sd_site_multi+X_multi[i,]*beta_multi, phi_multi);
    }
  }
}
 generated quantities{
 //vector[N_reef] y_predict_reef; //generate model predictions to assess fit
 //vector[N_ow] y_predict_ow; 
 //vector[N_iucn] y_predict_iucn; 
 //vector[N_aba] y_predict_aba; 
 //vector[N_multi] y_predict_multi; 
 vector[N_reef+N_ow+N_iucn+N_aba+N_multi] log_lik; //generate log likelihoods 

 //make separate y_predicts for each dataset to make them easy to compare to 
 //the raw data
 //for (i in 1:N_reef) y_predict_reef[i] = ordered_logistic_rng(a_yr_reef[year_id_reef[i]]+a_region[region_reef]*sd_region+a_site_reef[site_reef[i]]*sd_site_reef+X_reef[i,]*beta_reef,c_reef); 
 //for (z in 1:N_ow) y_predict_ow[z] = ordered_logistic_rng(a_yr_ow[year_id_ow[z]]+a_region[region_ow]*sd_region+a_site_ow[site_ow[z]]*sd_site_ow+a_dv_ow[diver_ow[z]]*sd_dv_ow+X_ow[z,]*beta_ow,c_ow); 
 //for (j in 1:N_iucn) y_predict_iucn[j] = neg_binomial_2_log_rng(a_yr_iucn[year_id_iucn[j]]+a_region[region_iucn]*sd_region+a_site_iucn[site_iucn[j]]*sd_site_iucn+a_source_iucn[source_iucn[j]]*sd_source_iucn+X_iucn[j,]*beta_iucn, phi_iucn);
 //for (k in 1:N_aba) y_predict_aba[k] = neg_binomial_2_log_rng(a_yr_aba[year_id_aba[k]]+a_region[region_aba]*sd_region+a_site_aba[site_aba[k]]*sd_site_aba+X_aba[k,]*beta_aba, phi_aba);
 //for (m in 1:N_multi) {
 //  if(bernoulli_logit_rng(theta_multi) == 1)
 //  y_predict_multi[m] = 0;
 //  else
 //  y_predict_multi[m] = neg_binomial_2_log_rng(a_yr_multi[year_id_multi[m]]+a_region[region_multi]*sd_region+a_site_multi[site_multi[m]]*sd_site_multi+X_multi[m,]*beta_multi, phi_multi);
 //}
 //and a single log likelihood for the whole model for calculating looic later
 //if necessary
 for (i in 1:N_reef) log_lik[i] = ordered_logistic_lpmf(y_reef[i]|a_yr_reef[year_id_reef[i]]+a_region_reef[year_id_reef[i],region_reef[i]]+a_site_reef[site_reef[i]]*sd_site_reef+X_reef[i,]*beta_reef,c_reef);
 for (z in 1:N_ow) log_lik[N_reef+z] = ordered_logistic_lpmf(y_ow[z]|a_yr_ow[year_id_ow[z]]+a_region_ow[year_id_ow[z],region_ow[z]]+a_site_ow[site_ow[z]]*sd_site_ow+a_dv_ow[diver_ow[z]]*sd_dv_ow+X_ow[z,]*beta_ow,c_ow);
 for (j in 1:N_iucn) log_lik[N_reef+N_ow+j] = neg_binomial_2_log_lpmf(y_iucn[j]|a_yr_iucn[year_id_iucn[j]]+a_region_iucn[year_id_iucn[j],region_iucn[j]]+a_site_iucn[site_iucn[j]]*sd_site_iucn+a_source_iucn[source_iucn[j]]*sd_source_iucn+X_iucn[j,]*beta_iucn, phi_iucn);
 for (k in 1:N_aba) log_lik[N_reef+N_ow+N_iucn+k] = neg_binomial_2_log_lpmf(y_aba[k]|a_yr_aba[year_id_aba[k]]+a_region_aba[year_id_aba[k],region_aba[k]]+a_site_aba[site_aba[k]]*sd_site_aba+X_aba[k,]*beta_aba, phi_aba);
  for (m in 1:N_multi){
   if(y_multi[m] == 0) {
     log_lik[N_reef+N_ow+N_iucn+N_aba+m] = log_sum_exp(bernoulli_logit_lpmf(1 | theta_multi),
                                                       bernoulli_logit_lpmf(0 | theta_multi)
                                                           + neg_binomial_2_log_lpmf(y_multi[m] | a_yr_multi[year_id_multi[m]]+a_region_multi[year_id_multi[m],region_multi[m]]+a_site_multi[site_multi[m]]*sd_site_multi+X_multi[m,]*beta_multi, phi_multi));
   }
   else {
     log_lik[N_reef+N_ow+N_iucn+N_aba+m] = bernoulli_logit_lpmf(0 | theta_multi)
                                               + neg_binomial_2_log_lpmf(y_multi[m] | a_yr_multi[year_id_multi[m]]+a_region_multi[year_id_multi[m],region_multi[m]]+a_site_multi[site_multi[m]]*sd_site_multi+X_multi[m,]*beta_multi, phi_multi);
   }
  }
 }


//leave an empty line at the end

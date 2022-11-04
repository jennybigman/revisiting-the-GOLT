data {
  int<lower=1> N; //the number of observations
  int<lower=1> J; //the number of groups
  int<lower=1,upper=J> sp[N]; //vector of group indices (group identifier)
  vector[N] LogGSAcm2; //the response variable model 1
  vector[N] LogCenterMassG; // preidctor model1 (mass)
  vector[J] GP; // response variable model 2
  vector[J] mean_GSA; // for estimating gill area index
  vector[J] mean_mass; // for estimating gill area index

  }

parameters {
  
  //level 1
  real global_int;
  real global_slope;
  vector[J] beta_int;
  vector[J] beta_slope;
  real<lower=0> sigma;
  
    
  //level 2
  real aGP_GAI;
  real bGP_GAI;
  real<lower=0> sigma_GP_GAI;
  
  }

transformed parameters {

  vector[N] mu_LogGSAcm2; //linear predictor model 1
  
  real beta_slope_ref; // intercept of first contrast
  vector[J] beta_slopes; // vector of intercepts extract from first model
  
  vector[J] GAI;
  vector[J] LogGAI;
  vector[J] LogGAI_std;

  vector[J] mu_GP_GAI;

  for(n in 1:N) {

  mu_LogGSAcm2[n] = ((global_int + beta_int[sp[n]]) + (global_slope + beta_slope[sp[n]]) * LogCenterMassG[n]);
  }
  
  beta_slope_ref = global_slope;
  beta_slopes[1:J] = beta_slope_ref + beta_slope[1:J];
  
  for(j in 1:J){
   GAI[j] = ((mean_GSA[j])/(mean_mass[j] ^ beta_slopes[j])); 
  }
  
  LogGAI = log10(GAI);
  
  LogGAI_std = ((LogGAI - mean(LogGAI)) / sd(LogGAI));
  
  mu_GP_GAI = aGP_GAI + bGP_GAI * LogGAI_std;
   
  }
  
model {
  
  global_int       ~ student_t(3, 0, 10);
  global_slope     ~ student_t(3, 0, 10);
  beta_int         ~ student_t(3, 0, 10);
  beta_slope       ~ student_t(3, 0, 10);
  sigma            ~ cauchy(0, 10); 
      
  LogGSAcm2 ~ normal(mu_LogGSAcm2,sigma);
  
  aGP_GAI      ~ student_t(3, 0, 10); // prior on intercept model 2
  bGP_GAI      ~ student_t(3, 0, 10); // prior on slope model 2
  sigma_GP_GAI  ~ cauchy(0, 10); // prior on variance model 2
  
  GP ~ normal(mu_GP_GAI, sigma_GP_GAI); 
  

  
  }

generated quantities {

   vector[J + N] log_lik;

  for (i in 1:N){

    log_lik[i] = normal_lpdf(LogGSAcm2[i]| mu_LogGSAcm2[i], sigma);

  }

  for(i in (N+1):(N+J)){

    log_lik[i] = normal_lpdf(GP[i-N] | mu_GP_GAI[i-N], sigma_GP_GAI);

  }}


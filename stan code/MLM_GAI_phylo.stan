data {
  int<lower=1> N; //the number of observations
  int<lower=1> J; //the number of groups
  int<lower=1,upper=J> sp[N]; //vector of group indices (group identifier)
  vector[N] LogGSAcm2; //the response variable model 1
  vector[N] LogCenterMassG; // preidctor model1 (mass)
  vector[J] GP; // response variable model 2
  vector[J] mean_GSA; // for estimating gill area index
  vector[J] mean_mass; // for estimating gill area index
  matrix[J, J] d_mat; // sigma matrix
  matrix[J, J] vcov_mat; // vcov matrix

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
  
  real<lower=0,upper=1>lambda; // phylogenetic signal

  
  }

transformed parameters {

  vector[N] mu_LogGSAcm2; //linear predictor model 1
  
  real beta_slope_ref; // intercept of first contrast
  vector[J] beta_slopes; // vector of intercepts extract from first model
  
  vector[J] GAI;
  vector[J] LogGAI;
  vector[J] LogGAI_std;

  vector[J] mu_GP_GAI;
  
  matrix[J, J] sigma_mat;
  matrix[J, J] sigma_total;


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
  
  sigma_mat = (1-lambda)*d_mat + lambda*vcov_mat;
  sigma_total = sigma_mat * sigma_GP_GAI;

   
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
  
  lambda  ~ uniform(0,1);
  
  }

generated quantities {

   real log_lik;

    log_lik = normal_lpdf(LogGSAcm2| mu_LogGSAcm2, sigma) +
              multi_normal_lpdf(GP | mu_GP_GAI, sigma_total);

  }


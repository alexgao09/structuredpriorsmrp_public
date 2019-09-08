data {
  int<lower=0> N; // number of samples
  int<lower=0> N_groups_age; // the number of groups for age
  int<lower=0> N_groups_income; // the number of groups for income
  int<lower=1,upper=N_groups_age> age[N]; // the column vector of design matrix X for age
  int<lower=1,upper=N_groups_income> income[N]; // the column vector of design matrix X for income
  
  int<lower=1,upper=51> state_index[N]; // the index (from 1 to 51) representing the state for every datapoint
  int<lower=1,upper=5> region_index[51]; // the index (from 1 to 5) representing the region the state is in for all 51 states
  real<lower=0,upper=1> state_vs[51]; // the 2004 Republican vote share for every state. This is hard coded for now.
  real<lower=0,upper=1> relig[51]; // the 2004 conservative religion percentage in every state
  
  int y[N]; // the response vector
}
parameters {
  vector[N_groups_age] U_age; // the random effect for age, not multiplied by sigma_age
  vector[N_groups_income] U_income; // the random effect for income, not multiplied by sigma_income
  
  vector[51] U_state; // the random effect for state, not multiplied by sigma_state
  vector[5] U_region; // the nested random effect for region, not multiplied by sigma_region
  
  real<lower=0> sigma_age; // sd of U_age (hyperparam). halfnormal prior put on this.
  real<lower=0> sigma_income; // sd of U_income (hyperparam). halfnormal prior put on this.
  
  real<lower=0> sigma_state; // sd of state
  real<lower=0> sigma_region; // sd of region
  real beta_state; // coeff. for Republican vote share in every state
  real beta_relig; // coeff. for conservative religion share in every state
  
  real intercept; // the intercept (global fixed effect)
  real<lower=0,upper=1> rho; // the autoregressive coefficient untransformed
}
transformed parameters { 
  vector[N_groups_age] U_age_transformed;
  vector[N_groups_income] U_income_transformed;
  
  vector[51] U_state_transformed;
  vector[5] U_region_transformed;
  
  vector[N] yhat;
  real<lower=-1,upper=1> rho_transformed;

  rho_transformed = (rho * 2) - 1; // the autoregressive coefficient

  U_age_transformed = sigma_age * U_age; // the random effect for age
  U_income_transformed = sigma_income * U_income; // the random effect for income
  U_region_transformed = sigma_region * U_region; // the random effect for region
  
  // noncentered parameterization for U_state_transformed
  for (j in 1:51) {
    U_state_transformed[j] = (U_state[j] * sigma_state) + ( U_region_transformed[region_index[j]] + (beta_state * state_vs[j]) + (beta_relig * relig[j]) );
  }

  for (i in 1:N) {
    yhat[i] = intercept + U_age_transformed[age[i]] + U_income_transformed[income[i]] + U_state_transformed[state_index[i]]; // the linear predictor at each point
  }
  
}
model {
  sigma_age ~ normal(0, 1); // sigma_A ~ halfnormal(0,1)
  sigma_income ~ normal(0, 1); // sigma_I ~ halfnormal(0,1)
  rho ~ beta(0.5, 0.5); // prior on autoregressive coefficient
  
  sigma_region ~ normal(0, 1);
  sigma_state ~ normal(0, 1);
  
  U_age[1] ~ normal(0, 1/sqrt(1-rho_transformed^2)); // before it was normal(0, 1) but this is wrong
  for (j in 2:N_groups_age) {
    U_age[j] ~normal(rho_transformed * U_age[j-1],1);
  }
  
  U_income ~ normal(0, 1); // random effect for income is normal
  
  U_region ~ normal(0, 1); 
    
  U_state ~ normal(0, 1);

  beta_state ~ normal(0, 1);
  beta_relig ~ normal(0, 1);
  intercept ~ normal(0, 1); 

  
  for (i in 1:N) {
    y[i] ~ bernoulli(inv_logit(yhat[i])); // the response
  }
  
}

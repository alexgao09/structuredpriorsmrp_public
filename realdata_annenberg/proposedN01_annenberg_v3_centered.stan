data {
  int<lower=0> N; // number of samples
  int<lower=0> N_groups_age; // the number of groups for age
  int<lower=0> N_groups_race; // the number of groups for race
  int<lower=0> N_groups_sex; // the number of groups for sex
  
  int<lower=0> N_groups_education; // the number of groups for education
  int<lower=0> N_groups_income; // the number of groups for income
  int<lower=0> N_groups_state; // the number of groups for state
  
  int<lower=1,upper=N_groups_age> age[N]; // the column vector of design matrix X for age
  int<lower=1,upper=N_groups_race> race[N]; // the column vector of design matrix X for race
  int<lower=0,upper=N_groups_sex-1> sex[N]; // the column vector of design matrix X for sex. Vector of 0/1
  
  int<lower=1,upper=N_groups_education> education[N]; // the column vector of design matrix X for education
  int<lower=1,upper=N_groups_income> income[N]; // the column vector of design matrix X for income
  int<lower=1,upper=N_groups_state> state[N]; // the column vector of design matrix X for state
  
  int<lower=1,upper=5> region_index[51]; // the index representing the region that a state is in according to jkastell 
  real<lower=0,upper=1> state_vs[51]; // the 2004 Republican vote share from jkastell
  real<lower=0,upper=1> relig[51]; // the 2004 conservative religion percentage
  
  int y[N]; // the response vector
}
parameters {
  vector[N_groups_age] U_age; // the random effect for age, not multiplied by sigma_age
  vector[N_groups_race] U_race; // the random effect for race, not multiplied by sigma_race
  vector[N_groups_income] U_income; 
  vector[N_groups_state] U_state; // we are using the centered parameterization for U_state now, so U_state_transformed = U_state now
  vector[N_groups_education] U_education;
  
  vector[5] U_region; // the nested random effect for region
  
  real<lower=0> sigma_age; // sd of U_age (hyperparam).
  real<lower=0> sigma_race; // sd of U_race (hyperparam).
  real<lower=0> sigma_income;
  real<lower=0> sigma_state;
  real<lower=0> sigma_education;
  
  real<lower=0> sigma_region;
  real beta_state; // coef for Republican vote share in every state
  real beta_relig; // coef for conservative religion share in every state
  
  real intercept; // the intercept (global fixed effect)
  real beta_sex;
  
  
}
transformed parameters { 
  vector[N_groups_age] U_age_transformed;
  vector[N_groups_race] U_race_transformed;
  vector[N_groups_income] U_income_transformed;
  vector[N_groups_state] U_state_transformed;
  vector[N_groups_education] U_education_transformed;
  
  vector[5] U_region_transformed;
  
  vector[N] yhat;

  U_age_transformed = sigma_age * U_age; // the random effect for age
  U_race_transformed = sigma_race * U_race; // the random effect for race
  U_income_transformed = sigma_income * U_income;
  //U_state_transformed = sigma_state * U_state;
  U_education_transformed = sigma_education * U_education;
  
  U_region_transformed = sigma_region * U_region;
  
  U_state_transformed = U_state;

  for (i in 1:N) {
    yhat[i] = intercept + U_age_transformed[age[i]] + U_race_transformed[race[i]] + U_income_transformed[income[i]] + U_state_transformed[state[i]] + U_education_transformed[education[i]] + beta_sex * sex[i]; // the linear predictor at each point
  }
  
}
model {
  sigma_age ~ normal(0,1); // sigma_A ~ lognormal(0,1). hyperparam.
  sigma_race ~ normal(0,1); // sigma_I ~ lognormal(0,1). hyperparam.
  sigma_income ~ normal(0,1);
  sigma_state ~ normal(0,1);
  sigma_education ~ normal(0,1);

  sigma_region ~ normal(0,1);

  for (j in 2:N_groups_age) {
    U_age[j] ~normal(U_age[j-1],1);
  }

  sum(U_age) ~ normal(0, 0.01 * N_groups_age); // constraint so we can write likelihood for rw(1).
  
  U_race ~ normal(0, 1); // random effect for race is normal
  U_income ~ normal(0, 1); 
  //U_state ~ normal(0, 1);
  U_education ~ normal(0, 1);
  
  U_region ~ normal(0,1);
  beta_state ~ normal(0,1);
  beta_relig ~ normal(0,1);
  
  for (j in 1:51) {
    U_state[j] ~ normal(U_region_transformed[region_index[j]] + (beta_state * state_vs[j]) + (beta_relig * relig[j]), sigma_state);
  }
  
  intercept ~ normal(0, 1); 
  beta_sex ~ normal(0, 1);
  
  
  for (i in 1:N) {
    y[i] ~ bernoulli(inv_logit(yhat[i])); // the response
  }
  
}

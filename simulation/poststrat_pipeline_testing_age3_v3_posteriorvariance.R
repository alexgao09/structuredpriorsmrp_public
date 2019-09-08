rm(list=ls()) # clear previous workspace
# posteriorvariance calculates posterior width for each poststratification cell

# v3 has 2004 Democratic vote share for each state as a predictor
# the v3 stan files accompany this file

options(bitmapType="cairo")

set.seed(21)

# this script allows us to dynamically change stuff 
# v2 gets a maximal poststrat matrix then reduces by age groups

# this creates the whole poststrat simulation process
# poststrat_pipeline_base.R uses sampling probabilities of each poststrat. cell as N_j / N

# load libraries --------------------------------------
library(rstan)
library(arm)
library(ggplot2)
library(ggridges)
library(reshape2)
library(stringr)
library(dplyr)
library(matrixStats)

options(mc.cores = parallel::detectCores()) # use multiple cores. default is 4.
rstan_options(auto_write = TRUE) # save compiled stan model to hard disk so no need to recompile

fig_width = 2400
fig_height = 2400

# -----------------------------------------------------
# Global variables ------------------------------------
# -----------------------------------------------------
save_ridgeplots = TRUE 
sample_size = 100
runs = 200
r = 1:9/10
store_stanobjects = FALSE # do you want to store fitted stan model after each run?
income_multiplier = 1 # partitions income into more categories
age_grouping_multiplier = 12 # how much we take the maximal poststratification. make sure this can divide 60

age = 21:80
age_delta = 60 # this equals length(age)
age_tr = (age-21)/age_delta # (transform to a unit scale for easy mapping)

# ----
# ----

# old truepref curves
#coef_age = dbeta(age_tr,2,2) # cap
#coef_age = 1-dbeta(age_tr,2,2) # cup
#coef_age = (1 - 2*exp(-age_tr/.2)) # increasing

# new truepref curves that are steeper
#coef_age = 2*dbeta(age_tr,2,2) # cap shape
#coef_age = 3-2*dbeta(age_tr,2,2) # true coefficients for d.g process, u shape
coef_age = (.7 - 3*exp(-age_tr/.2)) # increasing shape

coef_income = c(.1 * rep(1, income_multiplier),
                0 * rep(1, income_multiplier),
                -.2 * rep(1, income_multiplier),
                .2 * rep(1, income_multiplier))

intercept_term = 0 # for increasing curve only
#intercept_term = -1.5

# Load state-level data from https://scholar.princeton.edu/jkastellec/publications
StateLevel = foreign::read.dta("state_level_update.dta",convert.underscore = TRUE)
StateLevel = StateLevel[order(StateLevel$sstate.initnum),]

# Load state-level data from https://scholar.princeton.edu/jkastellec/publications
MarriageData = foreign::read.dta("gay_marriage_megapoll.dta", convert.underscore = TRUE) 

# check that row names of Statelevel are equal to sstate.initnum
sum(as.numeric(as.character(rownames(StateLevel))) == StateLevel$sstate.initnum)

stateregions = unique(MarriageData[,c("statename", "state.initnum", "region.cat", "region")])

stateregions_final = dplyr::left_join(x = StateLevel, 
                                      y = stateregions,
                                      by = c("sstate.initnum"="state.initnum"))[,c("sstate.initnum", 
                                                                                   "sstate",
                                                                                   "sstatename",
                                                                                   "region",
                                                                                   "region.cat",
                                                                                   "kerry.04",
                                                                                   "p.evang",
                                                                                   "p.mormon")]

# Alaska and Hawaii doesn't have a region so we  set it to west
stateregions_final[stateregions_final$sstatename == "Alaska", c("region")] = "west"
stateregions_final[stateregions_final$sstatename == "Hawaii", c("region")] = "west"
stateregions_final[stateregions_final$sstatename == "Alaska", c("region.cat")] = 4
stateregions_final[stateregions_final$sstatename == "Hawaii", c("region.cat")] = 4

# we will use the sstate.initnum order from stateregions_final
# coef_state comes from data used in http://www.princeton.edu/~jkastell/MRP_primer/mrp_primer.pdf
coef_state = (stateregions_final$kerry.04)/100 # value for 2004 Democratic vote share in every state. This is X_{\text{state-vs},j} for j \in [51]

# coef_relig comes from data used in http://www.princeton.edu/~jkastell/MRP_primer/mrp_primer.pdf. this is the percentage of conservative religions in every state
coef_relig = (stateregions_final$p.evang + stateregions_final$p.mormon)/100

beta_state = 0.5 # coefficient for state-level predictor democratic vote share. play around with this to get true pref on average age-cat in (0.1, 0.9)
beta_relig = -0.5

# ----
# ----

# below are subject to change
p_age = rep(1, length(age_tr)) / length(age_tr)
p_income = c(0.1 * rep(1, income_multiplier),
             0.2 * rep(1, income_multiplier),
             0.3 * rep(1, income_multiplier),
             0.4 * rep(1, income_multiplier)) # probability for landing in bucket for income
p_state = rep(1, 51)/51

# assume equal response rate
p_response_age = rep(1, length(p_age)) / length(p_age)
p_response_income = rep(1, length(p_income)) / length(p_income)
p_response_state = rep(1, 51)/51

total_pop = 1e8 # / (length(p_age) * length(p_income)) # the number of people in the overall population

response_binary = TRUE # this indicator determines whether response is binary or not
response_normal_sd = 100

stan_filename_baseline = "baselinemeanzeroN01_v3.stan"
modelname_baseline = "baselinemeanzero"

stan_filename_ar = "proposedarN01_v3.stan" # the stan file to compile
modelname_ar = "proposedarN01"

stan_filename_rw = "proposedN01_v3.stan"
modelname_rw = "proposedN01"

iterations = 2000 # the number of times to run markov chain (50% burn-in is default)
num_chains = 4 # number of markov chains to run

store_divs = FALSE # if this is true then we store the stanfit objects that have divergent transitions

# used to calculate posterior width
lowerquantile = 0.1
upperquantile = 0.9

# If save ridgeplots is FALSE, then we save ridgeplot dataframe as csv instead. 
# we do this because the jupiter server can't plot pretty pngs
# else we save pngs

# -----------------------------------------------------
# -----------------------------------------------------
# -----------------------------------------------------

# GET MAXIMAL POSTSTRAT MATRIX

poststrat = expand.grid(age = 1:length(p_age),
                        income = 1:length(p_income),
                        state = 1:length(p_state))

# get number of people in each group for the whole population
for (j in 1:(length(p_age) * length(p_income) * length(p_state))) { # this for loop goes through every row in the post-strat. df
  poststrat$N[j] = round(total_pop * 
                           p_age[poststrat[j, 1]] * 
                           p_income[poststrat[j, 2]] *
                           p_state[poststrat[j, 3]] )
}

# MAIN LOOP
if (store_stanobjects==TRUE){
  stanobjects_list_baseline = list() # stores the stan objects after each run
  postpred_sim_list_baseline = list() # a list containing the posterior linear predictors at each run
  
  stanobjects_list_ar = list() # stores the stan objects after each run
  postpred_sim_list_ar = list() # a list containing the posterior linear predictors at each run
  
  stanobjects_list_rw = list() # stores the stan objects after each run
  postpred_sim_list_rw = list() # a list containing the posterior linear predictors at each run
  runs = 1 # set this to 1 so we can store stuff in two lists above
}

poststrat_final_list = list() # stores the maximal poststrat_final matrices for each iteration 
poststrat_final_reduced_list = list() # stores the reduced poststrat matrices for each p 

# compile all three models
m_baseline = stan_model(file = stan_filename_baseline,
                        model_name = modelname_baseline)

m_ar = stan_model(file = stan_filename_ar,
                  model_name = modelname_ar)

m_rw = stan_model(file = stan_filename_rw,
                  model_name = modelname_rw)

# stores the vectors that contain the number of divergences for each p, for each model
num_divs_list_baseline = list()
num_divs_list_ar = list() 
num_divs_list_rw = list() 

# below three matrices store the num of divergences for each run, for each probability index
num_divs_counter_mat_baseline = matrix(0, length(r), runs) # should be 9 x 200
num_divs_counter_mat_ar = matrix(0, length(r), runs) # should be 9 x 200
num_divs_counter_mat_rw = matrix(0, length(r), runs) # should be 9 x 200
p_counter = 1
r_numericindex = 1:length(r)

for (p in r) { # setting p = 1 results in model fitting error
  
  p_response_age = rep(0, length(p_age))
  p_response_age[1:(length(p_age) * 2/3)] = (1-p)
  p_response_age[(length(p_age) * 2/3 + 1):length(p_age)] = p
  
  p_response = rep(NA, length(p_age) * length(p_income) * length(p_state)) # this contains the response rate for a group in post-strat. df. 
  for (j in 1:(length(p_age) * length(p_income) * length(p_state))) { # this for loop
    p_response[j] = p_response_age[poststrat[j,1]] *
      p_response_income[poststrat[j,2]] * 
      p_response_state[poststrat[j,3]]
  }
  
  
  # THIS CHUNK OF CODE NEEDS TO CHANGE
  # WE NEED A true_pref_age THEN DO LOOP BELOW
  true_pref = rep(NA, length(p_age) * length(p_income) * length(p_state)) # this contains the true preference probability of a group in the population
  for (j in 1:(length(p_age) * length(p_income) * length(p_state))) {
    true_pref[j] = invlogit(intercept_term + 
      coef_age[poststrat[j, 1]] + 
        coef_income[poststrat[j, 2]] +
        beta_state * coef_state[poststrat[j, 3]] +
        beta_relig * coef_relig[poststrat[j, 3]] )
  }
  
  # poststrat_final is the final postratification population matrix
  poststrat_final = cbind(poststrat,
                          p_response = p_response,
                          true_pref = true_pref) 
  
  age_group = data.frame(age_cat = rep(sapply(1:age_grouping_multiplier, toString),
                                       each = age_delta/age_grouping_multiplier),
                         age = age - min(age - 1))
  age_group$age_cat = factor(age_group$age_cat, levels = sapply(1:age_grouping_multiplier, toString)) # makes sure levels are ordered
  
  
  # maximal poststratification matrix
  poststrat_final = inner_join(poststrat_final,
                               age_group,
                               by = "age")
  
  poststrat_final_temp = poststrat_final %>% 
    group_by(age_cat) %>% 
    mutate(s_N = sum(N))
  
  # mrp_subestimates are the true preferences for each age group
  poststrat_final_temp$mrp = poststrat_final_temp$true_pref * poststrat_final_temp$N / poststrat_final_temp$s_N
  mrp_subestimates = poststrat_final_temp %>% group_by(age_cat) %>% summarise(mrp = sum(mrp))
  
  # get reduced poststratification matrix now
  poststrat_final_reduced_temp = poststrat_final %>%
    group_by(income, age_cat, state) %>%
    mutate(s_N = sum(N))
  
  # poststrat_final_reduced_temp contains MRP estimates for the subgroup group_by(income, age_cat, state)
  poststrat_final_reduced_temp$true_pref_grouped = poststrat_final_reduced_temp$true_pref * poststrat_final_reduced_temp$N/poststrat_final_reduced_temp$s_N
  
  # the reduced poststratification matrix used for model evaluation
  poststrat_final_reduced = data.frame(poststrat_final_reduced_temp %>% 
                                         group_by(income, age_cat, state) %>% 
                                         summarise(true_pref_grouped_final = sum(true_pref_grouped),
                                                   s_N = first(s_N)))
  
  poststrat_final_reduced$age_cat = as.numeric(as.character(poststrat_final_reduced$age_cat))
  
  # store maximal and reduced poststrat matrix in list
  poststrat_final_list[[toString(p)]] = poststrat_final 
  poststrat_final_reduced_list[[toString(p)]] = poststrat_final_reduced 
  
  
  # FOR LOOP ----------------------------------------------
  
  model_popn_pref_mat_baseline = rep(NA, runs) # stores the preference estimate from MRP
  model_popn_pref_sd_mat_baseline = rep(NA, runs) # stores sd of above estimate
  
  model_popn_pref_mat_ar = rep(NA, runs) # stores the preference estimate from MRP
  model_popn_pref_sd_mat_ar = rep(NA, runs) # stores sd of above estimate
  
  model_popn_pref_mat_rw = rep(NA, runs) # stores the preference estimate from MRP
  model_popn_pref_sd_mat_rw = rep(NA, runs) # stores sd of above estimate
  
  sample_mat = rep(NA, runs) # stores the preference estimate from mean of sample
  sample_sd_mat = rep(NA, runs) # stores sd of above estimate
  
  # means of posterior linear predictors for all runs k
  sample_cell_estimates_baseline = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  sample_cell_estimates_ar = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  sample_cell_estimates_rw = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  
  # medians of posterior linear predictors for all runs k
  sample_cell_estimates_median_baseline = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  sample_cell_estimates_median_ar = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  sample_cell_estimates_median_rw = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  
  # posterior width
  sample_cell_estimates_width_baseline = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  sample_cell_estimates_width_ar = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  sample_cell_estimates_width_rw = matrix(NA, runs, dim(poststrat_final_reduced)[1])
  
  # below four are lists of length 9, with each list element being a matrix of dimension runs x age_grouping_multiplier
  quantile10_list_baseline = matrix(NA, runs, age_grouping_multiplier) # every column contains 10 percentile for a age_cat
  quantile90_list_baseline = matrix(NA, runs, age_grouping_multiplier)
  median_quantile_list_baseline = matrix(NA, runs, age_grouping_multiplier)
  mean_list_baseline = matrix(NA, runs, age_grouping_multiplier)

  quantile10_list_ar = matrix(NA, runs, age_grouping_multiplier) # every column contains 10 percentile for a age_cat
  quantile90_list_ar = matrix(NA, runs, age_grouping_multiplier)
  median_quantile_list_ar = matrix(NA, runs, age_grouping_multiplier)
  mean_list_ar = matrix(NA, runs, age_grouping_multiplier)

  quantile10_list_rw = matrix(NA, runs, age_grouping_multiplier) # every column contains 10 percentile for a age_cat
  quantile90_list_rw = matrix(NA, runs, age_grouping_multiplier)
  median_quantile_list_rw = matrix(NA, runs, age_grouping_multiplier)
  mean_list_rw = matrix(NA, runs, age_grouping_multiplier)

  
  if (response_binary==TRUE) {
    response_tag = "binary"
  }else {
    response_tag = "cts"
  }
  
  num_divs_baseline = rep(0, runs) # stores the number of divergences for each k
  num_divs_ar = rep(0, runs) # stores the number of divergences for each k
  num_divs_rw = rep(0, runs) # stores the number of divergences for each k
  
  
  for (k in 1:runs) {
    print(paste("On run", k, ", on p_response_age[3] =", p))
    
    # Sample from population ------------------------------
    
    # sample from population in terms of group indices 
    # we sample from each group porportional to the response probability weighted by the subgroup size in population
    sample_ = sample(length(p_age) * length(p_income) * length(p_state),
                     sample_size,
                     replace = TRUE,
                     prob = (poststrat_final$p_response * poststrat_final$N)/sum((poststrat_final$p_response * poststrat_final$N)) )
    
    # get response of sample
    if (response_binary == TRUE) { # binary response
      y_sample_ = rbinom(sample_size, 1, poststrat_final$true_pref[sample_])
    }else{ # normal response
      y_sample_ = rnorm(sample_size, logit(poststrat_final$true_pref[sample_]), response_normal_sd)
    }
    
    # get covariates for every row of sample
    age_sample = poststrat_final[sample_, 1]
    income_sample = poststrat_final[sample_, 2]
    state_sample = poststrat_final[sample_, 3]
    
    sample_final = data.frame(pref = y_sample_,
                              age = age_sample,
                              income = income_sample,
                              state = state_sample,
                              region_index_sample = stateregions_final$region.cat[state_sample] ) # grab the region that every state in the sample corresponds to
    
    sample_final_ = inner_join(x = sample_final, y = age_group,
                               by = "age")
    sample_final_$age_cat = as.numeric(as.character(sample_final_$age_cat)) # stan needs numeric entries, not factors
    
    # get (# said yes)/(# in age subgroup)
    empirical_probs = data.frame(sample_final_ %>% group_by(age_cat) %>% summarise(empirical_prob=sum(pref)/n(), n=n()))
    
    
    # fit model -------------------------------------------
    
    # sample from model 
    fit_baseline = sampling(m_baseline, 
                            data = list(N = dim(sample_final_)[1],
                                        N_groups_age = age_grouping_multiplier,
                                        N_groups_income = 4 * income_multiplier, # there are 4 categories in income
                                        age = sample_final_$age_cat,
                                        income = sample_final_$income,
                                        
                                        state_index = sample_final_$state,
                                        state_vs = coef_state,
                                        relig = coef_relig,
                                        region_index = stateregions_final$region.cat,
                                        
                                          y = sample_final_$pref),
                            iter=iterations, chains=num_chains,
                            control=list(max_treedepth=15, adapt_delta=0.99),
                            seed = 21,
                            chain_id = num_chains*3*(k-1) + 1 + num_chains*3*(r_numericindex[p_counter] - 1)*(runs - 1) + num_chains*3*(r_numericindex[p_counter] - 1))
    
    fit_ar = sampling(m_ar, 
                      data = list(N = dim(sample_final_)[1],
                                  N_groups_age = age_grouping_multiplier,
                                  N_groups_income = 4 * income_multiplier, # there are 4 categories in income
                                  age = sample_final_$age_cat,
                                  income = sample_final_$income,
                                  
                                  state_index = sample_final_$state,
                                  state_vs = coef_state,
                                  relig = coef_relig,
                                  region_index = stateregions_final$region.cat,
                                  
                                  y = sample_final_$pref),
                      iter=iterations, chains=num_chains,
                      control=list(max_treedepth=15, adapt_delta=0.99),
                      seed = 21,
                      chain_id = num_chains*3*(k-1) + 5 + num_chains*3*(r_numericindex[p_counter] - 1)*(runs - 1) + num_chains*3*(r_numericindex[p_counter] - 1))
    
    fit_rw = sampling(m_rw, 
                      data = list(N = dim(sample_final_)[1],
                                  N_groups_age = age_grouping_multiplier,
                                  N_groups_income = 4 * income_multiplier, # there are 4 categories in income
                                  age = sample_final_$age_cat,
                                  income = sample_final_$income,
                                  
                                  state_index = sample_final_$state,
                                  state_vs = coef_state,
                                  relig = coef_relig,
                                  region_index = stateregions_final$region.cat,
                                  
                                  y = sample_final_$pref),
                      iter=iterations, chains=num_chains,
                      control=list(max_treedepth=15, adapt_delta=0.99),
                      seed = 21,
                      chain_id = num_chains*3*(k-1) + 9 + num_chains*3*(r_numericindex[p_counter] - 1)*(runs - 1) + num_chains*3*(r_numericindex[p_counter] - 1))
    
    print(paste("Number of divergence transitions:",
                sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[1]][,'divergent__']) +
                  sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[2]][,'divergent__']) +
                  sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[3]][,'divergent__']) +
                  sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[4]][,'divergent__'])
    ))
    
    print(paste("Number of divergence transitions:",
                sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[1]][,'divergent__']) +
                  sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[2]][,'divergent__']) +
                  sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[3]][,'divergent__']) +
                  sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[4]][,'divergent__'])
    ))
    
    print(paste("Number of divergence transitions:",
                sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[1]][,'divergent__']) +
                  sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[2]][,'divergent__']) +
                  sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[3]][,'divergent__']) +
                  sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[4]][,'divergent__'])
    ))
    
    num_divs_baseline[k] = sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[1]][,'divergent__']) +
      sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[2]][,'divergent__']) +
      sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[3]][,'divergent__']) +
      sum(get_sampler_params(fit_baseline, inc_warmup=FALSE)[[4]][,'divergent__'])
    
    num_divs_ar[k] = sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[1]][,'divergent__']) +
      sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[2]][,'divergent__']) +
      sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[3]][,'divergent__']) +
      sum(get_sampler_params(fit_ar, inc_warmup=FALSE)[[4]][,'divergent__'])
    
    num_divs_rw[k] = sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[1]][,'divergent__']) +
      sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[2]][,'divergent__']) +
      sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[3]][,'divergent__']) +
      sum(get_sampler_params(fit_rw, inc_warmup=FALSE)[[4]][,'divergent__'])
    
    # save divergent stan objects
    if (store_divs == TRUE) {
      if (num_divs_baseline[k] > 0) {
        num_divs_list_baseline[[paste(p, k)]] = fit_baseline
      }
      
      if (num_divs_ar[k] > 0) {
        num_divs_list_ar[[paste(p, k)]] = fit_ar
      }
      
      if (num_divs_rw[k] > 0) {
        num_divs_list_rw[[paste(p, k)]] = fit_rw
      }
    }
    
    fit_samples_baseline = as.data.frame(rstan::extract(fit_baseline, permuted=FALSE)) # extract posterior samples
    fit_samples_ar = as.data.frame(rstan::extract(fit_ar, permuted=FALSE)) # extract posterior samples
    fit_samples_rw = as.data.frame(rstan::extract(fit_rw, permuted=FALSE)) # extract posterior samples
    
    N_sub = rep(0, age_grouping_multiplier) # the N for each age subgroup
    for (g in 1:age_grouping_multiplier) {
      N_sub[g] = sum(poststrat_final_reduced[which(poststrat_final_reduced$age_cat==g),]$s_N)
    }
    
    # BASELINE MODEL ---------------------------
    # get posterior samples from stanfit object
    intercept_samples_baseline = tidyr::gather(fit_samples_baseline[,grepl("intercept", 
                                                                           names(fit_samples_baseline),
                                                                           fixed=TRUE)])$value # get post. samples for intercept
    
    # this matrix will store the posterior samples of the random effects for age
    U_age_samples_baseline = matrix(0, iterations * num_chains/2, age_grouping_multiplier)
    
    # extract U_age_samples
    for (j in 1:age_grouping_multiplier) {
      U_age_samples_baseline[,j] = tidyr::gather(fit_samples_baseline[,grepl(paste0("U_age_transformed[",
                                                                                    j
                                                                                    ,"]"), 
                                                                             names(fit_samples_baseline),
                                                                             fixed=TRUE)])$value
    }
    
    # this matrix will store the posterior samples of the random effects for income
    U_income_samples_baseline = matrix(0, iterations * num_chains/2, length(p_income))
    
    # extract U_income_samples
    for (j in 1:length(p_income)) {
      U_income_samples_baseline[,j] = tidyr::gather(fit_samples_baseline[,grepl(paste0("U_income_transformed[",
                                                                                       j
                                                                                       ,"]"), 
                                                                                names(fit_samples_baseline),
                                                                                fixed=TRUE)])$value
    }
    
    # this matrix will store the posterior samples for the random effects for state
    U_state_samples_baseline = matrix(0, iterations * num_chains/2, length(p_state))
    
    # extract U_state_samples
    for (j in 1:length(p_state)) {
      U_state_samples_baseline[,j] = tidyr::gather(fit_samples_baseline[,grepl(paste0("U_state_transformed[",
                                                                                      j
                                                                                      ,"]"), 
                                                                               names(fit_samples_baseline),
                                                                               fixed=TRUE)])$value
    }
    
    # this matrix stores the posterior linear predictors for each poststrat. cell
    postpred_sim_baseline = matrix(0, 
                                   dim(U_age_samples_baseline)[1], 
                                   dim(U_age_samples_baseline)[2] * dim(U_income_samples_baseline)[2] * dim(U_state_samples_baseline)[2])
    
    for (i in 1:(dim(U_age_samples_baseline)[2] * dim(U_income_samples_baseline)[2] * dim(U_state_samples_baseline)[2])) {
      postpred_sim_baseline[,i] = invlogit(intercept_samples_baseline + 
                                             U_age_samples_baseline[, as.numeric(poststrat_final_reduced[i,c("age_cat")]) ] +
                                             U_income_samples_baseline[, as.numeric(poststrat_final_reduced[i,c("income")]) ] +
                                             U_state_samples_baseline[, as.numeric(poststrat_final_reduced[i,c("state")])] )
    }
    
    mrp_samples_baseline = list() # this list stores posterior samples for subgroups
    for (g in 1:age_grouping_multiplier) {
      mrp_samples_baseline[[toString(g)]] = rep(0, dim(U_age_samples_baseline)[1]) # default is length 4000 
      
      for (i in as.numeric(rownames(poststrat_final_reduced[which(poststrat_final_reduced$age_cat==g),])) ) {
        mrp_samples_baseline[[toString(g)]] = mrp_samples_baseline[[toString(g)]] + 
          (postpred_sim_baseline[,i] * poststrat_final_reduced$s_N[i]/N_sub[g])
      }
    }

    # get quantiles for the 12 age_cat posteriors that come from the baseline model
    for (g in 1:age_grouping_multiplier) {
      quantile10_list_baseline[k, g] = quantile(mrp_samples_baseline[[toString(g)]], 0.1)
      quantile90_list_baseline[k, g] = quantile(mrp_samples_baseline[[toString(g)]], 0.9)
      median_quantile_list_baseline[k, g] = median(mrp_samples_baseline[[toString(g)]])
      mean_list_baseline[k, g] = mean(mrp_samples_baseline[[toString(g)]])
    }
    
    mrp_subsamples_allmodels = list() # this list stores matrices of posterior estimates for subgroups, where each matrix contains the three models
    for (g in 1:age_grouping_multiplier) {
      mrp_subsamples_allmodels[[toString(g)]] = mrp_samples_baseline[[toString(g)]]
    }

    # AUTOREGRESSIVE MODEL ------------------------------
    
    intercept_samples_ar = tidyr::gather(fit_samples_ar[,grepl("intercept", 
                                                               names(fit_samples_ar),
                                                               fixed=TRUE)])$value # get post. samples for intercept
    
    # this matrix will store the posterior samples of the random effects for age
    U_age_samples_ar = matrix(0, iterations * num_chains/2, age_grouping_multiplier)
    
    # extract U_age_samples
    for (j in 1:age_grouping_multiplier) {
      U_age_samples_ar[,j] = tidyr::gather(fit_samples_ar[,grepl(paste0("U_age_transformed[",
                                                                        j
                                                                        ,"]"), 
                                                                 names(fit_samples_ar),
                                                                 fixed=TRUE)])$value
    }
    
    # this matrix will store the posterior samples of the random effects for income
    U_income_samples_ar = matrix(0, iterations * num_chains/2, length(p_income))
    
    # extract U_income_samples
    for (j in 1:length(p_income)) {
      U_income_samples_ar[,j] = tidyr::gather(fit_samples_ar[,grepl(paste0("U_income_transformed[",
                                                                           j
                                                                           ,"]"), 
                                                                    names(fit_samples_ar),
                                                                    fixed=TRUE)])$value
    }
    
    # this matrix will store the posterior samples for the random effects for state
    U_state_samples_ar = matrix(0, iterations * num_chains/2, length(p_state))
    
    # extract U_state_samples
    for (j in 1:length(p_state)) {
      U_state_samples_ar[,j] = tidyr::gather(fit_samples_ar[,grepl(paste0("U_state_transformed[",
                                                                          j
                                                                          ,"]"), 
                                                                   names(fit_samples_ar),
                                                                   fixed=TRUE)])$value
    }
    
    # this matrix stores the posterior linear predictors for each poststrat. cell
    postpred_sim_ar = matrix(0, 
                             dim(U_age_samples_ar)[1], 
                             dim(U_age_samples_ar)[2] * dim(U_income_samples_ar)[2] * dim(U_state_samples_ar)[2])
    
    for (i in 1:(dim(U_age_samples_ar)[2] * dim(U_income_samples_ar)[2] * dim(U_state_samples_ar)[2])) {
      postpred_sim_ar[,i] = invlogit(intercept_samples_ar + 
                                       U_age_samples_ar[, as.numeric(poststrat_final_reduced[i,c("age_cat")]) ] +
                                       U_income_samples_ar[, as.numeric(poststrat_final_reduced[i,c("income")]) ] + 
                                       U_state_samples_ar[, as.numeric(poststrat_final_reduced[i,c("state")]) ] )
    }
    
    mrp_samples_ar = list() # this list stores posterior samples for subgroups
    for (g in 1:age_grouping_multiplier) {
      mrp_samples_ar[[toString(g)]] = rep(0, dim(U_age_samples_ar)[1]) # default is length 4000 
      
      for (i in as.numeric(rownames(poststrat_final_reduced[which(poststrat_final_reduced$age_cat==g),])) ) {
        mrp_samples_ar[[toString(g)]] = mrp_samples_ar[[toString(g)]] + 
          (postpred_sim_ar[,i] * poststrat_final_reduced$s_N[i]/N_sub[g])
      }
    }
    
    # get quantiles for the 12 age_cat posteriors that come from the ar model
    for (g in 1:age_grouping_multiplier) {
      quantile10_list_ar[k, g] = quantile(mrp_samples_ar[[toString(g)]], 0.1)
      quantile90_list_ar[k, g] = quantile(mrp_samples_ar[[toString(g)]], 0.9)
      median_quantile_list_ar[k, g] = median(mrp_samples_ar[[toString(g)]])
      mean_list_ar[k, g] = mean(mrp_samples_ar[[toString(g)]])
    }
    
    for (g in 1:age_grouping_multiplier) {
      mrp_subsamples_allmodels[[toString(g)]] = cbind(mrp_subsamples_allmodels[[toString(g)]], 
                                                      mrp_samples_ar[[toString(g)]])
    }
    
    # RANDOMWALK MODEL ----------------------------------
    
    intercept_samples_rw = tidyr::gather(fit_samples_rw[,grepl("intercept", 
                                                               names(fit_samples_rw),
                                                               fixed=TRUE)])$value # get post. samples for intercept
    
    # this matrix will store the posterior samples of the random effects for age
    U_age_samples_rw = matrix(0, iterations * num_chains/2, age_grouping_multiplier)
    
    # extract U_age_samples
    for (j in 1:age_grouping_multiplier) {
      U_age_samples_rw[,j] = tidyr::gather(fit_samples_rw[,grepl(paste0("U_age_transformed[",
                                                                        j
                                                                        ,"]"), 
                                                                 names(fit_samples_rw),
                                                                 fixed=TRUE)])$value
    }
    
    # this matrix will store the posterior samples of the random effects for income
    U_income_samples_rw = matrix(0, iterations * num_chains/2, length(p_income))
    
    # extract U_income_samples
    for (j in 1:length(p_income)) {
      U_income_samples_rw[,j] = tidyr::gather(fit_samples_rw[,grepl(paste0("U_income_transformed[",
                                                                           j
                                                                           ,"]"), 
                                                                    names(fit_samples_rw),
                                                                    fixed=TRUE)])$value
    }
    
    # this matrix will store the posterior samples of the random effects for state
    U_state_samples_rw = matrix(0,  iterations * num_chains/2, length(p_state))
    
    for (j in 1:length(p_state))  {
      U_state_samples_rw[,j] = tidyr::gather(fit_samples_rw[,grepl(paste0("U_state_transformed[",
                                                                          j
                                                                          ,"]"), 
                                                                   names(fit_samples_rw),
                                                                   fixed=TRUE)])$value
    }
    
    
    
    # this matrix stores the posterior linear predictors for each poststrat. cell
    postpred_sim_rw = matrix(0, 
                             dim(U_age_samples_rw)[1], 
                             dim(U_age_samples_rw)[2] * dim(U_income_samples_rw)[2] * dim(U_state_samples_rw)[2])
    
    for (i in 1:(dim(U_age_samples_rw)[2] * dim(U_income_samples_rw)[2] * dim(U_state_samples_rw)[2])) {
      postpred_sim_rw[,i] = invlogit(intercept_samples_rw + 
                                       U_age_samples_rw[, as.numeric(poststrat_final_reduced[i,c("age_cat")]) ] +
                                       U_income_samples_rw[, as.numeric(poststrat_final_reduced[i,c("income")]) ] + 
                                       U_state_samples_rw[, as.numeric(poststrat_final_reduced[i,c("state")]) ])
    }
    
    mrp_samples_rw = list() # this list stores posterior samples for subgroups
    for (g in 1:age_grouping_multiplier) {
      mrp_samples_rw[[toString(g)]] = rep(0, dim(U_age_samples_rw)[1]) # default is length 4000 
      
      for (i in as.numeric(rownames(poststrat_final_reduced[which(poststrat_final_reduced$age_cat==g),])) ) {
        mrp_samples_rw[[toString(g)]] = mrp_samples_rw[[toString(g)]] + 
          (postpred_sim_rw[,i] * poststrat_final_reduced$s_N[i]/N_sub[g])
      }
    }
    
    # get quantiles for the 12 age_cat posteriors that come from the rw model
    for (g in 1:age_grouping_multiplier) {
      quantile10_list_rw[k, g] = quantile(mrp_samples_rw[[toString(g)]], 0.1)
      quantile90_list_rw[k, g] = quantile(mrp_samples_rw[[toString(g)]], 0.9)
      median_quantile_list_rw[k, g] = median(mrp_samples_rw[[toString(g)]])
      mean_list_rw[k, g] = mean(mrp_samples_rw[[toString(g)]])
    }
    
    for (g in 1:age_grouping_multiplier) {
      mrp_subsamples_allmodels[[toString(g)]] = cbind(mrp_subsamples_allmodels[[toString(g)]], 
                                                      mrp_samples_rw[[toString(g)]])
    }
    
    
    # ---------------------------------------------------
    
    mrp_subsamples_melted_allmodels = list() # melted posterior samples for age subgroups for all three models
    for (g in 1:age_grouping_multiplier) {
      colnames(mrp_subsamples_allmodels[[toString(g)]]) = paste(toString(g), c("BL", "AR", "RW"))
      mrp_subsamples_melted_allmodels[[toString(g)]] = melt(mrp_subsamples_allmodels[[toString(g)]])
    }
    
    # combine all three dataframes
    mrp_subsamples_melted = c()
    for (g in 1:age_grouping_multiplier) {
      mrp_subsamples_melted = rbind(mrp_subsamples_melted, mrp_subsamples_melted_allmodels[[toString(g)]])
    }
    mrp_subsamples_melted$mod = rep(0, dim(mrp_subsamples_melted)[1])
    
    # name column indices
    colnames(mrp_subsamples_melted) = c("n", "Type", "Preference", "Model")
    
    mrp_subsamples_melted[grepl("RW", mrp_subsamples_melted$Type, fixed=TRUE),c("Model")] = "Random walk"
    mrp_subsamples_melted[grepl("AR", mrp_subsamples_melted$Type, fixed=TRUE),c("Model")] = "Autoregressive"
    mrp_subsamples_melted[grepl("BL", mrp_subsamples_melted$Type, fixed=TRUE),c("Model")] = "Baseline"
    
    
    
    points_df = data.frame(cbind(word(levels(mrp_subsamples_melted$Type), 1),
                                 levels(mrp_subsamples_melted$Type)))
    colnames(points_df) = c("age_cat", "Type")
    points_df = merge(mrp_subestimates, points_df, by="age_cat")
    points_df$age_cat = as.numeric(as.character(points_df$age_cat))
    
    # inner join to get empirical preferences on 
    points_df_final = cbind(left_join(points_df, empirical_probs, by = "age_cat"), 1)
    names(points_df_final)[dim(points_df_final)[2]] = "one"
    
    
    # save as csv/png
    if (k==1) { # create directory on run 1
      dir.create(paste0(runs,
                        "_",
                        age_grouping_multiplier,
                        "_",
                        income_multiplier,
                        "_",
                        sample_size,
                        "_",
                        response_tag,
                        "_oldest_",
                        p * 10))
    }
    
    # fix labelling of legend
    mrp_subsamples_melted$Model = factor(mrp_subsamples_melted$Model, 
                                         levels = c("Random walk",
                                                    "Autoregressive",
                                                    "Baseline"))
    
    
    
      # save png in directory created above
      png(paste0(runs,
                 "_",
                 age_grouping_multiplier,
                 "_",
                 income_multiplier,
                 "_",
                 sample_size,
                 "_",
                 response_tag,
                 "_oldest_",
                 p * 10,
                 "/",
                 toString(k),
                 "_", 
                 runs,
                 "_",
                 age_grouping_multiplier,
                 "_",
                 income_multiplier,
                 "_",
                 sample_size,
                 "_",
                 response_tag,
                 "_oldest_",
                 p * 10,
                 ".png"), 
          width = fig_width, height = fig_height)
      
      plot(ggplot(mrp_subsamples_melted, 
                  aes(x = Preference, y = Type, fill=Model)) + 
             geom_density_ridges2(alpha=0.5, 
                                  quantile_lines = TRUE, 
                                  quantiles = c(0.1, 0.5, 0.9), 
                                  vline_size = 0.5,
                                  vline_color = "black",
                                  scale = 1) + 
             xlab("Probability of voting yes") +
             xlim(0,1) +
             theme_bw() + 
             
             theme(plot.title = element_text(size = 50, face = "bold"),
                   axis.text=element_text(size=35),
                   axis.title=element_text(size=35, face="bold"),
                   legend.text = element_text(size=35),
                   legend.title = element_text(size=35, face="bold")) +
             
             ylab("Age group") +
             ggtitle(paste0("Sample size=", sample_size,
                            ", Num. age groups=", age_grouping_multiplier,
                            ", Response prob. of old=", p)) +
             geom_point(aes(x = mrp, y = Type), 
                        size=5,
                        colour = "black", 
                        show.legend = F,
                        data = points_df_final,
                        inherit.aes = F) +
             geom_point(aes(x = empirical_prob, y = Type), 
                        size=5,
                        colour = "#dd1c77", 
                        show.legend = F,
                        data = points_df_final,
                        inherit.aes = F) +
             geom_text(data = points_df_final,
                       aes(x = one, y = Type, label = n),
                       inherit.aes = F,
                       position=position_nudge(x = -0.02, y=0.35), colour="black", size=20)
      )
      
      dev.off()
    
    
    
    # store stuff if runs = 1 and we set store_stanobjects==TRUE
    if (store_stanobjects==TRUE) {
      postpred_sim_list_baseline[[toString(p * 10)]] = postpred_sim_baseline # store posterior linear predictor
      stanobjects_list_baseline[[toString(p * 10)]] = fit_baseline # store fitted stan object
      
      postpred_sim_list_ar[[toString(p * 10)]] = postpred_sim_ar # store posterior linear predictor
      stanobjects_list_ar[[toString(p * 10)]] = fit_ar # store fitted stan object
      
      postpred_sim_list_rw[[toString(p * 10)]] = postpred_sim_rw # store posterior linear predictor
      stanobjects_list_rw[[toString(p * 10)]] = fit_rw # store fitted stan object
    }
    
    # this gets the average preference of population based on poststratification on the sample
    
      # poststrat_sim is the model-based estimate. get it for all three models 
      poststrat_sim_baseline = postpred_sim_baseline %*% poststrat_final_reduced$s_N/sum(poststrat_final_reduced$s_N)
      model_popn_pref_baseline = c(round(mean(poststrat_sim_baseline),4), 
                                   round(sd(poststrat_sim_baseline),4))
      
      poststrat_sim_ar = postpred_sim_ar %*% poststrat_final_reduced$s_N/sum(poststrat_final_reduced$s_N)
      model_popn_pref_ar = c(round(mean(poststrat_sim_ar),4), 
                             round(sd(poststrat_sim_ar),4))
      
      poststrat_sim_rw = postpred_sim_rw %*% poststrat_final_reduced$s_N/sum(poststrat_final_reduced$s_N)
      model_popn_pref_rw = c(round(mean(poststrat_sim_rw),4), 
                             round(sd(poststrat_sim_rw),4))
    
    # ---
    
      print(paste0("population mean maximal poststrat. matrix : ",
                   sum(poststrat_final$true_pref * (poststrat_final$N)/sum(poststrat_final$N))))
      print(paste0("population mean reduced poststrat. matrix : ",
                   sum(poststrat_final_reduced$true_pref_grouped_final * (poststrat_final_reduced$s_N)/sum(poststrat_final_reduced$s_N))))
      
      print(paste0("poststrat posterior mean baseline: ",
                   model_popn_pref_baseline[1]))
      print(paste0("poststrat posterior mean ar: ",
                   model_popn_pref_ar[1]))
      print(paste0("poststrat posterior mean rw: ",
                   model_popn_pref_rw[1]))
      
      print("MRP subestimates :")
      print(mrp_subestimates)
      
      for (g in 1:age_grouping_multiplier) {
        print(paste("Posterior means for age group", g))
        print(colMeans(mrp_subsamples_allmodels[[toString(g)]]))
      }
      
      
      sample_pref = round(mean(sample_final$pref), 4)
      print(paste0("sample mean : ",
                   sample_pref))
      
      print("Printing warnings ...")
      warnings()
      
    model_popn_pref_mat_baseline[k] = model_popn_pref_baseline[1]
    model_popn_pref_sd_mat_baseline[k] = model_popn_pref_baseline[2]
    
    model_popn_pref_mat_ar[k] = model_popn_pref_ar[1]
    model_popn_pref_sd_mat_ar[k] = model_popn_pref_ar[2]
    
    model_popn_pref_mat_rw[k] = model_popn_pref_rw[1]
    model_popn_pref_sd_mat_rw[k] = model_popn_pref_rw[2]
    
    sample_mat[k] = round(mean(sample_final$pref), 4)
    sample_sd_mat[k] = round(sd(sample_final$pref), 4)
    
    sample_cell_estimates_baseline[k, ] = colMeans(postpred_sim_baseline)
    sample_cell_estimates_ar[k, ] = colMeans(postpred_sim_ar)
    sample_cell_estimates_rw[k, ] = colMeans(postpred_sim_rw)
    
    sample_cell_estimates_median_baseline[k, ] = colMedians(postpred_sim_baseline)
    sample_cell_estimates_median_ar[k, ] = colMedians(postpred_sim_ar)
    sample_cell_estimates_median_rw[k, ] = colMedians(postpred_sim_rw)
    
    # quantile width for baseline 
    sample_cell_estimates_width_baseline[k, ] = apply(X = postpred_sim_baseline,
                                                      MARGIN = 2,
                                                      FUN = quantile, 
                                                      probs = upperquantile) -
      apply(X = postpred_sim_baseline, 
            MARGIN = 2,
            FUN = quantile,
            probs = lowerquantile)
    
    # quantile width for ar
    sample_cell_estimates_width_ar[k, ] = apply(X = postpred_sim_ar,
                                                      MARGIN = 2,
                                                      FUN = quantile, 
                                                      probs = upperquantile) -
      apply(X = postpred_sim_ar, 
            MARGIN = 2,
            FUN = quantile,
            probs = lowerquantile)
    
    # quantile width for rw
    sample_cell_estimates_width_rw[k, ] = apply(X = postpred_sim_rw,
                                                      MARGIN = 2,
                                                      FUN = quantile, 
                                                      probs = upperquantile) -
      apply(X = postpred_sim_rw, 
            MARGIN = 2,
            FUN = quantile,
            probs = lowerquantile)
      
    # k MUST be an even number for the save to work
    # remove stanfit objects so we don't waste storage space
    rm(fit_baseline)
    rm(fit_ar)
    rm(fit_rw)
    rm(fit_samples_baseline)
    rm(fit_samples_ar)
    rm(fit_samples_rw)
    
    # save on every 10 runs. MAKE SURE K >= 10 an even number
    if ((k %% 10)==0) {    
      save.image(paste0(runs,
                        "_",
                        age_grouping_multiplier,
                        "_",
                        income_multiplier,
                        "_",
                        sample_size,
                        "_",
                        response_tag,
                        "_oldest_",
                        p * 10,
                        ".RData"))
    }
    
  }
  
  num_divs_counter_mat_baseline[p_counter,] = num_divs_baseline
  num_divs_counter_mat_ar[p_counter,] = num_divs_ar
  num_divs_counter_mat_rw[p_counter,] = num_divs_rw
  p_counter = p_counter + 1
  
}

print("Printing warnings ...")
print(warnings())

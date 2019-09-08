rm(list=ls())

#load into the ACS and make PS table.
options(bitmapType="cairo")

library(dplyr)
library(readr)
library(data.table)
library(rstan)
library(arm)
library(ggplot2)
library(ggridges)
library(reshape2)
library(stringr)
library(tidyr)

# ------------------------------------------------------------------------


# this is the maximal poststratification matrix

age_grouping_multiplier_vector = c(12, 48, 72)

iterations = 2000
num_chains = 4
mtreedepth = 15
ad = 0.99

use_state = TRUE # if this is true we use state as response

m_name_bl = "baseline_model_annenberg_centered" # name for baseline model
m_name_ar = "ar_model_annenberg_centered" # name for autoregressive model
m_name_rw = "rw_model_annenberg_centered" # name for random walk model

png_width = 2400
png_height = 2400

#age_cat_posteriorsummarystats_lineplot_list = list() # contains a list of ggplots
age_cat_posteriorsummarystats_list = list() # contains list of dataframes
empirical_pref_list = list() # contains a list of empirical preferences for the above list
counter = 1

covariates_list = c("CEc01_c", # FAVOR SAME-SEX MARRIAGE - THE RESPONSE
                    "WA01_c", # SEX
                    "WA02_c", # AGE
                    "WA03_c", # EDUCATION
                    "WC03_c", # RACE
                    "WA04_c", # INCOME QUESTION 1
                    "WA05_c", # INCOME QUESTION 2
                    "WB03_c", # LINE OF WORK
                    "WFc01_c", # STATE OF RESIDENCY
                    "WC06_c", # YEARS LIVED IN THE US
                    "WFa01_c", # ADULTS IN HOUSEHOLD
                    "WC01_c", # Are you of Hispanic or Latino origin or descent?
                    "WHb01_c") # SEXUAL ORIENTATION

posteriorsamples_staterandomeffect = list() # contains posterior samples for all 51 state random effects, for each age-category (12, 48, 72)
posteriorsamples_state = list() # contains posterior linear predictors for all 51 states, for each age-category (12, 48, 72)
posteriorsamples_overall = list() # contains posterior samples for the majority vote, for each age-category (12, 48, 72)

for (age_grouping_multiplier in age_grouping_multiplier_vector) {
  
  acs_ps = data.frame(readRDS("~/Desktop/annenbergdata/fiveyearacs/acs_ps.rds"))
  # "AGEP"        "SEX"         "race_x"      "education_x" "state_x"     "income_x"    "N" 
  colnames(acs_ps) = c("age", "sex", "race", "education", "state", "income", "N")

  # ------------------------------------------------------------------------
  
  # get mrp estimates based on reduced poststrat. matrix. do for all 3 models, all 2 age_cats
  
  acs_ps$age_cat = cut(acs_ps$age, age_grouping_multiplier)
  levels_age_cat = levels(acs_ps$age_cat)
  levels(acs_ps$age_cat) = 1:age_grouping_multiplier
  
  # --------------------------------------------------------------------------------------------
  # get reduced poststratification matrix ------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  #acs_ps$N[is.na(acs_ps$N)] = 0
  
  acs_reduced = data.frame(acs_ps %>% 
                             group_by(age_cat, sex, race, education, state, income) %>% 
                             summarise(sum(N)))
  
  acs_reduced$sex = factor(acs_reduced$sex)
  acs_reduced$race = factor(acs_reduced$race)
  acs_reduced$education = factor(acs_reduced$education)
  acs_reduced$state = factor(acs_reduced$state)
  acs_reduced$income = factor(acs_reduced$income)
  
  
  
  # format reduced poststratification matrix
  
  
  
  # get levels coming from ACS data
  acs_sex_levels = levels(acs_reduced$sex)
  acs_race_levels = levels(acs_reduced$race)
  acs_education_levels = levels(acs_reduced$education)
  acs_state_levels = levels(acs_reduced$state)
  acs_income_levels = levels(acs_reduced$income)
  
  # change ACS levels to the levels that I specified
  levels(acs_reduced$race) = c(1, 2, 3, 4, 5, 6) #  americanindian=1, asian=2, black=3, hisp=4, other=5, white=6
  levels(acs_reduced$sex) = c(0, 1) # men==0, women==1
  levels(acs_reduced$education) = c(1, 2, 3, 4, 5, 6) # fouryeardegree=1, fouryeardegreeplus=2, highschool=3, nohighschool=4, somecollege=5, twoyeardegree=6
  levels(acs_reduced$state) = 1:51
  
  # convert columns of acs_reduced to numbers
  acs_reduced$age_cat = as.numeric(as.character(acs_reduced$age_cat))
  acs_reduced$sex = as.numeric(as.character(acs_reduced$sex))
  acs_reduced$education = as.numeric(as.character(acs_reduced$education))
  acs_reduced$race = as.numeric(as.character(acs_reduced$race))
  acs_reduced$state = as.numeric(as.character(acs_reduced$state))
  acs_reduced$income = as.numeric(as.character(acs_reduced$income))
  
  N_sub = rep(0, age_grouping_multiplier) # the N for each age subgroup
  for (g in 1:age_grouping_multiplier) {
    N_sub[g] = sum(acs_reduced[which(acs_reduced$age_cat==g),]$sum.N.)
  }
  
  N_state = (acs_reduced %>% group_by(state) %>% summarize(N=sum(sum.N.)) %>% arrange(state))$N
  # --------------------------------------------------------------------------------------------
  # load fitted stan model for baseline --------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  
  fit_realdata = readRDS(paste0(m_name_bl,
              "_",
              age_grouping_multiplier,
              ".rds"))
  
  print(paste("Number of divergence transitions:",
              sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[1]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[2]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[3]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[4]][,'divergent__'])
  ))
  
  # extract samples
  fit_samples_realdata_bl = as.data.frame(rstan::extract(fit_realdata, permuted=FALSE))
  
  # remove the stanfit object because it's massive .. in the GBs
  rm(fit_realdata)
  
  # get postpred_sim. this is dim 4000 x (#poststrat cells)
  
  # extract posterior samples for intercept
  intercept_samples_bl = tidyr::gather(fit_samples_realdata_bl[,grepl("intercept", 
                                                                names(fit_samples_realdata_bl),
                                                                fixed=TRUE)])$value # get post. samples for intercept
  
  # 1. get random effects' posterior samples for age_cat
  U_age_samples_bl = matrix(0, iterations * num_chains/2, age_grouping_multiplier)
  
  for (j in 1:age_grouping_multiplier) {
  	U_age_samples_bl[,j] = tidyr::gather(fit_samples_realdata_bl[,grepl(paste0("U_age_transformed[", j, "]"), 
  		names(fit_samples_realdata_bl),fixed=TRUE)])$value
  }
  
  # 2. get random effects' posterior samples for race
  U_race_samples_bl = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$race)))
  
  for (j in 1:length(unique(acs_reduced$race))) {
  	U_race_samples_bl[,j] = tidyr::gather(fit_samples_realdata_bl[,grepl(paste0("U_race_transformed[", j, "]"), 
  		names(fit_samples_realdata_bl),fixed=TRUE)])$value
  }
  
  # 3. extract posterior samples for beta_gender coefficient
  beta_gender_samples_bl = tidyr::gather(fit_samples_realdata_bl[,grepl("beta_sex", 
  		names(fit_samples_realdata_bl),fixed=TRUE)])$value
  
  # 4. extract posterior samples for income
  U_income_samples_bl = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$income)))
  
  for (j in 1:length(unique(acs_reduced$income))) {
    U_income_samples_bl[,j] = tidyr::gather(fit_samples_realdata_bl[,grepl(paste0("U_income_transformed[", j, "]"), 
                                                                           names(fit_samples_realdata_bl),fixed=TRUE)])$value
  }
  
  # 5. extract posterior samples for state
  U_state_samples_bl = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$state)))
  
  for (j in 1:length(unique(acs_reduced$state))) {
    U_state_samples_bl[,j] = tidyr::gather(fit_samples_realdata_bl[,grepl(paste0("U_state_transformed[", j, "]"), 
                                                                          names(fit_samples_realdata_bl),fixed=TRUE)])$value   
  }
  
  # 6. extract posterior samples for education
  U_education_samples_bl = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$education)))
  
  for (j in 1:length(unique(acs_reduced$education))) {
    U_education_samples_bl[,j] = tidyr::gather(fit_samples_realdata_bl[,grepl(paste0("U_education_transformed[", j, "]"), 
                                                                              names(fit_samples_realdata_bl),fixed=TRUE)])$value   
  }
  
  rm(fit_samples_realdata_bl)
  gc()
  

  
  mrp_samples_bl = list() # this list stores posterior samples for age
  for (g in 1:age_grouping_multiplier) {
    mrp_samples_bl[[toString(g)]] = rep(0, dim(U_age_samples_bl)[1])
    
    for (i in as.numeric(rownames(acs_reduced[which(acs_reduced$age_cat==g),]))) {
      mrp_samples_bl[[toString(g)]] = mrp_samples_bl[[toString(g)]] + 
        (invlogit(intercept_samples_bl +
                    U_age_samples_bl[, as.numeric(acs_reduced[i, c("age_cat")])] +
                    U_race_samples_bl[, as.numeric(acs_reduced[i, c("race")])] + 
                    U_income_samples_bl[, as.numeric(acs_reduced[i, c("income")])] +
                    U_state_samples_bl[, as.numeric(acs_reduced[i, c("state")])] +
                    U_education_samples_bl[, as.numeric(acs_reduced[i, c("education")])] +
                    beta_gender_samples_bl * acs_reduced[i, c("sex")]
        ) * acs_reduced$sum.N.[i]/N_sub[g])  
    }
    
  }
  # ---
  
  mrp_subsamples_allmodels = list() # this list stores matrices of posterior estimates for subgroups, where each matrix contains the three models
  for (g in 1:age_grouping_multiplier) {
    mrp_subsamples_allmodels[[toString(g)]] = mrp_samples_bl[[toString(g)]]
  }
  
  #mrp_estimates_bl = postpred_sim_bl %*% acs_reduced$sum.N./sum(acs_reduced$sum.N.)
  
  # ---
  mrp_estimates_bl = rep(0, dim(U_age_samples_bl)[1])
  
  for (i in 1:dim(acs_reduced)[1]) {
    mrp_estimates_bl = mrp_estimates_bl + 
      (invlogit(intercept_samples_bl +
                  U_age_samples_bl[, as.numeric(acs_reduced[i, c("age_cat")])] +
                  U_race_samples_bl[, as.numeric(acs_reduced[i, c("race")])] + 
                  U_income_samples_bl[, as.numeric(acs_reduced[i, c("income")])] +
                  U_state_samples_bl[, as.numeric(acs_reduced[i, c("state")])] +
                  U_education_samples_bl[, as.numeric(acs_reduced[i, c("education")])] +
                  beta_gender_samples_bl * acs_reduced[i, c("sex")]
      ) * acs_reduced$sum.N.[i])  
  }
  
  mrp_estimates_bl = mrp_estimates_bl/sum(acs_reduced$sum.N.)
  # ---
  
  print(paste("mean:",
              mean(mrp_estimates_bl)))
  print(paste("sd",
              sd(mrp_estimates_bl)))
  
  # -----
  posteriorsamples_overall[[paste(age_grouping_multiplier, "bl")]] = mrp_estimates_bl # store posterior samples for majorityvote
  
  pld_state = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$state)))
  
  for (i in 1:length(unique(acs_reduced$state))) {
    
    for (j in as.numeric(rownames(acs_reduced[which(acs_reduced$state==i),]))) {
      pld_state[,i] = pld_state[,i] + (invlogit(intercept_samples_bl +
                                                  U_age_samples_bl[, as.numeric(acs_reduced[j, c("age_cat")])] +
                                                  U_race_samples_bl[, as.numeric(acs_reduced[j, c("race")])] + 
                                                  U_income_samples_bl[, as.numeric(acs_reduced[j, c("income")])] +
                                                  U_state_samples_bl[, as.numeric(acs_reduced[j, c("state")])] +
                                                  U_education_samples_bl[, as.numeric(acs_reduced[j, c("education")])] +
                                                  beta_gender_samples_bl * acs_reduced[j, c("sex")]
      ) * acs_reduced$sum.N.[j]/N_state[i])
    }
      
  }
  
  posteriorsamples_state[[paste(age_grouping_multiplier, "bl")]] = pld_state
  
  posteriorsamples_staterandomeffect[[paste(age_grouping_multiplier, "bl")]] = U_state_samples_bl # store posterior samples for each state
  # -----
  
  #rm(postpred_sim_bl) # remove this because its 6 gb
  gc() # IMPORTANT
  
  # --------------------------------------------------------------------------------------------
  # load fitted stan model for autoregressive --------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  fit_realdata = readRDS(paste0(m_name_ar,
              "_",
              age_grouping_multiplier,
              ".rds"))
  
  print(paste("Number of divergence transitions:",
              sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[1]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[2]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[3]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[4]][,'divergent__'])
  ))
  
  # extract samples
  fit_samples_realdata_ar = as.data.frame(rstan::extract(fit_realdata, permuted=FALSE))
  
  # remove the stanfit object because it's massive .. in the GBs
  rm(fit_realdata)
  
  # get postpred_sim. this is dim 4000 x (#poststrat cells)
  
  # extract posterior samples for intercept
  intercept_samples_ar = tidyr::gather(fit_samples_realdata_ar[,grepl("intercept", 
                                                                      names(fit_samples_realdata_ar),
                                                                      fixed=TRUE)])$value # get post. samples for intercept
  
  # 1. get random effects' posterior samples for age_cat
  U_age_samples_ar = matrix(0, iterations * num_chains/2, age_grouping_multiplier)
  
  for (j in 1:age_grouping_multiplier) {
    U_age_samples_ar[,j] = tidyr::gather(fit_samples_realdata_ar[,grepl(paste0("U_age_transformed[", j, "]"), 
                                                                        names(fit_samples_realdata_ar),fixed=TRUE)])$value
  }
  
  # 2. get random effects' posterior samples for race
  U_race_samples_ar = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$race)))
  
  for (j in 1:length(unique(acs_reduced$race))) {
    U_race_samples_ar[,j] = tidyr::gather(fit_samples_realdata_ar[,grepl(paste0("U_race_transformed[", j, "]"), 
                                                                         names(fit_samples_realdata_ar),fixed=TRUE)])$value
  }
  
  # 3. extract posterior samples for beta_gender coefficient
  beta_gender_samples_ar = tidyr::gather(fit_samples_realdata_ar[,grepl("beta_sex", 
                                                                        names(fit_samples_realdata_ar),fixed=TRUE)])$value
  
  # 4. extract posterior samples for income
  U_income_samples_ar = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$income)))
  
  for (j in 1:length(unique(acs_reduced$income))) {
    U_income_samples_ar[,j] = tidyr::gather(fit_samples_realdata_ar[,grepl(paste0("U_income_transformed[", j, "]"), 
                                                                           names(fit_samples_realdata_ar),fixed=TRUE)])$value
  }
  
  # 5. extract posterior samples for state
  U_state_samples_ar = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$state)))
  
  for (j in 1:length(unique(acs_reduced$state))) {
    U_state_samples_ar[,j] = tidyr::gather(fit_samples_realdata_ar[,grepl(paste0("U_state_transformed[", j, "]"), 
                                                                          names(fit_samples_realdata_ar),fixed=TRUE)])$value   
  }
  
  # 6. extract posterior samples for education
  U_education_samples_ar = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$education)))
  
  for (j in 1:length(unique(acs_reduced$education))) {
    U_education_samples_ar[,j] = tidyr::gather(fit_samples_realdata_ar[,grepl(paste0("U_education_transformed[", j, "]"), 
                                                                              names(fit_samples_realdata_ar),fixed=TRUE)])$value   
  }
  
  rm(fit_samples_realdata_ar)
  gc()
  
  
  mrp_samples_ar = list() # this list stores posterior samples for age
  for (g in 1:age_grouping_multiplier) {
    mrp_samples_ar[[toString(g)]] = rep(0, dim(U_age_samples_ar)[1])
    
    for (i in as.numeric(rownames(acs_reduced[which(acs_reduced$age_cat==g),]))) {
      mrp_samples_ar[[toString(g)]] = mrp_samples_ar[[toString(g)]] + 
        (invlogit(intercept_samples_ar +
                    U_age_samples_ar[, as.numeric(acs_reduced[i, c("age_cat")])] +
                    U_race_samples_ar[, as.numeric(acs_reduced[i, c("race")])] + 
                    U_income_samples_ar[, as.numeric(acs_reduced[i, c("income")])] +
                    U_state_samples_ar[, as.numeric(acs_reduced[i, c("state")])] +
                    U_education_samples_ar[, as.numeric(acs_reduced[i, c("education")])] +
                    beta_gender_samples_ar * acs_reduced[i, c("sex")]
        ) * acs_reduced$sum.N.[i]/N_sub[g])  
    }
    
  }
  
  for (g in 1:age_grouping_multiplier) {
    mrp_subsamples_allmodels[[toString(g)]] = cbind(mrp_subsamples_allmodels[[toString(g)]], 
                                                    mrp_samples_ar[[toString(g)]])
  }
  
  #mrp_estimates_ar = postpred_sim_ar %*% acs_reduced$sum.N./sum(acs_reduced$sum.N.)
  
  mrp_estimates_ar = rep(0, dim(U_age_samples_ar)[1])
  
  for (i in 1:dim(acs_reduced)[1]) {
    mrp_estimates_ar = mrp_estimates_ar + 
      (invlogit(intercept_samples_ar +
                  U_age_samples_ar[, as.numeric(acs_reduced[i, c("age_cat")])] +
                  U_race_samples_ar[, as.numeric(acs_reduced[i, c("race")])] + 
                  U_income_samples_ar[, as.numeric(acs_reduced[i, c("income")])] +
                  U_state_samples_ar[, as.numeric(acs_reduced[i, c("state")])] +
                  U_education_samples_ar[, as.numeric(acs_reduced[i, c("education")])] +
                  beta_gender_samples_ar * acs_reduced[i, c("sex")]
      ) * acs_reduced$sum.N.[i])  
  }
  
  mrp_estimates_ar = mrp_estimates_ar/sum(acs_reduced$sum.N.)
  
  print(paste("mean:",
              mean(mrp_estimates_ar)))
  print(paste("sd",
              sd(mrp_estimates_ar)))
  
  # -----
  posteriorsamples_overall[[paste(age_grouping_multiplier, "ar")]] = mrp_estimates_ar # store posterior samples for majorityvote
  
  pld_state = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$state)))
  
  for (i in 1:length(unique(acs_reduced$state))) {
    
    for (j in as.numeric(rownames(acs_reduced[which(acs_reduced$state==i),]))) {
      pld_state[,i] = pld_state[,i] + (invlogit(intercept_samples_ar +
                                                  U_age_samples_ar[, as.numeric(acs_reduced[j, c("age_cat")])] +
                                                  U_race_samples_ar[, as.numeric(acs_reduced[j, c("race")])] + 
                                                  U_income_samples_ar[, as.numeric(acs_reduced[j, c("income")])] +
                                                  U_state_samples_ar[, as.numeric(acs_reduced[j, c("state")])] +
                                                  U_education_samples_ar[, as.numeric(acs_reduced[j, c("education")])] +
                                                  beta_gender_samples_ar * acs_reduced[j, c("sex")]
      ) * acs_reduced$sum.N.[j]/N_state[i])
    }
    
  }
  
  posteriorsamples_state[[paste(age_grouping_multiplier, "ar")]] = pld_state
  
  posteriorsamples_staterandomeffect[[paste(age_grouping_multiplier, "ar")]] = U_state_samples_ar # store posterior samples for each state
  # -----
  
  #rm(postpred_sim_ar)
  gc()
  
  # --------------------------------------------------------------------------------------------
  # load fitted stan model for random walk -----------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  fit_realdata = readRDS(paste0(m_name_rw,
              "_",
              age_grouping_multiplier,
              ".rds"))
  
  print(paste("Number of divergence transitions:",
              sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[1]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[2]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[3]][,'divergent__']) +
                sum(get_sampler_params(fit_realdata, inc_warmup=FALSE)[[4]][,'divergent__'])
  ))
  
  # extract samples
  fit_samples_realdata_rw = as.data.frame(rstan::extract(fit_realdata, permuted=FALSE))
  
  # remove the stanfit object because it's massive .. in the GBs
  rm(fit_realdata)
  
  # get postpred_sim. this is dim 4000 x (#poststrat cells)
  
  # extract posterior samples for intercept
  intercept_samples_rw = tidyr::gather(fit_samples_realdata_rw[,grepl("intercept", 
                                                                      names(fit_samples_realdata_rw),
                                                                      fixed=TRUE)])$value # get post. samples for intercept
  
  # 1. get random effects' posterior samples for age_cat
  U_age_samples_rw = matrix(0, iterations * num_chains/2, age_grouping_multiplier)
  
  for (j in 1:age_grouping_multiplier) {
    U_age_samples_rw[,j] = tidyr::gather(fit_samples_realdata_rw[,grepl(paste0("U_age_transformed[", j, "]"), 
                                                                        names(fit_samples_realdata_rw),fixed=TRUE)])$value
  }
  
  # 2. get random effects' posterior samples for race
  U_race_samples_rw = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$race)))
  
  for (j in 1:length(unique(acs_reduced$race))) {
    U_race_samples_rw[,j] = tidyr::gather(fit_samples_realdata_rw[,grepl(paste0("U_race_transformed[", j, "]"), 
                                                                         names(fit_samples_realdata_rw),fixed=TRUE)])$value
  }
  
  # 3. extract posterior samples for beta_gender coefficient
  beta_gender_samples_rw = tidyr::gather(fit_samples_realdata_rw[,grepl("beta_sex", 
                                                                        names(fit_samples_realdata_rw),fixed=TRUE)])$value
  
  # 4. extract posterior samples for income
  U_income_samples_rw = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$income)))
  
  for (j in 1:length(unique(acs_reduced$income))) {
    U_income_samples_rw[,j] = tidyr::gather(fit_samples_realdata_rw[,grepl(paste0("U_income_transformed[", j, "]"), 
                                                                           names(fit_samples_realdata_rw),fixed=TRUE)])$value
  }
  
  # 5. extract posterior samples for state
  U_state_samples_rw = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$state)))
  
  for (j in 1:length(unique(acs_reduced$state))) {
    U_state_samples_rw[,j] = tidyr::gather(fit_samples_realdata_rw[,grepl(paste0("U_state_transformed[", j, "]"), 
                                                                          names(fit_samples_realdata_rw),fixed=TRUE)])$value   
  }
  
  # 6. extract posterior samples for education
  U_education_samples_rw = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$education)))
  
  for (j in 1:length(unique(acs_reduced$education))) {
    U_education_samples_rw[,j] = tidyr::gather(fit_samples_realdata_rw[,grepl(paste0("U_education_transformed[", j, "]"), 
                                                                              names(fit_samples_realdata_rw),fixed=TRUE)])$value   
  }
  
  rm(fit_samples_realdata_rw)
  gc()
  
  
  mrp_samples_rw = list() # this list stores posterior samples for age
  for (g in 1:age_grouping_multiplier) {
    mrp_samples_rw[[toString(g)]] = rep(0, dim(U_age_samples_rw)[1])
    
    for (i in as.numeric(rownames(acs_reduced[which(acs_reduced$age_cat==g),]))) {
      mrp_samples_rw[[toString(g)]] = mrp_samples_rw[[toString(g)]] + 
        (invlogit(intercept_samples_rw +
                    U_age_samples_rw[, as.numeric(acs_reduced[i, c("age_cat")])] +
                    U_race_samples_rw[, as.numeric(acs_reduced[i, c("race")])] + 
                    U_income_samples_rw[, as.numeric(acs_reduced[i, c("income")])] +
                    U_state_samples_rw[, as.numeric(acs_reduced[i, c("state")])] +
                    U_education_samples_rw[, as.numeric(acs_reduced[i, c("education")])] +
                    beta_gender_samples_rw * acs_reduced[i, c("sex")]
        ) * acs_reduced$sum.N.[i]/N_sub[g])  
    }
    
  }
  
  for (g in 1:age_grouping_multiplier) {
    mrp_subsamples_allmodels[[toString(g)]] = cbind(mrp_subsamples_allmodels[[toString(g)]], 
                                                    mrp_samples_rw[[toString(g)]])
  }
  
  #mrp_estimates_rw = postpred_sim_rw %*% acs_reduced$sum.N./sum(acs_reduced$sum.N.)
  
  mrp_estimates_rw = rep(0, dim(U_age_samples_rw)[1])
  
  for (i in 1:dim(acs_reduced)[1]) {
    mrp_estimates_rw = mrp_estimates_rw + 
      (invlogit(intercept_samples_rw +
                  U_age_samples_rw[, as.numeric(acs_reduced[i, c("age_cat")])] +
                  U_race_samples_rw[, as.numeric(acs_reduced[i, c("race")])] + 
                  U_income_samples_rw[, as.numeric(acs_reduced[i, c("income")])] +
                  U_state_samples_rw[, as.numeric(acs_reduced[i, c("state")])] +
                  U_education_samples_rw[, as.numeric(acs_reduced[i, c("education")])] +
                  beta_gender_samples_rw * acs_reduced[i, c("sex")]
      ) * acs_reduced$sum.N.[i])  
  }
  
  mrp_estimates_rw = mrp_estimates_rw/sum(acs_reduced$sum.N.)
  
  print(paste("mean:",
              mean(mrp_estimates_rw)))
  print(paste("sd",
              sd(mrp_estimates_rw)))
  
  # -----
  posteriorsamples_overall[[paste(age_grouping_multiplier, "rw")]] = mrp_estimates_rw # store posterior samples for majorityvote
  
  pld_state = matrix(0, iterations * num_chains/2, length(unique(acs_reduced$state)))
  
  for (i in 1:length(unique(acs_reduced$state))) {
    
    for (j in as.numeric(rownames(acs_reduced[which(acs_reduced$state==i),]))) {
      pld_state[,i] = pld_state[,i] + (invlogit(intercept_samples_rw +
                                                  U_age_samples_rw[, as.numeric(acs_reduced[j, c("age_cat")])] +
                                                  U_race_samples_rw[, as.numeric(acs_reduced[j, c("race")])] + 
                                                  U_income_samples_rw[, as.numeric(acs_reduced[j, c("income")])] +
                                                  U_state_samples_rw[, as.numeric(acs_reduced[j, c("state")])] +
                                                  U_education_samples_rw[, as.numeric(acs_reduced[j, c("education")])] +
                                                  beta_gender_samples_rw * acs_reduced[j, c("sex")]
      ) * acs_reduced$sum.N.[j]/N_state[i])
    }
    
  }
  
  posteriorsamples_state[[paste(age_grouping_multiplier, "rw")]] = pld_state
  
  posteriorsamples_staterandomeffect[[paste(age_grouping_multiplier, "rw")]] = U_state_samples_rw # store posterior samples for each state
  # -----
  
  #rm(postpred_sim_rw)
  gc()
  
  
  # --------------------------------------------------------------------------------------------
  # Do plotting for all three models now -------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
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
  
  
  # --------------------------------------------------------------------------------------------
  # Ridge plots --------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  

  
  survey_dat = readr::read_tsv("annenbergphone2008.txt", col_names = TRUE)
  
  survey_dat_selected = survey_dat %>% 
    dplyr::select(covariates_list) %>% 
    filter(WA02_c < 100, # Age less than 100
           CEc01_c %in% c(1,3), # 1 is vote yes for gay marriage, 3 is no
    ) 
  
  survey_dat_selected$CEc01_c[survey_dat_selected$CEc01_c==3] = 0 # change 3 to 0. 0 is no to gay marriage, 1 is yes.
  
  # Preprocess race
  survey_dat_selected = survey_dat_selected %>% filter(WC03_c < 10) # remove 998 and 999 from race
  
  survey_dat_selected$race_x = ifelse(survey_dat_selected$WC03_c==1 & survey_dat_selected$WC01_c==2, "white",
                                      ifelse(survey_dat_selected$WC03_c==2 & survey_dat_selected$WC01_c==2, "black",
                                             ifelse(survey_dat_selected$WC03_c==3 & survey_dat_selected$WC01_c==2, "asian",
                                                    ifelse(survey_dat_selected$WC03_c==4 & survey_dat_selected$WC01_c==2, "americanindian",
                                                           ifelse(survey_dat_selected$WC03_c==5 & survey_dat_selected$WC01_c==1, "hisp",
                                                                  "other"))
                                             )
                                      )
  )
  
  # Preprocess education
  survey_dat_selected = survey_dat_selected %>% filter(WA03_c < 100) # remove 998 and 999
  
  survey_dat_selected$education_x = ifelse(survey_dat_selected$WA03_c %in% c(1,2), "nohighschool",
                                           ifelse(survey_dat_selected$WA03_c %in% c(3,4), "highschool",
                                                  ifelse(survey_dat_selected$WA03_c==5, "somecollege",
                                                         ifelse(survey_dat_selected$WA03_c==6, "twoyeardegree",
                                                                ifelse(survey_dat_selected$WA03_c==7, "fouryeardegree",
                                                                       "fouryeardegreeplus")
                                                         )
                                                  )
                                           )
  )
  
  
  # Preprocess state
  survey_dat_selected$state_x = case_when(survey_dat_selected$WFc01_c == 1 ~ "alabama",
                                          survey_dat_selected$WFc01_c == 4 ~ "arizona",
                                          survey_dat_selected$WFc01_c == 5 ~ "arkansas",
                                          survey_dat_selected$WFc01_c == 6 ~ "california",
                                          survey_dat_selected$WFc01_c == 8 ~ "colorado",
                                          survey_dat_selected$WFc01_c == 9 ~ "connecticut",
                                          survey_dat_selected$WFc01_c == 10 ~ "delaware",
                                          survey_dat_selected$WFc01_c == 11 ~ "columbia",
                                          survey_dat_selected$WFc01_c == 12 ~ "florida",
                                          survey_dat_selected$WFc01_c == 13 ~ "georgia",
                                          survey_dat_selected$WFc01_c == 16 ~ "idaho",
                                          survey_dat_selected$WFc01_c == 17 ~ "illinois",
                                          survey_dat_selected$WFc01_c == 18 ~ "indiana",
                                          survey_dat_selected$WFc01_c == 19 ~ "iowa",
                                          survey_dat_selected$WFc01_c == 20 ~ "kansas",
                                          survey_dat_selected$WFc01_c == 21 ~ "kentucky",
                                          survey_dat_selected$WFc01_c == 22 ~ "louisiana",
                                          survey_dat_selected$WFc01_c == 23 ~ "maine",
                                          survey_dat_selected$WFc01_c == 24 ~ "maryland",
                                          survey_dat_selected$WFc01_c == 25 ~ "massachusetts",
                                          survey_dat_selected$WFc01_c == 26 ~ "michigan",
                                          survey_dat_selected$WFc01_c == 27 ~ "minnesota",
                                          survey_dat_selected$WFc01_c == 28 ~ "mississippi",
                                          survey_dat_selected$WFc01_c == 29 ~ "missouri",
                                          survey_dat_selected$WFc01_c == 30 ~ "montana",
                                          survey_dat_selected$WFc01_c == 31 ~ "nebraska",
                                          survey_dat_selected$WFc01_c == 32 ~ "nevada",
                                          survey_dat_selected$WFc01_c == 33 ~ "newhampshire",
                                          survey_dat_selected$WFc01_c == 34 ~ "newjersey",
                                          survey_dat_selected$WFc01_c == 35 ~ "newmexico",
                                          survey_dat_selected$WFc01_c == 36 ~ "newyork",
                                          survey_dat_selected$WFc01_c == 37 ~ "northcarolina",
                                          survey_dat_selected$WFc01_c == 38 ~ "northdakota",
                                          survey_dat_selected$WFc01_c == 39 ~ "ohio",
                                          survey_dat_selected$WFc01_c == 40 ~ "oklahoma",
                                          survey_dat_selected$WFc01_c == 41 ~ "oregon",
                                          survey_dat_selected$WFc01_c == 42 ~ "pennsylvania",
                                          survey_dat_selected$WFc01_c == 44 ~ "rhodeisland",
                                          survey_dat_selected$WFc01_c == 45 ~ "southcarolina",
                                          survey_dat_selected$WFc01_c == 46 ~ "southdakota",
                                          survey_dat_selected$WFc01_c == 47 ~ "tennessee",
                                          survey_dat_selected$WFc01_c == 48 ~ "texas",
                                          survey_dat_selected$WFc01_c == 49 ~ "utah",
                                          survey_dat_selected$WFc01_c == 50 ~ "vermont",
                                          survey_dat_selected$WFc01_c == 51 ~ "virginia",
                                          survey_dat_selected$WFc01_c == 53 ~ "washington",
                                          survey_dat_selected$WFc01_c == 54 ~ "westvirginia",
                                          survey_dat_selected$WFc01_c == 55 ~ "wisconsin",
                                          survey_dat_selected$WFc01_c == 56 ~ "wyoming",
                                          TRUE ~ as.character(survey_dat_selected$WFc01_c))
  
  
  # Preprocess income
  survey_dat_selected = survey_dat_selected %>% filter(WA04_c < 10) 
  survey_dat_selected = survey_dat_selected[,c("CEc01_c", # favors gay marriage
                                               "WA01_c", # sex
                                               "WA02_c", # age
                                               "education_x", # education
                                               "race_x", # race
                                               "WA04_c", # household income
                                               "state_x")] %>% data.frame() # state
  
  survey_dat_selected = survey_dat_selected[complete.cases(survey_dat_selected),]
  
  gaymarriagedata_rm_clean_mapped = survey_dat_selected %>% filter(WA02_c <= 95 & WA02_c >= 18)
  gaymarriagedata_rm_clean_mapped$age_cat = cut(gaymarriagedata_rm_clean_mapped$WA02_c,
                                                age_grouping_multiplier)
  
  levels(gaymarriagedata_rm_clean_mapped$age_cat) = 1:age_grouping_multiplier
  
  
  empirical_pref = data.frame(gaymarriagedata_rm_clean_mapped %>% 
                                group_by(age_cat) %>% 
                                summarise(empirical_pref = sum(CEc01_c)/n(), N=n()))
  
  
  points_df = data.frame(cbind(word(levels(mrp_subsamples_melted$Type), 1),
                               levels(mrp_subsamples_melted$Type)))
  colnames(points_df) = c("age_cat", "Type")
  
  points_df_final = merge(points_df, empirical_pref, by = "age_cat")
  points_df_final = cbind(points_df_final, 1)
  colnames(points_df_final)[dim(points_df_final)[2]] = "one"
  
  png(paste0(age_grouping_multiplier, ".png"), 
      width = png_width, height = png_height)
  
  # change order of levels so that legend of ridgeplots is in proper order
  mrp_subsamples_melted$Model = factor(mrp_subsamples_melted$Model, levels = c("Random walk", "Autoregressive", "Baseline"))
  
  plot(
    ggplot(mrp_subsamples_melted, 
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
      ggtitle(paste0("Num. age groups=", age_grouping_multiplier)) +
      
      geom_point(aes(x = empirical_pref, y = Type), 
                 size=5,
                 colour = "#dd1c77", 
                 show.legend = F,
                 data = points_df_final,
                 inherit.aes = F) + 
      geom_text(data = points_df_final,
                aes(x = one, y = Type, label = N),
                inherit.aes = F,
                position=position_nudge(x = -0.02, y=0.35), colour="black", size=20) 
  )
  
  dev.off()
  
  # --------------------------------------------------------------------------------------------
  # Overall est. -------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------
  
  # plot ridgeplot for overall MRP estimate
  
  mrp_estimates_bl_df = cbind(data.frame(mrp_estimates_bl), "Baseline")
  mrp_estimates_ar_df = cbind(data.frame(mrp_estimates_ar), "Autoregressive")
  mrp_estimates_rw_df = cbind(data.frame(mrp_estimates_rw), "Random walk")
  
  colnames(mrp_estimates_bl_df) = c("Preference", "Model")
  colnames(mrp_estimates_ar_df) = c("Preference", "Model")
  colnames(mrp_estimates_rw_df) = c("Preference", "Model")
  
  overall_estimates = rbind(mrp_estimates_bl_df, mrp_estimates_ar_df, mrp_estimates_rw_df)
  
  overall_pref = data.frame(cbind(c("Baseline", "Autoregressive", "Random walk"),
                                  mean(gaymarriagedata_rm_clean_mapped$CEc01_c)))
  colnames(overall_pref) = c("Model", "empirical_pref")
  overall_pref$empirical_pref = as.numeric(as.character(overall_pref$empirical_pref))
  
  png(paste0("overall_", age_grouping_multiplier, ".png"), 
        width = 3000, height = 2400)
  
  plot(
    ggplot(overall_estimates, 
           aes(x = Preference, y = Model, fill=Model)) + 
      geom_density_ridges2(alpha=0.5, 
                           quantile_lines = TRUE, 
                           quantiles = c(0.05, 0.5, 0.95), 
                           vline_size = 0.5,
                           vline_color = "black",
                           scale = 1) + 
      xlab("\n Probability of voting yes \n") +
      theme_bw() + 
      theme(plot.title = element_text(size = 50, face = "bold"),
            axis.text=element_text(size=35),
            axis.title=element_text(size=50, face="bold",margin=200),
            legend.text = element_text(size=50),
            legend.position  = "bottom",
            legend.key.size = unit(30,"line"),
            legend.key.height = unit(7,"line"),
            legend.key = element_rect(fill = "transparent",color = "transparent"),
            legend.title = element_text(size=50, face="bold"),
            strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
            strip.background = element_rect(fill="transparent",color="transparent")) +
      ylab("\n Model \n") +
      geom_point(aes(x = empirical_pref, y = Model), 
                 size=7,
                 colour = "#dd1c77", 
                 show.legend = F,
                 data = overall_pref,
                 inherit.aes = F) +
      scale_fill_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
      guides(color = guide_legend(override.aes = list(size=20)))
  )
  
  dev.off()
  
  # geom_ribbon plots --------------------------------------------------------------------
  
  age_cat_posteriorsummarystats = mrp_subsamples_melted %>% 
    group_by(Type) %>% 
    summarise(posterior_median = median(Preference), 
    upper = quantile(Preference, 0.95), # upper band is 90% quantile for age-category Posterior
    lower = quantile(Preference, 0.05), # upper band is 10% quantile for age-category Posterior
    posterior_mean = mean(Preference),
    age_cat = as.numeric("[["(strsplit(as.character(Type), " "), 1)[1]),
    Model = "[["(strsplit(as.character(Type), " "), 1)[2]
    ) %>% data.frame()
  
  age_cat_posteriorsummarystats$Model = factor(age_cat_posteriorsummarystats$Model)
  
  # rename levels
  levels(age_cat_posteriorsummarystats$Model)[levels(age_cat_posteriorsummarystats$Model)=="AR"]  = "Autoregressive"
  levels(age_cat_posteriorsummarystats$Model)[levels(age_cat_posteriorsummarystats$Model)=="BL"]  = "Baseline"
  levels(age_cat_posteriorsummarystats$Model)[levels(age_cat_posteriorsummarystats$Model)=="RW"]  = "Random walk"
  
  # change order of levels
  age_cat_posteriorsummarystats$Model = factor(age_cat_posteriorsummarystats$Model,
                                               levels = c("Baseline", "Autoregressive", "Random walk"))
  
  empirical_pref$age_cat = as.numeric(as.character(empirical_pref$age_cat))
  
  empirical_pref_list[[toString(age_grouping_multiplier)]] = empirical_pref
  
  samplesize_df_lineplot = ggplot(gaymarriagedata_rm_clean_mapped, aes(x = WA02_c)) + 
    geom_histogram(bins = 78) + # number of ages is 80
    ggtitle("Number of samples") +
    theme_bw() +
    theme(plot.title = element_text(size = 84, face = "bold"),
          axis.text=element_text(size=59),
          axis.title=element_text(size=59, face="bold"),
          legend.text = element_text(size=84),
          legend.position  = "bottom",
          legend.key.size = unit(50,"line"),
          legend.key.height = unit(12,"line"),
          legend.key = element_rect(fill = "transparent",color = "transparent"),
          legend.title = element_text(size=84, face="bold")) +
    xlab("\n Age \n") +
    ylab("\n Number of samples \n")
  
  age_cat_posteriorsummarystats_list[[toString(age_grouping_multiplier)]] = age_cat_posteriorsummarystats
  
  age_cat_posteriorsummarystats_lineplot = ggplot(age_cat_posteriorsummarystats, 
                                                  aes(y = posterior_median, x = age_cat)) + 
    geom_line(aes(color = Model), size = 2.5) + 
    scale_colour_manual(values = c("#4575b4", "#d73027","#fdae61")) +
    geom_ribbon(aes(ymin=lower, ymax=upper, x=age_cat, fill = Model), alpha = 0.2) + 
    scale_fill_manual(values = c("#4575b4", "#d73027","#fdae61")) +
    theme_bw() +
    theme(plot.title = element_text(size = 84, face = "bold"),
          axis.text=element_text(size=59),
          axis.title=element_text(size=59, face="bold"),
          legend.text = element_text(size=84),
          legend.position  = "bottom",
          legend.key.size = unit(50,"line"),
          legend.key.height = unit(12,"line"),
          legend.key = element_rect(fill = "transparent",color = "transparent"),
          legend.title = element_text(size=84, face="bold")) +
    xlab("\n  Age category \n") +
    ylab("\n Median of Posterior with 5-95 quantile bands \n") +
    geom_point(aes(x=rep(empirical_pref$age_cat, 3), y=rep(empirical_pref$empirical_pref, 3)), 
              size = 7.5,
              color="red")
  
  
  png(paste0("age_cat_ribbon_", age_grouping_multiplier ,".png"),
      width = 5000, height = 2400)
  
  # plot(
  #   age_cat_posteriorsummarystats_lineplot
  # )
  
  plot(cowplot::plot_grid(age_cat_posteriorsummarystats_lineplot,
                     NULL,
                     samplesize_df_lineplot,
             align = "v",
             ncol=1,
             rel_heights = c(3, 0.1, 1)))
  dev.off()
  
  #print("Saving image ...")
  #save.image(paste0(age_grouping_multiplier, ".RData"))
  
  
  #rm(list=ls()) # delete everything so we don't save RData which slows down next sessions
  counter = counter + 1
}

#this part is manual
p1 = cowplot::plot_grid(ggplot(age_cat_posteriorsummarystats_list[["12"]],
                               aes(y = posterior_median, x = age_cat)) +
                          geom_line(aes(color = Model), size = 2.5) +
                          scale_colour_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                          geom_ribbon(aes(ymin=lower, ymax=upper, x=age_cat, fill = Model), alpha = 0.2) +
                          scale_fill_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                          xlab("\n Age Category \n") +
                          theme_bw() +
                          theme(plot.title = element_text(size = 84, face = "bold"),
                                axis.text=element_text(size=84),
                                axis.title=element_text(size=84, face="bold"),
                                #axis.title.x=element_blank(),
                                #axis.text.x=element_blank(),
                                #axis.ticks.x=element_blank(),
                                axis.title.y=element_blank(),
                                legend.position = "none") +
                          ylim(0, 0.7) +
                          scale_x_continuous(breaks= (1:6)*2 ) +
                          geom_point(aes(x=rep(empirical_pref_list[["12"]]$age_cat, 3),
                                         y=rep(empirical_pref_list[["12"]]$empirical_pref, 3)),
                                     size = 7.5,
                                     color="red"),

                        ggplot(age_cat_posteriorsummarystats_list[["48"]],
                               aes(y = posterior_median, x = age_cat)) +
                          geom_line(aes(color = Model), size = 2.5) +
                          scale_colour_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                          geom_ribbon(aes(ymin=lower, ymax=upper, x=age_cat, fill = Model), alpha = 0.2) +
                          scale_fill_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                          theme_bw() +
                          xlab("\n Age Category \n") +
                          theme(plot.title = element_text(size = 84, face = "bold"),
                                axis.text=element_text(size=84),
                                #axis.title.x=element_blank(),
                                #axis.text.x=element_blank(),
                                #axis.ticks.x=element_blank(),
                                axis.title=element_text(size=84, face="bold"),
                                legend.position = "none") +
                          ylab("\n Posterior 5-50-95 quantiles \n") +
                          ylim(0, 0.7) +
                          scale_x_continuous(breaks= (1:4)*12 ) +
                          geom_point(aes(x=rep(empirical_pref_list[["48"]]$age_cat, 3),
                                         y=rep(empirical_pref_list[["48"]]$empirical_pref, 3)),
                                     size = 7.5,
                                     color="red"),

                        ggplot(age_cat_posteriorsummarystats_list[["72"]],
                               aes(y = posterior_median, x = age_cat)) +
                          geom_line(aes(color = Model), size = 2.5) +
                          scale_colour_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                          geom_ribbon(aes(ymin=lower, ymax=upper, x=age_cat, fill = Model), alpha = 0.2) +
                          scale_fill_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                          theme_bw() +
                          xlab("\n Age Category \n") +
                          theme(plot.title = element_text(size = 84, face = "bold"),
                                axis.text=element_text(size=84),
                                axis.title.y=element_blank(),
                                #axis.title.x=element_blank(),
                                #axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.title=element_text(size=84, face="bold"),
                                legend.position = "none") +
                          ylim(0, 0.7) +
                          scale_x_continuous(breaks= (1:6)*12 ) +
                          geom_point(aes(x=rep(empirical_pref_list[["72"]]$age_cat, 3),
                                         y=rep(empirical_pref_list[["72"]]$empirical_pref, 3)),
                                     size = 7.5,
                                     color="red"),
                        align = "v",
                        ncol=1)


legend_p1 = cowplot::get_legend(ggplot(age_cat_posteriorsummarystats_list[["48"]],
                                       aes(y = posterior_median, x = age_cat)) +
                                  geom_line(aes(color = Model), size = 2.5) +
                                  scale_colour_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                                  geom_ribbon(aes(ymin=lower, ymax=upper, x=age_cat, fill = Model), alpha = 0.2) +
                                  scale_fill_manual(values = c("#4575b4", "#d73027","#fdae61")) +
                                  theme_bw() +
                                  theme(plot.title = element_text(size = 84, face = "bold"),
                                        axis.text=element_text(size=59),
                                        axis.title.x=element_blank(),
                                        axis.title=element_text(size=59, face="bold"),
                                        legend.text = element_text(size=84),
                                        legend.position  = "bottom",
                                        legend.key.size = unit(50,"line"),
                                        legend.key.height = unit(12,"line"),
                                        legend.key = element_rect(fill = "transparent",color = "transparent"),
                                        legend.title = element_text(size=84, face="bold")) +
                                  ylab("\n Posterior 5-50-95 quantiles \n") +
                                  ylim(0, 0.7) +
                                  scale_x_continuous(breaks= scales::pretty_breaks()) +
                                  geom_point(aes(x=rep(empirical_pref_list[["48"]]$age_cat, 3),
                                                 y=rep(empirical_pref_list[["48"]]$empirical_pref, 3)),
                                             size = 7.5,
                                             color="red"))

annenberg_ages = cbind(gaymarriagedata_rm_clean_mapped$WA02_c, "Annenberg") %>% data.frame()
annenberg_ages$X1 = as.numeric(as.character(annenberg_ages$X1))

acs_ages = cbind(acs_ps$age, "ACS") %>% data.frame()
acs_ages$X1 = as.numeric(as.character(acs_ages$X1))

ages_combined = rbind(annenberg_ages, acs_ages)
colnames(ages_combined) = c("X1", "Survey")

colnames(acs_ages) = c("X1_", "X2_")


acs_sample = sample(acs_ps$age, 
                    dim(gaymarriagedata_rm_clean_mapped)[1],
                    replace = TRUE,
                    prob = acs_ps$N)


dat = data.frame(annenberg = gaymarriagedata_rm_clean_mapped$WA02_c,
                 acs = acs_sample)

samplesize_df_lineplot = ggplot(dat) +
  geom_histogram(mapping=aes(x=annenberg,y=..count..), bins = 78, fill="gray80")+
  geom_density(mapping = aes(x=acs, y=..count..), size=5, linetype = "dashed") + 
  #ggtitle("Ages for Annenberg survey (histogram) and ages for ACS (density) ") +
  theme_bw() +
  theme(plot.title = element_text(size = 84, face = "bold"),
        axis.text=element_text(size=84),
        axis.title=element_text(size=84, face="bold"),
        legend.text = element_text(size=84),
        legend.position  = "bottom",
        legend.key.size = unit(50,"line"),
        legend.key.height = unit(12,"line"),
        legend.key = element_rect(fill = "transparent",color = "transparent"),
        legend.title = element_text(size=84, face="bold")) +
  xlab("\n Age \n") +
  ylab("\n Count \n")


p2 = cowplot::plot_grid(legend_p1,
                        NULL,
                        p1,
                        NULL,
                        samplesize_df_lineplot,
                        ncol=1,
                        rel_heights = c(0.3,
                                        0.3,
                                        4,
                                        0.3,
                                        1.5))



png(paste0("combined_agecatribbons.png"),
    width = 4800, height = 4800)

plot(p2)

dev.off()

saveRDS(posteriorsamples_staterandomeffect, "U_state.rds")
saveRDS(posteriorsamples_state, "lp_state.rds")
saveRDS(posteriorsamples_overall, "overall.rds")

print("Done script.")
  
  
  

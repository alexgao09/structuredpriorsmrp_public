set.seed(21)


rm(list=ls())
gc()

library(dplyr)
library(readr)
library(data.table)
library(sp)
library(ggplot2)
library(INLA)

use_spatial_effect = TRUE # the default is TRUE. This will define the true_pref with coef_puma, the spatial effect

poverty_poststrat = readRDS("poverty_poststrat.rds")
us_pumas = readRDS("us_pumas.rds")
ma_adj_matrix = readRDS("ma_adj_matrix.rds")
ma_adj_matrix_sparse = readRDS("ma_adj_matrix_sparse.rds")
x = readRDS("smooth_x.rds")
  
num_education = length(unique(poverty_poststrat$education)) # number of strata for education
num_race = length(unique(poverty_poststrat$race_x)) # number of strata for race
num_puma = length(unique(poverty_poststrat$PUMA.x))

total_poverty = poverty_poststrat %>% group_by(PUMA.x) %>% summarise(total = sum(N),
                                                                     total_in_poverty = sum(y)) %>% arrange(PUMA.x) %>% as.data.frame()

total_poverty$proportion_in_poverty = total_poverty$total_in_poverty/total_poverty$total


ss = "massachusetts"
ss_code = 25 # use data dictionary to get ss code
continental_states = c("MA")

r = 1:9/10
runs = 200
sample_size = 500 # 500, 1000, 2000, 4000

puma_overundersample_index = c(3301, # "Boston City--Allston, Brighton & Fenway"
                               3305,  # "Boston City--Hyde Park, Jamaica Plain, Roslindale & West Roxbury"
                               3303, # "Boston City--Dorchester & South Boston"
                               3302, # "Boston City--Back Bay, Beacon Hill, Charlestown, East Boston, Central & South End" 
                               3304, # "Boston City--Mattapan & Roxbury" 
                               3306,
                               503,
                               504,
                               505,
                               506,
                               507,
                               508,
                               3500,
                               3601,
                               3602,
                               3603,
                               3400
)

# used to calculate posterior width
lowerquantile = 0.1
upperquantile = 0.9

temp_x = cbind(x, us_pumas$PUMACE10) %>% as.data.frame()
x_ordered = temp_x %>% arrange(V2) %>% as.data.frame()
coef_puma = x_ordered$V1 # coef_puma is in alphabetical order

coef_race = (1:num_race - mean(1:num_race))*0.1

coef_education = (1:num_education - mean(1:num_education) - 2)*0.15

intercept_term = 0

# http://www.paulamoraga.com/book-geospatial/sec-inla.html
# https://mc-stan.org/users/documentation/case-studies/icar_stan.html

icar_formula = n ~ 1 + 
  f(education_f,
    model="iid",
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.1)))
    ) +
  f(race_f,
    model="iid",
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.1)))
    ) +
  f(ID,
    model = "bym2",
    graph = ma_adj_matrix_sparse,
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.1))),
    adjust.for.con.comp = TRUE
  )

iid_formula = n ~ 1 + 
  f(education_f,
    model="iid",
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.1)))
    ) +
  f(race_f,
    model="iid",
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.1)))
    ) +
  f(ID_f,
    model = "iid",
    hyper = list(prec = list(prior="pc.prec", param=c(1, 0.1)))
  )

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

us_pumas_df = as.data.frame(us_pumas)
# ---

a = as.numeric(us_pumas$PUMACE10) # converts string to numeric
b = unique(poverty_poststrat$PUMA.x)

sort(a)
sort(b)

sort(a) - sort(b) # great

# https://www.census.gov/geographies/reference-maps/2010/geo/2010-pumas/massachusetts.html
us_pumas$PUMACE10 = as.numeric(us_pumas$PUMACE10) # just run once

us_pumas$proportion_in_poverty = total_poverty[base::match(us_pumas$PUMACE10, total_poverty$PUMA.x),c("proportion_in_poverty")]

# fit simple ICAR model to MA w/ response being proportion_in_poverty
us_pumas$ID = 1:nrow(us_pumas)

# -----------------------------------------------------------------------
# simulation pipeline ---------------------------------------------------
# -----------------------------------------------------------------------

num_education = length(unique(poverty_poststrat$education)) # number of strata for education
num_race = length(unique(poverty_poststrat$race_x)) # number of strata for race
num_puma = length(unique(poverty_poststrat$PUMA.x))

poststrat_final = expand.grid(puma_area = unique(poverty_poststrat$PUMA.x),
                              race = unique(poverty_poststrat$race_x),
                              education = unique(poverty_poststrat$education))

poststrat_final_joined = left_join(poststrat_final, poverty_poststrat, by=c("puma_area"="PUMA.x",
                                                                            "race"="race_x",
                                                                            "education"="education_x")) %>% as.data.frame()

# Order of factors for poststrat_final_joined
levels_education = unique(poverty_poststrat$education)
levels_race = unique(poverty_poststrat$race_x)
levels_puma = unique(poverty_poststrat$PUMA.x)


poststrat_final_join_numeric = expand.grid(puma_area = 1:length(unique(poverty_poststrat$PUMA.x)),
                                           race = 1:length(unique(poverty_poststrat$race_x)),
                                           education = 1:length(unique(poverty_poststrat$education)))

# https://stackoverflow.com/questions/18562680/replacing-nas-with-0s-in-r-dataframe
# replace all zeros in dataframe 
poststrat_final_joined[is.na(poststrat_final_joined)] = 0

true_pref = rep(NA, num_puma * num_race * num_education) # this contains the true preference probability of a group in the population
if (use_spatial_effect==TRUE) {
  
  for (j in 1:(num_education * num_race * num_puma)) {
    true_pref[j] = arm::invlogit(intercept_term +
                                   coef_puma[poststrat_final_join_numeric[j, 1]] +
                                   coef_race[poststrat_final_join_numeric[j, 2]] +
                                   coef_education[poststrat_final_join_numeric[j, 3]]
    )
  }
  
}else {
  
  for (j in 1:(num_education * num_race * num_puma)) {
    true_pref[j] = arm::invlogit(intercept_term +
                                   coef_race[poststrat_final_join_numeric[j, 2]] +
                                   coef_education[poststrat_final_join_numeric[j, 3]]
    )
  }
  
}

# y is number of people in poverty weighted by PWGTP, N is number of people in the PUMA area (weighted dby PWGTP)
poststrat_final_joined = cbind(poststrat_final_joined, true_pref)
poststrat_puma_pref = poststrat_final_joined %>% group_by(puma_area) %>% summarise(N_times_theta = sum(N*true_pref)) %>% as.data.frame()
poststrat_puma_pref = inner_join(poststrat_puma_pref,
                                 poststrat_final_joined %>% group_by(puma_area) %>% summarise(N_sub = sum(N)) %>% as.data.frame(),
                                 by = c("puma_area"="puma_area"))

# calculate poststratified true preference for every puma in MA
poststrat_puma_pref$puma_ps = poststrat_puma_pref$N_times_theta/poststrat_puma_pref$N_sub
summary(poststrat_puma_pref$puma_ps)

# true poststratified values based on dg process
us_pumas$puma_ps = poststrat_puma_pref$puma_ps[rank(us_pumas$PUMACE10)]

# Augment p in the simulation scheme
cbind(us_pumas$NAME10, us_pumas$PUMACE10)

p_response_race = rep(1, num_race)
p_response_education = rep(1, num_education)

N_sub = rep(0, num_puma) # the N for each PUMA based off poststrat matrix
g_counter = 1
for (g in levels_puma) {
  N_sub[g_counter] = sum(poststrat_final_joined[which(poststrat_final_joined$puma_area==g),]$N)
  g_counter = g_counter + 1
}


for (p in r) {
  
  # storage matrices for posterior estimates -----------------------------------------
  model_popn_pref_mat_icar = rep(NA, runs) # stores the preference estimate from MRP
  model_popn_pref_sd_mat_icar = rep(NA, runs) # stores sd of above estimate
  
  model_popn_pref_mat_iid = rep(NA, runs) # stores the preference estimate from MRP
  model_popn_pref_sd_mat_iid = rep(NA, runs) # stores sd of above estimate
  
  # medians of posterior linear predictors for all runs k
  sample_cell_estimates_median_icar = matrix(NA, runs, dim(poststrat_final_joined)[1])
  sample_cell_estimates_median_iid = matrix(NA, runs, dim(poststrat_final_joined)[1])
  
  # posterior width
  sample_cell_estimates_width_icar = matrix(NA, runs, dim(poststrat_final_joined)[1])
  sample_cell_estimates_width_iid = matrix(NA, runs, dim(poststrat_final_joined)[1])
  
  # the below lists stores puma posterior estimates. columns go by increasing order of puma code
  quantile10_puma_icar = matrix(NA, runs, num_puma)
  quantile90_puma_icar = matrix(NA, runs, num_puma)
  quantile10_puma_iid = matrix(NA, runs, num_puma)
  quantile90_puma_iid = matrix(NA, runs, num_puma)
  
  median_puma_icar = matrix(NA, runs, num_puma)
  median_puma_iid = matrix(NA, runs, num_puma)
  
  # we will under/oversample certain PUMA
  p_response_puma = rep(0, num_puma)
  p_response_puma[which((unique(poverty_poststrat$PUMA.x) %in% puma_overundersample_index) == TRUE)] = p
  p_response_puma[which((unique(poverty_poststrat$PUMA.x) %in% puma_overundersample_index) == FALSE)] = (1-p)
  
  p_response = rep(0, num_education * num_race * num_puma)
  for (j in 1:(num_education * num_race * num_puma)) {
    p_response[j] = p_response_puma[poststrat_final_join_numeric[j, 1]] *
      p_response_race[poststrat_final_join_numeric[j, 2]] *
      p_response_education[poststrat_final_join_numeric[j, 3]]
  }
  
  poststrat_final_joined_p = cbind(poststrat_final_joined, p_response)

  
 for (k in 1:runs) {
   print(paste("On run", k, ", p =", p))

    sample_ = sample(num_education * num_race * num_puma,
                     sample_size,
                     replace = TRUE,
                     prob = (poststrat_final_joined_p$p_response * poststrat_final_joined_p$N)/sum(poststrat_final_joined_p$p_response * poststrat_final_joined_p$N)
    )
    
    y_sample_ = rbinom(sample_size, 1, poststrat_final_joined_p$true_pref[sample_])
    
    # get covariates for every row of sample
    puma_sample = poststrat_final_joined_p[sample_,1]
    race_sample = poststrat_final_joined_p[sample_,2]
    education_sample = poststrat_final_joined_p[sample_,3]
    
    sample_final = data.frame(pref = y_sample_,
                              puma = puma_sample,
                              race = race_sample,
                              education = education_sample)
    
    # with all the factors
    sample_final$race_f = factor(sample_final$race, levels = levels_race)
    sample_final$education_f = factor(sample_final$education, levels = levels_education)
    
    
    us_pumas_df = as.data.frame(us_pumas)
    sample_final_df = left_join(sample_final, us_pumas_df, by=c("puma"="PUMACE10"))
    
    sample_final_df_binomial = sample_final_df %>% 
      group_by(puma, race_f,education_f) %>% 
      summarise(n=sum(pref),N=n())
    
    sample_final_df_binomial_df = left_join(sample_final_df_binomial, us_pumas_df, by=c("puma"="PUMACE10")) %>%
      as.data.frame()
    
    
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # prior specification inla 
    # https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html
    
    icar_model = inla(icar_formula, # pc prior: https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html#sec:pcpriors
                      data = sample_final_df_binomial_df, family ="binomial",
                      Ntrials = N,
                      control.predictor = list(compute = TRUE),
                      control.compute = list(dic = TRUE, waic = TRUE, config =TRUE),
                      control.fixed = list(prec.intercept=1) # makes precision of intercept fixed at 1
                      )
    
    # draw posterior samples for icar model 
    icar_samples = inla.posterior.sample(n=4000, icar_model)
    
    # intercept -----------------------
    icar_intercept_index_icar = grep("Intercept", rownames(icar_samples[[1]]$latent))
    
    icar_intercept_samples_f = function(x) { #x = samples[[i]]
      return(x$latent[icar_intercept_index_icar])
    }
    
    intercept_samples_icar = unlist(lapply(icar_samples, icar_intercept_samples_f))
    # ---------------------------------
    
    # race ----------------------------
    race_indices_icar = grep("race_f", rownames(icar_samples[[1]]$latent))
    
    icar_race_samples_f = function(x) {
      return(x$latent[race_indices_icar])
    }
    
    U_race_samples_icar = matrix(unlist(lapply(icar_samples, icar_race_samples_f)), 
                                 ncol = num_race,
                                 byrow = TRUE)
    # ---------------------------------
    
    # education -----------------------
    education_indices_icar = grep("education_f", rownames(icar_samples[[1]]$latent))
    
    icar_education_samples_f = function(x) {
      return(x$latent[education_indices_icar])
    }
    
    U_education_samples_icar = matrix(unlist(lapply(icar_samples, icar_education_samples_f)),
                                      ncol = num_education,
                                      byrow = TRUE)
    # ---------------------------------
    
    # puma area -----------------------
    puma_indices_icar = grep("ID:", rownames(icar_samples[[1]]$latent))[1:num_puma] # https://inla.r-inla-download.org/r-inla.org/doc/latent/bym2.pdf
    
    icar_puma_samples_f = function(x) {
      return(x$latent[puma_indices_icar])
    }
    
    U_puma_samples_icar = matrix(unlist(lapply(icar_samples, icar_puma_samples_f)),
                                 ncol = num_puma,
                                 byrow = TRUE)
    # ---------------------------------
    # get PUMA <-> ID match, where ID is used in f(ID,...)
    pumacode_id_equivalency = us_pumas_df[,c("PUMACE10","ID")]
    
    # the poststrat matrix 
    poststrat_final_joined_p_ID = left_join(poststrat_final_joined_p,
                                            pumacode_id_equivalency,
                                            by = c("puma_area"="PUMACE10"))
    
    # this matrix stores the posterior linear predictors for each poststrat cell in poststrat_final_joined_p
    postpred_sim_icar = matrix(0,
                               length(intercept_samples_icar), 
                               dim(poststrat_final_joined_p)[1])
    
    for (i in 1:(num_education * num_puma * num_race)) {
      postpred_sim_icar[,i] = arm::invlogit(intercept_samples_icar +
                                              U_puma_samples_icar[, poststrat_final_joined_p_ID[i,c("ID")] ] +
                                              U_race_samples_icar[, which(poststrat_final_joined_p_ID[i,c("race")]==levels_race) ] +
                                              U_education_samples_icar[, which(poststrat_final_joined_p_ID[i,c("education")]==levels_education) ]
      )
    }
    
    puma_mrp_samples_icar = matrix(0, 4000, num_puma) # this list stores posterior MRP samples for every puma
    
    g_counter = 1
    for (g in levels_puma) {
      #print(g)
      for (i in which(poststrat_final_joined_p_ID$puma_area==g)) {
        puma_mrp_samples_icar[,g_counter] = puma_mrp_samples_icar[,g_counter] + 
          (postpred_sim_icar[,i] * poststrat_final_joined_p_ID$N[i]/N_sub[g_counter])
      }
      g_counter = g_counter + 1
    }
    
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------------------------
    sample_final_df_binomial_df_iid = sample_final_df_binomial_df
    sample_final_df_binomial_df_iid$ID_f = factor(sample_final_df_binomial_df_iid$ID, levels = 1:num_puma)
    
    iid_model = inla(iid_formula, # pc prior: https://becarioprecario.bitbucket.io/inla-gitbook/ch-priors.html#sec:pcpriors
                     data = sample_final_df_binomial_df_iid, family ="binomial",
                     Ntrials = N,
                     control.predictor = list(compute = TRUE),
                     control.compute = list(dic = TRUE, waic = TRUE, config =TRUE),
                     control.fixed = list(prec.intercept=1) # makes precision of intercept fixed at 1
                     )
    
    
    # draw posterior samples for iid model 
    iid_samples = inla.posterior.sample(n=4000, iid_model)
    
    # intercept -----------------------
    iid_intercept_index_iid = grep("Intercept", rownames(iid_samples[[1]]$latent))
    
    iid_intercept_samples_f = function(x) { #x = samples[[i]]
      return(x$latent[iid_intercept_index_iid])
    }
    
    intercept_samples_iid = unlist(lapply(iid_samples, iid_intercept_samples_f))
    # ---------------------------------
    
    # race ----------------------------
    race_indices_iid = grep("race_f", rownames(iid_samples[[1]]$latent))
    
    iid_race_samples_f = function(x) {
      return(x$latent[race_indices_iid])
    }
    
    U_race_samples_iid = matrix(unlist(lapply(iid_samples, iid_race_samples_f)), 
                                ncol = num_race,
                                byrow = TRUE)
    # ---------------------------------
    
    # education -----------------------
    education_indices_iid = grep("education_f", rownames(iid_samples[[1]]$latent))
    
    iid_education_samples_f = function(x) {
      return(x$latent[education_indices_iid])
    }
    
    U_education_samples_iid = matrix(unlist(lapply(iid_samples, iid_education_samples_f)),
                                     ncol = num_education,
                                     byrow = TRUE)
    # ---------------------------------
    
    # puma area -----------------------
    puma_indices_iid = grep("ID_f:", rownames(iid_samples[[1]]$latent))
    
    iid_puma_samples_f = function(x) {
      return(x$latent[puma_indices_iid])
    }
    
    U_puma_samples_iid = matrix(unlist(lapply(iid_samples, iid_puma_samples_f)),
                                ncol = num_puma,
                                byrow = TRUE)
    # ---------------------------------
    
    # this matrix stores the posterior linear predictors for each poststrat cell in poststrat_final_joined_p
    postpred_sim_iid = matrix(0,
                              length(intercept_samples_iid), 
                              dim(poststrat_final_joined_p)[1])
    
    for (i in 1:(num_education * num_puma * num_race)) {
      postpred_sim_iid[,i] = arm::invlogit(intercept_samples_iid +
                                             U_puma_samples_iid[, poststrat_final_joined_p_ID[i,c("ID")] ] +
                                             U_race_samples_iid[, which(poststrat_final_joined_p_ID[i,c("race")]==levels_race) ] +
                                             U_education_samples_iid[, which(poststrat_final_joined_p_ID[i,c("education")]==levels_education) ]
      )
    }
    
    puma_mrp_samples_iid = matrix(0, 4000, num_puma) # this list stores posterior MRP samples for every puma
    
    g_counter = 1
    for (g in levels_puma) {
      #print(g)
      for (i in which(poststrat_final_joined_p_ID$puma_area==g)) {
        puma_mrp_samples_iid[,g_counter] = puma_mrp_samples_iid[,g_counter] + 
          (postpred_sim_iid[,i] * poststrat_final_joined_p_ID$N[i]/N_sub[g_counter])
      }
      g_counter = g_counter + 1
    }
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------------------------------
    true_popn_pref = poststrat_final_joined_p$true_pref %*% poststrat_final_joined_p$N / sum(poststrat_final_joined_p$N)
    
    model_popn_pref_mat_icar[k] = mean(postpred_sim_icar %*% poststrat_final_joined_p$N / sum(poststrat_final_joined_p$N))
    model_popn_pref_sd_mat_icar[k] = sd(postpred_sim_icar %*% poststrat_final_joined_p$N / sum(poststrat_final_joined_p$N))
    
    model_popn_pref_mat_iid[k] = mean(postpred_sim_iid %*% poststrat_final_joined_p$N / sum(poststrat_final_joined_p$N))
    model_popn_pref_sd_mat_iid[k] = sd(postpred_sim_iid %*% poststrat_final_joined_p$N / sum(poststrat_final_joined_p$N))
    
    # medians of posterior linear predictors for all runs k
    sample_cell_estimates_median_icar[k,] = matrixStats::colMedians(postpred_sim_icar)
    sample_cell_estimates_median_iid[k,] = matrixStats::colMedians(postpred_sim_iid)
    
    # posterior width
    sample_cell_estimates_width_icar[k,] = apply(X = postpred_sim_icar,
                                                 MARGIN = 2,
                                                 FUN = quantile, 
                                                 probs = upperquantile) -
      apply(X = postpred_sim_icar, 
            MARGIN = 2,
            FUN = quantile,
            probs = lowerquantile)
    
      
    sample_cell_estimates_width_iid[k,] = apply(X = postpred_sim_iid,
                                                MARGIN = 2,
                                                FUN = quantile, 
                                                probs = upperquantile) -
      apply(X = postpred_sim_iid, 
            MARGIN = 2,
            FUN = quantile,
            probs = lowerquantile)
    
    # # the below lists stores puma posterior estimates. columns go by increasing order of puma code

    median_puma_icar[k,] = matrixStats::colMedians(puma_mrp_samples_icar)
    median_puma_iid[k,] = matrixStats::colMedians(puma_mrp_samples_iid)
    
    quantile10_puma_icar[k,] = apply(X = puma_mrp_samples_icar,
                                     MARGIN = 2,
                                     FUN = quantile,
                                     probs = lowerquantile)
    
    quantile90_puma_icar[k,] = apply(X = puma_mrp_samples_icar,
                                     MARGIN = 2,
                                     FUN = quantile,
                                     probs = upperquantile)
    
    quantile10_puma_iid[k,] = apply(X = puma_mrp_samples_iid,
                                    MARGIN = 2,
                                    FUN = quantile,
                                    probs = lowerquantile)
    
    quantile90_puma_iid[k,] = apply(X = puma_mrp_samples_iid,
                                    MARGIN = 2,
                                    FUN = quantile,
                                    probs = upperquantile)
    
    
    if ((k %% 10)==0) { 
      print(k)
      save.image(paste0(runs,
                        "_",
                        sample_size,
                        "_spatial_",
                        p * 10,
                        ".RData"))
    }

  }
}









rm(list=ls())

library(rstan)
library(dplyr)

options(mc.cores = parallel::detectCores()) # use multiple cores. default is 4.
rstan_options(auto_write = TRUE) # save compiled stan model to hard disk so no need to recompile

# ---------------------------------------------------
age_grouping_multiplier = 12 # the number of strata for age

iterations = 2000
num_chains = 4
mtreedepth = 15
ad = 0.99

use_state = TRUE # if this is true we use state as response

# m_file = "baselinemeanzeroN01_annenberg_v3_centered.stan"
# m_name = "baseline_model_annenberg_centered"

m_file = "proposedarN01_annenberg_v3_centered.stan"
m_name = "ar_model_annenberg_centered"

# m_file = "proposedN01_annenberg_v3_centered.stan"
# m_name = "rw_model_annenberg_centered"

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


m = stan_model(file = m_file, model_name = m_name)

acs_ps = readRDS("acs_ps.rds") # load poststratification matrix
# "AGEP"        "SEX"         "race_x"      "education_x" "state_x"     "income_x"    "N"   
colnames(acs_ps) = c("age", "sex", "race", "education", "state", "income", "N")
# ---------------------------------------------------

survey_dat = readr::read_tsv("annenbergphone2008.txt", col_names = TRUE)

# Preprocess -------------------------------------------------------------------------

# Load state data from jkastell -----------------------------------------------------
# Load state-level data from http://www.princeton.edu/~jkastell/mrp_primer.html
StateLevel = foreign::read.dta("state_level_update.dta",convert.underscore = TRUE)
StateLevel = StateLevel[order(StateLevel$sstate.initnum),]

# Load state-level data from http://www.princeton.edu/~jkastell/mrp_primer.html
MarriageDataSP = foreign::read.dta("gay_marriage_megapoll.dta", convert.underscore = TRUE) 

# check that row names of Statelevel are equal to sstate.initnum
sum(as.numeric(as.character(rownames(StateLevel))) == StateLevel$sstate.initnum)

stateregions = unique(MarriageDataSP[,c("statename", "state.initnum", "region.cat", "region")])

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

stateregions_final$lowercasename = tolower(stringr::str_replace_all(stateregions_final$sstatename, " ", ""))

stateregions_final$lowercasename[stateregions_final$lowercasename=="d.c."]="columbia"

# ------------------------------------------------------------------------------------

survey_dat_selected = survey_dat %>% 
  select(covariates_list) %>% 
  filter(WA02_c < 100, # Age less than 100
         CEc01_c %in% c(1,3), # 1 is vote yes for gay marriage, 3 is no
  ) #%>%
  # mutate(gayfavor = as.numeric(CEc01_c),
  #        sex = as.numeric(WA01_c),
  #        age = as.numeric(WA02_c),
  #        education = as.numeric(WA03_c),
  #        race = as.numeric(WC03_c),
  #        income = as.numeric(WA04_c),
  #        state = as.numeric(WFc01_c),
  #        yearsinus = as.numeric(WC06_c),
  #        adultsinhouse = as.numeric(WFa01_c),
  #        orientation = as.numeric(WHb01_c)
  # ) # create new variable names

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

# no need to preprocess sex covariate

# no need to preprocess age covariate

# Select covariates that we're interested in now ----------------------------------------------------

gaymarriagedata = survey_dat_selected[,c("CEc01_c", # favors gay marriage
                                         "WA01_c", # sex
                                         "WA02_c", # age
                                         "education_x", # education
                                         "race_x", # race
                                         "WA04_c", # household income
                                         "state_x")] %>% data.frame() # state

# change to characters because levels in maximal poststrat matrix are in characters 
gaymarriagedata$CEc01_c = as.character(gaymarriagedata$CEc01_c)
gaymarriagedata$WA01_c = as.character(gaymarriagedata$WA01_c)
gaymarriagedata$WA04_c = as.character(gaymarriagedata$WA04_c)

gaymarriagedata_rm = gaymarriagedata[complete.cases(gaymarriagedata),] # remove NAs


acs_ps$race = factor(acs_ps$race)
acs_ps$sex = factor(acs_ps$sex)
acs_ps$education = factor(acs_ps$education)
acs_ps$income = factor(acs_ps$income)
acs_ps$state = factor(acs_ps$state)

# store original levels from maximal poststrat matrix
levels_race = levels(acs_ps$race)
levels_sex = levels(acs_ps$sex)
levels_education = levels(acs_ps$education)
levels_income = levels(acs_ps$income)
levels_state = levels(acs_ps$state)

# do some checks to see that we didn't make a typo. the levels from the maximal poststrat matrix should be the superset
unique(gaymarriagedata_rm$race_x) %in% levels_race
unique(gaymarriagedata_rm$WA01_c) %in% levels_sex
unique(gaymarriagedata_rm$education) %in% levels_education
unique(gaymarriagedata_rm$WA04_c) %in% levels_income
unique(gaymarriagedata_rm$state_x) %in% levels_state

# get state index from jkastell in this dataframe as well

# check that all state names in the annenberg survey are the same as in  jkastell's file
sum(unique(gaymarriagedata_rm$state_x) %in% unique(stateregions_final$lowercasename)) # should equal 49

# do left_join 
gaymarriagedata_rm = dplyr::left_join(x = gaymarriagedata_rm, y = stateregions_final[,c("sstate.initnum", "region.cat", "lowercasename")], 
                 by = c("state_x"="lowercasename"))


# use above levels to factor survey data
gaymarriagedata_rm$race_x = factor(gaymarriagedata_rm$race_x, levels = levels_race) # factor race
gaymarriagedata_rm$WA01_c = factor(gaymarriagedata_rm$WA01_c, levels = levels_sex) # factor sex
gaymarriagedata_rm$education_x = factor(gaymarriagedata_rm$education_x, levels = levels_education) # factor education
gaymarriagedata_rm$WA04_c = factor(gaymarriagedata_rm$WA04_c, levels = levels_income) # factor income
gaymarriagedata_rm$state_x = factor(gaymarriagedata_rm$state_x, levels = levels_state) # factor state


gaymarriagedata_rm_clean_mapped = gaymarriagedata_rm
# rename levels in our survey data to numerical values so we can feed it into a stan file
levels(gaymarriagedata_rm_clean_mapped$race_x) = c(1, 2, 3, 4, 5, 6) # americanindian=1, asian=2, black=3, hisp=4, other=5, white=6
levels(gaymarriagedata_rm_clean_mapped$WA01_c) = c(0, 1) # men=0, women=1
levels(gaymarriagedata_rm_clean_mapped$education_x) = c(1, 2, 3, 4, 5, 6) # fouryeardegree=1, fouryeardegreeplus=2, highschool=3, nohighschool=4, somecollege=5, twoyeardegree=6
levels(gaymarriagedata_rm_clean_mapped$state_x) = 1:51

#levels(gaymarriagedata_rm_clean_mapped$gayFavorFederalMarriage) = c(0, 1) # no==0, yes==1

summary(gaymarriagedata_rm_clean_mapped)
summary(gaymarriagedata_rm)
# create age-category covariate i.e stratify age
# 1. acs_ps$age ranges from 18 to 95 so we remove 96 and 97
gaymarriagedata_rm_clean_mapped = gaymarriagedata_rm_clean_mapped %>% filter(WA02_c <= 95 & WA02_c >= 18)
gaymarriagedata_rm_clean_mapped$age_cat = cut(gaymarriagedata_rm_clean_mapped$WA02_c,
                                              age_grouping_multiplier)

# store levels of age_cat
levels_age_cat = levels(gaymarriagedata_rm_clean_mapped$age_cat)

# change levels of age_cat to numeric
levels(gaymarriagedata_rm_clean_mapped$age_cat) = 1:age_grouping_multiplier


# -----------------------------------------------------------------------------

# rename columns 
print(colnames(gaymarriagedata_rm_clean_mapped))
# "CEc01_c"     "WA01_c"      "WA02_c"      "education_x" "race_x"      "WA04_c"      "state_x"     "age_cat"    
colnames(gaymarriagedata_rm_clean_mapped) = c("gayfavor", "sex", "age", "education", "race", "income", "state", "jk_state_index", "jk_region_cat", "age_cat")


# we do a left join because we're using the state index ordering that I defined
stateregions_final_ = dplyr::left_join(x = data.frame(levels_state, stringsAsFactors = FALSE), y = stateregions_final, by = c("levels_state"="lowercasename"))

# we will use the sstate.initnum order from stateregions_final
# coef_state comes from data used in http://www.princeton.edu/~jkastell/MRP_primer/mrp_primer.pdf
coef_state = (stateregions_final_$kerry.04)/100 # value for 2004 Democratic vote share in every state. This is X_{\text{state-vs},j} for j \in [51]

# coef_relig comes from data used in http://www.princeton.edu/~jkastell/MRP_primer/mrp_primer.pdf. this is the percentage of conservative religions in every state
coef_relig = (stateregions_final_$p.evang + stateregions_final_$p.mormon)/100

# unfactor stuff in gaymarriagedata_rm_clean_mapped
X = gaymarriagedata_rm_clean_mapped[,!(names(gaymarriagedata_rm_clean_mapped) %in% c("gayfavor", "age"))]
Y = as.numeric(as.character(gaymarriagedata_rm_clean_mapped[,names(gaymarriagedata_rm_clean_mapped)=="gayfavor"]))

# change covariates to numeric so we can feed into stan model
X$sex = as.numeric(as.character(X$sex))
X$age_cat = as.numeric(as.character(X$age_cat))
X$education = as.numeric(as.character(X$education))
X$race = as.numeric(as.character(X$race))
X$income = as.numeric(as.character(X$income))
X$state = as.numeric(as.character(X$state))

# note that states 2 and 12 (alaska and hawaii) are not present in the survey sample

fit_realdata = sampling(m, 
                        data = list(N = dim(X)[1],
                                    N_groups_age = age_grouping_multiplier, # dependent on number of age groups
                                    N_groups_race = length(unique(acs_ps$race)), # there are 4 categories in race
                                    N_groups_sex = length(unique(acs_ps$sex)), # there are 2 categories in gender
                                    N_groups_education = length(unique(acs_ps$education)),
                                    N_groups_income = length(unique(acs_ps$income)),
                                    
                                    state_vs = coef_state, # the 51 values for 2004 Democratic vote share based on the state index ordering I defined for levels_state
                                    relig = coef_relig, # the 51 values for conservative religion based on the state index ordering I defined for levels_state
                                    region_index =  stateregions_final_$region.cat, # the 51 values for state region based on the state index ordering I defined for levels_state
                                    
                                    N_groups_state = length(unique(acs_ps$state)),
                                    age = X$age_cat,
                                    race = X$race,
                                    sex = X$sex,
                                    education = X$education,
                                    income = X$income,
                                    state = X$state, # we are using the state index that I defined to fit the model
                                    y = Y),
                        iter=iterations, 
                        chains=num_chains,
                        seed=21,
                        control=list(max_treedepth=mtreedepth, adapt_delta=ad))

saveRDS(fit_realdata,
        paste0(m_name,
               "_",
               age_grouping_multiplier,
               ".rds"))

print("Done script.")

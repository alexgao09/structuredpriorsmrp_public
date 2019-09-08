rm(list=ls())

# for 5 YEAR ACS 2006 - 2010
#load into the ACS and make PS table.
library(dplyr)
library(readr)
library(data.table)

# data dictionary below
# https://www2.census.gov/programs-surveys/acs/tech_docs/pums/data_dict/PUMS_Data_Dictionary_2006-2010.pdf?#

# load data
acs1 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_pus/ss10pusa.csv", sep = ",", header =TRUE)
acs1 = acs1[acs1$AGEP>17,]
acs1_tr = acs1[,c("serialno",
                  "AGEP", # age
                  "SEX", # sex
                  "PWGTP", # person's weight
                  "HISP", # hispanic ethnicity
                  "RAC1P", # race
                  "ST", # state
                  "SCHL", # Educational attainment,
                  "RT" # record type
)]
rm(acs1)
gc() # free up space


acs2 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_pus/ss10pusb.csv", sep = ",", header =TRUE)
acs2 = acs2[acs2$AGEP>17,]
acs2_tr = acs2[,c("serialno",
                  "AGEP", # age
                  "SEX", # sex
                  "PWGTP", # person's weight
                  "HISP", # hispanic ethnicity: FOR RACE
                  "RAC1P", # race: FOR RACE
                  "ST", # state
                  "SCHL", # Educational attainment
                  "RT" # record type
)] 
rm(acs2)
gc() # free up space


acs3 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_pus/ss10pusc.csv", sep = ",", header =TRUE)
acs3 = acs3[acs3$AGEP>17,]
acs3_tr = acs3[,c("serialno",
                  "AGEP", # age
                  "SEX", # sex
                  "PWGTP", # person's weight
                  "HISP", # hispanic ethnicity
                  "RAC1P", # race
                  "ST", # state
                  "SCHL", # Educational attainment
                  "RT"
)]
rm(acs3)
gc()


acs4 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_pus/ss10pusd.csv", sep = ",", header =TRUE)
acs4 = acs4[acs4$AGEP>17,]
acs4_tr = acs4[,c("serialno",
                  "AGEP", # age
                  "SEX", # sex
                  "PWGTP", # person's weight
                  "HISP", # hispanic ethnicity: FOR RACE
                  "RAC1P", # race: FOR RACE
                  "ST", # state
                  "SCHL", # Educational attainment,
                  "RT"
)] 
rm(acs4)
gc()

acs_tr = rbind(acs1_tr, acs2_tr, acs3_tr, acs4_tr)



acsh1 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_hus/ss10husa.csv", sep = ",", header =TRUE)
acsh1_tr = acsh1[,c("serialno", "HINCP")]
rm(acsh1)
gc()


acsh2 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_hus/ss10husb.csv", sep = ",", header =TRUE)
acsh2_tr = acsh2[,c("serialno","HINCP")]
rm(acsh2)
gc()


acsh3 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_hus/ss10husc.csv", sep = ",", header =TRUE)
acsh3_tr = acsh3[,c("serialno","HINCP")]
rm(acsh3)
gc()


acsh4 = fread("~/Desktop/annenbergdata/fiveyearacs/csv_hus/ss10husd.csv", sep = ",", header =TRUE)
acsh4_tr = acsh4[,c("serialno","HINCP")]
rm(acsh4)
gc()



#sum(names(acsh1)==names(acsh2))
#sum(names(acsh2)==names(acsh3))
#sum(names(acsh3)==names(acsh4))


acsh_tr = rbind(acsh1_tr, acsh2_tr, acsh3_tr, acsh4_tr)

# check that there are no duplicate rows in acsh_tr
print(length(unique(acsh_tr$serialno)))
print(dim(acsh_tr))

acs_tr = dplyr::left_join(acs_tr, acsh_tr, by = c("serialno"="serialno"))

rm(acs1_tr, acs2_tr, acs3_tr, acs4_tr)
rm(acsh1_tr, acsh2_tr, acsh3_tr, acsh4_tr)

acs_tr = acs_tr[complete.cases(acs_tr),] # remove NAs
# PREPROCESS RACE

# Recoded detailed race code
# 1 .White alone
# 2 .Black or African American alone
# 3 .American Indian alone
# 4 .Alaska Native alone
# 5 .American Indian and Alaska Native tribes specified; or American
# .Indian or Alaska native, not specified and no other races
# 6 .Asian alone
# 7 .Native Hawaiian and Other Pacific Islander alone
# 8 .Some other race alone
# 9 .Two or more major race groups

acs_tr$race_x = ifelse(acs_tr$RAC1P==1 & acs_tr$HISP==1,"white",
                       ifelse(acs_tr$RAC1P==2 & acs_tr$HISP==1,"black",
                              ifelse(acs_tr$RAC1P==6 & acs_tr$HISP==1, "asian",
                                     ifelse(acs_tr$RAC1P==3 & acs_tr$HISP==1,"americanindian",
                                            ifelse(!(acs_tr$RAC1P %in% c(1,2,6,3)) & acs_tr$HISP %in% 2:24, "hisp",
                                                   "other")
                                     )
                              )
                       )
)


# PREPROCESS EDUCATION

# Educational attainment
# bb .N/A (less than 3 years old)
# 01 .No schooling completed
# 02 .Nursery school to grade 4
# 03 .Grade 5 or grade 6
# 04 .Grade 7 or grade 8
# 05 .Grade 9
# 06 .Grade 10
# 07 .Grade 11
# 08 .12th grade, no diploma
# 09 .High school graduate
# 10 .Some college, but less than 1 year
# 11 .One or more years of college, no degree
# 12 .Associate's degree
# 13 .Bachelor's degree
# 14 .Master's degree
# 15 .Professional school degree
# 16 .Doctorate degree

acs_tr$education_x = ifelse(acs_tr$SCHL %in% 1:8, "nohighschool",
                            ifelse(acs_tr$SCHL==9, "highschool",
                                   ifelse(acs_tr$SCHL %in% c(10,11), "somecollege",
                                          ifelse(acs_tr$SCHL==12, "twoyeardegree",
                                                 ifelse(acs_tr$SCHL==13, "fouryeardegree",
                                                        "fouryeardegreeplus")
                                          )
                                   )
                            )
)

# PREPROCESS STATE

# 01 .Alabama/AL
# 02 .Alaska/AK
# 04 .Arizona/AZ
# 05 .Arkansas/AR
# 06 .California/CA
# 08 .Colorado/CO
# 09 .Connecticut/CT
# 10 .Delaware/DE
# 11 .District of Columbia/DC
# 12 .Florida/FL
# 13 .Georgia/GA
# 15 .Hawaii/HI
# 16 .Idaho/ID
# 17 .Illinois/IL
# 18 .Indiana/IN
# 19 .Iowa/IA
# 20 .Kansas/KS
# 21 .Kentucky/KY
# 22 .Louisiana/LA
# 23 .Maine/ME
# 24 .Maryland/MD
# 25 .Massachusetts/MA
# 26 .Michigan/MI
# 27 .Minnesota/MN
# 28 .Mississippi/MS
# 29 .Missouri/MO
# 30 .Montana/MT
# 31 .Nebraska/NE
# 32 .Nevada/NV
# 33 .New Hampshire/NH
# 34 .New Jersey/NJ 
# 35 .New Mexico/NM
# 36 .New York/NY
# 37 .North Carolina/NC
# 38 .North Dakota/ND
# 39 .Ohio/OH
# 40 .Oklahoma/OK
# 41 .Oregon/OR
# 42 .Pennsylvania/PA
# 44 .Rhode Island/RI
# 45 .South Carolina/SC
# 46 .South Dakota/SD
# 47 .Tennessee/TN
# 48 .Texas/TX
# 49 .Utah/UT
# 50 .Vermont/VT
# 51 .Virginia/VA
# 53 .Washington/WA
# 54 .West Virginia/WV
# 55 .Wisconsin/WI
# 56 .Wyoming/WY
# 72 .Puerto Rico/PR

# DOUBLE CHECK BELOW
acs_tr$state_x = case_when(acs_tr$ST == 1 ~ "alabama",
                           acs_tr$ST == 2 ~ "alaska",
                           acs_tr$ST == 4 ~ "arizona",
                           acs_tr$ST == 5 ~ "arkansas",
                           acs_tr$ST == 6 ~ "california",
                           acs_tr$ST == 8 ~ "colorado",
                           acs_tr$ST == 9 ~ "connecticut",
                           acs_tr$ST == 10 ~ "delaware",
                           acs_tr$ST == 11 ~ "columbia",
                           acs_tr$ST == 12 ~ "florida",
                           acs_tr$ST == 13 ~ "georgia",
                           acs_tr$ST == 15 ~ "hawaii",
                           acs_tr$ST == 16 ~ "idaho",
                           acs_tr$ST == 17 ~ "illinois",
                           acs_tr$ST == 18 ~ "indiana",
                           acs_tr$ST == 19 ~ "iowa",
                           acs_tr$ST == 20 ~ "kansas",
                           acs_tr$ST == 21 ~ "kentucky",
                           acs_tr$ST == 22 ~ "louisiana",
                           acs_tr$ST == 23 ~ "maine",
                           acs_tr$ST == 24 ~ "maryland",
                           acs_tr$ST == 25 ~ "massachusetts",
                           acs_tr$ST == 26 ~ "michigan",
                           acs_tr$ST == 27 ~ "minnesota",
                           acs_tr$ST == 28 ~ "mississippi",
                           acs_tr$ST == 29 ~ "missouri",
                           acs_tr$ST == 30 ~ "montana",
                           acs_tr$ST == 31 ~ "nebraska",
                           acs_tr$ST == 32 ~ "nevada",
                           acs_tr$ST == 33 ~ "newhampshire",
                           acs_tr$ST == 34 ~ "newjersey",
                           acs_tr$ST == 35 ~ "newmexico",
                           acs_tr$ST == 36 ~ "newyork",
                           acs_tr$ST == 37 ~ "northcarolina",
                           acs_tr$ST == 38 ~ "northdakota",
                           acs_tr$ST == 39 ~ "ohio",
                           acs_tr$ST == 40 ~ "oklahoma",
                           acs_tr$ST == 41 ~ "oregon",
                           acs_tr$ST == 42 ~ "pennsylvania",
                           acs_tr$ST == 44 ~ "rhodeisland",
                           acs_tr$ST == 45 ~ "southcarolina",
                           acs_tr$ST == 46 ~ "southdakota",
                           acs_tr$ST == 47 ~ "tennessee",
                           acs_tr$ST == 48 ~ "texas",
                           acs_tr$ST == 49 ~ "utah",
                           acs_tr$ST == 50 ~ "vermont",
                           acs_tr$ST == 51 ~ "virginia",
                           acs_tr$ST == 53 ~ "washington",
                           acs_tr$ST == 54 ~ "westvirginia",
                           acs_tr$ST == 55 ~ "wisconsin",
                           acs_tr$ST == 56 ~ "wyoming",
                           TRUE ~ as.character(acs_tr$ST))


# PREPROCESS INCOME
# match annenberg coding
acs_tr$income_x = case_when(acs_tr$HINCP < 10000 ~ 1,
                            acs_tr$HINCP >= 10000 & acs_tr$HINCP < 15000 ~ 2,
                            acs_tr$HINCP >= 15000 & acs_tr$HINCP < 25000 ~ 3,
                            acs_tr$HINCP >= 25000 & acs_tr$HINCP < 35000 ~ 4,
                            acs_tr$HINCP >= 35000 & acs_tr$HINCP < 50000 ~ 5,
                            acs_tr$HINCP >= 50000 & acs_tr$HINCP < 75000 ~ 6,
                            acs_tr$HINCP >= 75000 & acs_tr$HINCP < 100000 ~ 7,
                            acs_tr$HINCP >= 100000 & acs_tr$HINCP < 150000 ~ 8,
                            acs_tr$HINCP >= 150000 ~ 9)

# PREPROCESS SEX
# 1 .Male
# 2 .Female

# PREPROCESS AGE
# no need to preprocess age

# construct maximal poststratification matrix --------------------------

# 1. do dplyr grouping
acs_ps = acs_tr %>%
  group_by(AGEP, SEX, race_x, education_x, state_x, income_x) %>%
  summarize(N = sum(PWGTP)) %>% # why sum over PWGTP?
  ungroup()

# check that you don't have all subpopulation cells
# prod(apply(acs_ps, 2, function(x) {length(unique(x))})[-7])

# 2. create empty maximal poststrat matrix
acs_ps_maximal_x = expand.grid(age = unique(acs_tr$AGEP),
                               sex = unique(acs_tr$SEX),
                               race = unique(acs_tr$race_x),
                               education = unique(acs_tr$education_x),
                               state = unique(acs_tr$state_x),
                               income = unique(acs_tr$income_x),
                               stringsAsFactors = FALSE
)

# 3. grab the groupings from acs_ps
acs_ps_maximal = dplyr::left_join(acs_ps_maximal_x, acs_ps, by = c("age"="AGEP",
                                                                   "sex"="SEX",
                                                                   "race"="race_x",
                                                                   "education"="education_x",
                                                                   "state"="state_x",
                                                                   "income"="income_x"))


# 4. Save poststratification matrix
saveRDS(data.frame(acs_ps),file='~/Desktop/annenbergdata/fiveyearacs/acs_ps.rds')

rm(list=ls())


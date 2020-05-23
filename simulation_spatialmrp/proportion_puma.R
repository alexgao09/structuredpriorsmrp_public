# this file calculates the proportion that structured priors outperform baseline priors in MRP for spatial simulation

rm(list=ls())
gc()


library(ggplot2)
library(dplyr)
library(sp)
library(reshape2)

sample_size = 1000
runs = 200
r = 1:9/10

poverty_poststrat = readRDS("poverty_poststrat.rds")
us_pumas = readRDS("us_pumas.rds")
ma_adj_matrix = readRDS("ma_adj_matrix.rds")
ma_adj_matrix_sparse = readRDS("ma_adj_matrix_sparse.rds")

num_education = length(unique(poverty_poststrat$education)) # number of strata for education
num_race = length(unique(poverty_poststrat$race_x)) # number of strata for race
num_puma = length(unique(poverty_poststrat$PUMA.x))


icar_iid_puma_bias_comparison = matrix(0, length(r), num_puma)
icar_iid_puma_sd_comparison = matrix(0, length(r), num_puma)


counter_ = 1
prob_sampling_old = rep(0, length(r))

for (p in r) {
  load(paste0(runs,
              "_",
              sample_size,
              "_spatial_",
              p * 10,
              ".RData"))
  
  # get the interpretable p which is the probability of sampling something in puma_overundersample_index
  poststrat_final_joined_p_interpretable = cbind(poststrat_final_joined_p,
                                                 (poststrat_final_joined_p$p_response * poststrat_final_joined_p$N)/sum(poststrat_final_joined_p$p_response * poststrat_final_joined_p$N)
  )
  # 
  colnames(poststrat_final_joined_p_interpretable)[length(colnames(poststrat_final_joined_p_interpretable))] = "prob_sampling_g1" 
  prob_sampling_old[counter_] = round(sum(poststrat_final_joined_p_interpretable[poststrat_final_joined_p_interpretable$puma_area %in% puma_overundersample_index,
                                                                                 c("prob_sampling_g1")]),
                                      2) # very important
  
  
  icar_iid_puma_bias_comparison[counter_,] = colSums((abs(sweep(median_puma_icar, 
                                                                2,
                                                                poststrat_puma_pref$puma_ps,
                                                                "-")) - 
                                                        abs(sweep(median_puma_iid,
                                                                  2, 
                                                                  poststrat_puma_pref$puma_ps,
                                                                  "-"))) <= 0)/dim(median_puma_icar)[1]
  
  
  icar_iid_puma_sd_comparison[counter_,] = colSums(((quantile90_puma_icar - quantile10_puma_icar) -
                                                      (quantile90_puma_iid - quantile10_puma_iid)) <= 0)/dim(median_puma_icar)[1]
  

  
  counter_ = counter_ + 1
}

colnames(icar_iid_puma_bias_comparison) = poststrat_puma_pref$puma_area
colnames(icar_iid_puma_sd_comparison) = poststrat_puma_pref$puma_area

icar_iid_puma_bias_comparison_melted = melt(icar_iid_puma_bias_comparison)
icar_iid_puma_sd_comparison_melted = melt(icar_iid_puma_sd_comparison)

colnames(icar_iid_puma_bias_comparison_melted) = c("p", "puma", "proportion")
colnames(icar_iid_puma_sd_comparison_melted) = c("p", "puma", "proportion")

icar_iid_puma_bias_comparison_melted$puma_group = ""
icar_iid_puma_sd_comparison_melted$puma_group = ""

icar_iid_puma_bias_comparison_melted[icar_iid_puma_bias_comparison_melted$puma %in% puma_overundersample_index,c("puma_group")] = "Near Boston"
icar_iid_puma_bias_comparison_melted[!(icar_iid_puma_bias_comparison_melted$puma %in% puma_overundersample_index),c("puma_group")] = "Away From Boston"

icar_iid_puma_sd_comparison_melted[icar_iid_puma_sd_comparison_melted$puma %in% puma_overundersample_index,c("puma_group")] = "Near Boston"
icar_iid_puma_sd_comparison_melted[!(icar_iid_puma_sd_comparison_melted$puma %in% puma_overundersample_index),c("puma_group")] = "Away From Boston"

icar_iid_puma_bias_comparison_melted$prob_sampling = rep(times = 52, x = prob_sampling_old)
icar_iid_puma_sd_comparison_melted$prob_sampling = rep(times = 52, x = prob_sampling_old)


# bias plot
saveRDS(icar_iid_puma_bias_comparison_melted,
        paste0("proportion_puma_", sample_size,".rds"))

ggplot(icar_iid_puma_bias_comparison_melted,
       aes(x=prob_sampling, y=proportion, group = puma)) + 
  geom_line(aes(col=as.factor(puma_group)),show.legend = TRUE, size = 0.5) +
  xlab("\n Probability of sampling for cluster of PUMA near Boston \n") +
  ylab("\n Improvement proportion for bias \n") +
  guides(col=guide_legend(title="Sampling Group",
                          override.aes = list(size = 10*1.6))) +
  geom_hline(yintercept=0.5, color="black", size=1, linetype = "dashed") + 
  scale_color_viridis_d(begin=0.25,end=.75)


# sd plot
saveRDS(icar_iid_puma_sd_comparison_melted,
        paste0("proportion_sd_puma_", sample_size,".rds"))

ggplot(icar_iid_puma_sd_comparison_melted,
       aes(x=prob_sampling, y=proportion, group = puma)) + 
  geom_line(aes(col=as.factor(puma_group)),show.legend = TRUE, size = 0.5) +
  xlab("\n Probability of sampling for cluster of PUMA near Boston \n") +
  ylab("\n Improvement proportion for sd \n") +
  guides(col=guide_legend(title="Sampling Group",
                          override.aes = list(size = 10*1.6))) +
  geom_hline(yintercept=0.5, color="black", size=1, linetype = "dashed") + 
  scale_color_viridis_d(begin=0.25,end=.75)


print(colMeans(icar_iid_puma_bias_comparison))
print(colMeans(icar_iid_puma_sd_comparison))

print(summary(colMeans(icar_iid_puma_bias_comparison)))
print(summary(colMeans(icar_iid_puma_sd_comparison)))


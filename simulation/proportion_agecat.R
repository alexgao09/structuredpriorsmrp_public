rm(list=ls()) 
gc()

options(bitmapType="cairo")

library(reshape2)
library(rstan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(knitr)
library(matrixStats)

# Global params -------------------------------------------

save_ridgeplots = TRUE 
sample_size = 1000
runs = 100
r = 1:9/10
income_multiplier = 1 # partitions income into more categories
age_grouping_multiplier = 12 # how much we take the maximal poststratification. make sure this can divide 60
number_of_states = 51

response_tag = "binary"

# ---------------------------------------------------------
rw_baseline_agecat_bias_comparison = matrix(0, length(r), age_grouping_multiplier)
ar_baseline_agecat_bias_comparison = matrix(0, length(r), age_grouping_multiplier)

rw_baseline_agecat_sd_comparison = matrix(0, length(r), age_grouping_multiplier)
ar_baseline_agecat_sd_comparison = matrix(0, length(r), age_grouping_multiplier)


counter = 1
prob_sampling_old = rep(0, length(r))

for (p in r) {
  
  load(paste0(
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
    ".RData"))
  
  points_df_final = unique(readRDS("points_df_final.rds")[,1:2]) %>% arrange(age_cat) # the true poststratified preferences for every age category
  
  # Get the interpretable p which is the probability of sampling someone in ages 61 - 80
  poststrat_final_p = cbind(poststrat_final, 
                            (poststrat_final$p_response * poststrat_final$N)/sum((poststrat_final$p_response * poststrat_final$N))) # this is the probability of sampling someone in ages  61 - 80
  colnames(poststrat_final_p)[length(colnames(poststrat_final_p))] = "prob_sampling_old" 
  prob_sampling_old[counter] = sum(poststrat_final_p[as.numeric(as.character(poststrat_final_p$age_cat)) >= (age_grouping_multiplier*2/3 + 1),
                                                     c("prob_sampling_old")]) # very important
  
  rw_baseline_agecat_bias_comparison[counter,] = colSums((abs(sweep(median_quantile_list_rw,
                                                                    2,
                                                                    points_df_final$mrp,
                                                                    "-")) - 
                                                            abs(sweep(median_quantile_list_baseline,
                                                                      2,
                                                                      points_df_final$mrp,
                                                                      "-"))) <= 0)/dim(median_quantile_list_baseline)[1]
  
  ar_baseline_agecat_bias_comparison[counter,] = colSums((abs(sweep(median_quantile_list_ar,
                                                                    2,
                                                                    points_df_final$mrp,
                                                                    "-")) - 
                                                            abs(sweep(median_quantile_list_baseline,
                                                                      2,
                                                                      points_df_final$mrp,
                                                                      "-"))) <= 0)/dim(median_quantile_list_baseline)[1]
  

  rw_baseline_agecat_sd_comparison[counter,] = colSums(((quantile90_list_rw - quantile10_list_rw) -
                                                          (quantile90_list_baseline - quantile10_list_baseline))<= 0)/dim(median_quantile_list_rw)[1]
  
  ar_baseline_agecat_sd_comparison[counter,] = colSums(((quantile90_list_ar - quantile10_list_ar) -
                                                          (quantile90_list_baseline - quantile10_list_baseline))<= 0)/dim(median_quantile_list_rw)[1]
  
  counter = counter + 1
}


# ar model
ar_baseline_agecat_bias_comparison_melted = melt(ar_baseline_agecat_bias_comparison)
colnames(ar_baseline_agecat_bias_comparison_melted) = c("p", "age_cat", "proportion")
ar_baseline_agecat_bias_comparison_melted$age_cat_grouped = ""

ar_baseline_agecat_bias_comparison_melted[ar_baseline_agecat_bias_comparison_melted$age_cat %in% 1:4, c("age_cat_grouped")] = "Age categories 1-4"
ar_baseline_agecat_bias_comparison_melted[ar_baseline_agecat_bias_comparison_melted$age_cat %in% 5:8, c("age_cat_grouped")] = "Age categories 5-8"
ar_baseline_agecat_bias_comparison_melted[ar_baseline_agecat_bias_comparison_melted$age_cat %in% 9:12, c("age_cat_grouped")] = "Age categories 9-12"


# rw model
rw_baseline_agecat_bias_comparison_melted = melt(rw_baseline_agecat_bias_comparison)
colnames(rw_baseline_agecat_bias_comparison_melted) = c("p", "age_cat", "proportion")
rw_baseline_agecat_bias_comparison_melted$age_cat_grouped = ""

rw_baseline_agecat_bias_comparison_melted[rw_baseline_agecat_bias_comparison_melted$age_cat %in% 1:4, c("age_cat_grouped")] = "Age categories 1-4"
rw_baseline_agecat_bias_comparison_melted[rw_baseline_agecat_bias_comparison_melted$age_cat %in% 5:8, c("age_cat_grouped")] = "Age categories 5-8"
rw_baseline_agecat_bias_comparison_melted[rw_baseline_agecat_bias_comparison_melted$age_cat %in% 9:12, c("age_cat_grouped")] = "Age categories 9-12"

rw_ar_baseline_agecat_bias_comparison_melted = rbind(cbind(ar_baseline_agecat_bias_comparison_melted, Model = "AR - Baseline priors comparison"),
                                                     cbind(rw_baseline_agecat_bias_comparison_melted, Model = "RW - Baseline priors comparison"))

mapdf_prob_sampling_old_bias = data.frame(old = sort(unique(rw_ar_baseline_agecat_bias_comparison_melted$p)), 
                                          new = prob_sampling_old)
rw_ar_baseline_agecat_bias_comparison_melted$prob_sampling_old = round(mapdf_prob_sampling_old_bias$new[match(rw_ar_baseline_agecat_bias_comparison_melted$p,
                                                                               mapdf_prob_sampling_old_bias$old)],
                                        digits = 2)


saveRDS(prob_sampling_old, 
        "prob_sampling_old.rds")

saveRDS(rw_ar_baseline_agecat_bias_comparison_melted, paste0("proportion_",
                                                             runs,
                                                             "_",
                                                             age_grouping_multiplier,
                                                             "_",
                                                             income_multiplier,
                                                             "_",
                                                             sample_size,".rds"))

png(filename = paste0("proportion_",
                      runs,
                      "_",
                      age_grouping_multiplier,
                      "_",
                      income_multiplier,
                      "_",
                      sample_size,"smaller.png"),
    width = 5400/5, height = 2700/5)

ggplot(rw_ar_baseline_agecat_bias_comparison_melted, aes(x=prob_sampling_old, y=proportion, group = age_cat)) + 
  geom_line(aes(col=as.factor(age_cat_grouped)), show.legend = TRUE, size = 1) + 
  scale_colour_manual(values=c("Age categories 1-4"="#fdcc8a",
                               "Age categories 5-8"="#b30000",
                               "Age categories 9-12"="#006d2c")) + 
  facet_wrap(.~Model,ncol=2,nrow=1,labeller = as_labeller(function(value){return(value)})) + 
  xlab("\n Probability of response for the elderly \n") +
  ylab("\n Proportion of the time that structured priors outperform \n") + 
  theme_bw() +
  theme(plot.title = element_text(size = 50/5, face = "bold"),
        axis.text=element_text(size=50*1.6/5),
        axis.title=element_text(size=50*1.6/5, face="bold",margin=200/5),
        legend.text = element_text(size=50*1.6/5),
        legend.position  = "bottom",
        legend.key.size = unit(30*1.6/5,"line"),
        legend.key.height = unit(7*1.6/5,"line"),
        legend.key = element_rect(fill = "transparent",color = "transparent"),
        legend.title = element_text(size=50*1.6/5, face="bold"),
        strip.text.x = element_text(size=50*1.6/5, face="bold", margin = margin(t=25/5,b=25/5) ),
        strip.background = element_rect(fill="transparent",color="transparent")
  ) + 
  ylim(c(.2,1)) +
  geom_hline(yintercept = 0.5,linetype=2,size=1) +
  guides(col=guide_legend(title="Age category",
                          override.aes = list(size = 10*1.6/5)))

dev.off()



# ar model sd
ar_baseline_agecat_sd_comparison_melted = melt(ar_baseline_agecat_sd_comparison)
colnames(ar_baseline_agecat_sd_comparison_melted) = c("p", "age_cat", "proportion")
ar_baseline_agecat_sd_comparison_melted$age_cat_grouped = ""

ar_baseline_agecat_sd_comparison_melted[ar_baseline_agecat_sd_comparison_melted$age_cat %in% 1:4, c("age_cat_grouped")] = "Age categories 1-4"
ar_baseline_agecat_sd_comparison_melted[ar_baseline_agecat_sd_comparison_melted$age_cat %in% 5:8, c("age_cat_grouped")] = "Age categories 5-8"
ar_baseline_agecat_sd_comparison_melted[ar_baseline_agecat_sd_comparison_melted$age_cat %in% 9:12, c("age_cat_grouped")] = "Age categories 9-12"


# rw model sd
rw_baseline_agecat_sd_comparison_melted = melt(rw_baseline_agecat_sd_comparison)
colnames(rw_baseline_agecat_sd_comparison_melted) = c("p", "age_cat", "proportion")
rw_baseline_agecat_sd_comparison_melted$age_cat_grouped = ""

rw_baseline_agecat_sd_comparison_melted[rw_baseline_agecat_sd_comparison_melted$age_cat %in% 1:4, c("age_cat_grouped")] = "Age categories 1-4"
rw_baseline_agecat_sd_comparison_melted[rw_baseline_agecat_sd_comparison_melted$age_cat %in% 5:8, c("age_cat_grouped")] = "Age categories 5-8"
rw_baseline_agecat_sd_comparison_melted[rw_baseline_agecat_sd_comparison_melted$age_cat %in% 9:12, c("age_cat_grouped")] = "Age categories 9-12"

rw_ar_baseline_agecat_sd_comparison_melted = rbind(cbind(ar_baseline_agecat_sd_comparison_melted, Model = "AR - Baseline priors comparison"),
                                                     cbind(rw_baseline_agecat_sd_comparison_melted, Model = "RW - Baseline priors comparison"))

rw_ar_baseline_agecat_sd_comparison_melted$prob_sampling_old = round(mapdf_prob_sampling_old_bias$new[match(rw_ar_baseline_agecat_sd_comparison_melted$p,
                                                                                                              mapdf_prob_sampling_old_bias$old)],
                                                                       digits = 2)


saveRDS(rw_ar_baseline_agecat_sd_comparison_melted, paste0("proportion_sd_",
                                                             runs,
                                                             "_",
                                                             age_grouping_multiplier,
                                                             "_",
                                                             income_multiplier,
                                                             "_",
                                                             sample_size,".rds"))

png(filename = paste0("proportion_sd_",
                      runs,
                      "_",
                      age_grouping_multiplier,
                      "_",
                      income_multiplier,
                      "_",
                      sample_size,"smaller.png"),
    width = 5400/5, height = 2700/5)

ggplot(rw_ar_baseline_agecat_sd_comparison_melted, aes(x=prob_sampling_old, y=proportion, group = age_cat)) + 
  geom_line(aes(col=as.factor(age_cat_grouped)), show.legend = TRUE, size = 1) + 
  scale_colour_manual(values=c("Age categories 1-4"="#fdcc8a",
                               "Age categories 5-8"="#b30000",
                               "Age categories 9-12"="#006d2c")) + 
  facet_wrap(.~Model,ncol=2,nrow=1,labeller = as_labeller(function(value){return(value)})) + 
  xlab("\n Probability of response for the elderly \n") +
  ylab("\n Proportion sd \n") + 
  theme_bw() +
  theme(plot.title = element_text(size = 50/5, face = "bold"),
        axis.text=element_text(size=50*1.6/5),
        axis.title=element_text(size=50*1.6/5, face="bold",margin=200/5),
        legend.text = element_text(size=50*1.6/5),
        legend.position  = "bottom",
        legend.key.size = unit(30*1.6/5,"line"),
        legend.key.height = unit(7*1.6/5,"line"),
        legend.key = element_rect(fill = "transparent",color = "transparent"),
        legend.title = element_text(size=50*1.6/5, face="bold"),
        strip.text.x = element_text(size=50*1.6/5, face="bold", margin = margin(t=25/5,b=25/5) ),
        strip.background = element_rect(fill="transparent",color="transparent")
  ) + 
  ylim(c(.1,1)) +
  geom_hline(yintercept = 0.5,linetype=2,size=1) +
  guides(col=guide_legend(title="Age category",
                          override.aes = list(size = 10*1.6/5)))

dev.off()

print(colMeans(rw_baseline_agecat_bias_comparison))
print(colMeans(ar_baseline_agecat_bias_comparison))

print(colMeans(rw_baseline_agecat_sd_comparison))
print(colMeans(ar_baseline_agecat_sd_comparison))









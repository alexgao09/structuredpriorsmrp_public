rm(list=ls()) 

options(bitmapType="cairo")
# bias plots for new pipeline

library(reshape2)
library(rstan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(knitr)
library(ggridges)
library(matrixStats)

# Global params -------------------------------------------

save_ridgeplots = TRUE 
sample_size = 100
runs = 200
r = 1:9/10
income_multiplier = 1 # partitions income into more categories
age_grouping_multiplier = 12 # how much we take the maximal poststratification. make sure this can divide 60
number_of_states = 51

response_tag = "binary"

# colour scheme for biasfacet
use_biascolourscheme = TRUE
biascolourscheme = c("1"="#fbb4b9",
                     "2"="#f768a1",
                     "3"="#c51b8a",
                     "4"="#7a0177",
                     "5"="#fdcc8a",
                     "6"="#fc8d59",
                     "7"="#e34a33",
                     "8"="#b30000",
                     "9"="#b2e2e2",
                     "10"="#66c2a4",
                     "11"="#2ca25f",
                     "12"="#006d2c")


# ---------------------------------------------------------
# df for 3 regression lines with their standard deviation bands : median of medians
regression_df_all200 = c()

# df for 3 regression lines with their standard deviation bands : median of means
regression_df_mean_all200 = c()

# bias for baseline model
avg_bias_mat_baseline = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the mean of bias cell-wise (mean of mean bias)
avg_bias_sd_mat_baseline = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the sd of the mean of bias cell-wise
avg_value_mat_baseline = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # posterior values of 12 cells, for each value of p

# bias for ar model
avg_bias_mat_ar = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the mean of bias cell-wise (mean of mean bias)
avg_bias_sd_mat_ar = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the sd of the mean of bias cell-wise
avg_value_mat_ar = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier)  * number_of_states) # posterior values of 12 cells, for each value of p

mean_of_medians_ar = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the mean of median bias cell-wise
mean_of_medians_sd_ar = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states)

# bias for rw model
avg_bias_mat_rw = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the mean of bias cell-wise (mean of mean bias)
avg_bias_sd_mat_rw = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the sd of the mean of bias cell-wise
avg_value_mat_rw = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # posterior values of 12 cells, for each value of p

mean_of_medians_rw = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states) # the mean of median bias cell-wise
mean_of_medians_sd_rw = matrix(0, length(r), age_grouping_multiplier * (4* income_multiplier) * number_of_states)

# df for 3 regression lines with their standard deviation bands : median of medians
regression_df = c()

# df for 3 regression lines with their standard deviation bands : median of means
regression_df_mean = c()

# df for majority vote for all three models
majorityvote_df = c()

lowerboundary = 0.05
upperboundary = 0.95

counter = 1
prob_sampling_old = rep(0, length(r))

# ---
cm_quantiles_mat_bl = matrix(0, length(r), age_grouping_multiplier)
cm_quantiles_mat_ar = matrix(0, length(r), age_grouping_multiplier)
cm_quantiles_mat_rw = matrix(0, length(r), age_grouping_multiplier)

regression_df_all200_postsd = c()
# ---

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
  
  # ---
  cm_quantiles_mat_bl[counter,] = colMeans(quantile90_list_baseline - quantile10_list_baseline)
  cm_quantiles_mat_ar[counter,] = colMeans(quantile90_list_ar - quantile10_list_ar)
  cm_quantiles_mat_rw[counter,] = colMeans(quantile90_list_rw - quantile10_list_rw)
  
  quantilediff_list_baseline_melted = melt(quantile90_list_baseline - quantile10_list_baseline)
  quantilediff_list_ar_melted = melt(quantile90_list_ar - quantile10_list_ar)
  quantilediff_list_rw_melted = melt(quantile90_list_rw - quantile10_list_rw)
  
  regression_df_all200_postsd = rbind(regression_df_all200_postsd,
                                      data.frame(cbind(quantilediff_list_baseline_melted,
                                                       type="Baseline",
                                                       p)),
                                      data.frame(cbind(quantilediff_list_ar_melted,
                                                       type="Autoregressive",
                                                       p)),
                                      data.frame(cbind(quantilediff_list_rw_melted,
                                                       type="Random walk",
                                                       p))
                                      )

  # ---
  
  # correct level order for the 4 ridgeplots below
  correct_levels = c()
  for (g in 1:age_grouping_multiplier) {
    correct_levels = c(correct_levels, 
                       seq(from = g, to = age_grouping_multiplier * 3, by = age_grouping_multiplier))
  }
  
  # Get the interpretable p which is the probability of sampling someone in ages 61 - 80
  poststrat_final_p = cbind(poststrat_final, 
                            (poststrat_final$p_response * poststrat_final$N)/sum((poststrat_final$p_response * poststrat_final$N))) # this is the probability of sampling someone in ages  61 - 80
  colnames(poststrat_final_p)[length(colnames(poststrat_final_p))] = "prob_sampling_old" 
  prob_sampling_old[counter] = sum(poststrat_final_p[as.numeric(as.character(poststrat_final_p$age_cat)) >= (age_grouping_multiplier*2/3 + 1),
                                                     c("prob_sampling_old")]) # very important
  
  # baseline   
  avg_bias_mat_baseline[counter,] = colMeans(sweep(sample_cell_estimates_median_baseline,
                                                   2,
                                                   poststrat_final_reduced$true_pref_grouped_final, 
                                                   "-"))
  
  avg_bias_sd_mat_baseline[counter,] = apply(sweep(sample_cell_estimates_median_baseline,
                                                   2,
                                                   poststrat_final_reduced$true_pref_grouped_final,
                                                   "-"),
                                             2,
                                             sd)
  
  avg_value_mat_baseline[counter,] = colMeans(sample_cell_estimates_median_baseline)
  
  # autoregressive
  avg_bias_mat_ar[counter,] = colMeans(sweep(sample_cell_estimates_median_ar,
                                             2,
                                             poststrat_final_reduced$true_pref_grouped_final, 
                                             "-"))
  
  avg_bias_sd_mat_ar[counter,] = apply(sweep(sample_cell_estimates_median_ar,
                                             2,
                                             poststrat_final_reduced$true_pref_grouped_final, "-"),
                                       2,
                                       sd)
  
  avg_value_mat_ar[counter,] = colMeans(sample_cell_estimates_median_ar)
  
  # random walk
  avg_bias_mat_rw[counter,] = colMeans(sweep(sample_cell_estimates_median_rw,
                                             2,
                                             poststrat_final_reduced$true_pref_grouped_final, 
                                             "-"))
  
  avg_bias_sd_mat_rw[counter,] = apply(sweep(sample_cell_estimates_median_rw,
                                             2,
                                             poststrat_final_reduced$true_pref_grouped_final, "-"),
                                       2,
                                       sd)
  
  avg_value_mat_rw[counter,] = colMeans(sample_cell_estimates_median_rw)
  
  # ---
  regression_df = rbind(regression_df, 
                        data.frame(cbind(1:age_grouping_multiplier, 
                              colMedians(median_quantile_list_baseline),
                              "Baseline",
                              p)),
                        data.frame(cbind(1:age_grouping_multiplier, 
                              colMedians(median_quantile_list_ar),
                              "Autoregressive",
                              p)),
                        data.frame(cbind(1:age_grouping_multiplier, 
                              colMedians(median_quantile_list_rw),
                              "Random walk",
                              p))
                        )
  # ---
  
  regression_df_mean = rbind(regression_df_mean, 
                        data.frame(cbind(1:age_grouping_multiplier, 
                                         colMedians(mean_list_baseline),
                                         "Baseline",
                                         p)),
                        data.frame(cbind(1:age_grouping_multiplier, 
                                         colMedians(mean_list_ar),
                                         "Autoregressive",
                                         p)),
                        data.frame(cbind(1:age_grouping_multiplier, 
                                         colMedians(mean_list_rw),
                                         "Random walk",
                                         p))
  )
  
  # ---
  
  median_quantile_list_baseline_melted = melt(median_quantile_list_baseline)
  median_quantile_list_ar_melted = melt(median_quantile_list_ar)
  median_quantile_list_rw_melted = melt(median_quantile_list_rw)
  
  regression_df_all200 = rbind(regression_df_all200, 
                               data.frame(cbind(median_quantile_list_baseline_melted,
                                                type="Baseline",  
                                                p)),
                               data.frame(cbind(median_quantile_list_ar_melted,
                                                type="Autoregressive",
                                                p)),
                               data.frame(cbind(median_quantile_list_rw_melted,
                                                type="Random walk",
                                                p))
  )
  
  # ---
  
  mean_list_baseline_melted = melt(mean_list_baseline)
  mean_list_ar_melted = melt(mean_list_ar)
  mean_list_rw_melted = melt(mean_list_rw)
  
  regression_df_mean_all200 = rbind(regression_df_mean_all200, 
                               data.frame(cbind(mean_list_baseline_melted,
                                                type="Baseline",  
                                                p)),
                               data.frame(cbind(mean_list_ar_melted,
                                                type="Autoregressive",
                                                p)),
                               data.frame(cbind(mean_list_rw_melted,
                                                type="Random walk",
                                                p))
  )
  
  # ---
  
  majorityvote_df = rbind(majorityvote_df, 
                          data.frame(majorityvote = model_popn_pref_mat_baseline, 
                                     type = "Baseline",
                                     p),
                          data.frame(majorityvote = model_popn_pref_mat_ar, 
                                     type = "Autoregressive",
                                     p),
                          data.frame(majorityvote = model_popn_pref_mat_rw, 
                                     type = "Random walk",
                                     p)
                          )
  
  counter = counter + 1
}

nonboundary_indices = (poststrat_final_reduced$true_pref_grouped_final > lowerboundary) & (poststrat_final_reduced$true_pref_grouped_final < upperboundary)

avg_value_mat_baseline = as.data.frame(avg_value_mat_baseline)
avg_value_mat_baseline = avg_value_mat_baseline[,nonboundary_indices]
df_val.melted_baseline = cbind( melt(avg_value_mat_baseline), rep(r, sum(nonboundary_indices)) ) 
colnames(df_val.melted_baseline) = c("cell", "value", "p_response_3")

avg_value_mat_ar = as.data.frame(avg_value_mat_ar)
avg_value_mat_ar = avg_value_mat_ar[,nonboundary_indices]
df_val.melted_ar = cbind(melt(avg_value_mat_ar), rep(r, sum(nonboundary_indices)) )
colnames(df_val.melted_ar) = c("cell", "value", "p_response_3")

avg_value_mat_rw = as.data.frame(avg_value_mat_rw)
avg_value_mat_rw = avg_value_mat_rw[,nonboundary_indices]
df_val.melted_rw = cbind(melt(avg_value_mat_rw), rep(r, sum(nonboundary_indices)) )
colnames(df_val.melted_rw) = c("cell", "value", "p_response_3")


poststrat_final_V = cbind(paste0("V", rownames(poststrat_final_reduced)),
                          poststrat_final_reduced,
                          stringsAsFactors=FALSE)
colnames(poststrat_final_V)[1] = "cell"

df_val.melted_baseline$cell = as.character(df_val.melted_baseline$cell)
df_val.melted_ar$cell = as.character(df_val.melted_ar$cell)
df_val.melted_rw$cell = as.character(df_val.melted_rw$cell)

df_val.melted_V_baseline = dplyr::inner_join(x = df_val.melted_baseline, y = poststrat_final_V, by="cell")
df_val.melted_V_unique_baseline = unique(df_val.melted_V_baseline[,c("cell", "true_pref_grouped_final")])
df_val.melted_V_baseline$true_pref_grouped_final = factor(df_val.melted_V_baseline$true_pref_grouped_final)

df_val.melted_V_ar = dplyr::inner_join(x = df_val.melted_ar, y = poststrat_final_V, by="cell")
df_val.melted_V_unique_ar = unique(df_val.melted_V_ar[,c("cell", "true_pref_grouped_final")])
df_val.melted_V_ar$true_pref_grouped_final = factor(df_val.melted_V_ar$true_pref_grouped_final)

df_val.melted_V_rw = dplyr::inner_join(x = df_val.melted_rw, y = poststrat_final_V, by="cell")
df_val.melted_V_unique_rw = unique(df_val.melted_V_rw[,c("cell", "true_pref_grouped_final")])
df_val.melted_V_rw$true_pref_grouped_final = factor(df_val.melted_V_rw$true_pref_grouped_final)


d_temp_baseline = ggplot(df_val.melted_V_baseline, aes(x=p_response_3, y=value, group=cell)) +
  geom_line(aes(col=as.factor(age_cat)), show.legend = TRUE, size = 1) +
  geom_hline(data = poststrat_final_V, mapping = aes(yintercept = true_pref_grouped_final, color = as.factor(age_cat)), linetype = "dashed", show.legend = FALSE) +
  xlab("Probability of response for oldest age group") +
  ylab("Baseline Avg. of avg. of PLD colour-coded by truth") + 
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position="top",
        legend.title = element_blank())

d_temp_ar = ggplot(df_val.melted_V_ar, aes(x=p_response_3, y=value, group=cell)) +
  geom_line(aes(col=as.factor(age_cat)), show.legend = TRUE, size = 1) +
  geom_hline(data = poststrat_final_V, mapping = aes(yintercept = true_pref_grouped_final, color = as.factor(age_cat)), linetype = "dashed", show.legend = FALSE) +
  xlab("Probability of response for oldest age group") +
  ylab("AR Avg. of avg. of PLD colour-coded by truth") + 
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position="top",
        legend.title = element_blank())

d_temp_rw = ggplot(df_val.melted_V_rw, aes(x=p_response_3, y=value, group=cell)) +
  geom_line(aes(col=as.factor(age_cat)), show.legend = TRUE, size = 1) +
  geom_hline(data = poststrat_final_V, mapping = aes(yintercept = true_pref_grouped_final, color = as.factor(age_cat)), linetype = "dashed", show.legend = FALSE) +
  xlab("Probability of response for oldest age group") +
  ylab("RW Avg. of avg. of PLD colour-coded by truth") + 
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position="top",
        legend.title = element_blank())



# BIAS -------

avg_bias_mat_baseline = as.data.frame(avg_bias_mat_baseline)
avg_bias_mat_baseline = avg_bias_mat_baseline[,nonboundary_indices]
df_bias.melted_baseline = cbind( melt(avg_bias_mat_baseline), rep(r, sum(nonboundary_indices)) ) 
colnames(df_bias.melted_baseline) = c("cell", "value", "p_response_3")

avg_bias_mat_ar = as.data.frame(avg_bias_mat_ar)
avg_bias_mat_ar = avg_bias_mat_ar[,nonboundary_indices]
df_bias.melted_ar = cbind(melt(avg_bias_mat_ar), rep(r, sum(nonboundary_indices)) )
colnames(df_bias.melted_ar) = c("cell", "value", "p_response_3")

avg_bias_mat_rw = as.data.frame(avg_bias_mat_rw)
avg_bias_mat_rw = avg_bias_mat_rw[,nonboundary_indices]
df_bias.melted_rw = cbind(melt(avg_bias_mat_rw), rep(r, sum(nonboundary_indices)) )
colnames(df_bias.melted_rw) = c("cell", "value", "p_response_3")

df_bias.melted_baseline$cell = as.character(df_bias.melted_baseline$cell)
df_bias.melted_ar$cell = as.character(df_bias.melted_ar$cell)
df_bias.melted_rw$cell = as.character(df_bias.melted_rw$cell)

df_bias.melted_V_baseline = dplyr::inner_join(x = df_bias.melted_baseline, y = poststrat_final_V, by="cell")
df_bias.melted_V_unique_baseline = unique(df_bias.melted_V_baseline[,c("cell", "true_pref_grouped_final")])
df_bias.melted_V_baseline$true_pref_grouped_final = factor(df_bias.melted_V_baseline$true_pref_grouped_final)

df_bias.melted_V_ar = dplyr::inner_join(x = df_bias.melted_ar, y = poststrat_final_V, by="cell")
df_bias.melted_V_unique_ar = unique(df_bias.melted_V_ar[,c("cell", "true_pref_grouped_final")])
df_bias.melted_V_ar$true_pref_grouped_final = factor(df_bias.melted_V_ar$true_pref_grouped_final)

df_bias.melted_V_rw = dplyr::inner_join(x = df_bias.melted_rw, y = poststrat_final_V, by="cell")
df_bias.melted_V_unique_rw = unique(df_bias.melted_V_rw[,c("cell", "true_pref_grouped_final")])
df_bias.melted_V_rw$true_pref_grouped_final = factor(df_bias.melted_V_rw$true_pref_grouped_final)


bias_facet_df = rbind(cbind(df_bias.melted_V_baseline, Model = "Baseline"),
                      cbind(df_bias.melted_V_ar, Model = "Autoregressive"),
                      cbind(df_bias.melted_V_rw, Model = "Random walk"))

mapdf_prob_sampling_old_bias = data.frame(old = sort(unique(bias_facet_df$p_response_3)), 
                                     new = prob_sampling_old)
bias_facet_df$prob_sampling_old = round(mapdf_prob_sampling_old_bias$new[match(bias_facet_df$p_response_3,
                                                                               mapdf_prob_sampling_old_bias$old)],
                                               digits = 2)

bias_facet_df$age_cat = factor(bias_facet_df$age_cat)

saveRDS(bias_facet_df, paste0("biasfacet_",
                              runs,
                              "_",
                              age_grouping_multiplier,
                              "_",
                              income_multiplier,
                              "_",
                              sample_size,".rds"))

png(filename = paste0("biasfacet_",
                      runs,
                      "_",
                      age_grouping_multiplier,
                      "_",
                      income_multiplier,
                      "_",
                      sample_size,".png"),
    width = 4800, height = 2400)

if (use_biascolourscheme==TRUE) {

    plot(
      ggplot(bias_facet_df, aes(x=prob_sampling_old, y=value, group=cell)) + 
        geom_line(aes(col=as.factor(age_cat)), show.legend = TRUE, size = 0.5) +
        facet_wrap(. ~ Model,ncol=3,nrow=1, labeller = as_labeller(
          function(value) {
            return(value)  # Lets you change the facet labels
          })
        ) +  
        xlab("\n Probability of response for the elderly \n") +
        ylab("\n Average Bias \n") + 
        scale_colour_manual(values=c("1"="#fbb4b9",
                                   "2"="#f768a1",
                                   "3"="#c51b8a",
                                   "4"="#7a0177",
                                   "5"="#fdcc8a",
                                   "6"="#fc8d59",
                                   "7"="#e34a33",
                                   "8"="#b30000",
                                   "9"="#b2e2e2",
                                   "10"="#66c2a4",
                                   "11"="#2ca25f",
                                   "12"="#006d2c")
                            ) +
        theme_bw() +
        theme(plot.title = element_text(size = 50, face = "bold"),
              axis.text=element_text(size=35*1.6),
              axis.title=element_text(size=50*1.6, face="bold",margin=200),
              legend.text = element_text(size=50*1.6),
              legend.position  = "bottom",
              legend.key.size = unit(30*1.6,"line"),
              legend.key.height = unit(7*1.6,"line"),
              legend.key = element_rect(fill = "transparent",color = "transparent"),
              legend.title = element_text(size=50*1.6, face="bold"),
              strip.text.x = element_text(size=50*1.6, face="bold", margin = margin(t=25,b=25) ),
              strip.background = element_rect(fill="transparent",color="transparent")
              ) + 
        guides(col=guide_legend(title="Age category",
                                override.aes = list(size = 10*1.6)))
  
    )
  
}else {
  plot(
    ggplot(bias_facet_df, aes(x=prob_sampling_old, y=value, group=cell)) + 
      geom_line(aes(col=as.factor(age_cat)), show.legend = TRUE, size = 0.5) +
      facet_wrap(. ~ Model,ncol=3,nrow=1, labeller = as_labeller(
        function(value) {
          return(value)  # Lets you change the facet labels
        })
      ) +  
      xlab("\n Probability of response for the elderly \n") +
      ylab("\n Average Bias \n") + 
      theme_bw() +
      theme(plot.title = element_text(size = 50, face = "bold"),
            axis.text=element_text(size=35*1.6),
            axis.title=element_text(size=50*1.6, face="bold",margin=200),
            legend.text = element_text(size=50*1.6),
            legend.position  = "bottom",
            legend.key.size = unit(30*1.6,"line"),
            legend.key.height = unit(7*1.6,"line"),
            legend.key = element_rect(fill = "transparent",color = "transparent"),
            legend.title = element_text(size=50*1.6, face="bold"),
            strip.text.x = element_text(size=50*1.6, face="bold", margin = margin(t=25,b=25) ),
            strip.background = element_rect(fill="transparent",color="transparent")) + 
      guides(col=guide_legend(title="Age category",
                              override.aes = list(size = 10*1.6)))
  )
}
  
dev.off()

  
  # all medians of posteriors
  
regression_df_all200$Var2 = factor(regression_df_all200$Var2)
colnames(regression_df_all200) = c("n", "age_cat", "posterior_medians", "Model", "p")

mapdf_prob_sampling_old = data.frame(old = sort(unique(regression_df_all200$p)), 
                                     new = prob_sampling_old)
regression_df_all200$prob_sampling_old = round(mapdf_prob_sampling_old$new[match(regression_df_all200$p,
                                                                                 mapdf_prob_sampling_old$old)],
                                               digits = 2)

regression_df_all200$age_cat = as.numeric(as.character(regression_df_all200$age_cat))

# save regression_df_all200
saveRDS(regression_df_all200, 
        paste0("allmedians_facet",
               "_",
               runs,
               "_",
               age_grouping_multiplier,
               "_",
               sample_size,
               "_",
               ".rds"))

saveRDS(points_df_final[,c("age_cat", "mrp", "Type")], 
        "points_df_final.rds")

png(paste0("allmedians_facet",
           "_",
           runs,
           "_",
           age_grouping_multiplier,
           "_",
           sample_size,
           "_",
           ".png"),
    width=3000, height=2400)


plot(
  ggplot(regression_df_all200, aes(x=age_cat, y=posterior_medians, color=Model)) + #geom_point(size=5)+
    geom_rect(regression_df_all200, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                  xmax = Inf,
                                                  ymin = -Inf,
                                                  ymax = Inf),
              alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
    geom_jitter(alpha=0.1) + 
    geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
    facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=3, labeller = as_labeller(
      function(value) {
        return(value)  # Lets you change the facet labels
      })
    ) +
    geom_point(aes(x = age_cat, y = mrp),
               size=7,
               colour = "black",
               fill="black",
               show.legend = F,
               data = points_df_final,
               inherit.aes = F) +
    xlab("\n Age Category \n") +
    ylab(paste("\n Median of", runs,"Posteriors \n")) +
    scale_x_continuous(breaks= scales::pretty_breaks()) + 
    theme_bw() +
    theme(plot.title = element_text(size = 50, face = "bold"),
          axis.text=element_text(size=50),
          axis.title=element_text(size=50, face="bold",margin=200),
          legend.text = element_text(size=50),
          legend.position  = "bottom",
          legend.key.size = unit(30,"line"),
          legend.key.height = unit(7,"line"),
          legend.key = element_rect(fill = "transparent",color = "transparent"),
          legend.title = element_text(size=50, face="bold"),
          strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
          strip.background = element_rect(fill="transparent",color="transparent")) +
    scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
    guides(color = guide_legend(override.aes = list(size=20)))
  
)

dev.off()
  

# plot of the posterior standard deviations 
regression_df_all200_postsd$Var2 = factor(regression_df_all200_postsd$Var2)
colnames(regression_df_all200_postsd) = c("n", "age_cat", "quantildiff90_10", "Model", "p")

#
mapdf_prob_sampling_old_postsd = data.frame(old = sort(unique(regression_df_all200_postsd$p)), 
                                     new = prob_sampling_old)
regression_df_all200_postsd$prob_sampling_old = round(mapdf_prob_sampling_old_postsd$new[match(regression_df_all200_postsd$p,
                                                                                               mapdf_prob_sampling_old_postsd$old)],
                                                      digits = 2)

regression_df_all200_postsd$age_cat = as.numeric(as.character(regression_df_all200_postsd$age_cat))
#

saveRDS(regression_df_all200_postsd,
        paste0("allquantilediff_facet",
               "_",
               runs,
               "_",
               age_grouping_multiplier,
               "_",
               sample_size,
               "_",
               ".rds"))

png(paste0("allquantilediff_facet",
           "_",
           runs,
           "_",
           age_grouping_multiplier,
           "_",
           sample_size,
           "_",
           ".png"),
    width=3000, height=2400)

plot(
  ggplot(regression_df_all200_postsd, aes(x=age_cat, y=quantildiff90_10, color=Model)) + #geom_point(size=5)+
    geom_rect(regression_df_all200_postsd, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                  xmax = Inf,
                                                  ymin = -Inf,
                                                  ymax = Inf),
              alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
    geom_jitter(alpha=0.1) + 
    geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
    facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=3, labeller = as_labeller(
      function(value) {
        return(value)  # Lets you change the facet labels
      })
    ) +
    xlab("\n Age Category \n") +
    ylab(paste("\n 90th - 10th quantile of", runs,"Posteriors \n")) +
    scale_x_continuous(breaks= scales::pretty_breaks()) + 
    theme_bw() +
    theme(plot.title = element_text(size = 50, face = "bold"),
          axis.text=element_text(size=50),
          axis.title=element_text(size=50, face="bold",margin=200),
          legend.text = element_text(size=50),
          legend.position  = "bottom",
          legend.key.size = unit(30,"line"),
          legend.key.height = unit(7,"line"),
          legend.key = element_rect(fill = "transparent",color = "transparent"),
          legend.title = element_text(size=50, face="bold"),
          strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
          strip.background = element_rect(fill="transparent",color="transparent")) +
    scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
    guides(color = guide_legend(override.aes = list(size=20)))
  
)

dev.off()

# --------------------------------------------------------------------------

print("Warnings:")
print(warnings())

    
    

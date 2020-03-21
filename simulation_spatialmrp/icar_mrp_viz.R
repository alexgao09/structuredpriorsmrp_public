# this script visualizes results from icar_mrp_simulation_bym2.R

rm(list=ls())
gc()

library(ggplot2)
library(dplyr)
library(sp)

sample_size = 500
runs = 200
r = 1:9/10

poverty_poststrat = readRDS("poverty_poststrat.rds")
us_pumas = readRDS("us_pumas.rds")
ma_adj_matrix = readRDS("ma_adj_matrix.rds")
ma_adj_matrix_sparse = readRDS("ma_adj_matrix_sparse.rds")

num_education = length(unique(poverty_poststrat$education)) # number of strata for education
num_race = length(unique(poverty_poststrat$race_x)) # number of strata for race
num_puma = length(unique(poverty_poststrat$PUMA.x))


# ----
avg_bias_mat_puma_icar = matrix(0, length(r), num_puma)
avg_bias_mat_puma_iid = matrix(0, length(r), num_puma)

# bias of all strata
avg_bias_mat_icar = matrix(0, length(r), num_education * num_race * num_puma)
avg_bias_mat_iid = matrix(0, length(r), num_education * num_race * num_puma)

# avg posterior width for 1872 poststratification cells
avg_posterior_cell_width_icar = matrix(0, length(r), num_education * num_race * num_puma)
avg_posterior_cell_width_iid = matrix(0, length(r), num_education * num_race * num_puma)

# avg posterior width for all 52 puma
avg_posterior_puma_width_icar = matrix(0, length(r), 52)
avg_posterior_puma_width_iid = matrix(0, length(r), 52)
# ----

# df for majoityvote for icar and iid prior
majorityvote_df = c()


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
                                      2)# very important
  
  avg_bias_mat_icar[counter_,] = colMeans(sweep(sample_cell_estimates_median_icar,
                                               2,
                                               poststrat_final_joined_p_ID$true_pref,
                                               "-"))
  
  avg_bias_mat_iid[counter_,] = colMeans(sweep(sample_cell_estimates_median_iid,
                                               2,
                                               poststrat_final_joined_p_ID$true_pref,
                                               "-"))
  

  avg_posterior_cell_width_icar[counter_,] = colMeans(sample_cell_estimates_width_icar)
  avg_posterior_cell_width_iid[counter_,] = colMeans(sample_cell_estimates_width_iid)
  
  avg_posterior_puma_width_icar[counter_,] = colMeans(quantile90_puma_icar - quantile10_puma_icar)
  avg_posterior_puma_width_iid[counter_,] = colMeans(quantile90_puma_iid - quantile10_puma_iid)
  
  
  avg_bias_mat_puma_icar[counter_,] = colMeans(sweep(median_puma_icar, 2, poststrat_puma_pref$puma_ps, "-"))
  avg_bias_mat_puma_iid[counter_,] = colMeans(sweep(median_puma_iid, 2, poststrat_puma_pref$puma_ps, "-"))
  
  # majorityvote -----------------------
  
  majorityvote_df = rbind(majorityvote_df, 
                          data.frame(majorityvote = model_popn_pref_mat_icar, 
                                     type = "ICAR",
                                     prob_sampling_old[counter_]),
                          data.frame(majorityvote = model_popn_pref_mat_iid, 
                                     type = "IID",
                                     prob_sampling_old[counter_])
  )
  
  counter_ = counter_ + 1
}

saveRDS(prob_sampling_old, "prob_sampling_old.rds")
# this being negative is good
summary(colMeans(abs(avg_bias_mat_icar) - 
                   abs(avg_bias_mat_iid))) * 100

# this being negative is good
summary(colMeans(abs(sweep(median_puma_icar,2,poststrat_puma_pref$puma_ps,"-")) - 
          abs(sweep(median_puma_iid,2,poststrat_puma_pref$puma_ps,"-")))) * 100

# this being negative is good
summary(colMeans(avg_posterior_cell_width_icar - avg_posterior_cell_width_iid))

# ---
poststrat_final_V = cbind(paste0("V", rownames(poststrat_final_joined)),
                          poststrat_final_joined,
                          stringsAsFactors=FALSE)
colnames(poststrat_final_V)[1] = "cell"
# ---


avg_bias_mat_icar = as.data.frame(avg_bias_mat_icar)
avg_bias_mat_iid = as.data.frame(avg_bias_mat_iid)

df_bias.melted_icar = cbind(reshape2::melt(avg_bias_mat_icar),
                            rep(r, dim(poststrat_final_joined_p_ID)[1]))
colnames(df_bias.melted_icar) = c("cell", "value", "p")
df_bias.melted_icar$cell = as.character(df_bias.melted_icar$cell)

df_bias.melted_V_icar = dplyr::inner_join(x = df_bias.melted_icar,
                                          y = poststrat_final_V,
                                          by ="cell")



df_bias.melted_iid = cbind(reshape2::melt(avg_bias_mat_iid),
                            rep(r, dim(poststrat_final_joined_p_ID)[1]))
colnames(df_bias.melted_iid) = c("cell", "value", "p")
df_bias.melted_iid$cell = as.character(df_bias.melted_iid$cell)

df_bias.melted_V_iid = dplyr::inner_join(x = df_bias.melted_iid,
                                          y = poststrat_final_V,
                                          by ="cell")

bias_facet_df = rbind(cbind(df_bias.melted_V_icar, Model = "ICAR"),
                      cbind(df_bias.melted_V_iid, Model = "IID"))

bias_facet_df$sampling_group = ifelse(bias_facet_df$puma_area %in% puma_overundersample_index, "Group 1", "Group 2")

png(filename = paste0("spatialbiasfacet_",
                      runs,
                      "_",
                      sample_size,
                      "_.png"),
    width = 4800, height = 2400)

ggplot(bias_facet_df,
       aes(x=p, y=value, group = cell)) + 
  geom_line(aes(col=as.factor(sampling_group)),show.legend = TRUE, size = 0.5) +
  facet_wrap(. ~ Model,ncol=2,nrow=1, labeller = as_labeller(
    function(value) {
      return(value)  # Lets you change the facet labels
    })
  ) +
  xlab("\n Probability of response for Group 1 \n") +
  ylab("\n Average Bias of Posterior Medians \n") + 
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
  guides(col=guide_legend(title="Sampling Group",
                          override.aes = list(size = 10*1.6)))

dev.off()

saveRDS(bias_facet_df, paste0("biasfacetspatial_",
                              runs,
                              "_",
                              sample_size,".rds"))

# Average bias of posterior medians for 52 PUMAs --------------------------------------------
avg_bias_mat_puma_icar = as.data.frame(avg_bias_mat_puma_icar)
avg_bias_mat_puma_iid = as.data.frame(avg_bias_mat_puma_iid)


colnames(avg_bias_mat_puma_icar) = levels_puma

df_bias_puma.melted_icar = cbind(reshape2::melt(avg_bias_mat_puma_icar),
                            rep(r, num_puma))
colnames(df_bias_puma.melted_icar) = c("puma_area", "value", "p")
df_bias_puma.melted_icar$puma_area = as.numeric(as.character(df_bias_puma.melted_icar$puma_area))

df_bias_puma.melted_icar$sampling_group = ifelse(df_bias_puma.melted_icar$puma_area %in% puma_overundersample_index, 
                                                 "Near Boston",
                                                 "Away from Boston")


colnames(avg_bias_mat_puma_iid) = levels_puma

df_bias_puma.melted_iid = cbind(reshape2::melt(avg_bias_mat_puma_iid),
                                 rep(r, num_puma))
colnames(df_bias_puma.melted_iid) = c("puma_area", "value", "p")
df_bias_puma.melted_iid$puma_area = as.numeric(as.character(df_bias_puma.melted_iid$puma_area))

df_bias_puma.melted_iid$sampling_group = ifelse(df_bias_puma.melted_iid$puma_area %in% puma_overundersample_index, 
                                                 "Near Boston",
                                                 "Away from Boston")


bias_facet_puma_df = rbind(cbind(df_bias_puma.melted_icar, Model = "ICAR"),
                      cbind(df_bias_puma.melted_iid, Model = "IID"))

bias_facet_puma_df$sampling_group = factor(bias_facet_puma_df$sampling_group,
                                           levels = c("Near Boston", "Away from Boston"))

png(filename = paste0("spatialbiasfacet_puma_",
                      runs,
                      "_",
                      sample_size,
                      "_.png"),
    width = 4800, height = 2400)

ggplot(bias_facet_puma_df,
       aes(x=p, y=value, group = puma_area)) + 
  geom_line(aes(col=as.factor(sampling_group)),show.legend = TRUE, size = 1.5) +
  facet_wrap(. ~ Model,ncol=2,nrow=1, labeller = as_labeller(
    function(value) {
      return(value)  # Lets you change the facet labels
    })
  ) +
  xlab("\n Probability of response for Group 1 \n") +
  ylab("\n Average Bias of Posterior Medians for PUMA \n") + 
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
  guides(col=guide_legend(title="Sampling Group",
                          override.aes = list(size = 10*1.6)))


dev.off()

saveRDS(bias_facet_puma_df, paste0("biasfacetspatialpuma_",
                              runs,
                              "_",
                              sample_size,".rds"))

# Average posterior sd width for 1872 poststrat cell posteriors ----------------------------

avg_posterior_cell_width_icar = as.data.frame(avg_posterior_cell_width_icar)
avg_posterior_cell_width_iid = as.data.frame(avg_posterior_cell_width_iid)


df_avg_posterior_cell_width_icar = cbind(reshape2::melt(avg_posterior_cell_width_icar),
                                         rep(r, dim(poststrat_final_joined_p_ID)[1]))

colnames(df_avg_posterior_cell_width_icar) = c("cell", "average_posterior_width", "p")
df_avg_posterior_cell_width_icar$cell = as.character(df_avg_posterior_cell_width_icar$cell)

df_avg_posterior_cell_width_icar_V = dplyr::inner_join(x = df_avg_posterior_cell_width_icar,
                                                       y = poststrat_final_V,
                                                       by = "cell")


df_avg_posterior_cell_width_iid = cbind(reshape2::melt(avg_posterior_cell_width_iid),
                                         rep(r, dim(poststrat_final_joined_p_ID)[1]))

colnames(df_avg_posterior_cell_width_iid) = c("cell", "average_posterior_width", "p")
df_avg_posterior_cell_width_iid$cell = as.character(df_avg_posterior_cell_width_iid$cell)

df_avg_posterior_cell_width_iid_V = dplyr::inner_join(x = df_avg_posterior_cell_width_iid,
                                                       y = poststrat_final_V,
                                                       by = "cell")
       
posterior_cell_width_facet_df = rbind(cbind(df_avg_posterior_cell_width_iid_V, Model = "ICAR"),
                                      cbind(df_avg_posterior_cell_width_iid_V, Model = "IID"))

posterior_cell_width_facet_df$sampling_group = ifelse(posterior_cell_width_facet_df$puma_area %in% puma_overundersample_index, "Group 1", "Group 2")

posterior_cell_width_facet_df$sampling_group[posterior_cell_width_facet_df$sampling_group=="Group 1"] = "Near Boston"
posterior_cell_width_facet_df$sampling_group[posterior_cell_width_facet_df$sampling_group=="Group 2"] = "Away from Boston"
posterior_cell_width_facet_df$sampling_group = factor(posterior_cell_width_facet_df$sampling_group,
                                                      levels = c("Near Boston", "Away from Boston"))

png(filename = paste0("spatialbiasfacet_width_",
                      runs,
                      "_",
                      sample_size,
                      "_.png"),
    width = 4800, height = 2400)

ggplot(posterior_cell_width_facet_df,
       aes(x=p, y=average_posterior_width, group = cell)) + 
  geom_line(aes(col=as.factor(sampling_group)),show.legend = TRUE, size = 0.5) +
  facet_wrap(. ~ Model,ncol=2,nrow=1, labeller = as_labeller(
    function(value) {
      return(value)  # Lets you change the facet labels
    })
  ) +
  xlab("\n Probability of sampling for cluster of PUMA near Boston \n") +
  ylab(paste0("\n 90th - 10th quantile of ", runs, " Posteriors \n")) + 
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
  guides(col=guide_legend(title="Sampling Group",
                          override.aes = list(size = 10*1.6)))

dev.off()

saveRDS(posterior_cell_width_facet_df, paste0("spatialbiasfacet_width_",
                                              runs,
                                              "_",
                                              sample_size,
                                              ".rds"))


# Average posterior sd of 52 PUMA -------------------------------------------------------------

avg_posterior_puma_width_icar = as.data.frame(avg_posterior_puma_width_icar)
avg_posterior_puma_width_iid = as.data.frame(avg_posterior_puma_width_iid)

colnames(avg_posterior_puma_width_icar) = levels_puma
colnames(avg_posterior_puma_width_iid) = levels_puma

df_avg_posterior_puma_width_icar = cbind(reshape2::melt(avg_posterior_puma_width_icar),
                                         rep(r, 52))

colnames(df_avg_posterior_puma_width_icar) = c("puma", "average_posterior_width", "p")
df_avg_posterior_puma_width_icar$puma = as.character(df_avg_posterior_puma_width_icar$puma)


df_avg_posterior_puma_width_iid = cbind(reshape2::melt(avg_posterior_puma_width_iid),
                                        rep(r, 52))

colnames(df_avg_posterior_puma_width_iid) = c("puma", "average_posterior_width", "p")
df_avg_posterior_puma_width_iid$puma = as.character(df_avg_posterior_puma_width_iid$puma)

posterior_puma_width_facet_df = rbind(cbind(df_avg_posterior_puma_width_icar, Model = "ICAR"),
                                      cbind(df_avg_posterior_puma_width_iid, Model = "IID"))

posterior_puma_width_facet_df$sampling_group = ifelse(posterior_puma_width_facet_df$puma %in% puma_overundersample_index, "Group 1", "Group 2")

posterior_puma_width_facet_df$sampling_group[posterior_puma_width_facet_df$sampling_group=="Group 1"] = "Near Boston"
posterior_puma_width_facet_df$sampling_group[posterior_puma_width_facet_df$sampling_group=="Group 2"] = "Away from Boston"
posterior_puma_width_facet_df$sampling_group = factor(posterior_puma_width_facet_df$sampling_group,
                                                      levels = c("Near Boston", "Away from Boston"))

png(filename = paste0("spatialbiasfacet_puma_width_",
                      runs,
                      "_",
                      sample_size,
                      "_.png"),
    width = 4800, height = 2400)

ggplot(posterior_puma_width_facet_df,
       aes(x=p, y=average_posterior_width, group = puma)) + 
  geom_line(aes(col=as.factor(sampling_group)),show.legend = TRUE, size = 1.5) +
  facet_wrap(. ~ Model,ncol=2,nrow=1, labeller = as_labeller(
    function(value) {
      return(value)  # Lets you change the facet labels
    })
  ) +
  xlab("\n Probability of sampling for cluster of PUMA near Boston \n") +
  ylab(paste0("\n 90th - 10th quantile of ", runs, " Posteriors \n")) + 
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
  guides(col=guide_legend(title="Sampling Group",
                          override.aes = list(size = 10*1.6)))

dev.off()

saveRDS(posterior_puma_width_facet_df, paste0("spatialbiasfacet_puma_width_",
                                              runs,
                                              "_",
                                              sample_size,
                                              ".rds"))


# heatmap plot of average posterior medians for 52 pumas, across 9 probability indices --------

us_pumas$ICAR_one = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                            y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[1,])),
                            by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_one = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                               y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[1,])),
                               by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_two = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                               y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[2,])),
                               by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_two = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                              y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[2,])),
                              by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_three = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                               y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[3,])),
                               by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_three = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                              y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[3,])),
                              by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_four = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                 y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[4,])),
                                 by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_four = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[4,])),
                                by = c("PUMACE10" = "puma_code"))[,2]
# ---
us_pumas$ICAR_five = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                            y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[5,])),
                            by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_five = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                           y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[5,])),
                           by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_six = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                 y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[6,])),
                                 by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_six = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[6,])),
                                by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_seven = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                               y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[7,])),
                               by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_seven = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                              y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[7,])),
                              by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_eight = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                 y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[8,])),
                                 by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_eight = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[8,])),
                                by = c("PUMACE10" = "puma_code"))[,2]

# ---
us_pumas$ICAR_nine = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                             y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_icar[9,])),
                             by = c("PUMACE10" = "puma_code"))[,2]

us_pumas$IID_nine = inner_join(x = data.frame(PUMACE10 = us_pumas$PUMACE10),
                                y = data.frame(puma_code = levels_puma, postmedian_bias = as.numeric(avg_bias_mat_puma_iid[9,])),
                                by = c("PUMACE10" = "puma_code"))[,2]

# https://edzer.github.io/sp/
png(filename = paste0("avg_bias_postmedian_puma",
                      runs,
                      "_",
                      sample_size,
                      "_.png"),
    width = 800, height = 2400)

spplot(us_pumas[,c("ICAR_one","IID_one",
                   "ICAR_two","IID_two",
                   "ICAR_three","IID_three",
                   "ICAR_four","IID_four",
                   "ICAR_five", "IID_five",
                   "ICAR_six","IID_six",
                   "ICAR_seven","IID_seven",
                   "ICAR_eight","IID_eight",
                   "ICAR_nine", "IID_nine")],
       
       names.attr = c("ICAR 0.05","IID 0.05",
                      "ICAR 0.10","IID 0.10",
                      "ICAR 0.17","IID 0.17",
                      "ICAR 0.24","IID 0.24",
                      "ICAR 0.32", "IID 0.32",
                      "ICAR 0.41","IID 0.41",
                      "ICAR 0.52","IID 0.52",
                      "ICAR 0.65","IID 0.65",
                      "ICAR 0.81", "IID 0.81"),
       
       #pretty = TRUE
       at = seq(-round(max(abs(as.data.frame(us_pumas[,c("ICAR_one","IID_one",
                                                    "ICAR_two","IID_two",
                                                    "ICAR_three","IID_three",
                                                    "ICAR_four","IID_four",
                                                    "ICAR_five","IID_five",
                                                    "ICAR_six","IID_six",
                                                    "ICAR_seven","IID_seven",
                                                    "ICAR_eight","IID_eight",
                                                    "ICAR_nine", "IID_nine")])))+0.01, 2),
                round(max(abs(as.data.frame(us_pumas[,c("ICAR_one","IID_one",
                                                    "ICAR_two","IID_two",
                                                    "ICAR_three","IID_three",
                                                    "ICAR_four","IID_four",
                                                    "ICAR_five","IID_five",
                                                    "ICAR_six","IID_six",
                                                    "ICAR_seven","IID_seven",
                                                    "ICAR_eight","IID_eight",
                                                    "ICAR_nine", "IID_nine")])))+0.01, 2), 
                length.out = 12),
       
      col.regions = RColorBrewer::brewer.pal(11,"RdBu")
) 
  

dev.off()

png(filename = paste0("avg_bias_postmedian_puma_forpaper",
                      runs,
                      "_",
                      sample_size,
                      "_.png"),
    width = 800, height = 800)

spplot(us_pumas[,c("ICAR_one","IID_one",
                   "ICAR_five", "IID_five",
                   "ICAR_nine", "IID_nine")],
       
       names.attr = c("BYM2 prior, probability of sampling near Boston = 0.05","IID prior, probability of sampling near Boston = 0.05",
                      "BYM2 prior, probability of sampling near Boston = 0.32", "IID prior, probability of sampling near Boston = 0.32",
                      "BYM2 prior, probability of sampling near Boston = 0.81", "IID prior, probability of sampling near Boston = 0.81"),
       
       #pretty = TRUE
       at = seq(-round(max(abs(as.data.frame(us_pumas[,c("ICAR_one","IID_one",
                                                         "ICAR_five","IID_five",
                                                         "ICAR_nine", "IID_nine")])))+0.01, 2),
                round(max(abs(as.data.frame(us_pumas[,c("ICAR_one","IID_one",
                                                        "ICAR_five","IID_five",
                                                        "ICAR_nine", "IID_nine")])))+0.01, 2), 
                length.out = 12),
       
       col.regions = RColorBrewer::brewer.pal(11,"RdBu"),
       
       colorkey = list(space="bottom")
) 


dev.off()

us_pumas$puma_ps_group1 = rep(0, 52)
us_pumas$puma_ps_group1[us_pumas$PUMACE10 %in% puma_overundersample_index] = us_pumas$puma_ps[us_pumas$PUMACE10 %in% puma_overundersample_index]

png(filename = paste0("truth_ps_spatial.png"),
    width = 1000, height = 1000)

spplot(us_pumas[,],c("puma_ps", "puma_ps_group1"),
       at = seq(0, max(as.data.frame(us_pumas[,c("puma_ps")]))+0.01, length.out=40),
       colorkey = list(space="bottom"),
       col="transparent",
       names.attr = c("Poststratified true preferences", "17 PUMA over/undersampled near Boston"))

dev.off()

# majorityvote plot -------------------------------------------------

colnames(majorityvote_df)[2] = "Model"
colnames(majorityvote_df)[3] = "p"


png(filename = paste0("majorityvote_",
                      runs,
                      "_",
                      sample_size,
                      ".png"),
    width = 4000, height = 2400)

  ggplot(majorityvote_df, aes(x = majorityvote, y = Model, fill = Model)) +
    ggridges::geom_density_ridges2(alpha=0.7, 
                         quantile_lines = TRUE, 
                         quantiles = c(0.1, 0.5, 0.9), 
                         vline_size = 0.5,
                         vline_color = "black",
                         scale = 1) + 
    xlab("\n Probability of voting yes \n") +
    ylab("\n Model \n") +
    facet_wrap(. ~ p,ncol=3,nrow=3, labeller = as_labeller(
      function(value) {
        return(value)  # Lets you change the facet labels
      })
    ) +
    geom_point(aes(x = as.numeric(true_popn_pref),
                   y = Model),
               size=7,
               colour = "black",
               show.legend = F,
               inherit.aes = F) +
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
    scale_fill_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
    guides(color = guide_legend(override.aes = list(size=20)))

  dev.off()






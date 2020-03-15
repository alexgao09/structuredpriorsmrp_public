# run this script after threemodelwriteup_v3_posteriorvariance.R to get the main plots in the paper
# this script corresponds to n=100, 500
rm(list=ls()) # clear previous workspace
gc()

library(ggplot2)

bias_facet_df_100 = readRDS("biasfacet_200_12_1_100.rds")
bias_facet_df_500 = readRDS("biasfacet_200_12_1_500.rds")

biasfacet_pngname = "biasfacet_200_12_100_500.png"
allmediansfacet_pngname = "allmediansfacet_200_12_100_500.png"
allquantilediff_pngname = "allquantilediff_facet_200_12_100_500.png"

runs = 200
age_grouping_multiplier = 12
  
# ---------------------------
bias_facet_df_100$age_cat_grouped = ifelse(bias_facet_df_100$age_cat %in% c(1,2,3,4), "Age categories 1-4",
                                           ifelse(bias_facet_df_100$age_cat %in% c(5,6,7,8), "Age categories 5-8",
                                                  "Age categories 9-12"))

bias_facet_df_500$age_cat_grouped = ifelse(bias_facet_df_500$age_cat %in% c(1,2,3,4), "Age categories 1-4",
                                           ifelse(bias_facet_df_500$age_cat %in% c(5,6,7,8), "Age categories 5-8",
                                                  "Age categories 9-12"))

p1 = cowplot::plot_grid(ggplot(bias_facet_df_100, aes(x=prob_sampling_old, y=value, group=cell)) + 
                         geom_line(aes(col=as.factor(age_cat_grouped)), show.legend = TRUE, size = 0.5) +
                         facet_wrap(. ~ Model,ncol=3,nrow=1, labeller = as_labeller(
                           function(value) {
                             return(value)  # Lets you change the facet labels
                           })
                         ) +  
                         xlab("\n Probability of response for the elderly \n") +
                         ylab("\n Average Bias \n") + 
                          ylim( -max(max(abs(bias_facet_df_500$value), abs(bias_facet_df_100$value))),
                                max(max(abs(bias_facet_df_500$value), abs(bias_facet_df_100$value))) ) +
                         # scale_colour_manual(values=c("1"="#fbb4b9",
                         #                              "2"="#f768a1",
                         #                              "3"="#c51b8a",
                         #                              "4"="#7a0177",
                         #                              "5"="#fdcc8a",
                         #                              "6"="#fc8d59",
                         #                              "7"="#e34a33",
                         #                              "8"="#b30000",
                         #                              "9"="#b2e2e2",
                         #                              "10"="#66c2a4",
                         #                              "11"="#2ca25f",
                         #                              "12"="#006d2c")
                         # ) +
                        scale_colour_manual(values=c("Age categories 1-4"="#fdcc8a",
                                                     "Age categories 5-8"="#b30000",
                                                     "Age categories 9-12"="#006d2c")
                        ) +
                         theme_bw() +
                         theme(plot.title = element_text(size = 50, face = "bold"),
                               axis.text=element_text(size=80),
                               axis.title=element_text(size=80, face="bold",margin=200),
                               strip.text.x = element_text(size=80, face="bold", margin = margin(t=25,b=25) ),
                               strip.background = element_rect(fill="transparent",color="transparent"),
                               legend.position = "none"
                         ) + 
                          geom_hline(yintercept=0, color="black", size=5, linetype = "dashed") +
                         guides(col=guide_legend(title="Age category",
                                                 override.aes = list(size = 10*1.6))),
                       
                       ggplot(bias_facet_df_500, aes(x=prob_sampling_old, y=value, group=cell)) + 
                         geom_line(aes(col=as.factor(age_cat_grouped)), show.legend = TRUE, size = 0.5) +
                         facet_wrap(. ~ Model,ncol=3,nrow=1, labeller = as_labeller(
                           function(value) {
                             return(value)  # Lets you change the facet labels
                           })
                         ) +  
                         xlab("\n Probability of response for the elderly \n") +
                         ylab("\n Average Bias \n") + 
                         ylim( -max(max(abs(bias_facet_df_500$value), abs(bias_facet_df_100$value))),
                               max(max(abs(bias_facet_df_500$value), abs(bias_facet_df_100$value))) ) +
                         # scale_colour_manual(values=c("1"="#fbb4b9",
                         #                              "2"="#f768a1",
                         #                              "3"="#c51b8a",
                         #                              "4"="#7a0177",
                         #                              "5"="#fdcc8a",
                         #                              "6"="#fc8d59",
                         #                              "7"="#e34a33",
                         #                              "8"="#b30000",
                         #                              "9"="#b2e2e2",
                         #                              "10"="#66c2a4",
                         #                              "11"="#2ca25f",
                         #                              "12"="#006d2c")
                         # ) +
                       scale_colour_manual(values=c("Age categories 1-4"="#fdcc8a",
                                                    "Age categories 5-8"="#b30000",
                                                    "Age categories 9-12"="#006d2c")
                                           ) +
                         theme_bw() +
                         theme(plot.title = element_text(size = 50, face = "bold"),
                               axis.text=element_text(size=80),
                               axis.title=element_text(size=80, face="bold",margin=200),
                               strip.text.x = element_text(size=80, face="bold", margin = margin(t=25,b=25) ),
                               strip.background = element_rect(fill="transparent",color="transparent"),
                               legend.position = "none"
                         ) + 
                         geom_hline(yintercept=0, color="black", size=5, linetype = "dashed") +
                         guides(col=guide_legend(title="Age category",
                                                 override.aes = list(size = 10*1.6)))
,
align="v",
nrow=2,

labels = c("n = 100", "n = 500"),

label_size = 85,

vjust = 1.25,
hjust = -0.25)

legend_p1 = cowplot::get_legend(ggplot(bias_facet_df_500, aes(x=prob_sampling_old, y=value, group=cell)) + 
                                  geom_line(aes(col=as.factor(age_cat_grouped)), show.legend = TRUE, size = 0.5) +
                                  facet_wrap(. ~ Model,ncol=3,nrow=1, labeller = as_labeller(
                                    function(value) {
                                      return(value)  # Lets you change the facet labels
                                    })
                                  ) +  
                                  xlab("\n Probability of response for the elderly \n") +
                                  ylab("\n Average Bias \n") + 
                                  # scale_colour_manual(values=c("1"="#fbb4b9",
                                  #                              "2"="#f768a1",
                                  #                              "3"="#c51b8a",
                                  #                              "4"="#7a0177",
                                  #                              "5"="#fdcc8a",
                                  #                              "6"="#fc8d59",
                                  #                              "7"="#e34a33",
                                  #                              "8"="#b30000",
                                  #                              "9"="#b2e2e2",
                                  #                              "10"="#66c2a4",
                                  #                              "11"="#2ca25f",
                                  #                              "12"="#006d2c")
                                  # ) +
                                scale_colour_manual(values=c("Age categories 1-4"="#fdcc8a",
                                                             "Age categories 5-8"="#b30000",
                                                             "Age categories 9-12"="#006d2c")
                                ) +
                                  theme_bw() +
                                  theme(plot.title = element_text(size = 50, face = "bold"),
                                        axis.text=element_text(size=50),
                                        axis.title=element_text(size=80, face="bold",margin=200),
                                        legend.text = element_text(size=80),
                                        legend.position  = "bottom",
                                        legend.key.size = unit(30*1.6,"line"),
                                        legend.key.height = unit(7*1.6,"line"),
                                        legend.key = element_rect(fill = "transparent",color = "transparent"),
                                        legend.title = element_text(size=80, face="bold"),
                                        strip.text.x = element_text(size=80, face="bold", margin = margin(t=25,b=25) ),
                                        strip.background = element_rect(fill="transparent",color="transparent")
                                  ) + 
                                  guides(col=guide_legend(title="Age category",
                                                          override.aes = list(size = 10*1.6))))

png(biasfacet_pngname,
    width = 5200, 
    height = 5000)

plot(cowplot::plot_grid(p1,
                        legend_p1,
                        ncol=1,
                        rel_heights = c(3, 0.3)))
dev.off()


# plot combined plots for allmedians_facet
allmedians_facet_df_100 = readRDS("allmedians_facet_200_12_100_.rds")
allmedians_facet_df_500 = readRDS("allmedians_facet_200_12_500_.rds")
points_df_final = readRDS("points_df_final.rds")

allmedians_facet_df_100_filtered = allmedians_facet_df_100[allmedians_facet_df_100$p %in% c(0.1, 0.5, 0.9),]

# labels for plots: Strongly undersampled, Representative sampling, Strongly oversampled
allmedians_facet_df_100_filtered[allmedians_facet_df_100_filtered$p==0.1, c("prob_sampling_old")] = "Strongly undersampling age categories 9-12"
allmedians_facet_df_100_filtered[allmedians_facet_df_100_filtered$p==0.5, c("prob_sampling_old")] = "Representatively sampling age categories 9-12"
allmedians_facet_df_100_filtered[allmedians_facet_df_100_filtered$p==0.9, c("prob_sampling_old")] = "Strongly oversampling age categories 9-12"

allmedians_facet_df_100_filtered$prob_sampling_old = factor(allmedians_facet_df_100_filtered$prob_sampling_old, levels = c("Strongly undersampling age categories 9-12",
                                                                                                                           "Representatively sampling age categories 9-12",
                                                                                                                           "Strongly oversampling age categories 9-12"))


allmedians_facet_df_500_filtered = allmedians_facet_df_500[allmedians_facet_df_500$p %in% c(0.1, 0.5, 0.9),]

# labels for plots: Strongly undersampled, Representative sampling, Strongly oversampled
allmedians_facet_df_500_filtered[allmedians_facet_df_500_filtered$p==0.1, c("prob_sampling_old")] = "Strongly undersampling age categories 9-12"
allmedians_facet_df_500_filtered[allmedians_facet_df_500_filtered$p==0.5, c("prob_sampling_old")] = "Representatively sampling age categories 9-12"
allmedians_facet_df_500_filtered[allmedians_facet_df_500_filtered$p==0.9, c("prob_sampling_old")] = "Strongly oversampling age categories 9-12"

allmedians_facet_df_500_filtered$prob_sampling_old = factor(allmedians_facet_df_500_filtered$prob_sampling_old, levels = c("Strongly undersampling age categories 9-12",
                                                                                                                           "Representatively sampling age categories 9-12",
                                                                                                                           "Strongly oversampling age categories 9-12"))


p2 = cowplot::plot_grid(ggplot(allmedians_facet_df_100_filtered, aes(x=age_cat, y=posterior_medians, color=Model, linetype=Model)) + #geom_point(size=5)+
                          geom_rect(allmedians_facet_df_100_filtered, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                                                    xmax = Inf,
                                                                                    ymin = -Inf,
                                                                                    ymax = Inf),
                                    alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
                          geom_jitter(alpha=0.1, show.legend=FALSE) + 
                          geom_point(aes(x = age_cat, y = mrp),
                                     size=7,
                                     colour = "black",
                                     fill="black",
                                     show.legend = F,
                                     data = points_df_final,
                                     inherit.aes = F) +
                          geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
                          facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=1, labeller = as_labeller(
                            function(value) {
                              return(value)  # Lets you change the facet labels
                            })
                          ) +
                          xlab("\n Age Category \n") +
                          ylab(paste("\n Median of", runs,"Posteriors \n")) +
                          scale_x_continuous(breaks= scales::pretty_breaks()) + 
                          coord_cartesian(ylim=c(0,1), expand=FALSE) +
                          theme_bw() +
                          theme(plot.title = element_text(size = 50, face = "bold"),
                                axis.text=element_text(size=50),
                                axis.title=element_text(size=50, face="bold",margin=200),
                                legend.text = element_text(size=50),
                                legend.position  = "none",
                                legend.key.size = unit(25,"line"),
                                legend.key.height = unit(7,"line"),
                                legend.key = element_rect(fill = "transparent",color = "transparent"),
                                legend.title = element_text(size=50, face="bold"),
                                strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
                                strip.background = element_rect(fill="transparent",color="transparent")) +
                          scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
                          guides(color = guide_legend(override.aes = list(size=20))),
                        
                        ggplot(allmedians_facet_df_500_filtered, aes(x=age_cat, y=posterior_medians, color=Model, linetype=Model)) + #geom_point(size=5)+
                          geom_rect(allmedians_facet_df_500_filtered, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                                                    xmax = Inf,
                                                                                    ymin = -Inf,
                                                                                    ymax = Inf),
                                    alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
                          geom_jitter(alpha=0.1, show.legend=FALSE) + 
                          geom_point(aes(x = age_cat, y = mrp),
                                     size=7,
                                     colour = "black",
                                     fill="black",
                                     show.legend = F,
                                     data = points_df_final,
                                     inherit.aes = F) +
                          geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
                          facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=3, labeller = as_labeller(
                            function(value) {
                              return(value)  # Lets you change the facet labels
                            })
                          ) +
                          xlab("\n Age Category \n") +
                          ylab(paste("\n Median of", runs,"Posteriors \n")) +
                          scale_x_continuous(breaks= scales::pretty_breaks()) + 
                          coord_cartesian(ylim=c(0,1), expand=FALSE) +
                          theme_bw() +
                          theme(plot.title = element_text(size = 50, face = "bold"),
                                axis.text=element_text(size=50),
                                axis.title=element_text(size=50, face="bold",margin=200),
                                legend.text = element_text(size=50),
                                legend.position  = "none",
                                legend.key.size = unit(25,"line"),
                                legend.key.height = unit(7,"line"),
                                legend.key = element_rect(fill = "transparent",color = "transparent"),
                                legend.title = element_text(size=50, face="bold"),
                                strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
                                strip.background = element_rect(fill="transparent",color="transparent")) +
                          scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
                          guides(color = guide_legend(override.aes = list(size=20))),
                        
                        ncol = 1,
                        
                        labels = c("n = 100", "n = 500"),
                        
                        label_size = 50,
                        
                        vjust = 1.5,
                        hjust = -0.25
)


legend_p2 = cowplot::get_legend(ggplot(allmedians_facet_df_500_filtered, aes(x=age_cat, y=posterior_medians, color=Model, linetype=Model)) + #geom_point(size=5)+
                                  geom_rect(allmedians_facet_df_500_filtered, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                                                            xmax = Inf,
                                                                                            ymin = -Inf,
                                                                                            ymax = Inf),
                                            alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
                                  geom_jitter(alpha=0.1, show.legend=FALSE) + 
                                  geom_point(aes(x = age_cat, y = mrp),
                                             size=7,
                                             colour = "black",
                                             fill="black",
                                             show.legend = F,
                                             data = points_df_final,
                                             inherit.aes = F) +
                                  geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
                                  facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=3, labeller = as_labeller(
                                    function(value) {
                                      return(value)  # Lets you change the facet labels
                                    })
                                  ) +
                                  xlab("\n Age Category \n") +
                                  ylab(paste("\n Median of", runs,"Posteriors \n")) +
                                  scale_x_continuous(breaks= scales::pretty_breaks()) + 
                                  theme_bw() +
                                  theme(plot.title = element_text(size = 50, face = "bold"),
                                        axis.text=element_text(size=35),
                                        axis.title=element_text(size=50, face="bold",margin=200),
                                        legend.text = element_text(size=50),
                                        legend.position  = "bottom",
                                        legend.key.size = unit(25,"line"),
                                        legend.key.height = unit(7,"line"),
                                        legend.key = element_rect(fill = "transparent",color = "transparent"),
                                        legend.title = element_text(size=50, face="bold"),
                                        strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
                                        strip.background = element_rect(fill="transparent",color="transparent")) +
                                  scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
                                  guides(color = guide_legend(override.aes = list(size=20))))

png(allmediansfacet_pngname,
    width = 3625, 
    height = 2200)

plot(cowplot::plot_grid(p2,
     legend_p2,
     ncol=1,
     rel_heights = c(3, 0.3)))
dev.off()


# plot combined plots for allquantilediff_facet
allquantilediff_facet_df_100 = readRDS("allquantilediff_facet_200_12_100_.rds")
allquantilediff_facet_df_500 = readRDS("allquantilediff_facet_200_12_500_.rds")

allquantilediff_facet_df_100_filtered = allquantilediff_facet_df_100[allquantilediff_facet_df_100$p %in% c(0.1, 0.5, 0.9),]

# labels for plots: Strongly undersampled, Representative sampling, Strongly oversampled
allquantilediff_facet_df_100_filtered[allquantilediff_facet_df_100_filtered$p==0.1, c("prob_sampling_old")] = "Strongly undersampling age categories 9-12"
allquantilediff_facet_df_100_filtered[allquantilediff_facet_df_100_filtered$p==0.5, c("prob_sampling_old")] = "Representatively sampling age categories 9-12"
allquantilediff_facet_df_100_filtered[allquantilediff_facet_df_100_filtered$p==0.9, c("prob_sampling_old")] = "Strongly oversampling age categories 9-12"

allquantilediff_facet_df_100_filtered$prob_sampling_old = factor(allquantilediff_facet_df_100_filtered$prob_sampling_old, levels = c("Strongly undersampling age categories 9-12",
                                                                                                                                     "Representatively sampling age categories 9-12",
                                                                                                                                     "Strongly oversampling age categories 9-12"))


allquantilediff_facet_df_500_filtered = allquantilediff_facet_df_500[allquantilediff_facet_df_500$p %in% c(0.1, 0.5, 0.9),]

# labels for plots: Strongly undersampled, Representative sampling, Strongly oversampled
allquantilediff_facet_df_500_filtered[allquantilediff_facet_df_500_filtered$p==0.1, c("prob_sampling_old")] = "Strongly undersampling age categories 9-12"
allquantilediff_facet_df_500_filtered[allquantilediff_facet_df_500_filtered$p==0.5, c("prob_sampling_old")] = "Representatively sampling age categories 9-12"
allquantilediff_facet_df_500_filtered[allquantilediff_facet_df_500_filtered$p==0.9, c("prob_sampling_old")] = "Strongly oversampling age categories 9-12"

allquantilediff_facet_df_500_filtered$prob_sampling_old = factor(allquantilediff_facet_df_500_filtered$prob_sampling_old, levels = c("Strongly undersampling age categories 9-12",
                                                                                                                                     "Representatively sampling age categories 9-12",
                                                                                                                                     "Strongly oversampling age categories 9-12"))



p3 = cowplot::plot_grid(ggplot(allquantilediff_facet_df_100_filtered, aes(x=age_cat, y=quantildiff90_10, color=Model, linetype=Model)) + #geom_point(size=5)+
                          geom_rect(allquantilediff_facet_df_100_filtered, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                                               xmax = Inf,
                                                                               ymin = -Inf,
                                                                               ymax = Inf),
                                    alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
                          geom_jitter(alpha=0.1, show.legend=FALSE) + 
                          geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
                          facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=3, labeller = as_labeller(
                            function(value) {
                              return(value)  # Lets you change the facet labels
                            })
                          ) +
                          xlab("\n Age Category \n") +
                          ylab(paste("\n 90th - 10th quantile of Posteriors \n")) +
                          scale_x_continuous(breaks= scales::pretty_breaks()) + 
                          coord_cartesian(ylim=c(0, max(allquantilediff_facet_df_100_filtered$quantildiff90_10)), expand=FALSE) +
                          theme_bw() +
                          theme(plot.title = element_text(size = 50, face = "bold"),
                                axis.text=element_text(size=50),
                                axis.title=element_text(size=50, face="bold",margin=200),
                                legend.text = element_text(size=50),
                                legend.position  = "none",
                                legend.key.size = unit(25,"line"),
                                legend.key.height = unit(7,"line"),
                                legend.key = element_rect(fill = "transparent",color = "transparent"),
                                legend.title = element_text(size=50, face="bold"),
                                strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
                                strip.background = element_rect(fill="transparent",color="transparent")) +
                          scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
                          guides(color = guide_legend(override.aes = list(size=20))),
                        
                        ggplot(allquantilediff_facet_df_500_filtered, aes(x=age_cat, y=quantildiff90_10, color=Model, linetype=Model)) + #geom_point(size=5)+
                          geom_rect(allquantilediff_facet_df_500_filtered, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                                                         xmax = Inf,
                                                                                         ymin = -Inf,
                                                                                         ymax = Inf),
                                    alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
                          geom_jitter(alpha=0.1, show.legend=FALSE) + 
                          geom_smooth(aes(group=Model), method="loess", size=2.5, se=FALSE) + 
                          facet_wrap(. ~ prob_sampling_old,ncol=3,nrow=3, labeller = as_labeller(
                            function(value) {
                              return(value)  # Lets you change the facet labels
                            })
                          ) +
                          xlab("\n Age Category \n") +
                          ylab(paste("\n 90th - 10th quantile of Posteriors \n")) +
                          scale_x_continuous(breaks= scales::pretty_breaks()) + 
                          coord_cartesian(ylim=c(0, max(allquantilediff_facet_df_100_filtered$quantildiff90_10)), expand=FALSE) +
                          theme_bw() +
                          theme(plot.title = element_text(size = 50, face = "bold"),
                                axis.text=element_text(size=50),
                                axis.title=element_text(size=50, face="bold",margin=200),
                                legend.text = element_text(size=50),
                                legend.position  = "none",
                                legend.key.size = unit(25,"line"),
                                legend.key.height = unit(7,"line"),
                                legend.key = element_rect(fill = "transparent",color = "transparent"),
                                legend.title = element_text(size=50, face="bold"),
                                strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
                                strip.background = element_rect(fill="transparent",color="transparent")) +
                          scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
                          guides(color = guide_legend(override.aes = list(size=20))),
                        
                        ncol=1,
                        
                        labels = c("n = 100", "n = 500"),
                        
                        label_size = 50,
                        
                        vjust = 1,
                        hjust = -0.25
                        )

legend_p3 = cowplot::get_legend(ggplot(allquantilediff_facet_df_500_filtered, aes(x=age_cat, y=quantildiff90_10, color=Model, linetype=Model)) + #geom_point(size=5)+
                                  geom_rect(allquantilediff_facet_df_500_filtered, mapping = aes(xmin = (2 * (age_grouping_multiplier * 2/3) + 1)/2,
                                                                                                 xmax = Inf,
                                                                                                 ymin = -Inf,
                                                                                                 ymax = Inf),
                                            alpha = 0.1,fill = "gray93", colour = NA, show.legend=FALSE) +
                                  geom_jitter(alpha=0.1, show.legend=FALSE) + 
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
                                        axis.text=element_text(size=35),
                                        axis.title=element_text(size=50, face="bold",margin=200),
                                        legend.text = element_text(size=50),
                                        legend.position  = "bottom",
                                        legend.key.size = unit(25,"line"),
                                        legend.key.height = unit(7,"line"),
                                        legend.key = element_rect(fill = "transparent",color = "transparent"),
                                        legend.title = element_text(size=50, face="bold"),
                                        strip.text.x = element_text(size=50, face="bold", margin = margin(t=25,b=25) ),
                                        strip.background = element_rect(fill="transparent",color="transparent")) +
                                  scale_color_manual(values=c("#4575b4", "#d73027","#fdae61")) + 
                                  guides(color = guide_legend(override.aes = list(size=20))))

png(allquantilediff_pngname,
    width = 3625, 
    height = 2200)

plot(cowplot::plot_grid(p3,
                        legend_p3,
                        ncol=1,
                        rel_heights = c(3, 0.3)))
dev.off()

print(warnings())

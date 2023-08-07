# Sensitivity analyses to select DNA QC metrics threshold to improve aTL-age correlation
# June 17, 2023

# Load packages
library(tidyverse)
library(correlation)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Generate the data matrix for plotting
#With the change of the threshold of DNA QC metrics, we recalculate the partial Pearson's correlation between aTL and age after accounting for sex.
dna_var <- c("X260.280", "X260.230", "A320", "DIN", "PUnfrag", "PSev", "Pico", "Nano")
dna_win_var <- paste0(dna_var, "_win")
dna_var_label <- c("260/280", "260/230", "A320", "DIN", "%Unfrag",
                   "%Severe Frag", "PicoGreen Conc (ng/uL)", "Nanodrop Conc (ng/uL)")

cohort <- unique(dataSW$Cohort)
tissue <- levels(dataSW$Tissue)
cor_res <- data.frame(NULL)
for (i in 1:length(dna_win_var)){
  for (j in 1:length(cohort)){
    for (k in 1:length(tissue)){
      data_sub <- dataSW %>%
        filter(Cohort == cohort[j] & Tissue == tissue[k])
      if (nrow(data_sub) == 0) {
        k <- k + 1
      } else {
        if (dna_win_var[i] %in% c("X260.280_win", "X260.230_win")){
          cor_res_1 <- data.frame(NULL)
          cor_res_2 <- data.frame(NULL)
          val_seq_1 <- seq(from = min(data_sub[,dna_win_var[i]],na.rm=T),
                           to = quantile(data_sub[,dna_win_var[i]],na.rm=T, prob = 0.75),
                           length.out = 50)
          val_seq_2 <- seq(from = max(data_sub[,dna_win_var[i]],na.rm=T),
                           to = quantile(data_sub[,dna_win_var[i]],na.rm=T, prob = 0.25),
                           length.out = 50)
          for (l in 1:length(val_seq_1)){
            data_temp <- data_sub %>%
              filter(!!as.name(dna_win_var[i]) >= val_seq_1[l]) %>%
              select(Age, Sex, aTL_win)
            cor_temp <- cor_test(data=data_temp, x="Age", y="aTL_win",
                                 method="pearson",include_factors=T, partial=T, 
                                 multilevel=F, use = "pairwise.complete.obs")
            cor_res_1 <- rbind(cor_res_1, 
                               as.data.frame(cor_temp) %>%
                                 mutate(Cohort = cohort[j],
                                        Tissue = tissue[k],
                                        DNAQCmetrics = dna_win_var[i],
                                        DNAQCmetrics_threshold_direction = ">=",
                                        DNAQCmetrics_threshold_value = val_seq_1[l]))
          }
          for (l in 1:length(val_seq_2)){
            data_temp <- data_sub %>%
              filter(!!as.name(dna_win_var[i]) <= val_seq_2[l]) %>%
              select(Age, Sex, aTL_win)
            cor_temp <- cor_test(data=data_temp, x="Age", y="aTL_win",
                                 method="pearson",include_factors=T, partial=T, 
                                 multilevel=F, use = "pairwise.complete.obs")
            cor_res_2 <- rbind(cor_res_2, 
                               as.data.frame(cor_temp) %>%
                                 mutate(Cohort = cohort[j],
                                        Tissue = tissue[k],
                                        DNAQCmetrics = dna_win_var[i],
                                        DNAQCmetrics_threshold_direction = "<=",
                                        DNAQCmetrics_threshold_value = val_seq_2[l]))
          }
          cor_res <- rbind(cor_res, cor_res_1, cor_res_2)
          
        } else if (dna_win_var[i] %in% c("DIN_win", "PUnfrag_win", "Pico_win", "Nano_win")){
          val_seq <- seq(from = min(data_sub[,dna_win_var[i]],na.rm=T),
                         to = quantile(data_sub[,dna_win_var[i]],na.rm=T, prob = 0.75),
                         length.out = 50)
          for (l in 1:length(val_seq)){
            data_temp <- data_sub %>%
              filter(!!as.name(dna_win_var[i]) >= val_seq[l]) %>%
              select(Age, Sex, aTL_win)
            cor_temp <- cor_test(data=data_temp, x="Age", y="aTL_win",
                                 method="pearson",include_factors=T, partial=T, 
                                 multilevel=F, use = "pairwise.complete.obs")
            cor_res <- rbind(cor_res, 
                             as.data.frame(cor_temp) %>%
                               mutate(Cohort = cohort[j],
                                      Tissue = tissue[k],
                                      DNAQCmetrics = dna_win_var[i],
                                      DNAQCmetrics_threshold_direction = ">=",
                                      DNAQCmetrics_threshold_value = val_seq[l]))
            
          }
        } else if (dna_win_var[i] %in% c("A320_win", "PSev_win")) {
          val_seq <- seq(from = max(data_sub[,dna_win_var[i]],na.rm=T),
                         to = quantile(data_sub[,dna_win_var[i]],na.rm=T, prob = 0.25),
                         length.out = 50)
          for (l in 1:length(val_seq)){
            data_temp <- data_sub %>%
              filter(!!as.name(dna_win_var[i]) <= val_seq[l]) %>%
              select(Age, Sex, aTL_win)
            cor_temp <- cor_test(data=data_temp, x="Age", y="aTL_win",
                                 method="pearson",include_factors=T, partial=T, 
                                 multilevel=F, use = "pairwise.complete.obs")
            cor_res <- rbind(cor_res, 
                             as.data.frame(cor_temp) %>%
                               mutate(Cohort = cohort[j],
                                      Tissue = tissue[k],
                                      DNAQCmetrics = dna_win_var[i],
                                      DNAQCmetrics_threshold_direction = "<=",
                                      DNAQCmetrics_threshold_value = val_seq[l]))
            
          }
        }
      }
    }
  }
}

head(cor_res)

# Generate line plots showing how the aTL-age correlation changes with the change of the threshold
#One plot panel per DNA QC metrics
cor_res <- cor_res %>%
  mutate(Cohort = factor(Cohort, levels = c("P50", "TL"), labels = c("Child", "Adult")),
         Cohort_Tissue = paste(Cohort, Tissue, sep = "_"),
         DNAQCmetrics = factor(DNAQCmetrics, levels = dna_win_var, labels = dna_var_label))

df_list <- cor_res %>%
  split(f = list(.$DNAQCmetrics_threshold_direction, .$DNAQCmetrics), drop = T)

length(df_list)

p_list <- purrr::pmap(
  .l = list(data = df_list),
  .f = function(data){
    ggplot(data = data, aes(x = DNAQCmetrics_threshold_value, y = r)) +
      geom_line(aes(group = Cohort_Tissue, color = Tissue, linetype = Cohort)) +
      labs(x = unique(data$DNAQCmetrics), 
           title = unique(data$DNAQCmetrics_threshold_direction)) +
      facet_wrap(facets = vars(DNAQCmetrics)) +
      scale_color_manual(values=c("red","blue","green4","purple","orange")) +
      scale_y_continuous(limits = c(-0.75, 0.60)) +
      theme_bw()
  }
)

pdf(file = "U01_Project/Results/DNAQC_threshold_selection.pdf", width = 8, height = 15)
ggstatsplot::combine_plots(
  plotlist = p_list,
  plotgrid.args = list(nrow = 5),
  annotation.args = list(tag_levels = "A")
) + labs(caption = "r: Paritial Pearson's correlation between aTL and age after accounting for sex;
         \"<=\" means keeping samples with DNA QC metrics of less than each value on the x axis;
         \">=\" means keeping samples with DNA QC metrics of greater than each value on the x axis.")

dev.off()

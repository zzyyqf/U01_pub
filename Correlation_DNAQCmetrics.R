# U01 Project
# Correlation among DNA QC metrics
# June 15, 2023

# Load packages
library(tidyverse)
library(psych)
library(ggpattern)
library(grid)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Rename the tissues
dataSW$Tissue <- factor(dataSW$Tissue, levels = c("Buccal","Saliva","DBS","Buffy","PBMC"),
                        labels = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"))

# Correlation between DNA QC metrics (winsorized) by cohort and tissue (adjusted for age and sex).# -----
cohort <- unique(dataSW$Cohort)
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230", 
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")

#generate the data matrix for plotting
cor_res <- data.frame(NULL)
for (i in 1:length(cohort)){
  for (j in 1:length(tissue)){
    data_temp <- dataSW %>% 
      filter(Cohort == cohort[i] & Tissue == tissue[j])
    if (nrow(data_temp) == 0){
      j <- j+1
    } else {
      for (k in 1:length(dna_qc_var)){
        for (l in 1:length(dna_qc_var)){
          if (k >= l){
            next
          } else {
            data_temp2 <- data_temp %>%
              select(all_of(c("Age", "Sex", dna_qc_var[k], dna_qc_var[l])))
            cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                              x=dna_qc_var[k], y=dna_qc_var[l],
                                              method="spearman",include_factors=T, partial=T, 
                                              multilevel=F)
            cor_res <- rbind(cor_res, 
                             as.data.frame(cor_temp) %>%
                               mutate(Cohort = cohort[i],
                                      Tissue = tissue[j]))
          }
        }
      }
    }
  }
}
cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Parameter2 <- factor(cor_res$Parameter2, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Cohort <- factor(cor_res$Cohort, levels = c("P50", "TL"),
                         labels = c("Child", "Adult"))
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#generate the plot
gg <- ggplot(data = cor_res, aes(x = Tissue, y = rho, group = Cohort)) +
  geom_col_pattern(aes(col = Tissue, fill = Tissue, pattern = Cohort, group = Cohort), 
                   pattern_color = "white", pattern_fill = "white",
                   position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.06, rho-0.12),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(rows=vars(Parameter1), cols=vars(Parameter2)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  #labs(caption = "Spearman's r among DNA QC metrics accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10))

#remove extra elements in the plot
grob <- ggplotGrob(gg)

#remove facets
grob_idx <- which(grob$layout$name %in% c("panel-2-1", 
                                          "panel-3-1", "panel-3-2",
                                          "panel-4-1", "panel-4-2", "panel-4-3",
                                          "panel-5-1", "panel-5-2", "panel-5-3", "panel-5-4",
                                          "panel-6-1", "panel-6-2", "panel-6-3", "panel-6-4", "panel-6-5",
                                          "panel-7-1", "panel-7-2", "panel-7-3", "panel-7-4", "panel-7-5", "panel-7-6",
                                          "panel-8-1", "panel-8-2", "panel-8-3", "panel-8-4", "panel-8-5", "panel-8-6", "panel-8-7"))

for(i in grob_idx) grob$grobs[[i]] <- grid::nullGrob()

#move the x axes up
x_idx <- which(grob$layout$name %in% paste0("axis-b-", 1:7))
grob$layout[x_idx, c("t", "b")] <- grob$layout[x_idx, c("t", "b")] - seq(14, 2, by = -2)

#move the y axes to the right
y_idx <- which(grob$layout$name %in% paste0("axis-l-", 2:8))
grob$layout[y_idx, c("l", "r")] <- grob$layout[y_idx, c("l", "r")] + seq(2, 14, by = 2)

#export the plot
pdf(file = "U01_Project/Results/Winsorized/Correlation_DNAQCmetrics_win.pdf", height = 11, width = 13)
grid::grid.draw(grob)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))

# Correlation between DNA QC metrics (not winsorized) by cohort and tissue (adjusted for age and sex).# -----
cohort <- unique(dataSW$Cohort)
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")

#generate the data matrix for plotting
cor_res <- data.frame(NULL)
for (i in 1:length(cohort)){
  for (j in 1:length(tissue)){
    data_temp <- dataSW %>% 
      filter(Cohort == cohort[i] & Tissue == tissue[j])
    if (nrow(data_temp) == 0){
      j <- j+1
    } else {
      for (k in 1:length(dna_qc_var)){
        for (l in 1:length(dna_qc_var)){
          if (k >= l){
            next
          } else {
            data_temp2 <- data_temp %>%
              select(all_of(c("Age", "Sex", dna_qc_var[k], dna_qc_var[l])))
            cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                              x=dna_qc_var[k], y=dna_qc_var[l],
                                              method="spearman",include_factors=T, partial=T, 
                                              multilevel=F)
            cor_res <- rbind(cor_res, 
                             as.data.frame(cor_temp) %>%
                               mutate(Cohort = cohort[i],
                                      Tissue = tissue[j]))
          }
        }
      }
    }
  }
}
cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Parameter2 <- factor(cor_res$Parameter2, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Cohort <- factor(cor_res$Cohort, levels = c("P50", "TL"),
                         labels = c("Child", "Adult"))
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#generate the plot
gg <- ggplot(data = cor_res, aes(x = Tissue, y = rho, group = Cohort)) +
  geom_col_pattern(aes(col = Tissue, fill = Tissue, pattern = Cohort, group = Cohort), 
                   pattern_color = "white", pattern_fill = "white",
                   position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.06, rho-0.12),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(rows=vars(Parameter1), cols=vars(Parameter2)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  #labs(caption = "Spearman's r among DNA QC metrics accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10))

#remove extra elements in the plot
grob <- ggplotGrob(gg)

#remove facets
grob_idx <- which(grob$layout$name %in% c("panel-2-1", 
                                          "panel-3-1", "panel-3-2",
                                          "panel-4-1", "panel-4-2", "panel-4-3",
                                          "panel-5-1", "panel-5-2", "panel-5-3", "panel-5-4",
                                          "panel-6-1", "panel-6-2", "panel-6-3", "panel-6-4", "panel-6-5",
                                          "panel-7-1", "panel-7-2", "panel-7-3", "panel-7-4", "panel-7-5", "panel-7-6",
                                          "panel-8-1", "panel-8-2", "panel-8-3", "panel-8-4", "panel-8-5", "panel-8-6", "panel-8-7"))

for(i in grob_idx) grob$grobs[[i]] <- grid::nullGrob()

#move the x axes up
x_idx <- which(grob$layout$name %in% paste0("axis-b-", 1:7))
grob$layout[x_idx, c("t", "b")] <- grob$layout[x_idx, c("t", "b")] - seq(14, 2, by = -2)

#move the y axes to the right
y_idx <- which(grob$layout$name %in% paste0("axis-l-", 2:8))
grob$layout[y_idx, c("l", "r")] <- grob$layout[y_idx, c("l", "r")] + seq(2, 14, by = 2)

#export the plot
pdf(file = "U01_Project/Results/Not_winsorized/Correlation_DNAQCmetrics_nowin.pdf", height = 11, width = 13)
grid::grid.draw(grob)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))

# Correlation between DNA QC metrics (winsorized) by tissue (adjusted for age and sex).# -----
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")

#generate the data matrix for plotting
cor_res <- data.frame(NULL)
for (j in 1:length(tissue)){
  data_temp <- dataSW %>% 
    filter(Tissue == tissue[j])
  if (nrow(data_temp) == 0){
    j <- j+1
  } else {
    for (k in 1:length(dna_qc_var)){
      for (l in 1:length(dna_qc_var)){
        if (k >= l){
          next
        } else {
          data_temp2 <- data_temp %>%
            select(all_of(c("Age", "Sex", dna_qc_var[k], dna_qc_var[l])))
          cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                            x=dna_qc_var[k], y=dna_qc_var[l],
                                            method="spearman",include_factors=T, partial=T, 
                                            multilevel=F)
          cor_res <- rbind(cor_res, 
                           as.data.frame(cor_temp) %>%
                             mutate(Tissue = tissue[j]))
        }
      }
    }
  }
}
cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Parameter2 <- factor(cor_res$Parameter2, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#export the correlation result
write_csv(cor_res, file = "U01_Project/Results/Winsorized/correlation_tab_DNAQCmetrics_combinedCohort_win.csv")

#generate the plot
gg <- ggplot(data = cor_res, aes(x = Tissue, y = rho)) +
  geom_col(aes(col = Tissue, fill = Tissue), 
           position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.06, rho-0.15),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(rows=vars(Parameter1), cols=vars(Parameter2)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  #labs(caption = "Spearman's r among DNA QC metrics accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10),
        legend.title = element_blank())

#remove extra elements in the plot
grob <- ggplotGrob(gg)

#remove facets
grob_idx <- which(grob$layout$name %in% c("panel-2-1", 
                                          "panel-3-1", "panel-3-2",
                                          "panel-4-1", "panel-4-2", "panel-4-3",
                                          "panel-5-1", "panel-5-2", "panel-5-3", "panel-5-4",
                                          "panel-6-1", "panel-6-2", "panel-6-3", "panel-6-4", "panel-6-5",
                                          "panel-7-1", "panel-7-2", "panel-7-3", "panel-7-4", "panel-7-5", "panel-7-6",
                                          "panel-8-1", "panel-8-2", "panel-8-3", "panel-8-4", "panel-8-5", "panel-8-6", "panel-8-7"))

for(i in grob_idx) grob$grobs[[i]] <- grid::nullGrob()

#move the x axes up
x_idx <- which(grob$layout$name %in% paste0("axis-b-", 1:7))
grob$layout[x_idx, c("t", "b")] <- grob$layout[x_idx, c("t", "b")] - seq(14, 2, by = -2)

#move the y axes to the right
y_idx <- which(grob$layout$name %in% paste0("axis-l-", 2:8))
grob$layout[y_idx, c("l", "r")] <- grob$layout[y_idx, c("l", "r")] + seq(2, 14, by = 2)

#export the plot
pdf(file = "U01_Project/Results/Winsorized/Correlation_DNAQCmetrics_win_combinedCohorts.pdf", height = 10, width = 11)
grid::grid.draw(grob)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))



# Correlation between DNA QC metrics (not winsorized) by tissue (adjusted for age and sex).# -----
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")

#generate the data matrix for plotting
cor_res <- data.frame(NULL)
for (j in 1:length(tissue)){
  data_temp <- dataSW %>% 
    filter(Tissue == tissue[j])
  if (nrow(data_temp) == 0){
    j <- j+1
  } else {
    for (k in 1:length(dna_qc_var)){
      for (l in 1:length(dna_qc_var)){
        if (k >= l){
          next
        } else {
          data_temp2 <- data_temp %>%
            select(all_of(c("Age", "Sex", dna_qc_var[k], dna_qc_var[l])))
          cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                            x=dna_qc_var[k], y=dna_qc_var[l],
                                            method="spearman",include_factors=T, partial=T, 
                                            multilevel=F)
          cor_res <- rbind(cor_res, 
                           as.data.frame(cor_temp) %>%
                             mutate(Tissue = tissue[j]))
        }
      }
    }
  }
}
cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Parameter2 <- factor(cor_res$Parameter2, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#export the correlation result
write_csv(cor_res, file = "U01_Project/Results/Not_winsorized/correlation_tab_DNAQCmetrics_combinedCohort_nowin.csv")

#generate the plot
gg <- ggplot(data = cor_res, aes(x = Tissue, y = rho)) +
  geom_col(aes(col = Tissue, fill = Tissue), 
           position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.06, rho-0.15),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(rows=vars(Parameter1), cols=vars(Parameter2)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  #labs(caption = "Spearman's r among DNA QC metrics accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10),
        legend.title = element_blank())

#remove extra elements in the plot
grob <- ggplotGrob(gg)

#remove facets
grob_idx <- which(grob$layout$name %in% c("panel-2-1", 
                                          "panel-3-1", "panel-3-2",
                                          "panel-4-1", "panel-4-2", "panel-4-3",
                                          "panel-5-1", "panel-5-2", "panel-5-3", "panel-5-4",
                                          "panel-6-1", "panel-6-2", "panel-6-3", "panel-6-4", "panel-6-5",
                                          "panel-7-1", "panel-7-2", "panel-7-3", "panel-7-4", "panel-7-5", "panel-7-6",
                                          "panel-8-1", "panel-8-2", "panel-8-3", "panel-8-4", "panel-8-5", "panel-8-6", "panel-8-7"))

for(i in grob_idx) grob$grobs[[i]] <- grid::nullGrob()

#move the x axes up
x_idx <- which(grob$layout$name %in% paste0("axis-b-", 1:7))
grob$layout[x_idx, c("t", "b")] <- grob$layout[x_idx, c("t", "b")] - seq(14, 2, by = -2)

#move the y axes to the right
y_idx <- which(grob$layout$name %in% paste0("axis-l-", 2:8))
grob$layout[y_idx, c("l", "r")] <- grob$layout[y_idx, c("l", "r")] + seq(2, 14, by = 2)

#export the plot
pdf(file = "U01_Project/Results/Not_winsorized/Correlation_DNAQCmetrics_nowin_combinedCohorts.pdf", height = 10, width = 11)
grid::grid.draw(grob)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))



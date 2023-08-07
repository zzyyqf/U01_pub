# U01 Project
# Correlation between telomere length and DNA QC metrics
# June 15, 2023

# Load packages
library(tidyverse)
library(psych)
library(correlation)
library(ggpattern)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Rename the tissues
dataSW$Tissue <- factor(dataSW$Tissue, levels = c("Buccal","Saliva","DBS","Buffy","PBMC"),
                        labels = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"))

# Correlation between aTL and DNA QC metrics (winsorized) by cohort and tissue # -----
#generate the matrix for plotting
cohort <- unique(dataSW$Cohort)
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                    "X260.280", "X260.230",
                    "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")

cor_res <- data.frame(NULL)
for (i in 1:length(cohort)){
  for (j in 1:length(tissue)){
    data_temp1 <- dataSW %>%
      filter(Cohort == cohort[i] & Tissue == tissue[j])
    if (nrow(data_temp1) == 0){
      j <- j+1
    } else {
      for (k in 1:length(dna_qc_var)){
        data_temp2 <- data_temp1 %>%
          select(all_of(c("Age", "Sex", "aTL_win", dna_qc_var[k])))
        cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                          x=dna_qc_var[k], y="aTL_win",
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

cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                                 labels = dna_qc_var_label)
cor_res$Cohort <- factor(cor_res$Cohort, levels = c("P50", "TL"),
                             labels = c("Child", "Adult"))
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#generate the plot and export
pdf(file = "U01_Project/Results/Winsorized/Correlation_aTL_DNAQCmetrics_win.pdf", width = 12, height = 4)
p1 <- ggplot(data = cor_res, aes(x = Tissue, y = rho, group = Cohort)) +
  geom_col_pattern(aes(col = Tissue, fill = Tissue, pattern = Cohort, group = Cohort), 
                   pattern_color = "white", pattern_fill = "white",
                   position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.03, rho-0.04),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(cols=vars(Parameter1)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  ylab("Partial Spearman's rho with aTL") +
  #labs(caption = "Partial Spearman's r accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10),
        legend.position = "bottom")
cowplot::plot_grid(p1)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))

# Correlation between aTL and DNA QC metrics (not winsorized) by cohort and tissue # -----
#generate the matrix for plotting
cohort <- unique(dataSW$Cohort)
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")
cor_res <- data.frame(NULL)
for (i in 1:length(cohort)){
  for (j in 1:length(tissue)){
    data_temp1 <- dataSW %>%
      filter(Cohort == cohort[i] & Tissue == tissue[j])
    if (nrow(data_temp1) == 0){
      j <- j+1
    } else {
      for (k in 1:length(dna_qc_var)){
        data_temp2 <- data_temp1 %>%
          select(all_of(c("Age", "Sex", "aTL", dna_qc_var[k])))
        cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                          x=dna_qc_var[k], y="aTL",
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

cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Cohort <- factor(cor_res$Cohort, levels = c("P50", "TL"),
                         labels = c("Child", "Adult"))
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#generate the plot and export
pdf(file = "U01_Project/Results/Not_winsorized/Correlation_aTL_DNAQCmetrics_nowin.pdf", width = 12, height = 4)
p1 <- ggplot(data = cor_res, aes(x = Tissue, y = rho, group = Cohort)) +
  geom_col_pattern(aes(col = Tissue, fill = Tissue, pattern = Cohort, group = Cohort), 
                   pattern_color = "white", pattern_fill = "white",
                   position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.03, rho-0.04),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(cols=vars(Parameter1)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  ylab("Partial Spearman's rho with aTL") +
  #labs(caption = "Partial Spearman's r accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10),
        legend.position = "bottom")
cowplot:: plot_grid(p1)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))

# Correlation between aTL and DNA QC metrics (winsorized) by tissue # -----
#generate the matrix for plotting
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")
cor_res <- data.frame(NULL)
for (j in 1:length(tissue)){
  data_temp1 <- dataSW %>%
    filter(Tissue == tissue[j])
  if (nrow(data_temp1) == 0){
    j <- j+1
  } else {
    for (k in 1:length(dna_qc_var)){
      data_temp2 <- data_temp1 %>%
        select(all_of(c("Age", "Sex", "aTL_win", dna_qc_var[k])))
      cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                        x=dna_qc_var[k], y="aTL_win",
                                        method="spearman",include_factors=T, partial=T, 
                                        multilevel=F)
      cor_res <- rbind(cor_res, 
                       as.data.frame(cor_temp) %>%
                         mutate(Tissue = tissue[j]))
    }
  }
}

cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#export the correlation result
write_csv(cor_res, file = "U01_Project/Results/Winsorized/correlation_tab_aTLandDNAQCmetrics_combinedCohort_win.csv")

#generate the plot and export
pdf(file = "U01_Project/Results/Winsorized/Correlation_aTL_DNAQCmetrics_win_combinedCohorts.pdf", width = 12, height = 4)
p1 <- ggplot(data = cor_res, aes(x = Tissue, y = rho)) +
  geom_col(aes(col = Tissue, fill = Tissue),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.03, rho-0.04),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(cols=vars(Parameter1)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  ylab("Partial Spearman's rho with aTL") +
  #labs(caption = "Partial Spearman's r accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom")
cowplot:: plot_grid(p1)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))

# Correlation between aTL and DNA QC metrics (not winsorized) by tissue # -----
#generate the matrix for plotting
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfragmented\n(>3000 bp)", "%Highly\nFragmented\n(250 to 3000 bp)", "%Severely\nFragmented\n(<250 bp)",
                      "A260/A280", "A260/A230",
                      "Nanodrop\nConcentration\n(ng/uL)", "PicoGreen\nConcentration\n(ng/uL)", "TapeStation\nConcentration\n(ng/uL)")
cor_res <- data.frame(NULL)
for (j in 1:length(tissue)){
  data_temp1 <- dataSW %>%
    filter(Tissue == tissue[j])
  if (nrow(data_temp1) == 0){
    j <- j+1
  } else {
    for (k in 1:length(dna_qc_var)){
      data_temp2 <- data_temp1 %>%
        select(all_of(c("Age", "Sex", "aTL", dna_qc_var[k])))
      cor_temp <- correlation::cor_test(data=data_temp2 %>% na.omit(), 
                                        x=dna_qc_var[k], y="aTL",
                                        method="spearman",include_factors=T, partial=T, 
                                        multilevel=F)
      cor_res <- rbind(cor_res, 
                       as.data.frame(cor_temp) %>%
                         mutate(Tissue = tissue[j]))
    }
  }
}

cor_res
cor_res$Parameter1 <- factor(cor_res$Parameter1, levels = dna_qc_var,
                             labels = dna_qc_var_label)
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))
cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#export the correlation result
write_csv(cor_res, file = "U01_Project/Results/Not_winsorized/correlation_tab_aTLandDNAQCmetrics_combinedCohort_nowin.csv")

#generate the plot and export
pdf(file = "U01_Project/Results/Not_winsorized/Correlation_aTL_DNAQCmetrics_nowin_combinedCohorts.pdf", width = 12, height = 4)
p1 <- ggplot(data = cor_res, aes(x = Tissue, y = rho)) +
  geom_col(aes(col = Tissue, fill = Tissue),
           position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_text(aes(y = ifelse(rho>=0, rho+0.03, rho-0.04),
                label = ifelse(p.adj <0.01, "*", "")),
            position = position_dodge(0.9)) +
  facet_grid(cols=vars(Parameter1)) + 
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  ylab("Partial Spearman's rho with aTL") +
  #labs(caption = "Partial Spearman's r accounting for age and sex\n*Statistical significant controlling FDR at <0.01") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom")
cowplot:: plot_grid(p1)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide")))


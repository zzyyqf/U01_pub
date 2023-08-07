# U01 Project
# DNA QC metrics cross-tissue comparison
# June 17, 2023

# Load packages
library(tidyverse)
library(ggstatsplot)
library(statsExpressions)
library(ggsignif)
library(ggpubr)
library(ggpattern)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Change the labels of the cohorts
dataSW <- dataSW %>%
  mutate(Cohort = factor(Cohort, levels = c("P50", "TL"), labels = c("Child", "Adult")))

# Rename the tissues
dataSW$Tissue <- factor(dataSW$Tissue, levels = c("Buccal","Saliva","DBS","Buffy","PBMC"),
                        labels = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"))

# Violin plots of DNA QC metrics cross-tissue comparison (winsorized, combined cohorts) # -----
cohort <- levels(dataSW$Cohort)
dna_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
             "X260.280", "X260.230", 
             "Nano", "Pico", "Tape")
dna_var <- paste0(dna_var, "_win")
dna_var_label <- c("DIN", "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
                   "A260/A280", "A260/A230",
                   "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")

for (i in 1:length(dna_var)){
  p <- ggplot(data = dataSW, aes(x = Tissue, y = .data[[dna_var[i]]])) +
    geom_violin(aes(x = Tissue, y = .data[[dna_var[i]]]), scale = "width") +
    geom_boxplot(aes(x = Tissue, y = .data[[dna_var[i]]], fill = Tissue),
                         width = 0.3, fatten = 3, alpha = 0.7,
                         position = position_dodge(width = 0.9, preserve = "total")) +
    scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
    scale_color_manual(values=c("red","blue","green4","purple","orange")) +
    labs(y = dna_var_label[i],
         title = dna_var_label[i]) +
    theme_bw() +
    theme(legend.title = element_blank())
  if(i == 1){
    p_list <- list(p)
  } else{
    p_list <- c(p_list, list(p))
  }
}

pdf(file = "U01_Project/Results/Winsorized/DNAQC_win_cross_tissue_comparison_combinedCohorts.pdf", width = 12, height = 10)
ggarrange(plotlist = p_list,
          ncol = 3,
          nrow = 3,
          labels = "AUTO",
          legend = "bottom",
          common.legend = T)
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plots of DNA QC metrics cross-tissue comparison (winsorized, split by cohort) # -----
cohort <- levels(dataSW$Cohort)
dna_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
             "X260.280", "X260.230", 
             "Nano", "Pico", "Tape")
dna_var <- paste0(dna_var, "_win")
dna_var_label <- c("DIN", "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
                   "A260/A280", "A260/A230",
                   "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")

for (i in 1:length(dna_var)){
  p <- ggplot(data = dataSW, aes(x = Tissue, y = .data[[dna_var[i]]])) +
    geom_violin_pattern(aes(x = Tissue, y = .data[[dna_var[i]]], pattern = Cohort),
                        scale = "width") +
    geom_boxplot_pattern(aes(x = Tissue, y = .data[[dna_var[i]]],
                             fill = Tissue, pattern = Cohort),
                         pattern_fill = "white",
                         width = 0.3, fatten = 3, alpha = 0.7,
                         position = position_dodge(width = 0.9, preserve = "total")) +
    scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
    scale_color_manual(values=c("red","blue","green4","purple","orange")) +
    labs(y = dna_var_label[i],
         title = dna_var_label[i]) +
    theme_bw() +
    theme(legend.title = element_blank())
  if(i == 1){
    p_list <- list(p)
  } else{
    p_list <- c(p_list, list(p))
  }
}

pdf(file = "U01_Project/Results/Winsorized/DNAQC_win_cross_tissue_comparison.pdf", width = 14, height = 10)
ggarrange(plotlist = p_list,
          ncol = 3,
          nrow = 3,
          labels = "AUTO",
          legend = "bottom",
          common.legend = T)
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plots of DNA QC metrics cross-tissue comparison (not winsorized, combined cohorts) # -----
cohort <- levels(dataSW$Cohort)
dna_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
             "X260.280", "X260.230", 
             "Nano", "Pico", "Tape")
dna_var_label <- c("DIN", "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
                   "A260/A280", "A260/A230",
                   "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")

for (i in 1:length(dna_var)){
  p <- ggplot(data = dataSW, aes(x = Tissue, y = .data[[dna_var[i]]])) +
    geom_violin(aes(x = Tissue, y = .data[[dna_var[i]]]), scale = "width") +
    geom_boxplot(aes(x = Tissue, y = .data[[dna_var[i]]], fill = Tissue),
                 width = 0.3, fatten = 3, alpha = 0.7,
                 position = position_dodge(width = 0.9, preserve = "total")) +
    scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
    scale_color_manual(values=c("red","blue","green4","purple","orange")) +
    labs(y = dna_var_label[i],
         title = dna_var_label[i]) +
    theme_bw() +
    theme(legend.title = element_blank())
  if(i == 1){
    p_list <- list(p)
  } else{
    p_list <- c(p_list, list(p))
  }
}

pdf(file = "U01_Project/Results/Not_winsorized/DNAQC_nowin_cross_tissue_comparison_combinedCohorts.pdf", width = 12, height = 10)
ggarrange(plotlist = p_list,
          ncol = 3,
          nrow = 3,
          labels = "AUTO",
          legend = "bottom",
          common.legend = T)
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plots of DNA QC metrics cross-tissue comparison (not winsorized, split by cohort) # -----
cohort <- levels(dataSW$Cohort)
dna_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
             "X260.280", "X260.230", 
             "Nano", "Pico", "Tape")
dna_var_label <- c("DIN", "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
                   "A260/A280", "A260/A230",
                   "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")

for (i in 1:length(dna_var)){
  p <- ggplot(data = dataSW, aes(x = Tissue, y = .data[[dna_var[i]]])) +
    geom_violin_pattern(aes(x = Tissue, y = .data[[dna_var[i]]], pattern = Cohort),
                        scale = "width") +
    geom_boxplot_pattern(aes(x = Tissue, y = .data[[dna_var[i]]],
                             fill = Tissue, pattern = Cohort),
                         pattern_fill = "white",
                         width = 0.3, fatten = 3, alpha = 0.7,
                         position = position_dodge(width = 0.9, preserve = "total")) +
    scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
    scale_color_manual(values=c("red","blue","green4","purple","orange")) +
    labs(y = dna_var_label[i],
         title = dna_var_label[i]) +
    theme_bw() +
    theme(legend.title = element_blank())
  if(i == 1){
    p_list <- list(p)
  } else{
    p_list <- c(p_list, list(p))
  }
}

pdf(file = "U01_Project/Results/Not_winsorized/DNAQC_nowin_cross_tissue_comparison.pdf", width = 14, height = 10)
ggarrange(plotlist = p_list,
          ncol = 3,
          nrow = 3,
          labels = "AUTO",
          legend = "bottom",
          common.legend = T)
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))




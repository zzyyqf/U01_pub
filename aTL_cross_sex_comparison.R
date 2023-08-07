# U01 Project
# aTL cross-sex comparison
# Aug 4, 2023

# Load packages
library(tidyverse)
library(ggpattern)
library(ggstatsplot)
library(statsExpressions)
library(ggsignif)
library(ggpubr)
library(ggmagnify)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Rename the tissues
dataSW$Tissue <- factor(dataSW$Tissue, levels = c("Buccal","Saliva","DBS","Buffy","PBMC"),
                        labels = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"))

# Violin plot of aTL across tissue (winsorized, no stats, split by sex and tissue) # -----
p <- ggplot(data = dataSW, aes(x = Sex, y = aTL_win)) +
  geom_violin(aes(x = Sex, y = aTL_win), scale = "width") +
  geom_boxplot(aes(x = Sex, y = aTL_win, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "aTL (kb)") + 
  facet_wrap(~Tissue, nrow = 1) +
  theme_bw() +
  theme(legend.title = element_blank())
  

pdf(file = "U01_Project/Results/Winsorized/aTLwin_cross_sex_tissue_comparison.pdf", width = 8, height = 4)
p
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plot of aTL across tissue (not winsorized, no stats, split by sex and tissue) # -----
p <- ggplot(data = dataSW, aes(x = Sex, y = aTL)) +
  geom_violin(aes(x = Sex, y = aTL), scale = "width") +
  geom_boxplot(aes(x = Sex, y = aTL, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "aTL (kb)") + 
  facet_wrap(~Tissue, nrow = 1) +
  theme_bw() +
  theme(legend.title = element_blank())


pdf(file = "U01_Project/Results/Not_winsorized/aTLnowin_cross_sex_tissue_comparison.pdf", width = 8, height = 4)
p
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))




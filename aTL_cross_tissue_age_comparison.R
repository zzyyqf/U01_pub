# U01 Project
# aTL cross-tissue-age comparison
# June 16, 2023

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

# Make the Cohort as a factor
dataSW$Cohort <- factor(dataSW$Cohort, levels = c("P50", "TL"),
                        labels = c("Child", "Adult"))
summarytools::freq(dataSW$Cohort)

# Rename the tissues
dataSW$Tissue <- factor(dataSW$Tissue, levels = c("Buccal","Saliva","DBS","Buffy","PBMC"),
                        labels = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"))

# Violin plot and scatter plot of aTL across tissue and age (winsorized, no stats, combined cohorts) # -----
p_scatter <- ggplot(data = dataSW, aes(x=Age, y=aTL_win, color=Tissue, shape=Tissue, fill = Tissue)) + 
  geom_point(na.rm=T) +
  geom_smooth(method="lm", na.rm=T) +
  theme_bw() + 
  labs(x="Age (years)",y="aTL (kb)") + 
  theme(legend.title = element_blank()) + 
  theme(legend.key.size = unit(0.5,'cm')) +
  scale_fill_manual(breaks = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(breaks = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),values=c("red","blue","green4","purple","orange")) +
  scale_shape_manual(breaks = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),values=c(0,1,2,3,4,5)) +
  scale_x_continuous(limits = c(5, 100), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(limits = c(0, 30)) +
  geom_magnify(from = c(8, 16, 0, 25), to = c(60, 100, 15, 30), shadow = T,
               recompute = T, expand = 0, aspect = "free", axes = "xy")

p_violin <- ggplot(data = dataSW, aes(x = Tissue, y = aTL_win)) +
  geom_violin(aes(x = Tissue, y = aTL_win), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = aTL_win, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(name = "", values=c("red","blue","green4","purple","orange")) + 
  labs(y = "aTL (kb)") + 
  theme_bw()

pdf(file = "U01_Project/Results/Winsorized/aTLwin_cross_tissue_age_comparison_combinedCohorts.pdf", width = 12, height = 5)
ggarrange(p_scatter, p_violin,
          ncol = 2,
          nrow = 1,
          labels = "AUTO",
          legend = "right",
          common.legend = T,
          widths = c(1.5, 1))
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plot of aTL across tissue (winsorized, no stats, split by cohort) # -----
p <- ggplot(data = dataSW, aes(x = Tissue, y = aTL_win)) +
  geom_violin_pattern(aes(x = Tissue, y = aTL_win, pattern = Cohort),
                      scale = "width") +
  geom_boxplot_pattern(aes(x = Tissue, y = aTL_win,
                           fill = Tissue, pattern = Cohort),
                       pattern_fill = "white",
                       width = 0.3, fatten = 3, alpha = 0.7,
                       position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "aTL (kb)") + 
  theme_bw()

pdf(file = "U01_Project/Results/Winsorized/aTLwin_cross_tissue_comparison.pdf", width = 7, height = 5)
cowplot::plot_grid(p)
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plot and scatter plot of aTL across tissue and age (not winsorized, no stats, combined cohorts) # -----
p_scatter <- ggplot(data = dataSW, aes(x=Age, y=aTL, color=Tissue, shape=Tissue, fill = Tissue)) + 
  geom_point(na.rm=T) +
  geom_smooth(method="lm", na.rm=T) +
  theme_bw() + 
  labs(x="Age (years)",y="aTL (kb)") + 
  theme(legend.title = element_blank()) + 
  theme(legend.key.size = unit(0.5,'cm')) +
  scale_fill_manual(breaks = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(breaks = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),values=c("red","blue","green4","purple","orange")) +
  scale_shape_manual(breaks = c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),values=c(0,1,2,3,4,5)) +
  scale_x_continuous(limits = c(5, 80), breaks = seq(0, 80, by = 10)) +
  scale_y_continuous(limits = c(0, 38)) +
  geom_magnify(from = c(8, 16, 0, 38), to = c(40, 80, 22, 38), shadow = T,
               recompute = T, expand = 0, aspect = "free", axes = "xy")


p_violin <- ggplot(data = dataSW, aes(x = Tissue, y = aTL)) +
  geom_violin(aes(x = Tissue, y = aTL), scale = "width") +
  geom_boxplot(aes(x = Tissue, y = aTL, fill = Tissue),
               width = 0.3, fatten = 3, alpha = 0.7,
               position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(name = "", values=c("red","blue","green4","purple","orange")) + 
  labs(y = "aTL (kb)") + 
  theme_bw()

pdf(file = "U01_Project/Results/Not_winsorized/aTLnowin_cross_tissue_age_comparison_combinedCohorts.pdf", width = 12, height = 5)
ggarrange(p_scatter, p_violin,
          ncol = 2,
          nrow = 1,
          labels = "AUTO",
          legend = "right",
          common.legend = T,
          widths = c(1.5, 1))
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Violin plot of aTL across tissue (not winsorized, no stats, split by cohort) # -----
p <- ggplot(data = dataSW, aes(x = Tissue, y = aTL)) +
  geom_violin_pattern(aes(x = Tissue, y = aTL, pattern = Cohort),
                      scale = "width") +
  geom_boxplot_pattern(aes(x = Tissue, y = aTL,
                           fill = Tissue, pattern = Cohort),
                       pattern_fill = "white",
                       width = 0.3, fatten = 3, alpha = 0.7,
                       position = position_dodge(width = 0.9, preserve = "total")) +
  scale_fill_manual(values=c("red","blue","green4","purple","orange")) + 
  scale_color_manual(values=c("red","blue","green4","purple","orange")) +
  labs(y = "aTL (kb)") + 
  theme_bw()

pdf(file = "U01_Project/Results/Not_winsorized/aTLnowin_cross_tissue_comparison.pdf", width = 7, height = 5)
cowplot::plot_grid(p)
dev.off()

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))


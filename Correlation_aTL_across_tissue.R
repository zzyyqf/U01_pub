# U01 Project
# Correlation of aTL across tissues
# June 16, 2023

# Load packages
library(tidyverse)
library(correlation)
library(psych)
library(corrplot)
source("U01_Project/Code/Rfunctions/corrplot_update.R") # updated corrplot() based on https://stackoverflow.com/questions/49695113/changing-the-position-of-the-signifigance-pch-symbols-in-corrplot

# Read in data
load("U01_Project/Data/u01_data.RData")

# Correlation of aTL (winsorized) across tissues (separate cohorts) #-----
#generate correlation matrix for plotting (child)
dataSW_wide_child <- dataSW_wide %>% filter(Cohort == "P50")
tl_win_child_var <- c("aTL_win_Buccal", "aTL_win_Saliva", "aTL_win_DBS", "aTL_win_Buffy")
cor_res_win_child <- data.frame(NULL)
for (i in 1:length(tl_win_child_var)) {
  for(j in 1:length(tl_win_child_var)) {
    if (i == j){
      next
    } else {
      data_temp <- dataSW_wide_child %>%
        select(all_of(c("Age", "Sex", tl_win_child_var[i], tl_win_child_var[j])))
      
      cor_temp <- correlation::cor_test(data = data_temp %>% na.omit(), x = tl_win_child_var[i], y = tl_win_child_var[j],
                                        method = "spearman", include_factors = T, partial = T,
                                        multilevel = F)
      
      cor_res_win_child <- rbind(cor_res_win_child,
                                 as.data.frame(cor_temp) %>%
                                   mutate(Cohort = "P50"))
    }
    
  }
}

#generate correlation matrix for plotting (adult)
dataSW_wide_adult <- dataSW_wide %>% filter(Cohort == "TL")
tl_win_adult_var <- c("aTL_win_Buccal", "aTL_win_Saliva", "aTL_win_DBS", "aTL_win_PBMC")
cor_res_win_adult <- data.frame(NULL)
for (i in 1:length(tl_win_adult_var)) {
  for(j in 1:length(tl_win_adult_var)) {
    if (i == j){
      next
    } else {
      data_temp <- dataSW_wide_adult %>%
        select(all_of(c("Age", "Sex", tl_win_adult_var[i], tl_win_adult_var[j])))
      
      cor_temp <- correlation::cor_test(data = data_temp %>% na.omit(), x = tl_win_adult_var[i], y = tl_win_adult_var[j],
                                        method = "spearman", include_factors = T, partial = T,
                                        multilevel = F)
      
      cor_res_win_adult <- rbind(cor_res_win_adult,
                                 as.data.frame(cor_temp) %>%
                                   mutate(Cohort = "TL"))
    }
    
  }
}

#adjust p values for multiple testing
cor_res_win <- rbind(cor_res_win_child, cor_res_win_adult)
cor_res_win$p.adj <- p.adjust(cor_res_win$p, method = "BH")
cor_res_win_child <- cor_res_win %>% filter(Cohort == "P50")
cor_res_win_adult <- cor_res_win %>% filter(Cohort == "TL")

#format the results into a correlation matrix (child)
colnames(cor_res_win_child)[c(1,2)] <- c("var1", "var2")
cor_res_win_child_rho_matrix <- as.matrix(rstatix::cor_spread(cor_res_win_child, value = "rho"))
rownames(cor_res_win_child_rho_matrix) <- cor_res_win_child_rho_matrix[, 1]
cor_res_win_child_rho_matrix <- cor_res_win_child_rho_matrix[, -1]
cor_res_win_child_rho_matrix <- cor_res_win_child_rho_matrix[, rownames(cor_res_win_child_rho_matrix)]
cor_res_win_child_rho_matrix <- apply(cor_res_win_child_rho_matrix, 2, as.numeric)
colnames(cor_res_win_child_rho_matrix) <- str_split_i(colnames(cor_res_win_child_rho_matrix), pattern = "_", i = -1)
colnames(cor_res_win_child_rho_matrix)[which(colnames(cor_res_win_child_rho_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_win_child_rho_matrix) <- colnames(cor_res_win_child_rho_matrix)

cor_res_win_child_p_matrix <- as.matrix(rstatix::cor_spread(cor_res_win_child, value = "p.adj"))
rownames(cor_res_win_child_p_matrix) <- cor_res_win_child_p_matrix[, 1]
cor_res_win_child_p_matrix <- cor_res_win_child_p_matrix[, -1]
cor_res_win_child_p_matrix <- cor_res_win_child_p_matrix[, rownames(cor_res_win_child_p_matrix)]
cor_res_win_child_p_matrix <- apply(cor_res_win_child_p_matrix, 2, as.numeric)
colnames(cor_res_win_child_p_matrix) <- str_split_i(colnames(cor_res_win_child_p_matrix), pattern = "_", i = -1)
colnames(cor_res_win_child_p_matrix)[which(colnames(cor_res_win_child_p_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_win_child_p_matrix) <- colnames(cor_res_win_child_p_matrix)

#format the results into a correlation matrix (adult)
colnames(cor_res_win_adult)[c(1,2)] <- c("var1", "var2")
cor_res_win_adult_rho_matrix <- as.matrix(rstatix::cor_spread(cor_res_win_adult, value = "rho"))
rownames(cor_res_win_adult_rho_matrix) <- cor_res_win_adult_rho_matrix[, 1]
cor_res_win_adult_rho_matrix <- cor_res_win_adult_rho_matrix[, -1]
cor_res_win_adult_rho_matrix <- cor_res_win_adult_rho_matrix[, rownames(cor_res_win_adult_rho_matrix)]
cor_res_win_adult_rho_matrix <- apply(cor_res_win_adult_rho_matrix, 2, as.numeric)
colnames(cor_res_win_adult_rho_matrix) <- str_split_i(colnames(cor_res_win_adult_rho_matrix), pattern = "_", i = -1)
rownames(cor_res_win_adult_rho_matrix) <- colnames(cor_res_win_adult_rho_matrix)

cor_res_win_adult_p_matrix <- as.matrix(rstatix::cor_spread(cor_res_win_adult, value = "p.adj"))
rownames(cor_res_win_adult_p_matrix) <- cor_res_win_adult_p_matrix[, 1]
cor_res_win_adult_p_matrix <- cor_res_win_adult_p_matrix[, -1]
cor_res_win_adult_p_matrix <- cor_res_win_adult_p_matrix[, rownames(cor_res_win_adult_p_matrix)]
cor_res_win_adult_p_matrix <- apply(cor_res_win_adult_p_matrix, 2, as.numeric)
colnames(cor_res_win_adult_p_matrix) <- str_split_i(colnames(cor_res_win_adult_p_matrix), pattern = "_", i = -1)
rownames(cor_res_win_adult_p_matrix) <- colnames(cor_res_win_adult_p_matrix)

# Combine the child and adult plots together
pdf("U01_project/Results/Winsorized/Correlation_aTLwin_across_tissue.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
par(xpd=TRUE)
corrplot_update(cor_res_win_child_rho_matrix, diag=F,type="lower",tl.col="black", number.cex = 1,tl.srt=0,
                method="ellipse", col=colorRampPalette(c("steelblue","white","tomato3"))(200), 
                outline=T, tl.offset = 1, cl.ratio = 0.3,
                tl.cex=1, cl.cex=1, addCoef.col = "black", mar = c(0,0,2,0),
                p.mat = cor_res_win_child_p_matrix, insig = "label_sig", sig.level = 0.01,
                title = "Child Cohort")

corrplot_update(cor_res_win_adult_rho_matrix, diag=F,type="lower",tl.col="black", number.cex = 1,tl.srt=0,
                method="ellipse", col=colorRampPalette(c("steelblue","white","tomato3"))(200), 
                outline=T, tl.offset = 1, cl.ratio = 0.3,
                tl.cex=1, cl.cex=1, addCoef.col = "black", mar = c(0,0,2,0),
                p.mat = cor_res_win_adult_p_matrix, insig = "label_sig", sig.level = 0.01,
                title = "Adult Cohort")
# title(sub = "Partial Spearman's r accounting for age and sex.\n*Statistical significant controlling FDR at <0.01", adj = 0.75, cex.sub = 0.8)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "corrplot_update", "draw_grid",
                          "draw_method_color", "draw_method_square")))

# Correlation of aTL (winsorized) across tissues (combined cohorts) # ------
#generate correlation matrix for plotting
tl_var <- c("aTL_win_Buccal", "aTL_win_Saliva", "aTL_win_DBS", "aTL_win_Buffy", "aTL_win_PBMC")
cor_res <- data.frame(NULL)
for (i in 1:length(tl_var)) {
  for(j in 1:length(tl_var)) {
    if (i == j){
      next
    } else if (tl_var[i] == "aTL_win_Buffy" & tl_var[j] == "aTL_win_PBMC"){
      next
    } else if (tl_var[i] == "aTL_win_PBMC" & tl_var[j] == "aTL_win_Buffy") {
      break
    } else {
      data_temp <- dataSW_wide %>%
        select(all_of(c("Age", "Sex", tl_var[i], tl_var[j])))
      
      cor_temp <- correlation::cor_test(data = data_temp %>% na.omit(), x = tl_var[i], y = tl_var[j],
                                        method = "spearman", include_factors = T, partial = T,
                                        multilevel = F)
      
      cor_res <- rbind(cor_res, as.data.frame(cor_temp))
    }
    
  }
}

#adjust p values for multiple testing
cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#format the results into a correlation matrix for plotting
colnames(cor_res)[c(1,2)] <- c("var1", "var2")
cor_res_rho_matrix <- as.matrix(rstatix::cor_spread(cor_res, value = "rho"))
rownames(cor_res_rho_matrix) <- cor_res_rho_matrix[, 1]
cor_res_rho_matrix <- cor_res_rho_matrix[, -1]
cor_res_rho_matrix <- cor_res_rho_matrix[, rownames(cor_res_rho_matrix)]
cor_res_rho_matrix <- apply(cor_res_rho_matrix, 2, as.numeric)
colnames(cor_res_rho_matrix) <- str_split_i(colnames(cor_res_rho_matrix), pattern = "_", i = -1)
colnames(cor_res_rho_matrix)[which(colnames(cor_res_rho_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_rho_matrix) <- colnames(cor_res_rho_matrix)

cor_res_p_matrix <- as.matrix(rstatix::cor_spread(cor_res, value = "p.adj"))
rownames(cor_res_p_matrix) <- cor_res_p_matrix[, 1]
cor_res_p_matrix <- cor_res_p_matrix[, -1]
cor_res_p_matrix <- cor_res_p_matrix[, rownames(cor_res_p_matrix)]
cor_res_p_matrix <- apply(cor_res_p_matrix, 2, as.numeric)
colnames(cor_res_p_matrix) <- str_split_i(colnames(cor_res_p_matrix), pattern = "_", i = -1)
colnames(cor_res_p_matrix)[which(colnames(cor_res_p_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_p_matrix) <- colnames(cor_res_p_matrix)

#make correlation plot
pdf("U01_project/Results/Winsorized/Correlation_aTLwin_across_tissue_combinedCohorts.pdf", width = 6, height = 6)
par(xpd=TRUE)
corrplot_update(cor_res_rho_matrix[,-(4:5)], diag=F,type="lower",tl.col="black", number.cex = 1, tl.srt=0,
                method="ellipse", col=colorRampPalette(c("steelblue","white","tomato3"))(200), 
                outline=T, tl.offset = 1, cl.ratio = 0.3,
                tl.cex=1, cl.cex=0.8, addCoef.col = "black", mar = c(0,0,2,0),
                p.mat = cor_res_p_matrix[,-(4:5)], insig = "label_sig", sig.level = 0.01)
# title(sub = "Partial Spearman's r accounting for age and sex.\n*Statistical significant controlling FDR at <0.01", adj = 0.75, cex.sub = 0.8, outer = F)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "corrplot_update", "draw_grid",
                          "draw_method_color", "draw_method_square")))

# Correlation of aTL (not winsorized) across tissues (separate cohorts) #-----
#generate correlation matrix for plotting (child)
dataSW_wide_child <- dataSW_wide %>% filter(Cohort == "P50")
tl_nowin_child_var <- c("aTL_Buccal", "aTL_Saliva", "aTL_DBS", "aTL_Buffy")
cor_res_nowin_child <- data.frame(NULL)
for (i in 1:length(tl_nowin_child_var)) {
  for(j in 1:length(tl_nowin_child_var)) {
    if (i == j){
      next
    } else {
      data_temp <- dataSW_wide_child %>%
        select(all_of(c("Age", "Sex", tl_nowin_child_var[i], tl_nowin_child_var[j])))
      
      cor_temp <- correlation::cor_test(data = data_temp %>% na.omit(), x = tl_nowin_child_var[i], y = tl_nowin_child_var[j],
                                        method = "spearman", include_factors = T, partial = T,
                                        multilevel = F)
      
      cor_res_nowin_child <- rbind(cor_res_nowin_child,
                                 as.data.frame(cor_temp) %>%
                                   mutate(Cohort = "P50"))
    }
    
  }
}

#generate correlation matrix for plotting (adult)
dataSW_wide_adult <- dataSW_wide %>% filter(Cohort == "TL")
tl_nowin_adult_var <- c("aTL_Buccal", "aTL_Saliva", "aTL_DBS", "aTL_PBMC")
cor_res_nowin_adult <- data.frame(NULL)
for (i in 1:length(tl_nowin_adult_var)) {
  for(j in 1:length(tl_nowin_adult_var)) {
    if (i == j){
      j <- j + 1
    } else {
      data_temp <- dataSW_wide_adult %>%
        select(all_of(c("Age", "Sex", tl_nowin_adult_var[i], tl_nowin_adult_var[j])))
      
      cor_temp <- correlation::cor_test(data = data_temp %>% na.omit(), x = tl_nowin_adult_var[i], y = tl_nowin_adult_var[j],
                                        method = "spearman", include_factors = T, partial = T,
                                        multilevel = F)
      
      cor_res_nowin_adult <- rbind(cor_res_nowin_adult,
                                   as.data.frame(cor_temp) %>%
                                     mutate(Cohort = "TL"))
    }
    
  }
}

#adjust p values for multiple testing
cor_res_nowin <- rbind(cor_res_nowin_child, cor_res_nowin_adult)
cor_res_nowin$p.adj <- p.adjust(cor_res_nowin$p, method = "BH")
cor_res_nowin_child <- cor_res_nowin %>% filter(Cohort == "P50")
cor_res_nowin_adult <- cor_res_nowin %>% filter(Cohort == "TL")

#format the results into a correlation matrix (child)
colnames(cor_res_nowin_child)[c(1,2)] <- c("var1", "var2")
cor_res_nowin_child_rho_matrix <- as.matrix(rstatix::cor_spread(cor_res_nowin_child, value = "rho"))
rownames(cor_res_nowin_child_rho_matrix) <- cor_res_nowin_child_rho_matrix[, 1]
cor_res_nowin_child_rho_matrix <- cor_res_nowin_child_rho_matrix[, -1]
cor_res_nowin_child_rho_matrix <- cor_res_nowin_child_rho_matrix[, rownames(cor_res_nowin_child_rho_matrix)]
cor_res_nowin_child_rho_matrix <- apply(cor_res_nowin_child_rho_matrix, 2, as.numeric)
colnames(cor_res_nowin_child_rho_matrix) <- str_split_i(colnames(cor_res_nowin_child_rho_matrix), pattern = "_", i = -1)
colnames(cor_res_nowin_child_rho_matrix)[which(colnames(cor_res_nowin_child_rho_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_nowin_child_rho_matrix) <- colnames(cor_res_nowin_child_rho_matrix)

cor_res_nowin_child_p_matrix <- as.matrix(rstatix::cor_spread(cor_res_nowin_child, value = "p.adj"))
rownames(cor_res_nowin_child_p_matrix) <- cor_res_nowin_child_p_matrix[, 1]
cor_res_nowin_child_p_matrix <- cor_res_nowin_child_p_matrix[, -1]
cor_res_nowin_child_p_matrix <- cor_res_nowin_child_p_matrix[, rownames(cor_res_nowin_child_p_matrix)]
cor_res_nowin_child_p_matrix <- apply(cor_res_nowin_child_p_matrix, 2, as.numeric)
colnames(cor_res_nowin_child_p_matrix) <- str_split_i(colnames(cor_res_nowin_child_p_matrix), pattern = "_", i = -1)
colnames(cor_res_nowin_child_p_matrix)[which(colnames(cor_res_nowin_child_p_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_nowin_child_p_matrix) <- colnames(cor_res_nowin_child_p_matrix)

#format the results into a correlation matrix (adult)
colnames(cor_res_nowin_adult)[c(1,2)] <- c("var1", "var2")
cor_res_nowin_adult_rho_matrix <- as.matrix(rstatix::cor_spread(cor_res_nowin_adult, value = "rho"))
rownames(cor_res_nowin_adult_rho_matrix) <- cor_res_nowin_adult_rho_matrix[, 1]
cor_res_nowin_adult_rho_matrix <- cor_res_nowin_adult_rho_matrix[, -1]
cor_res_nowin_adult_rho_matrix <- cor_res_nowin_adult_rho_matrix[, rownames(cor_res_nowin_adult_rho_matrix)]
cor_res_nowin_adult_rho_matrix <- apply(cor_res_nowin_adult_rho_matrix, 2, as.numeric)
colnames(cor_res_nowin_adult_rho_matrix) <- str_split_i(colnames(cor_res_nowin_adult_rho_matrix), pattern = "_", i = -1)
rownames(cor_res_nowin_adult_rho_matrix) <- colnames(cor_res_nowin_adult_rho_matrix)

cor_res_nowin_adult_p_matrix <- as.matrix(rstatix::cor_spread(cor_res_nowin_adult, value = "p.adj"))
rownames(cor_res_nowin_adult_p_matrix) <- cor_res_nowin_adult_p_matrix[, 1]
cor_res_nowin_adult_p_matrix <- cor_res_nowin_adult_p_matrix[, -1]
cor_res_nowin_adult_p_matrix <- cor_res_nowin_adult_p_matrix[, rownames(cor_res_nowin_adult_p_matrix)]
cor_res_nowin_adult_p_matrix <- apply(cor_res_nowin_adult_p_matrix, 2, as.numeric)
colnames(cor_res_nowin_adult_p_matrix) <- str_split_i(colnames(cor_res_nowin_adult_p_matrix), pattern = "_", i = -1)
rownames(cor_res_nowin_adult_p_matrix) <- colnames(cor_res_nowin_adult_p_matrix)

# Combine the child and adult plots together
pdf("U01_project/Results/Not_winsorized/Correlation_aTLnowin_across_tissue.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
par(xpd=TRUE)
corrplot_update(cor_res_nowin_child_rho_matrix, diag=F,type="lower",tl.col="black", number.cex = .8,tl.srt=0,
                method="ellipse", col=colorRampPalette(c("steelblue","white","tomato3"))(200), 
                outline=T, tl.offset = 1, cl.ratio = 0.3,
                tl.cex=1, cl.cex=1, addCoef.col = "black", mar = c(0,0,2,0),
                p.mat = cor_res_nowin_child_p_matrix, insig = "label_sig", sig.level = 0.01,
                title = "Child Cohort")

corrplot_update(cor_res_nowin_adult_rho_matrix, diag=F,type="lower",tl.col="black", number.cex = .8,tl.srt=0,
                method="ellipse", col=colorRampPalette(c("steelblue","white","tomato3"))(200), 
                outline=T, tl.offset = 1, cl.ratio = 0.3,
                tl.cex=1, cl.cex=1, addCoef.col = "black", mar = c(0,0,2,0),
                p.mat = cor_res_nowin_adult_p_matrix, insig = "label_sig", sig.level = 0.01,
                title = "Adult Cohort")
# title(sub = "Partial Spearman's r accounting for age and sex.\n*Statistical significant controlling FDR at <0.01", adj = 0.75, cex.sub = 0.8)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "corrplot_update", "draw_grid",
                          "draw_method_color", "draw_method_square")))

# Correlation of aTL (not winsorized) across tissues (combined cohorts) # ------
#generate correlation matrix for plotting
tl_var <- c("aTL_Buccal", "aTL_Saliva", "aTL_DBS", "aTL_Buffy", "aTL_PBMC")
cor_res <- data.frame(NULL)
for (i in 1:length(tl_var)) {
  for(j in 1:length(tl_var)) {
    if (i == j){
      next
    } else if (tl_var[i] == "aTL_Buffy" & tl_var[j] == "aTL_PBMC"){
      next
    } else if (tl_var[i] == "aTL_PBMC" & tl_var[j] == "aTL_Buffy") {
      break
    } else {
      data_temp <- dataSW_wide %>%
        select(all_of(c("Age", "Sex", tl_var[i], tl_var[j])))
      
      cor_temp <- correlation::cor_test(data = data_temp %>% na.omit(), x = tl_var[i], y = tl_var[j],
                                        method = "spearman", include_factors = T, partial = T,
                                        multilevel = F)
      
      cor_res <- rbind(cor_res, as.data.frame(cor_temp))
    }
    
  }
}

#adjust p values for multiple testing
cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")

#format the results into a correlation matrix
colnames(cor_res)[c(1,2)] <- c("var1", "var2")
cor_res_rho_matrix <- as.matrix(rstatix::cor_spread(cor_res, value = "rho"))
rownames(cor_res_rho_matrix) <- cor_res_rho_matrix[, 1]
cor_res_rho_matrix <- cor_res_rho_matrix[, -1]
cor_res_rho_matrix <- cor_res_rho_matrix[, rownames(cor_res_rho_matrix)]
cor_res_rho_matrix <- apply(cor_res_rho_matrix, 2, as.numeric)
colnames(cor_res_rho_matrix) <- str_split_i(colnames(cor_res_rho_matrix), pattern = "_", i = -1)
colnames(cor_res_rho_matrix)[which(colnames(cor_res_rho_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_rho_matrix) <- colnames(cor_res_rho_matrix)

cor_res_p_matrix <- as.matrix(rstatix::cor_spread(cor_res, value = "p.adj"))
rownames(cor_res_p_matrix) <- cor_res_p_matrix[, 1]
cor_res_p_matrix <- cor_res_p_matrix[, -1]
cor_res_p_matrix <- cor_res_p_matrix[, rownames(cor_res_p_matrix)]
cor_res_p_matrix <- apply(cor_res_p_matrix, 2, as.numeric)
colnames(cor_res_p_matrix) <- str_split_i(colnames(cor_res_p_matrix), pattern = "_", i = -1)
colnames(cor_res_p_matrix)[which(colnames(cor_res_p_matrix) == "Buffy")] <- "Buffy Coat"
rownames(cor_res_p_matrix) <- colnames(cor_res_p_matrix)

#make correlation plot
pdf("U01_project/Results/Not_winsorized/Correlation_aTLnowin_across_tissue_combinedCohorts.pdf", width = 6, height = 6)
par(xpd=TRUE)
corrplot_update(cor_res_rho_matrix[,-(4:5)], diag=F,type="lower",tl.col="black", number.cex = 1, tl.srt=0,
                method="ellipse", col=colorRampPalette(c("steelblue","white","tomato3"))(200), 
                outline=T, tl.offset = 1, cl.ratio = 0.3,
                tl.cex=1, cl.cex=0.8, addCoef.col = "black", mar = c(0,0,2,0),
                p.mat = cor_res_p_matrix[,-(4:5)], insig = "label_sig", sig.level = 0.01)
# title(sub = "Partial Spearman's r accounting for age and sex.\n*Statistical significant controlling FDR at <0.01", adj = 0.75, cex.sub = 0.8, outer = F)
dev.off()

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "corrplot_update", "draw_grid",
                          "draw_method_color", "draw_method_square")))


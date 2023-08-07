# U01 Project
# Linear Mixed Effect Models for aTL and DNA QC metrics variation by different factors
# June 19, 2023

# Load packages
library(tidyverse)
library(nlme)
library(lme4)
library(emmeans)
library(openxlsx)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Make the Cohort variable as a factor
dataSW$Cohort <- factor(dataSW$Cohort, levels = c("P50", "TL"), labels = c("Child", "Adult"))
summarytools::freq(dataSW$Cohort)

# Change the name of the EthCorr to Race
summarytools::freq(dataSW$EthCorr)
dataSW <- dataSW %>%
  mutate(Race = EthCorr)
summarytools::freq(dataSW$Race)

# Modeling with the combined cohorts (winsorized; 5 tissues).# -----
dv <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
        "X260.280", "X260.230",
        "Nano", "Pico", "Tape")
dv <- paste0(dv, "_win")
dv_label <- c("TL (kb)", "DIN", 
              "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
              "A260/A280", "A260/230",
              "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")
m.data <- dataSW
post_hoc_res <- list()
for (i in 1:length(dv)){
  m.formula <- as.formula(paste0(dv[i], "~ Age + Sex + Tissue + Race + Age:Tissue"))
  mod <- lme(m.formula, data=m.data, random=~1|ID, na.action=na.omit)
  p <- plot(mod, main = dv_label[i])
  m_coef_temp <- as_tibble(coef(summary(mod)), rownames = "Var_d") %>%
    mutate(Value = as.character(round(Value, 3)),
           Std.Error = round(Std.Error, 3),
           Pvalue = `p-value`) %>%
    select(-all_of(c("DF", "t-value", "p-value")))
  
  m_coef_temp <- m_coef_temp %>%
    add_row(Var_d = paste0("Sex", levels(m.data$Sex)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Sex", levels(m.data$Sex)[2]))) %>%
    add_row(Var_d = paste0("Tissue", levels(m.data$Tissue)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Tissue", levels(m.data$Tissue)[2]))+1) %>%
    add_row(Var_d = paste0("Race", levels(m.data$Race)[1]),
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Race", levels(m.data$Race)[2]))+2) %>%
    mutate(dv = dv[i])

  m_anova_temp <- as_tibble(anova(mod), rownames = "Var_c") %>% 
    transmute(Var_c = Var_c,
              anova_numDF = numDF,
              anova_denDF = denDF,
              anova_Fvalue = round(`F-value`, 3),
              anova_Pvalue = `p-value`,
              dv = dv[i])
  
  # post-hoc analyses
  ph_race <- as.data.frame(emmeans(mod, pairwise~Race, adjust = "none")$contrast) %>% 
    mutate(dv = dv[i])

  ph_tissue <- as.data.frame(emmeans(mod, pairwise~Tissue|Age, adjust = "none")$contrasts) %>% 
    mutate(dv = dv[i])

  ph_ageTrend_tissue <- as.data.frame(emtrends(mod, pairwise~Tissue, var = "Age")$emtrends) %>%
    mutate(dv = dv[i])
  
  if(i == 1){
    p_list <- list(p)
    m_coef_tab <- m_coef_temp
    m_anova_tab <- m_anova_temp
    post_hoc_res[["diff_across_race"]] <- ph_race
    post_hoc_res[["diff_across_tissue"]] <- ph_tissue
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- ph_ageTrend_tissue
    
  } else {
    p_list <- c(p_list, list(p))
   
    m_coef_tab <- rbind(m_coef_tab, m_coef_temp)
    
    m_anova_tab <- rbind(m_anova_tab, m_anova_temp)
    
    post_hoc_res[["diff_across_race"]] <- rbind(post_hoc_res[["diff_across_race"]],
                                                ph_race)
    post_hoc_res[["diff_across_tissue"]] <- rbind(post_hoc_res[["diff_across_tissue"]],
                                                  ph_tissue)
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- rbind(post_hoc_res[["ageTrendDiff_across_tissue"]],
                                                          ph_ageTrend_tissue)
    
    
  }
}

# Adjust the p value independently for the coef tab and the anova tab
m_coef_tab$p.adj <- p.adjust(m_coef_tab$Pvalue, method = "BH")
m_anova_tab$p.adj <- p.adjust(m_anova_tab$anova_Pvalue, method = "BH")
m_coef_tab <- m_coef_tab %>%
  mutate(Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue, 3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj)
m_coef_tab_wide <- m_coef_tab %>%
  pivot_wider(names_from = dv, values_from = c(Value, Std.Error, Pvalue)) %>%
  select(all_of(c("Var_d", paste0(rep(c("Value_", "Std.Error_", "Pvalue_"), 3),
                                  rep(dv, each = 3)))))

m_anova_tab <- m_anova_tab %>%
  mutate(anova_Pvalue = ifelse(anova_Pvalue<0.001, "<0.001", round(anova_Pvalue, 3))) %>%
  mutate(anova_Pvalue = ifelse(p.adj<0.01, paste0(anova_Pvalue, "*"), anova_Pvalue)) %>%
  select(-p.adj)

m_anova_tab_wide <- m_anova_tab %>%
  mutate(anova_df = paste(anova_numDF, anova_denDF, sep = ", ")) %>%
  select(-c(anova_numDF, anova_denDF)) %>%
  pivot_wider(names_from = dv, values_from = c(anova_df, anova_Fvalue, anova_Pvalue)) %>%
  select(all_of(c("Var_c", paste0(rep(c("anova_df_", "anova_Fvalue_", "anova_Pvalue_"), 3),
                                  rep(dv, each = 3)))))

# Pick out the aTL coef and anova tables
tl_var <- "aTL_win"
tl_coef_tab <- m_coef_tab_wide %>%
  select(Var_d, all_of(paste(c("Value", "Std.Error", "Pvalue"), tl_var, sep = "_"))) %>%
  mutate(Var_c = apply(str_extract_all(Var_d, 
                                       pattern = "\\(Intercept\\)|Age|Cohort|Sex|Race|Tissue|:", 
                                       simplify = T),
                       1, function(x) paste0(x, collapse = "")))

tl_anova_tab <- m_anova_tab_wide %>%
  select(Var_c, all_of(paste(c("anova_df", "anova_Fvalue", "anova_Pvalue"), tl_var, sep = "_")))

tl_tab <- tl_coef_tab %>%
  left_join(tl_anova_tab, by = "Var_c") %>%
  select(Var_c, Var_d, all_of(paste(c("Value", "Std.Error", "Pvalue"), tl_var, sep = "_")), 
         all_of(paste(c("anova_df", "anova_Fvalue", "anova_Pvalue"), tl_var, sep = "_")))

# Adjust the p value for post-hoc analyses
post_hoc_res[["diff_across_race"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_race"]])$p.value, method = "BH")
post_hoc_res[["diff_across_race"]] <- post_hoc_res[["diff_across_race"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))
post_hoc_res[["diff_across_tissue"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_tissue"]])$p.value, method = "BH")
post_hoc_res[["diff_across_tissue"]] <- post_hoc_res[["diff_across_tissue"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))

# Plot the residuals of each model
pdf(file = "U01_Project/Results/Winsorized/LMEM_aTL_DNAQC_win_residual_plot.pdf", width = 15, height = 10)
cowplot::plot_grid(plotlist = p_list, ncol = 4)
dev.off()

# Output the results
out_list <- list("tl_tab" = tl_tab,
                 "coef_tab" = m_coef_tab_wide, 
                 "anova_tab" = m_anova_tab_wide,
                 "post_hoc_race" = post_hoc_res[["diff_across_race"]],
                 "post_hoc_tissue" = post_hoc_res[["diff_across_tissue"]],
                 "post_hoc_ageTrend_tissue" = post_hoc_res[["ageTrendDiff_across_tissue"]])

write.xlsx(out_list, file = "U01_Project/Results/Winsorized/LMEM_aTL_DNAQC_win_result.xlsx")

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))

# Modeling with the combined cohorts (not winsorized; 5 tissues).# -----
dv <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
        "X260.280", "X260.230",
        "Nano", "Pico", "Tape")
dv_label <- c("TL (kb)", "DIN", 
              "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
              "A260/A280", "A260/230",
              "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")

m.data <- dataSW
post_hoc_res <- list()
for (i in 1:length(dv)){
  m.formula <- as.formula(paste0(dv[i], "~ Age + Sex + Tissue + Race + Age:Tissue"))
  mod <- lme(m.formula, data=m.data, random=~1|ID, na.action=na.omit)
  p <- plot(mod, main = dv_label[i])
  m_coef_temp <- as_tibble(coef(summary(mod)), rownames = "Var_d") %>%
    mutate(Value = as.character(round(Value, 3)),
           Std.Error = round(Std.Error, 3),
           Pvalue = `p-value`) %>%
    select(-all_of(c("DF", "t-value", "p-value")))
  
  m_coef_temp <- m_coef_temp %>%
    add_row(Var_d = paste0("Sex", levels(m.data$Sex)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Sex", levels(m.data$Sex)[2]))) %>%
    add_row(Var_d = paste0("Tissue", levels(m.data$Tissue)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Tissue", levels(m.data$Tissue)[2]))+1) %>%
    add_row(Var_d = paste0("Race", levels(m.data$Race)[1]),
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Race", levels(m.data$Race)[2]))+2) %>%
    mutate(dv = dv[i])
  
  m_anova_temp <- as_tibble(anova(mod), rownames = "Var_c") %>% 
    transmute(Var_c = Var_c,
              anova_numDF = numDF,
              anova_denDF = denDF,
              anova_Fvalue = round(`F-value`, 3),
              anova_Pvalue = `p-value`,
              dv = dv[i])
  
  # post-hoc analyses
  ph_race <- as.data.frame(emmeans(mod, pairwise~Race, adjust = "none")$contrast) %>% 
    mutate(dv = dv[i])
  
  ph_tissue <- as.data.frame(emmeans(mod, pairwise~Tissue|Age, adjust = "none")$contrasts) %>% 
    mutate(dv = dv[i])
  
  ph_ageTrend_tissue <- as.data.frame(emtrends(mod, pairwise~Tissue, var = "Age")$emtrends) %>%
    mutate(dv = dv[i])
  
  if(i == 1){
    p_list <- list(p)
    m_coef_tab <- m_coef_temp
    m_anova_tab <- m_anova_temp
    post_hoc_res[["diff_across_race"]] <- ph_race
    post_hoc_res[["diff_across_tissue"]] <- ph_tissue
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- ph_ageTrend_tissue
    
  } else {
    p_list <- c(p_list, list(p))
    
    m_coef_tab <- rbind(m_coef_tab, m_coef_temp)
    
    m_anova_tab <- rbind(m_anova_tab, m_anova_temp)
    
    post_hoc_res[["diff_across_race"]] <- rbind(post_hoc_res[["diff_across_race"]],
                                                ph_race)
    post_hoc_res[["diff_across_tissue"]] <- rbind(post_hoc_res[["diff_across_tissue"]],
                                                  ph_tissue)
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- rbind(post_hoc_res[["ageTrendDiff_across_tissue"]],
                                                          ph_ageTrend_tissue)
    
    
  }
}

# Adjust the p value independently for the coef tab and the anova tab
m_coef_tab$p.adj <- p.adjust(m_coef_tab$Pvalue, method = "BH")
m_anova_tab$p.adj <- p.adjust(m_anova_tab$anova_Pvalue, method = "BH")
m_coef_tab <- m_coef_tab %>%
  mutate(Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue, 3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj)
m_coef_tab_wide <- m_coef_tab %>%
  pivot_wider(names_from = dv, values_from = c(Value, Std.Error, Pvalue)) %>%
  select(all_of(c("Var_d", paste0(rep(c("Value_", "Std.Error_", "Pvalue_"), 3),
                                  rep(dv, each = 3)))))

m_anova_tab <- m_anova_tab %>%
  mutate(anova_Pvalue = ifelse(anova_Pvalue<0.001, "<0.001", round(anova_Pvalue, 3))) %>%
  mutate(anova_Pvalue = ifelse(p.adj<0.01, paste0(anova_Pvalue, "*"), anova_Pvalue)) %>%
  select(-p.adj)

m_anova_tab_wide <- m_anova_tab %>%
  mutate(anova_df = paste(anova_numDF, anova_denDF, sep = ", ")) %>%
  select(-c(anova_numDF, anova_denDF)) %>%
  pivot_wider(names_from = dv, values_from = c(anova_df, anova_Fvalue, anova_Pvalue)) %>%
  select(all_of(c("Var_c", paste0(rep(c("anova_df_", "anova_Fvalue_", "anova_Pvalue_"), 3),
                                  rep(dv, each = 3)))))

# Pick out the aTL coef and anova tables
tl_var <- "aTL"
tl_coef_tab <- m_coef_tab_wide %>%
  select(Var_d, all_of(paste(c("Value", "Std.Error", "Pvalue"), tl_var, sep = "_"))) %>%
  mutate(Var_c = apply(str_extract_all(Var_d, 
                                       pattern = "\\(Intercept\\)|Age|Cohort|Sex|Race|Tissue|:", 
                                       simplify = T),
                       1, function(x) paste0(x, collapse = "")))

tl_anova_tab <- m_anova_tab_wide %>%
  select(Var_c, all_of(paste(c("anova_df", "anova_Fvalue", "anova_Pvalue"), tl_var, sep = "_")))

tl_tab <- tl_coef_tab %>%
  left_join(tl_anova_tab, by = "Var_c") %>%
  select(Var_c, Var_d, all_of(paste(c("Value", "Std.Error", "Pvalue"), tl_var, sep = "_")), 
         all_of(paste(c("anova_df", "anova_Fvalue", "anova_Pvalue"), tl_var, sep = "_")))

# Adjust the p value for post-hoc analyses
post_hoc_res[["diff_across_race"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_race"]])$p.value, method = "BH")
post_hoc_res[["diff_across_race"]] <- post_hoc_res[["diff_across_race"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))
post_hoc_res[["diff_across_tissue"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_tissue"]])$p.value, method = "BH")
post_hoc_res[["diff_across_tissue"]] <- post_hoc_res[["diff_across_tissue"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))

# Plot the residuals of each model
pdf(file = "U01_Project/Results/Not_winsorized/LMEM_aTL_DNAQC_nowin_residual_plot.pdf", width = 15, height = 10)
cowplot::plot_grid(plotlist = p_list, ncol = 4)
dev.off()

# Output the results
out_list <- list("tl_tab" = tl_tab,
                 "coef_tab" = m_coef_tab_wide, 
                 "anova_tab" = m_anova_tab_wide,
                 "post_hoc_race" = post_hoc_res[["diff_across_race"]],
                 "post_hoc_tissue" = post_hoc_res[["diff_across_tissue"]],
                 "post_hoc_ageTrend_tissue" = post_hoc_res[["ageTrendDiff_across_tissue"]])

write.xlsx(out_list, file = "U01_Project/Results/Not_winsorized/LMEM_aTL_DNAQC_nowin_result.xlsx")

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))


# Modeling with the combined cohorts (winsorized; 5 tissues; add DNA QC metrics as covariates).# -----
dv <- c("aTL")
dv <- paste0(dv, "_win")
covar <- c("DIN", "PUnfrag", "PHigh", "PSev",
           "X260.280", "X260.230",
           "Nano", "Pico", "Tape")
covar <- paste0(covar, "_win")
dv_label <- c("TL (kb)")
m.data <- dataSW
post_hoc_res <- list()
for (i in 1:length(dv)){
  m.formula <- as.formula(paste0(dv[i], "~ Age + Sex + Tissue + Race + Age:Tissue",
                                 "+", paste0(covar, collapse = "+")))
  mod <- lme(m.formula, data=m.data, random=~1|ID, na.action=na.omit)
  p <- plot(mod, main = dv_label[i])
  m_coef_temp <- as_tibble(coef(summary(mod)), rownames = "Var_d") %>%
    mutate(Value = as.character(round(Value, 3)),
           Std.Error = round(Std.Error, 3),
           Pvalue = `p-value`) %>%
    select(-all_of(c("DF", "t-value", "p-value")))
  
  m_coef_temp <- m_coef_temp %>%
    add_row(Var_d = paste0("Sex", levels(m.data$Sex)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Sex", levels(m.data$Sex)[2]))) %>%
    add_row(Var_d = paste0("Tissue", levels(m.data$Tissue)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Tissue", levels(m.data$Tissue)[2]))+1) %>%
    add_row(Var_d = paste0("Race", levels(m.data$Race)[1]),
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Race", levels(m.data$Race)[2]))+2) %>%
    mutate(dv = dv[i])
  
  m_anova_temp <- as_tibble(anova(mod), rownames = "Var_c") %>% 
    transmute(Var_c = Var_c,
              anova_numDF = numDF,
              anova_denDF = denDF,
              anova_Fvalue = round(`F-value`, 3),
              anova_Pvalue = `p-value`,
              dv = dv[i])
  
  # post-hoc analyses
  ph_race <- as.data.frame(emmeans(mod, pairwise~Race, adjust = "none")$contrast) %>% 
    mutate(dv = dv[i])
  
  ph_tissue <- as.data.frame(emmeans(mod, pairwise~Tissue|Age, adjust = "none")$contrasts) %>% 
    mutate(dv = dv[i])
  
  ph_ageTrend_tissue <- as.data.frame(emtrends(mod, pairwise~Tissue, var = "Age")$emtrends) %>%
    mutate(dv = dv[i])
  
  if(i == 1){
    p_list <- list(p)
    m_coef_tab <- m_coef_temp
    m_anova_tab <- m_anova_temp
    post_hoc_res[["diff_across_race"]] <- ph_race
    post_hoc_res[["diff_across_tissue"]] <- ph_tissue
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- ph_ageTrend_tissue
    
  } else {
    p_list <- c(p_list, list(p))
    
    m_coef_tab <- rbind(m_coef_tab, m_coef_temp)
    
    m_anova_tab <- rbind(m_anova_tab, m_anova_temp)
    
    post_hoc_res[["diff_across_race"]] <- rbind(post_hoc_res[["diff_across_race"]],
                                                ph_race)
    post_hoc_res[["diff_across_tissue"]] <- rbind(post_hoc_res[["diff_across_tissue"]],
                                                  ph_tissue)
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- rbind(post_hoc_res[["ageTrendDiff_across_tissue"]],
                                                          ph_ageTrend_tissue)
    
    
  }
}

# Adjust the p value independently for the coef tab and the anova tab
m_coef_tab$p.adj <- p.adjust(m_coef_tab$Pvalue, method = "BH")
m_anova_tab$p.adj <- p.adjust(m_anova_tab$anova_Pvalue, method = "BH")
m_coef_tab <- m_coef_tab %>%
  mutate(Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue, 3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj)
m_coef_tab_wide <- m_coef_tab %>%
  pivot_wider(names_from = dv, values_from = c(Value, Std.Error, Pvalue)) %>%
  select(all_of(c("Var_d", paste0(rep(c("Value_", "Std.Error_", "Pvalue_"), 3),
                                  rep(dv, each = 3)))))

m_anova_tab <- m_anova_tab %>%
  mutate(anova_Pvalue = ifelse(anova_Pvalue<0.001, "<0.001", round(anova_Pvalue, 3))) %>%
  mutate(anova_Pvalue = ifelse(p.adj<0.01, paste0(anova_Pvalue, "*"), anova_Pvalue)) %>%
  select(-p.adj)

m_anova_tab_wide <- m_anova_tab %>%
  mutate(anova_df = paste(anova_numDF, anova_denDF, sep = ", ")) %>%
  select(-c(anova_numDF, anova_denDF)) %>%
  pivot_wider(names_from = dv, values_from = c(anova_df, anova_Fvalue, anova_Pvalue)) %>%
  select(all_of(c("Var_c", paste0(rep(c("anova_df_", "anova_Fvalue_", "anova_Pvalue_"), 3),
                                  rep(dv, each = 3)))))

# Adjust the p value for post-hoc analyses
post_hoc_res[["diff_across_race"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_race"]])$p.value, method = "BH")
post_hoc_res[["diff_across_race"]] <- post_hoc_res[["diff_across_race"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))
post_hoc_res[["diff_across_tissue"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_tissue"]])$p.value, method = "BH")
post_hoc_res[["diff_across_tissue"]] <- post_hoc_res[["diff_across_tissue"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))

# Output the results
out_list <- list("coef_tab" = m_coef_tab_wide, 
                 "anova_tab" = m_anova_tab_wide,
                 "post_hoc_race" = post_hoc_res[["diff_across_race"]],
                 "post_hoc_tissue" = post_hoc_res[["diff_across_tissue"]],
                 "post_hoc_ageTrend_tissue" = post_hoc_res[["ageTrendDiff_across_tissue"]])

write.xlsx(out_list, file = "U01_Project/Results/Winsorized/LMEM_aTL_DNAQCasCovar_win_result.xlsx")

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))


# Modeling with the combined cohorts (not winsorized; 5 tissues; add DNA QC metrics as covariates).# -----
dv <- c("aTL")
covar <- c("DIN", "PUnfrag", "PHigh", "PSev",
           "X260.280", "X260.230",
           "Nano", "Pico", "Tape")
dv_label <- c("TL (kb)")
m.data <- dataSW
post_hoc_res <- list()
for (i in 1:length(dv)){
  m.formula <- as.formula(paste0(dv[i], "~ Age + Sex + Tissue + Race + Age:Tissue",
                                 "+", paste0(covar, collapse = "+")))
  mod <- lme(m.formula, data=m.data, random=~1|ID, na.action=na.omit)
  p <- plot(mod, main = dv_label[i])
  m_coef_temp <- as_tibble(coef(summary(mod)), rownames = "Var_d") %>%
    mutate(Value = as.character(round(Value, 3)),
           Std.Error = round(Std.Error, 3),
           Pvalue = `p-value`) %>%
    select(-all_of(c("DF", "t-value", "p-value")))
  
  m_coef_temp <- m_coef_temp %>%
    add_row(Var_d = paste0("Sex", levels(m.data$Sex)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Sex", levels(m.data$Sex)[2]))) %>%
    add_row(Var_d = paste0("Tissue", levels(m.data$Tissue)[1]), 
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Tissue", levels(m.data$Tissue)[2]))+1) %>%
    add_row(Var_d = paste0("Race", levels(m.data$Race)[1]),
            Value = "Reference", .before = which(m_coef_temp$Var_d == paste0("Race", levels(m.data$Race)[2]))+2) %>%
    mutate(dv = dv[i])
  
  m_anova_temp <- as_tibble(anova(mod), rownames = "Var_c") %>% 
    transmute(Var_c = Var_c,
              anova_numDF = numDF,
              anova_denDF = denDF,
              anova_Fvalue = round(`F-value`, 3),
              anova_Pvalue = `p-value`,
              dv = dv[i])
  
  # post-hoc analyses
  ph_race <- as.data.frame(emmeans(mod, pairwise~Race, adjust = "none")$contrast) %>% 
    mutate(dv = dv[i])
  
  ph_tissue <- as.data.frame(emmeans(mod, pairwise~Tissue|Age, adjust = "none")$contrasts) %>% 
    mutate(dv = dv[i])
  
  ph_ageTrend_tissue <- as.data.frame(emtrends(mod, pairwise~Tissue, var = "Age")$emtrends) %>%
    mutate(dv = dv[i])
  
  if(i == 1){
    p_list <- list(p)
    m_coef_tab <- m_coef_temp
    m_anova_tab <- m_anova_temp
    post_hoc_res[["diff_across_race"]] <- ph_race
    post_hoc_res[["diff_across_tissue"]] <- ph_tissue
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- ph_ageTrend_tissue
    
  } else {
    p_list <- c(p_list, list(p))
    
    m_coef_tab <- rbind(m_coef_tab, m_coef_temp)
    
    m_anova_tab <- rbind(m_anova_tab, m_anova_temp)
    
    post_hoc_res[["diff_across_race"]] <- rbind(post_hoc_res[["diff_across_race"]],
                                                ph_race)
    post_hoc_res[["diff_across_tissue"]] <- rbind(post_hoc_res[["diff_across_tissue"]],
                                                  ph_tissue)
    post_hoc_res[["ageTrendDiff_across_tissue"]] <- rbind(post_hoc_res[["ageTrendDiff_across_tissue"]],
                                                          ph_ageTrend_tissue)
    
    
  }
}

# Adjust the p value independently for the coef tab and the anova tab
m_coef_tab$p.adj <- p.adjust(m_coef_tab$Pvalue, method = "BH")
m_anova_tab$p.adj <- p.adjust(m_anova_tab$anova_Pvalue, method = "BH")
m_coef_tab <- m_coef_tab %>%
  mutate(Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue, 3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj)
m_coef_tab_wide <- m_coef_tab %>%
  pivot_wider(names_from = dv, values_from = c(Value, Std.Error, Pvalue)) %>%
  select(all_of(c("Var_d", paste0(rep(c("Value_", "Std.Error_", "Pvalue_"), 3),
                                  rep(dv, each = 3)))))

m_anova_tab <- m_anova_tab %>%
  mutate(anova_Pvalue = ifelse(anova_Pvalue<0.001, "<0.001", round(anova_Pvalue, 3))) %>%
  mutate(anova_Pvalue = ifelse(p.adj<0.01, paste0(anova_Pvalue, "*"), anova_Pvalue)) %>%
  select(-p.adj)

m_anova_tab_wide <- m_anova_tab %>%
  mutate(anova_df = paste(anova_numDF, anova_denDF, sep = ", ")) %>%
  select(-c(anova_numDF, anova_denDF)) %>%
  pivot_wider(names_from = dv, values_from = c(anova_df, anova_Fvalue, anova_Pvalue)) %>%
  select(all_of(c("Var_c", paste0(rep(c("anova_df_", "anova_Fvalue_", "anova_Pvalue_"), 3),
                                  rep(dv, each = 3)))))

# Adjust the p value for post-hoc analyses
post_hoc_res[["diff_across_race"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_race"]])$p.value, method = "BH")
post_hoc_res[["diff_across_race"]] <- post_hoc_res[["diff_across_race"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))
post_hoc_res[["diff_across_tissue"]][,"p.adj"] <- p.adjust((post_hoc_res[["diff_across_tissue"]])$p.value, method = "BH")
post_hoc_res[["diff_across_tissue"]] <- post_hoc_res[["diff_across_tissue"]] %>%
  mutate(`sig.p (FDR<0.01)` = ifelse(p.adj < 0.01, "*", ""))

# Output the results
out_list <- list("coef_tab" = m_coef_tab_wide, 
                 "anova_tab" = m_anova_tab_wide,
                 "post_hoc_race" = post_hoc_res[["diff_across_race"]],
                 "post_hoc_tissue" = post_hoc_res[["diff_across_tissue"]],
                 "post_hoc_ageTrend_tissue" = post_hoc_res[["ageTrendDiff_across_tissue"]])

write.xlsx(out_list, file = "U01_Project/Results/Not_winsorized/LMEM_aTL_DNAQCasCovar_nowin_result.xlsx")

rm(list=setdiff(ls(), c("dataSW", "dataSW_wide")))






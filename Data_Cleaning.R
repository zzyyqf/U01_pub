# U01 Project
# Data Cleaning (Winsorizing)
# June 14, 2023

# Load packages
library(tidyverse)
library(psych)
library(summarytools)
library(DescTools)
library(openxlsx)

# Read in data
data <- read.csv("U01_Project/Data/U01 DNA and Telomere Data.csv", header=T, stringsAsFactors=FALSE, fileEncoding="latin1")
add_data_320 <- read.csv("U01_Project/Data/tl_nanodrop_data.csv")

# Drop one extra row in the "add_data_320"
add_data_320 <- add_data_320 %>%
  filter(!(Origin == 723 & Sample_ID == "TL12_saliva_Tube1.2"))

# Subset "add_data_320"
add_data_320 <- add_data_320 %>%
  select(ID, Baseline_Corr_Abs_320, Tissue.x) %>%
  mutate(Tissue = Tissue.x) %>%
  select(-Tissue.x)

# Subset "data"
dataS <- data[,c(1:8,10:12,15:19,21,24:28,30)]

# Join the the A320nm data into the dataS
dataS <- dataS %>%
  left_join(add_data_320, by = c("ID", "Tissue"))

names(dataS) <- c("Origin","ID","Cohort","Age","Sex","EthO","EthCorr","Tissue","Nano","X260.280","X260.230","DIN",
                  "Tape","PUnfrag","PHigh","PSev","Pico","aTL","CV_Telo","rep_Telo","CV_IFNB1","rep_IFNB1","Rerun", "A320")

# Winsorizing aTL variable and DNA QC metrics within each cohort/sex/tissue subset
# Outliers are defined as values outside the range of (Q1-1.5*IQR) to (Q3+1.5*IQR).
# Q1 = lower quartile; Q3 = upper quartile; IQR = inter-quartile range.
# Outliers were winsorized to have value of Q1-1.5*IQR or Q3+1.5*IQR.
cohort <- unique(dataS$Cohort)
tissue <- unique(dataS$Tissue)
win_var <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev", 
             "X260.280", "X260.230",
             "Nano", "Pico", "Tape")
win_var_label <- c("aTL", "DIN", 
                   "%Unfragmented (>3000 bp)", "%Highly Fragmented (250 to 3000 bp)", "%Severely Fragmented (<250 bp)",
                   "A260/A280", "A260/A230",
                   "Nanodrop Concentration (ng/uL)", "PicoGreen Concentration (ng/uL)", "TapeStation Concentration (ng/uL)")

dataSW <- data.frame(NULL)
for (i in (1:length(cohort))) {
  for (j in (1:length(tissue))){
    data_temp <- dataS %>%
      filter(Cohort == cohort[i] & Tissue == tissue[j])
    for (k in 1:length(win_var)) {
      var_median <- median(data_temp[, win_var[k]], na.rm = T)
      var_iqr <- IQR(data_temp[, win_var[k]], na.rm = T)
      var_q1 <- quantile(data_temp[, win_var[k]], na.rm = T, probs = 0.25)
      var_q3 <- quantile(data_temp[, win_var[k]], na.rm = T, probs = 0.75)
      data_temp <- data_temp %>%
        dplyr::mutate(!!as.name(paste0(win_var[k],"_win")) := Winsorize(data_temp[, win_var[k]], 
                                                                 minval = var_q1-1.5*var_iqr,
                                                                 maxval = var_q3+1.5*var_iqr,
                                                                 probs = NULL,
                                                                 na.rm = T))
      
    }
    dataSW <- rbind(dataSW, data_temp)
  }
}

# Look at how many samples were winsorized.
data_sum <- data.frame(NULL)
for (i in 1:length(win_var)) {
  data_sum_temp <- dataSW %>%
    filter(!!as.name(win_var[i]) != !!as.name(paste0(win_var[i],"_win"))) %>%
    group_by(Cohort, Tissue) %>%
    dplyr::count() %>%
    mutate(Variable = win_var[i])
  
  data_sum <- rbind(data_sum, data_sum_temp)
}

View(data_sum)
write_csv(data_sum, file = "U01_Project/Results/number_of_samples_winsorized.csv")

# Look at the summary stats of variables before and after winsorizing.
des_tab_stratified <- describeBy(dataSW[,c(win_var, paste0(win_var, "_win"))], group = list(dataSW$Cohort, dataSW$Tissue), mat = T)
des_tab_stratified <- des_tab_stratified %>% as_tibble(rownames = "var")
des_tab_combined <- describe(dataSW[,c(win_var, paste0(win_var, "_win"))])
des_tab_combined <- des_tab_combined %>% as_tibble(rownames = "var")
outlist <- list("stratified" = des_tab_stratified,
                "combined" = des_tab_combined)

write.xlsx(outlist, file = "U01_Project/Results/skewness_kurtosis_tab.xlsx")

# Look at the distribution of the variables before and after winsorizing.
pdf(file = "U01_Project/Results/Distribution_aTL_and_DNAQC_before_and_after_winsorization.pdf", width = 11, height = 8)
par(mfrow=c(3, 4))
for (i in (1:length(win_var))) {
  {hist(dataSW[,win_var[i]], col = 'gray', freq = T,
        main = win_var_label[i], xlab = win_var_label[i])
    hist(dataSW[,paste0(win_var[i], "_win")], col = 'blue', 
         freq = T, add = T)}
}
dev.off()

# Change the class of the Tissue variable to factor.
dataSW$Tissue <- factor(dataSW$Tissue, levels=c("Buccal","Saliva","DBS","Buffy Coat","PBMC"),
                        labels = c("Buccal","Saliva","DBS","Buffy","PBMC"))
freq(dataSW$Tissue)

# Check the empty value of sex
class(dataSW$Sex)
freq(dataSW$Sex)
dataSW[dataSW$Sex == "", "ID"]
dataSW[dataSW$ID == 201681,"Sex"] <- dataSW[dataSW$ID == 201681 & dataSW$Tissue == "Buffy","Sex"]
dataSW$Sex <- factor(dataSW$Sex, levels = c("Female", "Male"))

# Check the empty value of ethnicity
class(dataSW$EthCorr)
freq(dataSW$EthCorr)
dataSW[dataSW$EthCorr == "", "ID"]
dataSW$EthCorr <- factor(dataSW$EthCorr, levels = c("White", "Black", "Other"))

# Pivot the data from the long format to the wide format (one participant per row).
win_var <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev", 
             "X260.280", "X260.230", "Nano", "Pico", "Tape")

target_var <- c(win_var, paste0(win_var, "_win"),
                "CV_Telo", "rep_Telo", "CV_IFNB1", "rep_IFNB1",
                "Rerun")

dataSW_wide <- dataSW %>%
  select(-c(Origin, A320)) %>%
  pivot_wider(names_from = "Tissue",
              values_from = all_of(target_var))

# Save the winsorized data
save(dataSW, dataSW_wide, file = "U01_Project/Data/u01_data.RData")

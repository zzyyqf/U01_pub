# U01 Project
# Descriptive Analyses
# June 14, 2023

# Load packages
library(tidyverse)
library(psych)
library(summarytools)
library(Rmisc)
library(crosstable)
library(labelled)
library(openxlsx)

# Read in data
load("U01_Project/Data/u01_data.RData")

# Make a table summarizing demographics by cohort. # -----
var_label(dataSW_wide$EthCorr) <- "Race"
demo_tab <- crosstable(data = dataSW_wide,
                       cols = c("Sex", "Age", "EthCorr"),
                       by = c("Cohort"),
                       showNA = "ifany",
                       funs = c("Mean (SD)" = meansd, 
                                "Min-Max" = function(x) paste0(round(range(x, na.rm=T),2), collapse = "-"), 
                                "N (NA)" = nna),
                       funs_arg = list(na.rm = T),
                       test = F,
                       num_digits = 2,
                       percent_digits = 1, percent_pattern = "{n} ({p_col})")

demo_tab %>% as_flextable()

write_csv(demo_tab, file = "U01_Project/Results/Demographics_table.csv")

# Make a table summarizing aTL and DNA QC metrics (winsorized) by cohort and tissue. # ------
var_label(dataSW$aTL_win) <- "aTL (kb)"
var_label(dataSW$DIN_win) <- "DIN"
var_label(dataSW$PUnfrag_win) <- "% Unfrag DNA"
var_label(dataSW$PHigh_win) <- "% Highly Degraded DNA"
var_label(dataSW$PSev_win) <- "% Severely Degraded DNA"
var_label(dataSW$X260.280_win) <- "A260/A280"
var_label(dataSW$X260.230_win) <- "A260/A230"
var_label(dataSW$Nano_win) <- "Nanodrop Conc (ng/uL)"
var_label(dataSW$Pico_win) <- "PicoGreen Conc (ng/uL)"
var_label(dataSW$Tape_win) <- "TapeStation Conc (ng/uL)"

metrics_win_tab <- crosstable(data = dataSW,
                             cols = c("aTL_win", "DIN_win", "PUnfrag_win", 
                                      "PHigh_win", "PSev_win",
                                      "X260.280_win", "X260.230_win",
                                      "Nano_win", "Pico_win", "Tape_win"),
                             by = c("Cohort", "Tissue"),
                             showNA = "ifany",
                             funs = c("Mean (SD)" = meansd),
                             funs_arg = list(na.rm = T),
                             test = F,
                             num_digits = 2)

metrics_win_tab %>% as_flextable()

write_csv(metrics_win_tab, file = "U01_Project/Results/Winsorized/aTL_DNAQC_win_table.csv")

# Make a table summarizing aTL and DNA QC metrics (not winsorized) by cohort and tissue. # -----
var_label(dataSW$aTL) <- "aTL (kb)"
var_label(dataSW$DIN) <- "DIN"
var_label(dataSW$PUnfrag) <- "% Unfrag DNA"
var_label(dataSW$PHigh) <- "% Highly Degraded DNA"
var_label(dataSW$PSev) <- "% Severely Degraded DNA"
var_label(dataSW$X260.280) <- "A260/A280"
var_label(dataSW$X260.230) <- "A260/A230"
var_label(dataSW$Nano) <- "Nanodrop Conc (ng/uL)"
var_label(dataSW$Pico) <- "PicoGreen Conc (ng/uL)"
var_label(dataSW$Tape) <- "TapeStation Conc (ng/uL)"

metrics_nowin_tab <- crosstable(data = dataSW,
                               cols = c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
                                        "X260.280", "X260.230",
                                        "Nano", "Pico", "Tape"),
                               by = c("Tissue", "Cohort"),
                               showNA = "ifany",
                               funs = c("Mean (SD)" = meansd),
                               funs_arg = list(na.rm = T),
                               test = F,
                               num_digits = 2)

metrics_nowin_tab %>% as_flextable()

write_csv(metrics_nowin_tab, file = "U01_Project/Results/Not_winsorized/aTL_DNAQC_nowin_table.csv")





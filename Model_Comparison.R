# U01 Project
# Incorpration of DNA metrics into models of predicting aTL
# June 17, 2023

# Load packages
library(tidyverse)
library(psych)
library(MuMIn)
library(lme4)
library(lmerTest)

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

# Make function that runs glm model comparison and spits out AIC values and also beta estimates and p-values for age
hfit <- function(variables,set){
  
  ## set vars
  vars=variables
  
  ## make matrix to hold results
  mvars=c("mod","k","AICc","delta","wi","R2m","R2c")
  mat=matrix(NA,ncol=length(mvars),nrow=length(vars))
  
  ## make matrix data frame and give names
  mat=data.frame(mat)
  names(mat)=mvars
  
  ## loop through
  for(i in 1:length(vars)){
    
    ## code model
    mod=glm(formula(vars[i]),data=set,family="gaussian")
    
    ## AICc
    mat[i,"AICc"]=round(AICc(mod),2)
    
    ## model formulas
    mat[i,"mod"]=vars[i]
    
    ## k
    mat[i,"k"]=length(coef(mod))
    
    ## r2 for mixed models
    #R2=r.squared(mod)
    mat[i,"R2m"]=round(r.squaredGLMM(mod)[1,1],3)
    mat[i,"R2c"]=round(r.squaredGLMM(mod)[1,2],3)
  }
  
  ## calculate delta AICc and weights
  rank=mat[with(mat,order(AICc)),]
  rank$delta=rank$AICc-rank$AICc[1]
  rank$wi=Weights(rank$AICc)
  rank$wi=round(rank$wi,3)
  
  ## clean 
  rownames(rank)=NULL
  rank
}

# Make function that runs lmer model comparison and spits out AIC values and also beta estimates and p-values for age
hfit.lmer <- function(variables,set){
  
  ## set vars
  vars=variables
  
  ## make matrix to hold results
  mvars=c("mod","k","AICc","delta","wi","R2m","R2c")
  mat=matrix(NA,ncol=length(mvars),nrow=length(vars))
  
  ## make matrix data frame and give names
  mat=data.frame(mat)
  names(mat)=mvars
  
  ## loop through
  for(i in 1:length(vars)){
    
    ## code model
    mod=lmer(formula(vars[i]),data=set)
    
    ## AICc
    mat[i,"AICc"]=round(AICc(mod),2)
    
    ## model formulas
    mat[i,"mod"]=vars[i]
    
    ## k
    mat[i,"k"]=length(fixef(mod))+2
    
    ## r2 for mixed models
    #R2=r.squared(mod)
    mat[i,"R2m"]=round(r.squaredGLMM(mod)[1,1],3)
    mat[i,"R2c"]=round(r.squaredGLMM(mod)[1,2],3)
  }
  
  ## calculate delta AICc and weights
  rank=mat[with(mat,order(AICc)),]
  rank$delta=rank$AICc-rank$AICc[1]
  rank$wi=Weights(rank$AICc)
  rank$wi=round(rank$wi,3)
  
  ## clean 
  rownames(rank)=NULL
  rank
}

## Winsorized variables (separate models for each tissue) ##----
# Calculate the correlation among DNA QC metrics (by tissue).
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")

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
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer", "cor_res",
                          "dna_qc_var")))

# Configure the subsetting expression term for dredge
tissue <- levels(dataSW$Tissue)
tape_var <- c("DIN", "PUnfrag", "PHigh", "PSev", "Tape")
tape_var <- paste0(tape_var, "_win")
s_expr <- list()
for (i in 1:length(tissue)){
  s_expr[[tissue[i]]] <- cor_res %>% 
    filter(Tissue == tissue[i])
  if (tissue[i] == "Buffy") {
    s_expr[[tissue[i]]] <- s_expr[[tissue[i]]] %>%
      filter((!(Parameter1 %in% tape_var))&(!(Parameter2 %in% tape_var)))
  } else{
    s_expr[[tissue[i]]] <- s_expr[[tissue[i]]]
  }
  s_expr[[tissue[i]]] <- s_expr[[tissue[i]]] %>%
    mutate(p1_p2 = ifelse(abs(rho) >=0.4, paste0("(", Parameter1, " && ", Parameter2, ")"), NA)) %>%
    select(p1_p2) %>%
    na.omit() %>%
    unlist() %>%
    paste0(collapse = " || ") %>%
    paste0("!", "(", ., ")")  
}

# Run model comparison across all the tissues and extract the top ranked models and averaged coefficients
dv <- "aTL_win"
fit_list <- list()
Fsum_list <- list()
Fcoef_list <- list()
for (i in 1:length(tissue)){
  fixed_iv <- c("Age", "Sex", "Race")
  m_data <- dataSW %>%
    filter(Tissue == tissue[i]) %>%
    select(all_of(c("ID", fixed_iv, dv, 
                    if (tissue[i] == "Buffy") {setdiff(dna_qc_var, tape_var)} else {dna_qc_var}))) %>%
    na.omit()
  m_formula <- as.formula(paste0(dv, "~", paste0(c(fixed_iv, if (tissue[i] == "Buffy") {setdiff(dna_qc_var, tape_var)} else {dna_qc_var}), collapse = "+")))
  lmglobal <- glm(formula=m_formula, data=m_data, na.action=na.fail, family="gaussian")
  mset <- dredge(lmglobal, 
                 m.lim=c(0,round(nrow(m_data)/10,0)), # limited the number of parameters per candidate model to approximately 1 per 10 observations
                 #fixed=fixed_iv,
                 subset=rlang::parse_expr(s_expr[[tissue[i]]]))
  # get model formula
  models <- sapply(get.models(mset,subset=T),function(x) formula(x))
  models <- paste(models)
  rmodels <- models
  null_mod <- paste0(dv, " ~ 1")
  age_mod <- paste0(dv, " ~ Age + 1")
  if (null_mod %in% rmodels){
    rmodels
  } else{
    rmodels <- c(rmodels, null_mod)
  }
  
  if (age_mod %in% rmodels){
    rmodels
  } else{
    rmodels <- c(rmodels, age_mod)
  }
  
  # get the model fit criteria
  fit <- hfit(variables=rmodels,set=m_data)
  fit <- fit[fit$delta<=2 | fit$mod == null_mod | fit$mod == age_mod,]
  
  # replace the model variable names with more reader-friendly names
  mod_var <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
               "X260.280", "X260.230", 
               "Nano", "Pico", "Tape")
  mod_var <- paste0(mod_var, "_win")
  mod_var_label <-  c("aTL", "DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")
  
  names(mod_var_label) <- mod_var
  
  fit_mod_origin <- fit$mod
  fit$mod <- str_replace_all(fit$mod, pattern = mod_var_label)
  
  # get averaged model coefficients across top ranked models
  if (length(fit_mod_origin) == 3) {
    avgs <- glm(formula = fit_mod_origin[1], data = m_data)
    Fsum <- summary(avgs)
    Fcoef <- as_tibble(coef(Fsum), rownames = "Predictor") %>%
      transmute(Predictor = Predictor,
                Estimate = Estimate,
                SE = `Std. Error`,
                Pvalue = `Pr(>|t|)`,
                Tissue = tissue[i])
  } else{
    avgs <- model.avg(mset,subset=delta<2, family=gaussian, revised.var=T, beta="none")
    Fsum <- summary(avgs)
    Fcoef <- as_tibble(Fsum$coefmat.subset, rownames = "Predictor") %>%
      transmute(Predictor = Predictor,
                Estimate = Estimate,
                SE = `Adjusted SE`,
                Pvalue = `Pr(>|z|)`,
                Tissue = tissue[i])
  }
  
  # put wanted results into a list
  fit_list[[tissue[i]]] <- fit
  Fsum_list[[tissue[i]]] <- Fsum
  Fcoef_list[[tissue[i]]] <- Fcoef
}

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "cor_res", "hfit", "hfit.lmer",
                          "fit_list", "Fsum_list", "Fcoef_list")))

# Format the model coefficient table
coef_tab <- rbind(Fcoef_list[[1]], Fcoef_list[[2]], Fcoef_list[[3]], 
                  Fcoef_list[[4]], Fcoef_list[[5]])
coef_tab$p.adj <- p.adjust(coef_tab$Pvalue, method = "BH")
coef_tab <- coef_tab %>%
  mutate(Estimate = round(Estimate,3),
         SE = round(SE,3),
         Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue,3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj) %>%
  pivot_wider(names_from = "Tissue", values_from = c("Estimate", "SE", "Pvalue")) %>%
  select(all_of(c("Predictor", 
                  paste0(rep(c("Estimate_", "SE_", "Pvalue_"), 3), 
                         rep(c("Buccal", "Saliva", "DBS", "Buffy", "PBMC"), each = 3)))))

dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")
iv_fixed <- c("Age", "Sex", "Race")
var_raw <- sapply(iv_fixed, function(x){paste0(x,levels(dataSW[[x]]))}) %>% unlist()
var_trans <- data.frame(Predictor = c("(Intercept)", var_raw, dna_qc_var),
                        Predictor_d = c("(Intercept)", var_raw, dna_qc_var_label))

coef_tab <- var_trans %>%
  full_join(coef_tab, by = "Predictor") %>% select(-Predictor)

# Export the results in one Excel file.
out_list <- list("coef_all_tissue" = coef_tab,
                 "Buccal_top_models" = fit_list[["Buccal"]], 
                 "Saliva_top_models" = fit_list[["Saliva"]],
                 "DBS_top_models" = fit_list[["DBS"]],
                 "Buffy_top_models" = fit_list[["Buffy"]],
                 "PBMC_top_models" = fit_list[["PBMC"]])

write.xlsx(out_list, file = "U01_Project/Results/Winsorized/Model_comparison_results_win.xlsx")
save(fit_list, Fsum_list, Fcoef_list, file = "U01_Project/Results/Winsorized/Model_comparison_unformattedResults_win.RData")

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer")))

## Not winsorized variables (separate models for each tissue) ##----
# Calculate the correlation among DNA QC metrics (by tissue).
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")

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
cor_res$Tissue <- factor(cor_res$Tissue, levels = levels(dataSW$Tissue))

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer", "cor_res",
                          "dna_qc_var")))

# Configure the subsetting expression term for dredge
tissue <- levels(dataSW$Tissue)
tape_var <- c("DIN", "PUnfrag", "PHigh", "PSev", "Tape")
s_expr <- list()
for (i in 1:length(tissue)){
  s_expr[[tissue[i]]] <- cor_res %>% 
    filter(Tissue == tissue[i])
  if (tissue[i] == "Buffy") {
    s_expr[[tissue[i]]] <- s_expr[[tissue[i]]] %>%
      filter((!(Parameter1 %in% tape_var))&(!(Parameter2 %in% tape_var)))
  } else{
    s_expr[[tissue[i]]] <- s_expr[[tissue[i]]]
  }
  s_expr[[tissue[i]]] <- s_expr[[tissue[i]]] %>%
    mutate(p1_p2 = ifelse(abs(rho) >=0.4, paste0("(", Parameter1, " && ", Parameter2, ")"), NA)) %>%
    select(p1_p2) %>%
    na.omit() %>%
    unlist() %>%
    paste0(collapse = " || ") %>%
    paste0("!", "(", ., ")")  
}

# Run model comparison across all the tissues and extract the top ranked models and averaged coefficients
dv <- "aTL"
fit_list <- list()
Fsum_list <- list()
Fcoef_list <- list()
for (i in 1:length(tissue)){
  fixed_iv <- c("Age", "Sex", "Race")
  m_data <- dataSW %>%
    filter(Tissue == tissue[i]) %>%
    select(all_of(c("ID", fixed_iv, dv, 
                    if (tissue[i] == "Buffy") {setdiff(dna_qc_var, tape_var)} else {dna_qc_var}))) %>%
    na.omit()
  m_formula <- as.formula(paste0(dv, "~", paste0(c(fixed_iv, if (tissue[i] == "Buffy") {setdiff(dna_qc_var, tape_var)} else {dna_qc_var}), collapse = "+")))
  lmglobal <- glm(formula=m_formula, data=m_data, na.action=na.fail, family="gaussian")
  mset <- dredge(lmglobal, 
                 m.lim=c(0,round(nrow(m_data)/10,0)), # limited the number of parameters per candidate model to approximately 1 per 10 observations
                 #fixed=fixed_iv,
                 subset=rlang::parse_expr(s_expr[[tissue[i]]]))
  # get model formula
  models <- sapply(get.models(mset,subset=T),function(x) formula(x))
  models <- paste(models)
  rmodels <- models
  null_mod <- paste0(dv, " ~ 1")
  age_mod <- paste0(dv, " ~ Age + 1")
  if (null_mod %in% rmodels){
    rmodels
  } else{
    rmodels <- c(rmodels, null_mod)
  }
  
  if (age_mod %in% rmodels){
    rmodels
  } else{
    rmodels <- c(rmodels, age_mod)
  }
  
  # get the model fit criteria
  fit <- hfit(variables=rmodels,set=m_data)
  fit <- fit[fit$delta<=2 | fit$mod == null_mod | fit$mod == age_mod,]
  
  # replace the model variable names with more reader-friendly names
  mod_var <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
               "X260.280", "X260.230", 
               "Nano", "Pico", "Tape")
  mod_var_label <-  c("aTL", "DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")
  
  names(mod_var_label) <- mod_var
  
  fit_mod_origin <- fit$mod
  fit$mod <- str_replace_all(fit$mod, pattern = mod_var_label)
  
  # get averaged model coefficients across top ranked models
  if (length(fit_mod_origin) == 3) {
    avgs <- glm(formula = fit_mod_origin[1], data = m_data)
    Fsum <- summary(avgs)
    Fcoef <- as_tibble(coef(Fsum), rownames = "Predictor") %>%
      transmute(Predictor = Predictor,
                Estimate = Estimate,
                SE = `Std. Error`,
                Pvalue = `Pr(>|t|)`,
                Tissue = tissue[i])
  } else{
    avgs <- model.avg(mset,subset=delta<2, family=gaussian, revised.var=T, beta="none")
    Fsum <- summary(avgs)
    Fcoef <- as_tibble(Fsum$coefmat.subset, rownames = "Predictor") %>%
      transmute(Predictor = Predictor,
                Estimate = Estimate,
                SE = `Adjusted SE`,
                Pvalue = `Pr(>|z|)`,
                Tissue = tissue[i])
  }
  
  # put wanted results into a list
  fit_list[[tissue[i]]] <- fit
  Fsum_list[[tissue[i]]] <- Fsum
  Fcoef_list[[tissue[i]]] <- Fcoef
}

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "cor_res", "hfit", "hfit.lmer",
                          "fit_list", "Fsum_list", "Fcoef_list")))

# Format the model coefficient table
coef_tab <- rbind(Fcoef_list[[1]], Fcoef_list[[2]], Fcoef_list[[3]], 
                  Fcoef_list[[4]], Fcoef_list[[5]])
coef_tab$p.adj <- p.adjust(coef_tab$Pvalue, method = "BH")
coef_tab <- coef_tab %>%
  mutate(Estimate = round(Estimate,3),
         SE = round(SE,3),
         Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue,3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj) %>%
  pivot_wider(names_from = "Tissue", values_from = c("Estimate", "SE", "Pvalue")) %>%
  select(all_of(c("Predictor", 
                  paste0(rep(c("Estimate_", "SE_", "Pvalue_"), 3), 
                         rep(c("Buccal", "Saliva", "DBS", "Buffy", "PBMC"), each = 3)))))

dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")
iv_fixed <- c("Age", "Sex", "Race")
var_raw <- sapply(iv_fixed, function(x){paste0(x,levels(dataSW[[x]]))}) %>% unlist()
var_trans <- data.frame(Predictor = c("(Intercept)", var_raw, dna_qc_var),
                        Predictor_d = c("(Intercept)", var_raw, dna_qc_var_label))

coef_tab <- var_trans %>%
  full_join(coef_tab, by = "Predictor") %>% select(-Predictor)

# Export the results in one Excel file.
out_list <- list("coef_all_tissue" = coef_tab,
                 "Buccal_top_models" = fit_list[["Buccal"]], 
                 "Saliva_top_models" = fit_list[["Saliva"]],
                 "DBS_top_models" = fit_list[["DBS"]],
                 "Buffy_top_models" = fit_list[["Buffy"]],
                 "PBMC_top_models" = fit_list[["PBMC"]])

write.xlsx(out_list, file = "U01_Project/Results/Not_winsorized/Model_comparison_results_nowin.xlsx")
save(fit_list, Fsum_list, Fcoef_list, file = "U01_Project/Results/Not_winsorized/Model_comparison_unformattedResults_nowin.RData")

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer")))


## Winsorized variables (combined models for all the tissues) ##----
# Calculate the correlation among DNA QC metrics (by tissue).
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")

cor_res <- data.frame(NULL)
data_temp <- dataSW
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
                       as.data.frame(cor_temp))
    }
  }
}

cor_res

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer", "cor_res",
                          "dna_qc_var")))

# Configure the subsetting expression term for dredge
s_expr <- cor_res %>%
  mutate(p1_p2 = ifelse(abs(rho) >=0.4, paste0("(", Parameter1, " && ", Parameter2, ")"), NA)) %>%
  select(p1_p2) %>%
  na.omit() %>%
  unlist() %>%
  paste0(collapse = " || ") %>%
  paste0("!", "(", ., ")")  


# Run model comparison across all the tissues and extract the top ranked models and averaged coefficients
dv <- "aTL_win"
fixed_iv <- c("Age", "Sex", "Race", "Tissue")
m_data <- dataSW %>%
  select(all_of(c("ID", fixed_iv, dv, dna_qc_var))) %>%
  na.omit()
m_formula <- as.formula(paste0(dv, "~", "(1|ID)+", paste0(c(fixed_iv, dna_qc_var), collapse = "+")))
lmglobal <- lmer(formula=m_formula, data=m_data, na.action=na.fail)
mset <- dredge(lmglobal, 
               m.lim=c(0,round(nrow(m_data)/10,0)), # limited the number of parameters per candidate model to approximately 1 per 10 observations
               #fixed=fixed_iv,
               subset=rlang::parse_expr(s_expr))
# get model formula
models <- sapply(get.models(mset,subset=T),function(x) formula(x))
models <- paste(models)
rmodels <- models
null_mod <- paste0(dv, " ~ (1 | ID)")
age_mod <- paste0(dv, " ~ Age + (1 | ID)")
if (null_mod %in% rmodels){
  rmodels <- rmodels
} else{
  rmodels <- c(rmodels, null_mod)
}

if (age_mod %in% rmodels){
  rmodels <- rmodels
} else{
  rmodels <- c(rmodels, age_mod)
}

# get the model fit criteria
fit <- hfit.lmer(variables=rmodels,set=m_data)
fit <- fit[fit$delta<=2 | fit$mod == null_mod | fit$mod == age_mod,]

# replace the model variable names with more reader-friendly names
mod_var <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
             "X260.280", "X260.230", 
             "Nano", "Pico", "Tape")
mod_var <- paste0(mod_var, "_win")
mod_var_label <-  c("aTL", "DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                    "A260/A280", "A260/A230",
                    "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")

names(mod_var_label) <- mod_var

fit_mod_origin <- fit$mod
fit$mod <- str_replace_all(fit$mod, pattern = mod_var_label)

# get averaged model coefficients across top ranked models
if (length(fit_mod_origin) == 3) {
  avgs <- lmer(formula = fit_mod_origin[1], data = m_data)
  Fsum <- summary(avgs)
  Fcoef <- as_tibble(coef(Fsum), rownames = "Predictor") %>%
    transmute(Predictor = Predictor,
              Estimate = Estimate,
              SE = `Std. Error`,
              Pvalue = `Pr(>|t|)`)
} else{
  avgs <- model.avg(mset,subset=delta<2, revised.var=T, beta="none")
  Fsum <- summary(avgs)
  Fcoef <- as_tibble(Fsum$coefmat.subset, rownames = "Predictor") %>%
    transmute(Predictor = Predictor,
              Estimate = Estimate,
              SE = `Adjusted SE`,
              Pvalue = `Pr(>|z|)`)
}

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "cor_res", "hfit", "hfit.lmer",
                          "fit", "Fsum", "Fcoef")))

# Format the model coefficient table
coef_tab <- Fcoef
coef_tab$p.adj <- p.adjust(coef_tab$Pvalue, method = "BH")
coef_tab <- coef_tab %>%
  mutate(Estimate = round(Estimate,3),
         SE = round(SE,3),
         Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue,3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj)

dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var <- paste0(dna_qc_var, "_win")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")
iv_fixed <- c("Age", "Sex", "Race", "Tissue")
var_raw <- sapply(iv_fixed, function(x){paste0(x,levels(dataSW[[x]]))}) %>% unlist()
var_trans <- data.frame(Predictor = c("(Intercept)", var_raw, dna_qc_var),
                        Predictor_d = c("(Intercept)", var_raw, dna_qc_var_label))

coef_tab <- var_trans %>%
  full_join(coef_tab, by = "Predictor") %>% select(-Predictor)

# Export the results in one Excel file.
out_list <- list("coef_average" = coef_tab,
                 "top_models" = fit)

write.xlsx(out_list, file = "U01_Project/Results/Winsorized/Model_comparison_results_combTissue_win.xlsx")
save(fit, Fsum, Fcoef, file = "U01_Project/Results/Winsorized/Model_comparison_unformattedResults_combTissue_win.RData")

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer")))



## Winsorized variables (combined models for all the tissues) ##----
# Calculate the correlation among DNA QC metrics (by tissue).
tissue <- levels(dataSW$Tissue)
dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")

cor_res <- data.frame(NULL)
data_temp <- dataSW
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
                       as.data.frame(cor_temp))
    }
  }
}

cor_res

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer", "cor_res",
                          "dna_qc_var")))

# Configure the subsetting expression term for dredge
s_expr <- cor_res %>%
  mutate(p1_p2 = ifelse(abs(rho) >=0.4, paste0("(", Parameter1, " && ", Parameter2, ")"), NA)) %>%
  select(p1_p2) %>%
  na.omit() %>%
  unlist() %>%
  paste0(collapse = " || ") %>%
  paste0("!", "(", ., ")")  


# Run model comparison across all the tissues and extract the top ranked models and averaged coefficients
dv <- "aTL_win"
fixed_iv <- c("Age", "Sex", "Race", "Tissue")
m_data <- dataSW %>%
  select(all_of(c("ID", fixed_iv, dv, dna_qc_var))) %>%
  na.omit()
m_formula <- as.formula(paste0(dv, "~", "(1|ID)+", paste0(c(fixed_iv, dna_qc_var), collapse = "+")))
lmglobal <- lmer(formula=m_formula, data=m_data, na.action=na.fail)
mset <- dredge(lmglobal, 
               m.lim=c(0,round(nrow(m_data)/10,0)), # limited the number of parameters per candidate model to approximately 1 per 10 observations
               #fixed=fixed_iv,
               subset=rlang::parse_expr(s_expr))
# get model formula
models <- sapply(get.models(mset,subset=T),function(x) formula(x))
models <- paste(models)
rmodels <- models
null_mod <- paste0(dv, " ~ (1 | ID)")
age_mod <- paste0(dv, " ~ Age + (1 | ID)")
if (null_mod %in% rmodels){
  rmodels <- rmodels
} else{
  rmodels <- c(rmodels, null_mod)
}

if (age_mod %in% rmodels){
  rmodels <- rmodels
} else{
  rmodels <- c(rmodels, age_mod)
}

# get the model fit criteria
fit <- hfit.lmer(variables=rmodels,set=m_data)
fit <- fit[fit$delta<=2 | fit$mod == null_mod | fit$mod == age_mod,]

# replace the model variable names with more reader-friendly names
mod_var <- c("aTL", "DIN", "PUnfrag", "PHigh", "PSev",
             "X260.280", "X260.230", 
             "Nano", "Pico", "Tape")
mod_var_label <-  c("aTL", "DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                    "A260/A280", "A260/A230",
                    "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")

names(mod_var_label) <- mod_var

fit_mod_origin <- fit$mod
fit$mod <- str_replace_all(fit$mod, pattern = mod_var_label)

# get averaged model coefficients across top ranked models
if (length(fit_mod_origin) == 3) {
  avgs <- lmer(formula = fit_mod_origin[1], data = m_data)
  Fsum <- summary(avgs)
  Fcoef <- as_tibble(coef(Fsum), rownames = "Predictor") %>%
    transmute(Predictor = Predictor,
              Estimate = Estimate,
              SE = `Std. Error`,
              Pvalue = `Pr(>|t|)`)
} else{
  avgs <- model.avg(mset,subset=delta<2, revised.var=T, beta="none")
  Fsum <- summary(avgs)
  Fcoef <- as_tibble(Fsum$coefmat.subset, rownames = "Predictor") %>%
    transmute(Predictor = Predictor,
              Estimate = Estimate,
              SE = `Adjusted SE`,
              Pvalue = `Pr(>|z|)`)
}

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "cor_res", "hfit", "hfit.lmer",
                          "fit", "Fsum", "Fcoef")))

# Format the model coefficient table
coef_tab <- Fcoef
coef_tab$p.adj <- p.adjust(coef_tab$Pvalue, method = "BH")
coef_tab <- coef_tab %>%
  mutate(Estimate = round(Estimate,3),
         SE = round(SE,3),
         Pvalue = ifelse(Pvalue<0.001, "<0.001", round(Pvalue,3))) %>%
  mutate(Pvalue = ifelse(p.adj<0.01, paste0(Pvalue, "*"), Pvalue)) %>%
  select(-p.adj)

dna_qc_var <- c("DIN", "PUnfrag", "PHigh", "PSev",
                "X260.280", "X260.230",
                "Nano", "Pico", "Tape")
dna_qc_var_label <- c("DIN", "%Unfrag", "%High Frag", "%Severe Frag",
                      "A260/A280", "A260/A230",
                      "Conc(Nano)", "Conc(Pico)", "Conc(Tape)")
iv_fixed <- c("Age", "Sex", "Race", "Tissue")
var_raw <- sapply(iv_fixed, function(x){paste0(x,levels(dataSW[[x]]))}) %>% unlist()
var_trans <- data.frame(Predictor = c("(Intercept)", var_raw, dna_qc_var),
                        Predictor_d = c("(Intercept)", var_raw, dna_qc_var_label))

coef_tab <- var_trans %>%
  full_join(coef_tab, by = "Predictor") %>% select(-Predictor)

# Export the results in one Excel file.
out_list <- list("coef_average" = coef_tab,
                 "top_models" = fit)

write.xlsx(out_list, file = "U01_Project/Results/Not_winsorized/Model_comparison_results_combTissue_nowin.xlsx")
save(fit, Fsum, Fcoef, file = "U01_Project/Results/Not_winsorized/Model_comparison_unformattedResults_combTissue_nowin.RData")

rm(list = setdiff(ls(), c("dataSW", "dataSW_wide", "hfit", "hfit.lmer")))



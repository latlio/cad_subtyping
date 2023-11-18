# Prepare UKB phenotype for PRSet
# What this does: linearizes binary phenotype after adjusting for covariates and creates training/test data
# Author: Lathan Liou
# Contact: lathan.liou@icahn.mssm.edu
# Inputs: CSV generated from SQL extraction script
# Outputs:

library(tidyverse)
library(data.table) #for fast reading in of data
library(here)
library(caret)
library(janitor)
# library(mice)
# library(doRNG)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

standardize <- function(x){
  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

normalize_x <- function(x) {
  (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

convert_idk_to_na <- function(x) {
  case_when(x %in% c(-1,-3) ~ NA_integer_,
            # is.factor(x) ~ NA_character_,
            # is.character(x) ~ NA_character_,
            TRUE ~ x)
}

convert_neg_10 <- function(x) {
  case_when(x == -10 ~ 0.5,
            TRUE ~ x)
}

pass_through <- function(data, fun) {fun(data); data}

impute_data <- function(df, iter = 5, seed = 47) {
  # ref: https://rpubs.com/kaz_yos/mice-exclude
  all_vars <- names(df)
  (miss_vars <- names(df)[colSums(is.na(df)) > 0])
  
  imputed_matrix <- imputer_matrix <- predictor_matrix <- matrix(0, ncol = length(all_vars), nrow = length(all_vars))
  rownames(predictor_matrix) <- all_vars
  colnames(predictor_matrix) <- all_vars
  rownames(imputer_matrix) <- all_vars
  colnames(imputer_matrix) <- all_vars
  rownames(imputed_matrix) <- all_vars
  colnames(imputed_matrix) <- all_vars
  
  imputer_vars <- c("age", "sex", "race", "ses", "centre", "cad_bin", "diabetes",
                    paste0("pc", 1:15), "ses_std", "bmi_std", "sys_bp_std", "smoke",
                    "pa", "diet", "sleep", "egfr", "trigly_std", "ldl_std", "hdl_std",
                    "chol_std", "lpa_std", "glucose_std", "crp_std", "apob_std")
  imputer_matrix[,imputer_vars] <- 1
  imputed_matrix[miss_vars,] <- 1
  predictor_matrix <- imputer_matrix * imputed_matrix
  diag(predictor_matrix) <- 0
  
  dry_run <- mice(data = df, m = 1, predictorMatrix = predictor_matrix, maxit = 0)
  dry_run$method <- str_replace_all(dry_run$method, "pmm", "rf")
  dry_run$method <- str_replace_all(dry_run$method, "polyreg", "rf")
  names(dry_run$method) <- all_vars
  
  set.seed(seed)
  
  pheno_mice <- foreach(i = seq_len(iter), .combine = ibind) %dorng% {
    cat("### Started iteration", i, "\n")
    miceout <- mice(data = df, m = 1, print = TRUE,
                    predictorMatrix = predictor_matrix, method = dry_run$method,
                    MaxNWts = 2000)
    cat("### Completed iteration", i, "\n")
    ## Make sure to return the output
    miceout
  }
  pheno_mice
}

prepare_phenotype <- function(input_path,
                              output_path,
                              output_root_file_name,
                              sql_csv,
                              qc_file,
                              dropout_file,
                              covariate_file,
                              outcome_var,
                              ses = FALSE,
                              impute = FALSE,
                              keep_base_only = TRUE,
                              european_only = TRUE) {
  browser()
  
  #outcome var: 
  # Setup ----
  if(file.exists(output_path)) {
    out_base_path <- paste0(output_path, output_root_file_name)
  } else {
    dir.create(file.path(output_path))
    out_base_path <- paste0(output_path, output_root_file_name)
  }
  
  qc <- fread(here(input_path, qc_file))
  dropout <- fread(here(input_path, dropout_file), header = F)
  cov <- fread(here(input_path, covariate_file))
  
  # Initial cleaning ----
  pheno_df <- fread(here(input_path, sql_csv)) %>%
    filter(!IID %in% dropout$V1,
           IID %in% qc$V2) %>%
    left_join(cov, by = c("FID", "IID")) %>%
    janitor::clean_names() %>%
    # select(-fid) %>%
    #sometimes prset has a weird error "Malformed pheno file"
    mutate(
      centre = as.factor(centre),
      race = as.factor(case_when(
        race %in% c(1, 1001, 1002, 1003) ~ 0,
        race %in% c(3001, 3002, 3003) ~ 1,
        race %in% c(3, 3004, 5) ~ 2,
        race %in% c(4, 4001, 4002, 4003) ~ 3,
        race %in% c(2, 2001, 2002, 2003, 2004) ~ 4,
        TRUE ~ 5
      ))
    ) %>%
    select(-c(dias_bp_auto, dx_date, dad_fam_hx, mom_fam_hx, sib_fam_hx, statin_follow)) %>%
    mutate(across(diabetes:bp_med2, convert_idk_to_na)) %>%
    mutate(across(diabetes:bp_med2, convert_neg_10)) %>%
    mutate(bp_med = coalesce(bp_med1, bp_med2),
           bp_med_bin = ifelse(bp_med == 2, 1, 0), #antihypertensives
           # premature_cad_bin = case_when(
           #0 = female, 1 = male
           # age is age at diagnosis
           #   sex == 1 & age < 40 & cad_bin == 1 ~ 1,
           #   sex == 1 & age >= 40 & cad_bin == 1 ~ 0,
           #   sex == 0 & age < 50 & cad_bin == 1 ~ 1,
           #   sex == 0 & age >= 50 & cad_bin == 1 ~ 0,
           #   TRUE ~ 0
           # ),
           ses_std = standardize(ses),
           bmi_std = standardize(bmi),
           sys_bp_std = standardize(sys_bp_auto),
           smoke = case_when(
             current_smoke == 1 | past_smoke == 1 ~ 1,
             TRUE ~ 0
           ),
           duration_moderate_pa = case_when(
             n_days_moderate_pa == 0 ~ 0,
             TRUE ~ duration_moderate_pa
           ),
           duration_vigorous_pa = case_when(
             n_days_vigorous_pa == 0 ~ 0,
             TRUE ~ duration_vigorous_pa
           ),
           pa = case_when(
             n_days_moderate_pa > 3 | 
               duration_moderate_pa > 150 | 
               n_days_vigorous_pa > 3 | 
               duration_vigorous_pa > 150 ~ 1,
             TRUE ~ 0
           ),
           fruit = normalize_x(fresh_fruit + dried_fruit/5),
           veg = normalize_x(cooked_veg/3 + raw_veg/3),
           wholegrain = case_when(
             bread_type == 3 & !cereal_type %in% c(1,3,4) ~ normalize_x(bread/7),
             bread_type != 3 & cereal_type %in% c(1,3,4) ~ normalize_x(cereal/7),
             bread_type == 3 & cereal_type %in% c(1,3,4) ~ normalize_x(bread/7 + cereal/7),
             TRUE ~ 0),
           refined_grain = 1-wholegrain, #not sure about this
           fish = normalize_x(oily_fish + non_oily_fish),
           dairy = normalize_x(cheese),
           milk_cat = ifelse(milk %in% 1:5, 1, 0),
           spread_cat = ifelse(spread %in% 1:3, 1, 0),
           vegetable_oils = normalize_x(milk_cat + spread_cat),
           processed_meat_norm = 1 - normalize_x(processed_meat),
           unprocessed_meat = 1 - normalize_x(poultry + beef + lamb + pork),
           sugar_cat = case_when(
             added_sugar == 4 ~ 1, #healthy
             TRUE ~ 0
           )) %>% # have to break up mutate here so that rowSums works
    mutate(diet = rowSums(select(., "fruit",
                                 "veg", 
                                 "wholegrain",
                                 "refined_grain",
                                 "fish",
                                 "dairy",
                                 "vegetable_oils",
                                 "processed_meat_norm",
                                 "unprocessed_meat",
                                 "sugar_cat"), na.rm = TRUE),
           sleep_duration_cat = case_when(
             sleep_duration >= 7 ~ 1, #healthy
             TRUE ~ 0
           ),
           insomnia_cat = case_when(
             insomnia == 1 ~ 1, #healthy
             TRUE ~ 0
           ),
           snoring_cat = case_when(
             snoring == 0 ~ 1, #healthy
             TRUE ~ 0
           ),
           narcolepsy_cat = case_when(
             narcolepsy == 0 ~ 1, #healthy
             TRUE ~ 0
           ),
           sleep = sleep_duration_cat + insomnia_cat + snoring_cat + narcolepsy_cat,
           creatinine_std = standardize(creatinine),
           egfr = case_when(
             sex == 0 ~ 142 * pmin(creatinine_std/0.9, 1)^-0.302 *
               pmax(creatinine_std/0.9, 1)^-1.200 *
               0.9938^age,
             sex == 1 ~ 142 * pmin(creatinine_std/0.7, 1)^-0.241 *
               pmax(creatinine_std/0.7, 1)^-1.200 *
               0.9938^age*1.012),
           trigly_std = standardize(trigly),
           ldl_chol = ldl_chol*18.0182,
           ldl_bin = ifelse(ldl_chol >= 70, 1, 0),
           ldl_std = standardize(ldl_chol),
           hdl_std = standardize(hdl_chol),
           chol_std = standardize(chol),
           lpa_bin = ifelse(lpa >= 150, 1, 0),
           lpa_std = standardize(lpa),
           glucose_std = standardize(glucose),
           crp_std = standardize(crp),
           apoa_std = standardize(apoa),
           apob_std = standardize(apob)) %>%
    select(-c(bp_med, bp_med1, bp_med2, current_smoke, past_smoke,
              n_days_moderate_pa, duration_moderate_pa,
              n_days_vigorous_pa, duration_vigorous_pa, fruit, fresh_fruit, 
              dried_fruit, veg, cooked_veg, raw_veg, wholegrain, bread, cereal,
              bread_type, cereal_type, refined_grain, fish, oily_fish, 
              non_oily_fish, dairy, cheese, milk_cat, milk, spread_cat,
              spread, vegetable_oils, processed_meat_norm, processed_meat, 
              unprocessed_meat, poultry, beef, lamb, pork, added_sugar, 
              sugar_cat, sleep_duration_cat, sleep_duration, insomnia_cat,
              insomnia, snoring_cat, snoring, narcolepsy_cat,
              narcolepsy, birth_country, bmi, trigly, 
              lpa, glucose, creatinine, crp, apoa, apob, creatinine_std)) %>%
    select(-c(pc16:pc40)) %>%
    mutate(cad_icd = str_replace_all(cad_icd, "\"", ""),
           age2 = age^2) %>%
    mutate(occlusive_cad_bin = case_when(
      str_detect(cad_icd, "I2111|I2121|I2129|I213|I220|I228|I240|
                        I2510|I257|I2581|I2582|I2583|I2584|41090|41181|4140|4142|4143|4144") |
        str_detect(operation_pheno, "K401|K402|K403|K404|K411|K412|K413|K414|K451|K452|K453|K454|K455|K491|K492|
                        K498|K499|K502|K751|K752|K753|K754|K758|K759") ~ 1,
      TRUE ~ 0),
      unstable_cad_bin = case_when(
        # 1 = ACS, 0 = stable
        str_detect(cad_icd, "I201|I202|I203|I208|I209|I251|I259") ~ 0,
        !str_detect(cad_icd, "I201|I202|I203|I208|I209|I251|I259") & cad_bin == 1 ~ 1,
        TRUE ~ NA))
  
  write_csv(pheno_df, "cad_subtype/output_data/all_clean_original_ukb.csv")
  
  if(european_only) {
    pheno_df <- pheno_df %>%
      filter(race == 0)
    write_csv(pheno_df, "cad_subtype/output_data/clean_original_ukb.csv")
  }
  
  #NOTE: UKB SQL extracted ICD codes do not have "." in them
  
  # Define models ----
  base_vars <- "sex+age+age2+centre+batch+"
  sociodemographic_vars <- "bmi_std+ses_std+sys_bp_std+diabetes+"
  lifestyle_vars <- "smoke+pa+sleep+diet+"
  biomarker_vars <- "chol_std+ldl_std+hdl_std+glucose_std+crp_std+trigly_std+lpa_std+egfr+apob_std+"
  
  # Model 1: base
  mod1 <- paste("pc", 1:15, sep = "", collapse = "+") %>%
    paste0(outcome_var, "~", base_vars, .) %>%
    as.formula()
  
  # Model 2: base + sociodemographic
  mod2 <- paste("pc", 1:15, sep = "", collapse = "+") %>%
    paste0(outcome_var, "~", base_vars, sociodemographic_vars, .) %>%
    as.formula()
  
  # Model 3: base + lifestyle
  # note: smoke, physical activity, sleep, and diet are derived variables
  mod3 <- paste("pc", 1:15, sep = "", collapse = "+") %>%
    paste0(outcome_var, "~", base_vars, lifestyle_vars, .) %>%
    as.formula()
  
  # Model 4: base + biomarker
  # note: egfr is a derived variable from creatinine
  mod4 <- paste("pc", 1:15, sep = "", collapse = "+") %>%
    paste0(outcome_var, "~", base_vars, biomarker_vars, .) %>%
    as.formula()
  
  # Model 5: base + sociodemographic + lifestyle
  mod5 <- paste("pc", 1:15, sep = "", collapse = "+") %>%
    paste0(outcome_var, "~", base_vars, sociodemographic_vars, lifestyle_vars, .) %>%
    as.formula()
  
  # Model 6: base + sociodemographic + lifestyle + biomarker
  mod6 <- paste("pc", 1:15, sep = "", collapse = "+") %>%
    paste0(outcome_var, "~", base_vars, sociodemographic_vars,
           lifestyle_vars, biomarker_vars, .) %>%
    as.formula()
  
  # Prepare subtype dataframes ----
  switch(outcome_var,
         cad_bin = {
           pheno_df <- pheno_df %>%
             filter(statin_baseline == 0 | is.na(statin_baseline)) %>%
             as.data.frame()
         },
         ldl_bin = {
           pheno_df <- pheno_df %>%
             as.data.frame() %>% 
             filter(cad_bin == 1,
                    statin_baseline == 0 | is.na(statin_baseline),
                    !is.na(!!ensym(outcome_var))) %>%
             janitor::remove_empty(which = "cols")
         },
         lpa_bin = {
           pheno_df <- pheno_df %>%
             as.data.frame() %>% 
             filter(cad_bin == 1,
                    statin_baseline == 0 | is.na(statin_baseline),
                    !is.na(!!ensym(outcome_var))) %>%
             janitor::remove_empty(which = "cols")
         },
         premature_cad_bin = {
           pheno_df <- pheno_df %>%
             as.data.frame() %>% 
             filter(cad_bin == 1,
                    statin_baseline == 0 | is.na(statin_baseline),
                    !is.na(!!ensym(outcome_var))) %>%
             janitor::remove_empty(which = "cols")
         },
         occlusive_cad_bin = {
           pheno_df <- pheno_df %>%
             as.data.frame() %>% 
             filter(cad_bin == 1,
                    statin_baseline == 0 | is.na(statin_baseline),
                    !is.na(!!ensym(outcome_var))) %>%
             janitor::remove_empty(which = "cols")
         },
         unstable_cad_bin = {
           pheno_df <- pheno_df %>%
             as.data.frame() %>% 
             filter(cad_bin == 1,
                    statin_baseline == 0 | is.na(statin_baseline),
                    !is.na(!!ensym(outcome_var))) %>%
             janitor::remove_empty(which = "cols")
         },
         ascvd_cad_bin = {
           pheno_df <- pheno_df %>%
             as.data.frame() %>% 
             drop_na(age, chol, hdl_chol, sys_bp_auto, smoke, diabetes) %>% 
             mutate(ascvd = case_when(
               # untreated male
               sex == 0 & bp_med_bin == 0 ~ 1 - 0.9144^exp((12.344*log(age) + 11.853*log(chol) - 2.664*log(age)*log(chol) - 7.99*log(hdl_chol) +
                                                              1.769*log(age)*log(hdl_chol) + 1.764*log(sys_bp_auto) + 7.837*smoke - 1.795*log(age)*smoke + 0.658*diabetes) - 61.18),
               # treated male
               sex == 0 & bp_med_bin == 1 ~ 1 - 0.9144^exp((12.344*log(age) + 11.853*log(chol) - 2.664*log(age)*log(chol) - 7.99*log(hdl_chol) +
                                                              1.769*log(age)*log(hdl_chol) + 1.797*log(sys_bp_auto) + 7.837*smoke - 1.795*log(age)*smoke + 0.658*diabetes) - 61.18),
               # untreated female
               sex == 1 & bp_med_bin == 0 ~ 1 - 0.9665^exp((-29.799*log(age) + 4.884*log(age)^2 + 13.540*log(chol) - 3.114*log(age)*log(chol) - 13.578*log(hdl_chol) +
                                                              3.149*log(age)*log(hdl_chol) + 1.957*log(sys_bp_auto) + 7.574*smoke - 1.665*log(age)*smoke + 0.661*diabetes) + 29.18),
               # treated female
               sex == 1 & bp_med_bin == 1 ~ 1 - 0.9665^exp((-29.799*log(age) + 4.884*log(age)^2 + 13.540*log(chol) - 3.114*log(age)*log(chol) - 13.578*log(hdl_chol) +
                                                              3.149*log(age)*log(hdl_chol) + 2.019*log(sys_bp_auto) + 7.574*smoke - 1.665*log(age)*smoke + 0.661*diabetes) + 29.18)
             ),
             ascvd_cad_bin = ifelse(ascvd > 0.075, 1, 0)) %>% 
             filter(cad_bin == 1,
                    statin_baseline == 0 | is.na(statin_baseline),
                    !is.na(!!ensym(outcome_var))) %>%
             janitor::remove_empty(which = "cols")
         })
  
  if(impute) {
    cat("Running imputation")
    
    df_clean <- impute_data(pheno_df)
    df_clean <- complete(df_clean)
    
    df_resid_mod1 <- df_clean %>%
      mutate(pheno = glm(mod1, ., family = binomial) %>%
               resid())
    
    df_resid_mod2 <- df_clean %>%
      mutate(pheno = glm(mod2, ., family = binomial) %>%
               resid())
    
    df_resid_mod3 <- df_clean %>%
      mutate(pheno = glm(mod3, ., family = binomial) %>%
               resid())
    
    df_resid_mod4 <- df_clean %>%
      mutate(pheno = glm(mod4, ., family = binomial) %>%
               resid())
    
    df_resid_mod5 <- df_clean %>%
      mutate(pheno = glm(mod5, ., family = binomial) %>%
               resid())
    
    df_resid_mod6 <- df_clean %>%
      mutate(pheno = glm(mod6, ., family = binomial) %>%
               resid())
    
    fwrite(df_resid_mod1, paste0(out_base_path, "_mod1.csv"), sep = "\t")
    fwrite(df_resid_mod2, paste0(out_base_path, "_mod2.csv"), sep = "\t")
    fwrite(df_resid_mod3, paste0(out_base_path, "_mod3.csv"), sep = "\t")
    fwrite(df_resid_mod4, paste0(out_base_path, "_mod4.csv"), sep = "\t")
    fwrite(df_resid_mod5, paste0(out_base_path, "_mod5.csv"), sep = "\t")
    fwrite(df_resid_mod6, paste0(out_base_path, "_mod6.csv"), sep = "\t")
    
    return(list(mod1 = df_resid_mod1,
                mod2 = df_resid_mod2,
                mod3 = df_resid_mod3,
                mod4 = df_resid_mod4,
                mod5 = df_resid_mod5,
                mod6 = df_resid_mod6))
    
  } else {
    if(keep_base_only) {
      
      cat("Fitting models...")
      
      df_resid_mod1 <- pheno_df %>%
        mutate(pheno = glm(mod1, ., family = binomial, na.action = na.exclude) %>%
                 resid()) %>%
        filter(!is.na(pheno)) %>%
        pass_through(function(x) print(table(x %>% select(!!enquo(outcome_var))))) %>%
        select(c(iid, fid, pheno, age, age2, sex, centre, batch, pc1:pc15, !!enquo(outcome_var)))
      
      fwrite(df_resid_mod1, paste0(out_base_path, "_mod1.csv"), sep = "\t")
      
      return(list(mod1 = df_resid_mod1))
      
    } else {
      df_clean <- pheno_df %>%
        na.omit()
      
      cat("Fitting models...")
      
      df_resid_mod1 <- df_clean %>%
        mutate(pheno = glm(mod1, ., family = binomial) %>%
                 resid())
      
      df_resid_mod2 <- df_clean %>%
        mutate(pheno = glm(mod2, ., family = binomial) %>%
                 resid())
      
      df_resid_mod3 <- df_clean %>%
        mutate(pheno = glm(mod3, ., family = binomial) %>%
                 resid())
      
      df_resid_mod4 <- df_clean %>%
        mutate(pheno = glm(mod4, ., family = binomial) %>%
                 resid())
      
      df_resid_mod5 <- df_clean %>%
        mutate(pheno = glm(mod5, ., family = binomial) %>%
                 resid())
      
      df_resid_mod6 <- df_clean %>%
        mutate(pheno = glm(mod6, ., family = binomial) %>%
                 resid())
      
      fwrite(df_resid_mod1, paste0(out_base_path, "_mod1.csv"), sep = "\t")
      fwrite(df_resid_mod2, paste0(out_base_path, "_mod2.csv"), sep = "\t")
      fwrite(df_resid_mod3, paste0(out_base_path, "_mod3.csv"), sep = "\t")
      fwrite(df_resid_mod4, paste0(out_base_path, "_mod4.csv"), sep = "\t")
      fwrite(df_resid_mod5, paste0(out_base_path, "_mod5.csv"), sep = "\t")
      fwrite(df_resid_mod6, paste0(out_base_path, "_mod6.csv"), sep = "\t")
      
      return(list(mod1 = df_resid_mod1,
                  mod2 = df_resid_mod2,
                  mod3 = df_resid_mod3,
                  mod4 = df_resid_mod4,
                  mod5 = df_resid_mod5,
                  mod6 = df_resid_mod6))
    }
  }
}

# European-only ----
cad_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                               output_path = "cad_subtype/output_data/",
                               output_root_file_name = "cad_bin",
                               sql_csv = "cad_subtypes_v3.csv",
                               qc_file = "ukb18177-v2-qc.fam",
                               dropout_file = "w18177_2023-04-25.csv",
                               covariate_file = "ukb18177-v2.covar",
                               outcome_var = "cad_bin",
                               ses = FALSE,
                               impute = FALSE,
                               keep_base_only = TRUE,
                               european_only = TRUE)

ldl_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                               output_path = "cad_subtype/output_data/",
                               output_root_file_name = "ldl_bin",
                               sql_csv = "cad_subtypes_v3.csv",
                               qc_file = "ukb18177-v2-qc.fam",
                               dropout_file = "w18177_2023-04-25.csv",
                               covariate_file = "ukb18177-v2.covar",
                               outcome_var = "ldl_bin",
                               ses = FALSE,
                               impute = FALSE,
                               keep_base_only = TRUE,
                               european_only = TRUE)

lpa_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                               output_path = "cad_subtype/output_data/",
                               output_root_file_name = "lpa_bin",
                               sql_csv = "cad_subtypes_v3.csv",
                               qc_file = "ukb18177-v2-qc.fam",
                               dropout_file = "w18177_2023-04-25.csv",
                               covariate_file = "ukb18177-v2.covar",
                               outcome_var = "lpa_bin",
                               ses = FALSE,
                               impute = FALSE,
                               keep_base_only = TRUE)

# premature_cad_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
#                                          output_path = "cad_subtype/output_data/",
#                                          output_root_file_name = "premature_cad_bin",
#                                          sql_csv = "cad_subtypes_v3.csv",
#                                          qc_file = "ukb18177-v2-qc.fam",
#                                          dropout_file = "w18177_2023-04-25.csv",
#                                          covariate_file = "ukb18177-v2.covar",
#                                          outcome_var = "premature_cad_bin",
#                                          ses = FALSE,
#                                          impute = FALSE,
#                                          keep_base_only = TRUE)

occlusive_cad_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                         output_path = "cad_subtype/output_data/",
                                         output_root_file_name = "occlusive_cad_bin",
                                         sql_csv = "cad_subtypes_v3.csv",
                                         qc_file = "ukb18177-v2-qc.fam",
                                         dropout_file = "w18177_2023-04-25.csv",
                                         covariate_file = "ukb18177-v2.covar",
                                         outcome_var = "occlusive_cad_bin",
                                         ses = FALSE,
                                         impute = FALSE,
                                         keep_base_only = TRUE)

unstable_cad_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                        output_path = "cad_subtype/output_data/",
                                        output_root_file_name = "unstable_cad_bin",
                                        sql_csv = "cad_subtypes_v3.csv",
                                        qc_file = "ukb18177-v2-qc.fam",
                                        dropout_file = "w18177_2023-04-25.csv",
                                        covariate_file = "ukb18177-v2.covar",
                                        outcome_var = "unstable_cad_bin",
                                        ses = FALSE,
                                        impute = FALSE,
                                        keep_base_only = TRUE)

ascvd_cad_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                     output_path = "cad_subtype/output_data/",
                                     output_root_file_name = "ascvd_cad_bin",
                                     sql_csv = "cad_subtypes_v3.csv",
                                     qc_file = "ukb18177-v2-qc.fam",
                                     dropout_file = "w18177_2023-04-25.csv",
                                     covariate_file = "ukb18177-v2.covar",
                                     outcome_var = "ascvd_cad_bin",
                                     ses = FALSE,
                                     impute = FALSE,
                                     keep_base_only = TRUE)

# All ancestry ----
cad_pheno_all <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                   output_path = "cad_subtype/output_data/",
                                   output_root_file_name = "all_cad_bin",
                                   sql_csv = "cad_subtypes_v3.csv",
                                   qc_file = "ukb18177-allpop-qc.fam",
                                   dropout_file = "w18177_2023-04-25.csv",
                                   covariate_file = "ukb18177-v2.covar",
                                   outcome_var = "cad_bin",
                                   ses = FALSE,
                                   impute = FALSE,
                                   keep_base_only = TRUE,
                                   european_only = FALSE)

ldl_pheno_all <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                   output_path = "cad_subtype/output_data/",
                                   output_root_file_name = "all_ldl_bin",
                                   sql_csv = "cad_subtypes_v3.csv",
                                   qc_file = "ukb18177-allpop-qc.fam",
                                   dropout_file = "w18177_2023-04-25.csv",
                                   covariate_file = "ukb18177-v2.covar",
                                   outcome_var = "ldl_bin",
                                   ses = FALSE,
                                   impute = FALSE,
                                   keep_base_only = TRUE,
                                   european_only = FALSE)

lpa_pheno_all <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                   output_path = "cad_subtype/output_data/",
                                   output_root_file_name = "all_lpa_bin",
                                   sql_csv = "cad_subtypes_v3.csv",
                                   qc_file = "ukb18177-allpop-qc.fam",
                                   dropout_file = "w18177_2023-04-25.csv",
                                   covariate_file = "ukb18177-v2.covar",
                                   outcome_var = "lpa_bin",
                                   ses = FALSE,
                                   impute = FALSE,
                                   keep_base_only = TRUE,
                                   european_only = FALSE)

# premature_cad_pheno <- prepare_phenotype(input_path = "cad_subtype/input_data/",
#                                          output_path = "cad_subtype/output_data/",
#                                          output_root_file_name = "premature_cad_bin",
#                                          sql_csv = "cad_subtypes_v3.csv",
#                                          qc_file = "ukb18177-v2-qc.fam",
#                                          dropout_file = "w18177_2023-04-25.csv",
#                                          covariate_file = "ukb18177-v2.covar",
#                                          outcome_var = "premature_cad_bin",
#                                          ses = FALSE,
#                                          impute = FALSE,
#                                          keep_base_only = TRUE)

occlusive_cad_pheno_all <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                             output_path = "cad_subtype/output_data/",
                                             output_root_file_name = "all_occlusive_cad_bin",
                                             sql_csv = "cad_subtypes_v3.csv",
                                             qc_file = "ukb18177-allpop-qc.fam",
                                             dropout_file = "w18177_2023-04-25.csv",
                                             covariate_file = "ukb18177-v2.covar",
                                             outcome_var = "occlusive_cad_bin",
                                             ses = FALSE,
                                             impute = FALSE,
                                             keep_base_only = TRUE,
                                             european_only = FALSE)

unstable_cad_pheno_all <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                            output_path = "cad_subtype/output_data/",
                                            output_root_file_name = "all_unstable_cad_bin",
                                            sql_csv = "cad_subtypes_v3.csv",
                                            qc_file = "ukb18177-allpop-qc.fam",
                                            dropout_file = "w18177_2023-04-25.csv",
                                            covariate_file = "ukb18177-v2.covar",
                                            outcome_var = "unstable_cad_bin",
                                            ses = FALSE,
                                            impute = FALSE,
                                            keep_base_only = TRUE,
                                            european_only = FALSE)

ascvd_cad_pheno_all <- prepare_phenotype(input_path = "cad_subtype/input_data/",
                                         output_path = "cad_subtype/output_data/",
                                         output_root_file_name = "all_ascvd_cad_bin",
                                         sql_csv = "cad_subtypes_v3.csv",
                                         qc_file = "ukb18177-allpop-qc.fam",
                                         dropout_file = "w18177_2023-04-25.csv",
                                         covariate_file = "ukb18177-v2.covar",
                                         outcome_var = "ascvd_cad_bin",
                                         ses = FALSE,
                                         impute = FALSE,
                                         keep_base_only = TRUE,
                                         european_only = FALSE)

# Scratchpad ----
cad_pheno$adj %>%
  group_by(cad_bin) %>%
  count()

# acs cases
acs_pheno$adj %>%
  group_by(acs_bin) %>%
  count()

# angina subtypes cases
angina_subtype_pheno$adj %>%
  group_by(angina_subtype) %>%
  count()

# number of angina and MI cases
final_pheno$adj %>%
  group_by(clinical_subtype) %>%
  count()

final_pheno$adj %>%
  mutate(race = as.numeric(race),
         race = case_when(
           race %in% c(1, 1001, 1002, 1003) ~ 0,
           race %in% c(3001, 3002, 3003) ~ 1,
           race %in% c(3, 3004, 5) ~ 2,
           race %in% c(4, 4001, 4002, 4003) ~ 3,
           race %in% c(2, 2001, 2002, 2003, 2004) ~ 4,
           TRUE ~ 5
         )) %>%
  group_by(race) %>%
  count()

final_pheno$adj %>%
  filter(birth_country %in% 1:6) %>% 
  mutate(birth_bin = case_when(
    birth_country %in% 1:5 ~ 0,
    birth_country == 6 ~ 1
  )) %>%
  group_by(birth_bin,
           clinical_subtype) %>%
  count()

#   filter(race != 0) %>%
#   group_by(clinical_subtype) %>%
#   summarize(n = n(),
#             mean_age = mean(age, na.rm = TRUE),
#             sd_age = sd(age, na.rm = TRUE))

# train/test race
final_pheno$adj %>% 
  mutate(race = as.numeric(race),
         race = case_when(
           race %in% c(1, 1001, 1002, 1003) ~ 0,
           race %in% c(3001, 3002, 3003) ~ 1,
           race %in% c(3, 3004, 5) ~ 2,
           race %in% c(4, 4001, 4002, 4003) ~ 3,
           race %in% c(2, 2001, 2002, 2003, 2004) ~ 4,
           TRUE ~ 5
         )) %>%
  filter(race == 0) %>%
  group_by(clinical_subtype) %>%
  summarize(n = n(),
            mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE),
            prop_male = sum(sex == 1)/n)

final_pheno$adj %>%
  filter(birth_country %in% 1:6) %>% 
  mutate(birth_bin = case_when(
    birth_country %in% 1:5 ~ 0,
    birth_country == 6 ~ 1
  )) %>%
  filter(birth_bin == 0) %>%
  group_by(clinical_subtype) %>%
  summarize(n = n(),
            mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE),
            prop_male = sum(sex == 1)/n)

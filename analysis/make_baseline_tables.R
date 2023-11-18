# Create baseline tables for CAD subtypes cohorts for UKB 
# Author: Lathan Liou
# Contact: lathan.liou@icahn.mssm.edu

library(tidyverse)
library(table1)
library(data.table)
library(janitor)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(stats::chisq.test)
conflicts_prefer(table1::`label<-`)

raw_ukb <- fread("cad_subtype/input_data/cad_subtypes_v3.csv")
dropout <- fread("cad_subtype/input_data/w18177_2023-04-25.csv", header = F)
qc <- fread("cad_subtype/input_data/ukb18177-v2-qc.fam")
qc_all <- fread("cad_subtype/input_data/ukb18177-allpop-qc.fam")
cov <- fread("cad_subtype/input_data/ukb18177-v2.covar")

# Functions ----
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

standardize <- function(x){
  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}

# Data Processing ----
processed_ukb <- raw_ukb %>%
  filter(!IID %in% dropout$V1,
         IID %in% qc_all$V2) %>%
  janitor::clean_names() %>%
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
  select(-c(dx_date, dad_fam_hx, mom_fam_hx, sib_fam_hx, statin_follow)) %>%
  mutate(across(diabetes:bp_med2, convert_idk_to_na)) %>%
  mutate(across(diabetes:bp_med2, convert_neg_10)) %>%
  mutate(age2 = age^2,
         bp_med = coalesce(bp_med1, bp_med2),
         bp_med_bin = ifelse(bp_med == 2, 1, 0), #antihypertensives
         premature_cad_bin = case_when(
           #0 = female, 1 = male
           # age is age at diagnosis
           sex == 1 & age < 40 & cad_bin == 1 ~ 1,
           sex == 1 & age >= 40 & cad_bin == 1 ~ 0,
           sex == 0 & age < 50 & cad_bin == 1 ~ 1,
           sex == 0 & age >= 50 & cad_bin == 1 ~ 0,
           TRUE ~ 0
         ),
         smoke = case_when(
           current_smoke == 1 | past_smoke == 1 ~ 1,
           TRUE ~ 0
         ),
         creatinine_std = standardize(creatinine),
         egfr = case_when(
           sex == 0 ~ 142 * pmin(creatinine_std/0.9, 1, na.rm = TRUE)^-0.302 *
             pmax(creatinine_std/0.9, 1, na.rm = TRUE)^-1.200 *
             0.9938^age,
           sex == 1 ~ 142 * pmin(creatinine_std/0.7, 1, na.rm = TRUE)^-0.241 *
             pmax(creatinine_std/0.7, 1, na.rm = TRUE)^-1.200 *
             0.9938^age*1.012),
         ldl_chol = ldl_chol*18.0182,
         ldl_bin = ifelse(ldl_chol >= 70, 1, 0),
         lpa_bin = ifelse(lpa >= 150, 1, 0),
         cad_icd = str_replace_all(cad_icd, "\"", ""),
         unstable_cad_bin = case_when(
           # 1 = ACS, 0 = stable
           str_detect(cad_icd, "I201|I202|I203|I208|I209|I251|I259") ~ 0,
           !str_detect(cad_icd, "I201|I202|I203|I208|I209|I251|I259") & cad_bin == 1 ~ 1,
           TRUE ~ NA),
         occlusive_cad_bin = case_when(
           str_detect(cad_icd, "I2111|I2121|I2129|I213|I220|I228|I240|
                        I2510|I257|I2581|I2582|I2583|I2584|41090|41181|4140|4142|4143|4144") |
             str_detect(operation_pheno, "K401|K402|K403|K404|K411|K412|K413|K414|K451|K452|K453|K454|K455|K491|K492|
                        K498|K499|K502|K751|K752|K753|K754|K758|K759") ~ 1,
           TRUE ~ 0)) %>% 
  filter(statin_baseline == 0 | is.na(statin_baseline))

processed_ukb2 <- processed_ukb %>%
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
  drop_na(ascvd_cad_bin)

# european counts
processed_ukb2 %>%
  filter(race == 0,
         cad_bin == 1) %>%
  group_by(ascvd_cad_bin) %>%
  count()

processed_ukb_renamed <- processed_ukb 
processed_ukb_renamed2 <- processed_ukb2

tstat <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$statistic
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$statistic
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", format(round(p, 1), nsmall = 1))
}

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

custom_render <- function(x, ...) {
  if(is.numeric(x)) {
    parse.abbrev.render.code(c("", "Mean (SD)"))(x)
  } else {
    return(render.categorical.default(x))
  }
}

# Baseline Table CAD ----
processed_ukb_renamed$sex <- factor(processed_ukb_renamed$sex,
                                    levels = 0:1,
                                    labels = c("Male",
                                               "Female"))

label(processed_ukb_renamed$sex) <- "Sex"

label(processed_ukb_renamed$age) <- "Age"

processed_ukb_renamed$race <- factor(processed_ukb_renamed$race,
                                     levels = 0:5,
                                     labels = c("White",
                                                "Indian",
                                                "East Asian",
                                                "African",
                                                "Mixed",
                                                "Unknown"))

label(processed_ukb_renamed$race) <- "Ancestry"

processed_ukb_renamed %>%
  group_by(cad_bin) %>%
  count(race)

processed_ukb_renamed$diabetes <- factor(processed_ukb_renamed$diabetes,
                                         levels = 0:1,
                                         labels = c("No",
                                                    "Yes"))
label(processed_ukb_renamed$diabetes) <- "Diabetes"

label(processed_ukb_renamed$hdl_chol) <- "HDL Cholesterol (mmol/L)"

label(processed_ukb_renamed$ldl_chol) <- "LDL Cholesterol (mmol/L)"

label(processed_ukb_renamed$chol) <- "Total Cholesterol (mmol/L)"

label(processed_ukb_renamed$trigly) <- "Triglycerides (mmol/L)"

label(processed_ukb_renamed$lpa) <- "Lipoprotein A (nmol/L)"

label(processed_ukb_renamed$apob) <- "Apolipoprotein B (g/L)"

label(processed_ukb_renamed$crp) <- "C-reactive Protein (mg/L)"

label(processed_ukb_renamed$egfr) <- "Estimated GFR (umol/L)"

processed_ukb_renamed$smoke <- factor(processed_ukb_renamed$smoke,
                                      levels = 0:1,
                                      labels = c("No",
                                                 "Yes"))
label(processed_ukb_renamed$smoke) <- "Smoked ever"

label(processed_ukb_renamed$sys_bp_auto) <- "Systolic Blood Pressure (mmHg)"

label(processed_ukb_renamed$dias_bp_auto) <- "Diastolic Blood Pressure (mmHg)"

label(processed_ukb_renamed$bmi) <- "BMI (kg/m2)"

processed_ukb_renamed$cad_bin <- factor(processed_ukb_renamed$cad_bin,
                                        levels = 0:1,
                                        labels = c("No CAD", "CAD"))

table1(~ age + sex + race + diabetes + smoke + hdl_chol + ldl_chol + chol + trigly +
         lpa + apob + crp + egfr + sys_bp_auto + dias_bp_auto + bmi | cad_bin,
       data = processed_ukb_renamed,
       overall = FALSE,
       extra.col=list(`Statistic`=tstat,
                      `P-value`=pvalue))

# Baseline Table Stable Subtype ----
processed_ukb_renamed$unstable_cad_bin <- factor(processed_ukb_renamed$unstable_cad_bin,
                                                  levels = 0:1,
                                                  labels = c("Stable Angina", 
                                                             "Acute Coronary Syndrome"))
table1(~ age + sex + race + diabetes + smoke + hdl_chol + ldl_chol + chol + trigly +
         lpa + apob + crp + egfr + sys_bp_auto + dias_bp_auto + bmi | unstable_cad_bin,
       data = processed_ukb_renamed %>%
         filter(cad_bin == "CAD",
                !is.na(unstable_cad_bin)),
       overall = FALSE,
       extra.col=list(`Statistic`=tstat,
                      `P-value`=pvalue))

# Baseline Table Occlusive Subtype ----
processed_ukb_renamed$occlusive_cad_bin <- factor(processed_ukb_renamed$occlusive_cad_bin,
                                                  levels = 0:1,
                                                  labels = c("Non-occlusive CAD", 
                                                             "Occlusive CAD"))

table1(~ age + sex + race + diabetes + smoke + hdl_chol + ldl_chol + chol + trigly +
         lpa + apob + crp + egfr + sys_bp_auto + dias_bp_auto + bmi | occlusive_cad_bin,
       data = processed_ukb_renamed %>%
         filter(cad_bin == "CAD",
                !is.na(occlusive_cad_bin)),
       overall = FALSE,
       extra.col=list(`Statistic`=tstat,
                      `P-value`=pvalue))

# Baseline Table LDL Subtype ----
processed_ukb_renamed$ldl_bin <- factor(processed_ukb_renamed$ldl_bin,
                                        levels = 0:1,
                                        labels = c("Normal LDL", "High LDL (>3.88 mmol/L)"))
table1(~ age + sex + race + diabetes + smoke + hdl_chol + chol + trigly +
         lpa + apob + crp + egfr + sys_bp_auto + dias_bp_auto + bmi | ldl_bin,
       data = processed_ukb_renamed %>%
         filter(cad_bin == "CAD",
                !is.na(ldl_bin)),
       overall = FALSE,
       extra.col=list(`Statistic`=tstat,
                      `P-value`=pvalue))

# Baseline Table Lpa Subtype ----
processed_ukb_renamed$lpa_bin <- factor(processed_ukb_renamed$lpa_bin,
                                        levels = 0:1,
                                        labels = c("Normal Lpa", "High Lpa (>150 nmol/L)"))
table1(~ age + sex + race + diabetes + smoke + ldl_chol + hdl_chol + chol + trigly +
         apob + crp + egfr + sys_bp_auto + dias_bp_auto + bmi | lpa_bin,
       data = processed_ukb_renamed %>%
         filter(cad_bin == "CAD",
                !is.na(lpa_bin)),
       overall = FALSE,
       extra.col=list(`Statistic`=tstat,
                      `P-value`=pvalue))

# Baseline Table ASCVD Subtype ----
processed_ukb_renamed2$sex <- factor(processed_ukb_renamed2$sex,
                                    levels = 0:1,
                                    labels = c("Male",
                                               "Female"))

label(processed_ukb_renamed2$sex) <- "Sex"

label(processed_ukb_renamed2$age) <- "Age"

processed_ukb_renamed2$race <- factor(processed_ukb_renamed2$race,
                                     levels = 0:5,
                                     labels = c("White",
                                                "Indian",
                                                "East Asian",
                                                "African",
                                                "Mixed",
                                                "Unknown"))

label(processed_ukb_renamed2$race) <- "Ancestry"

processed_ukb_renamed2$diabetes <- factor(processed_ukb_renamed2$diabetes,
                                         levels = 0:1,
                                         labels = c("No",
                                                    "Yes"))
label(processed_ukb_renamed2$diabetes) <- "Diabetes"

label(processed_ukb_renamed2$hdl_chol) <- "HDL Cholesterol (mmol/L)"

label(processed_ukb_renamed2$ldl_chol) <- "LDL Cholesterol (mmol/L)"

label(processed_ukb_renamed2$chol) <- "Total Cholesterol (mmol/L)"

label(processed_ukb_renamed2$trigly) <- "Triglycerides (mmol/L)"

label(processed_ukb_renamed2$lpa) <- "Lipoprotein A (nmol/L)"

label(processed_ukb_renamed2$apob) <- "Apolipoprotein B (g/L)"

label(processed_ukb_renamed2$crp) <- "C-reactive Protein (mg/L)"

label(processed_ukb_renamed2$egfr) <- "Estimated GFR (umol/L)"

processed_ukb_renamed2$smoke <- factor(processed_ukb_renamed2$smoke,
                                      levels = 0:1,
                                      labels = c("No",
                                                 "Yes"))
label(processed_ukb_renamed2$smoke) <- "Smoked ever"

label(processed_ukb_renamed2$sys_bp_auto) <- "Systolic Blood Pressure (mmHg)"

label(processed_ukb_renamed2$dias_bp_auto) <- "Diastolic Blood Pressure (mmHg)"

label(processed_ukb_renamed2$bmi) <- "BMI (kg/m2)"

processed_ukb_renamed2$cad_bin <- factor(processed_ukb_renamed2$cad_bin,
                                        levels = 0:1,
                                        labels = c("No CAD", "CAD"))

processed_ukb_renamed2$ascvd_cad_bin <- factor(processed_ukb_renamed2$ascvd_cad_bin,
                                        levels = 0:1,
                                        labels = c("Low ASCVD Risk", "High ASCVD Risk"))
table1(~ age + sex + race + diabetes + smoke + ldl_chol + hdl_chol + chol + trigly +
         apob + crp + egfr + sys_bp_auto + dias_bp_auto + bmi | ascvd_cad_bin,
       data = processed_ukb_renamed2 %>%
         filter(cad_bin == "CAD",
                !is.na(ascvd_cad_bin)),
       overall = FALSE,
       extra.col=list(`Statistic`=tstat,
                      `P-value`=pvalue))

# Baseline Table BioMe ----
# biome_df_renamed$sex <- factor(biome_df_renamed$sex,
#                                levels = 1:2,
#                                labels = c("Male", "Female"))
# label(biome_df_renamed$sex) <- "Sex"
# 
# label(biome_df_renamed$age) <- "Age"
# 
# label(biome_df_renamed$ancestry) <- "Ancestry"
# 
# label(biome_df_renamed$chol) <- "Total Cholesterol (mmol/L)"
# 
# label(biome_df_renamed$hdl) <- "HDL Cholesterol (mmol/L)"
# 
# label(biome_df_renamed$ldl) <- "LDL Cholesterol (mmol/L)"
# 
# label(biome_df_renamed$sys_bp) <- "Systolic Blood Pressure (mmHg)"
# 
# label(biome_df_renamed$pers_hx_smoking) <- "Smoked Ever"
# 
# biome_df_renamed$clinical_subtype <- factor(biome_df_renamed$clinical_subtype,
#                                             levels = 0:1,
#                                             labels = c("Angina", "MI"))
# 
# table1(~ age + sex + ancestry + pers_hx_smoking + hdl + ldl + chol + sys_bp | clinical_subtype,
#        data = biome_df_renamed,
#        overall = FALSE)

# Prepare base GWAS data
# What this does: processes/standardizes base GWAS data 
# Author: Lathan Liou
# Contact: lathan.liou@icahn.mssm.edu
# Inputs: base GWAS data
# Outputs:

library(data.table)
library(tidyverse)
library(here)
library(janitor)
input_path <- "cad_subtype/input_data/"
out_path <- "cad_subtype/output_data/"

sumstat <- fread(here(input_path, "GCST90132314_buildGRCh37_WITHrsID.tsv")) %>%
  janitor::clean_names()
  mutate(a2 = case_when(
    a1 == "A" ~ "T",
    a1 == "T" ~ "A",
    a1 == "C" ~ "G",
    a1 == "G" ~ "C"
  ),
  lb = as.numeric(str_extract(x95_percent_ci_text, "(?<=\\[)(.*)(?=\\-)")),
  se = (lb - or_or_beta)/(-1.96)) 
  
is_beta <- TRUE
data <- "ukb"

# Process base GWAS data ----

if(is_beta & data == "ukb") {
  
  sumstat_dt <- data.table(sumstat)
  setnames(sumstat_dt,c("snp", "a1", "a2", "p_value", "or_or_beta", "se"),
           c("snp", "a1", "a2", "p", "beta", "se"))
  sumstat_clean_dt <- sumstat_dt[, .(chromosome = chr_id[which.min(p)],
                                     base_pair_location = chr_pos[which.min(p)],
                                     a1 = a1[which.min(p)],
                                     a2 = a2[which.min(p)],
                                     beta = beta[which.min(p)], p = min(p), se = se[which.min(p)], n = .N),
                                 by = snp][order(p)]
} else if(!is_beta & data == "ukb") {
  sumstat_clean <- sumstat %>%
    rename(snp = markername,
           a1 = effect_allele,
           a2 = noneffect_allele,
           p = p_dgc,
           se = se_dgc) %>%
    select(snp, a1, a2, beta, p, se) %>%
    mutate(beta = log(beta),
           n = n())
}

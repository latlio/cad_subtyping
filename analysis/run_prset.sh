#!/bin/bash
# Perform analyses for CAD

### All Pop
# 1. run prset (AUC, binomial)
./PRSice --base /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad-std.sumstat.txt \
        --target /sc/arion/projects/data-ark/ukb/application/ukb18177/genotyped/ukb18177 \
        --out cad_bin/c4d_prset_all_ukb_cad_bin \
        --keep /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad_bin/all_cad_bin_mod1.csv \
        --extract /sc/arion/projects/data-ark/ukb/application/ukb18177/shared_pheno/ukb18177-allpop-qc.snplist \
        --snp snp \
        --pvalue p \
        --beta \
        --stat beta \
        --a1 a1 \
        --a2 a2 \
        --binary-target T \
        --pheno /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad_bin/all_cad_bin_mod1.csv \
        --pheno-col cad_bin \
        --cov /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad_bin/all_cad_bin_mod1.csv \
        --cov-col age,age2,sex,centre,batch,@PC[1-15] \
        --cov-factor centre,batch \
        --gtf /sc/arion/projects/psychgen/forLathan/Homo_sapiens.GRCh37.75.gtf.gz \
        --msigdb ebi_cvd_exponly.gmt,gene_sets/reactome_clean.gmt,gene_sets/pid_clean.gmt,gene_sets/kegg_clean.gmt,gene_sets/GO_clean.gmt,gene_sets/biocarta_clean.gmt,gene_sets/MGI_clean.gmt \
        --print-snp \
        --ultra \
        --wind-3 5 \
        --wind-5 10 \
        --thread 10 \
        --set-perm 10000

# 2. run prsice (AUC, bin)

./PRSice --base /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad-std.sumstat.txt \
        --target /sc/arion/projects/data-ark/ukb/application/ukb18177/genotyped/ukb18177 \
        --out cad_bin/c4d_prsice_all_ukb_cad_bin \
        --keep /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad_bin/all_cad_bin_mod1.csv \
        --extract /sc/arion/projects/data-ark/ukb/application/ukb18177/shared_pheno/ukb18177-allpop-qc.snplist \
        --snp snp \
        --pvalue p \
        --beta \
        --stat beta \
        --a1 a1 \
        --a2 a2 \
        --binary-target T \
        --pheno /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad_bin/all_cad_bin_mod1.csv \
        --pheno-col cad_bin \
        --cov /sc/arion/projects/psychgen/projects/prs/cad_subtype/cad_bin/all_cad_bin_mod1.csv \
        --cov-col age,age2,sex,centre,batch,@PC[1-15] \
        --cov-factor centre,batch \
        --print-snp \
        --ultra 

#3. build biomarker models + test (AUC), requires ssh -X to forward graphical display

Rscript analyze_univariate_models.R --original_ukb_path "all_clean_original_ukb.csv" \
                        --pheno_df_path "cad_bin/all_cad_bin_mod1.csv" \
                        --results_path "cad_bin" \
                        --train_prop 0.8 \
                        --run_name "c4d_all_ukb_cad_bin" \
                        --prsice_run_name "c4d_prsice_all_ukb_cad_bin" \
                        --prset_run_name "c4d_prset_all_ukb_cad_bin" \
                        --outcome_var "cad_bin" \
                        --seed 47 \
                        --threshold 0.5
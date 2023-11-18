# Code for "Prediction of Coronary Artery Disease Subtypes Using Clinical and Genetic Risk Factors"

**Lathan Liou**
Judit Garcia Gonzalez
Beatrice Wu
Zhe Wang
Shing Wan Choi
Clive J Hoggart
Amy R Kontorovich
Jason C Kovacic
Paul O'Reilly

This paper has been submitted for publication in Nature Medicine

![](/figures/main_figure.png)

*Comparative performance of single variable, genome-wide polygenic risk score, and pathway polygenic risk scores in subtype classification*


## Abstract

Coronary Artery Disease (CAD) is a complex, heterogeneous disease with distinct underlying etiological mechanisms. These different etiologies may give rise to multiple subtypes of CAD that would benefit from alternative therapeutic approaches. However, whether CAD subtypes can be successfully predicted using clinical and genetic factors is unknown. Here we developed and optimized statistical models incorporating clinical and genetic factors to predict CAD subtypes in 20,935 CAD patients in the UK Biobank, with external validation in the US-based All of Us (9,962 CAD patients). Subtypes were defined as high vs. normal low-density lipoprotein (LDL) levels, high vs. normal lipoprotein A (Lpa) levels, occlusive vs. non-occlusive CAD, stable vs. unstable, high vs. low 10-year atherosclerotic cardiovascular disease (ASCVD risk) score. Clinical risk factors included levels of apolipoprotein A, apolipoprotein B, high-density lipoprotein, triglycerides, and C-reactive protein. Genetic factors were genome-wide and pathway-based polygenic risk scores (PRS).  Results show models that combine clinical and genetic factors lead to greater predictive power than clinical-only and genetics-only models. Pathway-based PRS typically increased discriminatory power more than genome-wide PRS, while providing insights into subtype etiology. For example, pathways most predictive of the high-Lpa subtype included fibronectin binding (P=6x10-260), apolipoprotein binding (P=6.2x10-252), and serine-type endopeptidase activity (P=1.7x10-203). Overall, we demonstrate that prediction models can distinguish CAD subtypes and that pathway PRSs hold particular promise for increased biological understanding and tailored treatment of CAD subtypes.


## Software implementation

This folder is organized into the `preprocessing`, `analysis`, `data`, and `figures` directories. Scripts were run within RStudio. For All of Us, scripts were copied and pasted into hosted Jupyter Notebooks. The Notebooks cannot be shared due to confidentiality agreements signed with All of Us.

* Note: `data` does not contain every data file that was used in our analysis (see [Not Included Here](#not-incluided-here)). It contains the drug targettor database we used for our SNP-gene-drug annotation analysis as well as the input pathway gene sets we used for PRSet.


## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://github.com/latlio/cad_subtyping.git) repository:

    git clone https://github.com/latlio/cad_subtyping.git

or [download a zip archive](https://github.com/latlio/cad_subtyping/archive/refs/heads/main.zip).

## Dependencies

    install.packages(c("nnet", "pROC", "tidyverse", "conflicted", "data.table", "dplyr", "glmnet", "janitor", "rmarkdown", "here", "knitr", "lubridate", "readxl", "broom", "readr", "stringr", "testthat", "tibble", "devtools", "PRSmix", "caret", "GO.db", "org.Hs.eg.db", "biomaRt", "GOfuncR", "networkD3", "qs", "DT", "ggVennDiagram", "GOSemSim", "kableExtra", "patchwork", "cluster", "magrittr", "doParallel", "boot", "ggrepel", "ggsci", "glue", "optparse", "PredictABEL", "magick", "AnnotationDbi", "base", "scales", "tidytext", "stats", "table1"))

And see the [PRSice documentation](https://github.com/choishingwan/PRSice) for how to install the PRSice software that is used for PRSice and PRSet implementation.

## Navigating the Code Directory
1. `preprocessing/00_prep_base_data.R` preps the base GWAS data Nikpay 2015 (https://doi.org/10.1038/ng.3396)
2. `preprocessing/01_prep_phenotype.R` preps the sample data (UK Biobank)
3. `analysis/make_baseline_table.R` contains code that creates summary baseline tables
4. `analysis/run_prset.sh` contains code that is meant to be run on the command line (e.g. `./analysis/run_prset.sh` on a Mac/Linux). This code has 3 functions: running PRSice, PRSet, and analyze_univariate_models.R
4. `analysis/analyze_univariate_models.R` contains the main code that runs most of the regressions and generates most of the figures. This can be run either interactively or via the command line (see `analysis/run_prset.sh`).
5. `analysis/make_dotplots.R` contains code that creates the dotplots (Figure 4) of pathways
6. `analysis/analyze_drug_targets.R` contains code that links our pathways to a gene drug database   
7. `analysis/common_prset_fcts` contains general purpose functions used throughout the codebase

Make sure you change paths to the relevant ones on your machine/coding environment.

## Not Included Here

* Script extracting data from our copy of UK Biobank as data access is provisioned upon application
* Jupyter notebooks containing processing and analysis code for All of Us, although they were essentially lifted from code made available here
* Our raw pathway PRS files as they are very large. We can potentially provide them upon reasonable request.


## License

All source code is made available under a BSD 3-clause license. You can freely
use and modify the code, without warranty, so long as you provide attribution
to the authors. See `LICENSE.md` for the full license text.

The manuscript text is not open source. The authors reserve the rights to the
article content, which is currently submitted for publication in the
Nature Medicine.
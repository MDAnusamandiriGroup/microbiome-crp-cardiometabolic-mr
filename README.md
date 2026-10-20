# Microbiome CRP and cardiometabolic Mendelian randomization

This repository packages the current hybrid R pipeline for microbiome MR using CRP as the main trait, with replication checks and portability analyses for type 2 diabetes and systolic blood pressure. It also includes an explainable ranking layer for integrating MR evidence.

The ranking model reweights features derived from the MR results. It is not an independent learner of causal truth. The code separates primary analyses from replication and supportive comparisons.

## Data

Input data are not included. The pipeline reads local, previously prepared CSV files and makes no GWAS downloads by itself. Keep them outside this prepared upload folder; a browser upload can include any file you add to the selected folder, even if it matches a .gitignore rule. See **data/README.md** for the expected file names and study IDs. Before sharing input files, check the terms of the original GWAS data providers.

## Requirements

Use R and install the following packages:

    install.packages(c("data.table", "xgboost", "ggplot2", "dplyr", "tidyr"))

Install **TwoSampleMR** using the current instructions in its official documentation: https://github.com/MRCIEU/TwoSampleMR

## Run

From the repository root, place the input CSVs in **data/private/**. Then run in R:

    Sys.setenv(AIXBIO_DATA_DIR = "data/private", AIXBIO_OUTPUT_DIR = "results")
    source("pipeline.R")

The default run uses ten repeats of five-fold cross-validation, 300 model-refit bootstraps, and 300 fixed-score bootstrap samples. A full run can take time. For a smoke test, first reduce those settings near the top of **pipeline.R**, then restore them for the full analysis.

The output folder contains run manifests, audit tables, rankings, and figures. This repository does not include generated outcome tables because they require the study input files and a fresh run.

The MIT license applies to the original code and documentation in this repository. It does not apply to third-party data.

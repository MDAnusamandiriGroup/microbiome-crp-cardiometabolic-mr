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

## Proxy sensitivity and colocalisation

The repository also includes separate scripts for proxy sensitivity and targeted colocalisation. They use prepared local input files and do not download GWAS data or retrieve proxy variants. The required proxy map, regional summary statistics and signed LD matrix are not included here.

These are new implementations based on the methods described in the manuscript, not copies of the original analysis scripts. They have not yet been checked against the manuscript's reported results. See [`scripts/ANALYSIS_EXTENSIONS.md`](scripts/ANALYSIS_EXTENSIONS.md) for the required input formats and assumptions.

Once the input files are prepared, run:

```sh
Rscript --vanilla scripts/run_proxy_sensitivity.R
Rscript --vanilla scripts/run_colocalisation.R

The proxy analysis requires TwoSampleMR. Colocalisation requires coloc; susieR is optional for the multi-signal analysis. Results are written under results/proxy_sensitivity/runs/ and results/colocalisation/runs/.

The MIT license applies to the original code and documentation in this repository. It does not apply to third-party data.

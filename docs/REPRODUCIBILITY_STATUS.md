# Reproducibility status

The repository contains the main R pipeline for the exact-match CRP analysis, CRP comparison, T2D and SBP portability analyses, and exploratory XGBoost/SHAP ranking. It also now includes separate local-input workflows for proxy sensitivity and targeted colocalisation in `scripts/`. See [`scripts/ANALYSIS_EXTENSIONS.md`](../scripts/ANALYSIS_EXTENSIONS.md) for required inputs, assumptions, and commands.

The extension scripts are new, documented implementations based on the methods described in the author's manuscript revision. They are not the original unpublished scripts. The required proxy map, regional summary statistics, and signed LD matrix are not included in this repository, so the analyses cannot yet be run against the original study inputs or checked against the manuscript's reported numerical results.

Other limitations to keep in mind:

- The repository does not include the full per-taxon result tables or underlying GWAS inputs. Input files belong under `data/private/` and should only be shared if their source terms allow redistribution.
- The main pipeline currently names `ieu-b-35` as its CRP comparison input. The manuscript summary reports external CRP results from the UK Household Longitudinal Study (`ebi-a-GCST005067`, N=9,541). These are different GWAS sources; the current pipeline should not be described as reproducing that UKHLS analysis without an input and code update.
- The reported results in `results/REPORTED_EMPIRICAL_SUMMARY.md` and `results/REPORTED_CANDIDATE_RESULTS.csv` were transcribed from the author's manuscript revision dated 12 October 2026. They are a record of reported findings, not output regenerated from the public repository.

Before describing the repository as a full reproduction, supply and document the original inputs and exact selection rules, verify the extension scripts against the original workflow, and resolve the difference between the pipeline's CRP comparison dataset and the manuscript's UKHLS dataset.

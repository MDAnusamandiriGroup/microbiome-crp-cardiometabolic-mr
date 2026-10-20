# Reported empirical summary

This is a short record of the results reported in the author's manuscript revision, **“Reliability of microbiome MR signals for C-reactive protein” (12 October 2026)**. It is not a fresh run of the public R pipeline. The underlying study inputs and full result tables are not included in this repository.

## Study scope

The manuscript reports a screen of 468 shotgun-metagenomic gut microbial taxa. The frozen exposure set contained 4,775 taxon–SNP rows and 4,135 unique variants. CRP was the main outcome; T2D and systolic blood pressure (SBP) were secondary portability outcomes.

## Main findings as reported

- In the exact-match primary CRP analysis, four taxa passed FDR correction. Two of the strongest estimates, for *Collinsella* and *Faecalicatena torques*, each used one instrument and were treated as exploratory.
- In the external UK Household Longitudinal Study CRP comparison, all four candidates had direction-concordant estimates, but none met the candidate-level replication threshold (two-sided P < 0.0125).
- The proxy-augmented CRP sensitivity analysis returned estimates for 461 taxa. It found 31 nominal associations but no FDR-significant associations. Among 363 taxa with results in both branches, 281 (77.4%) retained the same direction. Exact-versus-proxy ranking agreement was modest (Spearman ρ = 0.326; top-20 Jaccard = 0.212).
- Targeted colocalisation supported shared primary CRP–microbiome locus signals for *Collinsella* and *F. torques* at the rs635634 region (maximum PP.H4 0.993 and 0.988). It did not provide convincing support for the tested *Comamonas* or CAG-390 loci, and external CRP colocalisation was not convincing. Shared-locus evidence does not establish that microbiome changes cause CRP changes.
- No taxa passed FDR correction in the primary T2D or SBP portability screens.
- The XGBoost/SHAP analysis was exploratory and reconstructed an MR-derived evidence grade; it was not treated as independent evidence of causality.

## Reading the candidate table

`REPORTED_CANDIDATE_RESULTS.csv` gives the four exact-primary CRP candidates, their exact-match estimates, external UKHLS estimates, proxy-branch estimates and reported colocalisation summaries. Blank values mean the manuscript table did not report that quantity for the candidate.

## Important reproducibility note

The repository now includes new proxy-sensitivity and colocalisation scripts, but the code and regional data needed to rerun the manuscript analyses are still incomplete. The scripts require the original proxy map, dense regional summary statistics and signed LD matrix, which are not included here, and have not yet been checked against the reported estimates. The pipeline also uses `ieu-b-35` as its CRP comparison input, while the manuscript's external results above use UKHLS `ebi-a-GCST005067`. These summary values should therefore be cited as manuscript-reported findings, not as outputs regenerated from the files currently in the repository. `ieu-b-35` is a CHARGE CRP dataset, whereas `ebi-a-GCST005067` is the UKHLS dataset reported in the manuscript; the studies are not interchangeable. See `docs/REPRODUCIBILITY_STATUS.md` and `scripts/ANALYSIS_EXTENSIONS.md`.

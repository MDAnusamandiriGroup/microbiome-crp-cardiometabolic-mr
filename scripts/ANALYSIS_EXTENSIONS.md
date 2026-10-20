# Proxy sensitivity and colocalisation extensions

These scripts add transparent, local-input workflows for the analyses described in the microbiome MR manuscript. They do not download summary statistics or query a proxy service. The required inputs must be obtained and checked against the source GWAS terms before running.

## Packages

The proxy workflow needs `data.table` and `TwoSampleMR`. The ABF colocalisation workflow needs `data.table` and `coloc`. Add `susieR` if a region manifest includes an LD matrix for multi-signal colocalisation.

```r
install.packages(c("data.table", "coloc", "susieR"))
```

Install `TwoSampleMR` using the package's current official instructions: <https://github.com/MRCIEU/TwoSampleMR>.

## Proxy sensitivity inputs

Put four CSV files in `data/private/proxy_sensitivity/`:

1. `exposure_instruments.csv` has `taxon`, `SNP`, `beta.exposure`, `se.exposure`, `effect_allele.exposure`, `other_allele.exposure`, `eaf.exposure`, and `pval.exposure`. The script also accepts `taxa` instead of `taxon` for the historical input file.
2. `crp_exact.csv` contains exact-match rows from the primary CRP GWAS, using `SNP`, `beta.outcome`, `se.outcome`, `effect_allele.outcome`, `other_allele.outcome`, `eaf.outcome`, and `pval.outcome`.
3. `crp_proxy.csv` contains outcome statistics for physical proxy variants, using the same outcome columns. Its `SNP` values must be proxy variant IDs, not target IDs.
4. `proxy_pairs.csv` maps each `target_SNP` to a `proxy_SNP` and contains `r2`, `signed_ld`, `target_allele`, `proxy_allele`, and `maf_eligible`.

`r2` must be at least 0.80 for a proxy to be eligible. `signed_ld` is the signed correlation between the allele named in `target_allele` and the allele named in `proxy_allele`, from the LD reference used for the analysis. The script checks that `r2` is approximately `signed_ld^2` and uses the sign to align the proxy outcome beta to the exposure target allele. It does not scale the beta by the magnitude of the LD correlation.

`maf_eligible` must be an explicit TRUE/FALSE flag computed using the exact minor-allele-frequency rule intended for the analysis. The manuscript describes a 0.30 MAF threshold, but the public repository does not preserve the original proxy-query script that operationalised that threshold. Record the exact rule used when preparing this flag; do not guess from the column name. If several eligible proxies are available for one target, the script selects the highest-r2 proxy with outcome data and records the alternatives in its audit file. This deterministic selection rule may differ from the original analysis if its query returned a different proxy.

Exact outcome matches take priority. The script reports exact-primary and exact-plus-proxy estimates separately, applies Benjamini–Hochberg FDR within each branch, and keeps all estimators returned by `TwoSampleMR::mr()`. Reused physical proxy variants within a taxon are flagged for review.

Run from the repository root:

```sh
Rscript --vanilla scripts/run_proxy_sensitivity.R
```

The script checkpoints one taxon at a time under `results/proxy_sensitivity/runs/`. It reuses valid taxon checkpoints, writes a completion manifest only after validation, and never overwrites a completed run. The outputs include primary estimates, all estimators, a proxy-selection audit, and a branch summary. The run folder is ignored by Git except for its README.

## Colocalisation inputs

Create `data/private/colocalisation/region_manifest.csv`, one row per taxon–locus–outcome comparison, with these columns:

`region_id,taxon,chromosome,lead_position,window_bp,trait1_file,trait2_file,n1,n2,sdY1,sdY2,genome_build,ld_rds`

`trait1_file` and `trait2_file` are paths relative to `data/private/colocalisation/`. Each regional summary-statistic file must contain dense, unfiltered rows for one trait and one region with columns:

`SNP,chr,pos,beta,se,effect_allele,other_allele,eaf,genome_build`

Use the same genome build and forward-strand allele coding in both files. The script retains shared variants with matching alleles, flips the second beta when effect and other alleles are reversed, and writes an audit of included and excluded variants. It stops if shared SNP positions or chromosomes disagree. It does not select variants based on P value; coloc needs dense regional coverage.

For quantitative traits, supply `n1` and `n2`. If `sdY1` or `sdY2` is unavailable, the corresponding file must provide EAF for every shared variant so the script can supply MAF and N to `coloc`. The `coloc` package requires beta and variance of beta; the script uses `se^2` for `varbeta`.

The script runs `coloc.abf` with p1=p2=1e-4 and p12 values 1e-6, 1e-5, and 1e-4. Set `window_bp` explicitly for each region. The manuscript used ±500 kb for most regions and ±1 Mb for two sparse CAG-390 regions.

To request optional multi-signal SuSiE colocalisation, set `ld_rds` to an RDS file containing a named, signed correlation matrix. Its row and column names must be the SNP IDs used in both summary-statistic files. The matrix must be aligned to the effect-allele coding in trait 1 after the allele harmonisation above, and the same reference panel should be appropriate for both traits. For the manuscript's rs635634 analysis this was a signed European-ancestry LD matrix. The script runs SuSiE once per region and reuses that fit for all p12 values.

Run from the repository root:

```sh
Rscript --vanilla scripts/run_colocalisation.R
```

The script checkpoints one region at a time under `results/colocalisation/runs/`. Outputs include ABF and optional SuSiE summaries, a variant alignment audit, and region-level input counts. Coloc's ABF method assumes at most one causal variant per trait in a region; SuSiE relaxes that assumption but depends on a well-aligned LD matrix. See the official `coloc` documentation for its [data requirements](https://chr1swallace.github.io/coloc/articles/a02_data.html) and [SuSiE workflow](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html).

## What this does and does not reproduce

These are new, explicit implementations of the methods; they are not recovered copies of the original unpublished scripts. The public repository does not contain the dense UKHLS/primary regional files, the original proxy map, or the LD matrix, so it cannot reproduce the manuscript's numerical results until those inputs and selection details are supplied. The current main pipeline also uses `ieu-b-35` for its CRP comparison, while the manuscript's external results use UKHLS `ebi-a-GCST005067`. Resolve that dataset difference before treating the pipeline as a full reproduction.

Use `Rscript --vanilla tests/test_analysis_extensions.R` for the cheap synthetic checks of proxy orientation and colocalisation allele harmonisation. These checks do not validate an empirical analysis.

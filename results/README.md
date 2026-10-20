# Results

The main pipeline writes generated tables and figures here when run with the study inputs. The proxy sensitivity and colocalisation scripts write versioned outputs under `proxy_sensitivity/runs/` and `colocalisation/runs/`. These generated files are not included in the repository.

The repository also contains a concise summary of findings reported in the author's manuscript revision:

- [`REPORTED_EMPIRICAL_SUMMARY.md`](REPORTED_EMPIRICAL_SUMMARY.md) gives the main findings and their limits.
- [`REPORTED_CANDIDATE_RESULTS.csv`](REPORTED_CANDIDATE_RESULTS.csv) gives the four exact-primary CRP candidates.

Those values were transcribed from the manuscript and are not outputs regenerated from the public code and data. The proxy and colocalisation workflows are new implementations that still need the original inputs and validation. See [`../docs/REPRODUCIBILITY_STATUS.md`](../docs/REPRODUCIBILITY_STATUS.md) and [`../scripts/ANALYSIS_EXTENSIONS.md`](../scripts/ANALYSIS_EXTENSIONS.md).

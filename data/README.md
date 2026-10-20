# Local input data

For analysis, place the prepared CSVs directly in **data/private/**. Keep them outside the prepared upload folder. A browser upload can include any file you add to the selected folder, even if it matches a .gitignore rule.

Required exposure file:

- **qin_IV_adaptive_combined.csv**

Required outcome files:

- **CRP_ebi-a-GCST90025959_selected_snps.csv** — CRP primary, GWAS Catalog ID **ebi-a-GCST90025959**
- **ieu-b-35_selected_snps.csv** — CRP replication
- **ebi-a-GCST90018926_selected_snps.csv** — T2D primary
- **finn-b-T2D_selected_snps.csv** — T2D replication
- **ieu-b-38_selected_snps.csv** — systolic blood pressure primary
- **ukb-b-20175_selected_snps.csv** — systolic blood pressure replication

Optional sensitivity file:

- **ieu-b-4817_selected_snps.csv** — within-family systolic blood pressure

The pipeline also searches the nested directory layouts used in the earlier analysis. The flat layout above is the simplest. Obtain data from their original sources and check each source's terms before redistribution.

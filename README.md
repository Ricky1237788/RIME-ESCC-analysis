# RIME_ESCC_code

R analysis scripts for the manuscript on plasma lncRNA RIME and PD-L1 in advanced ESCC treated with first-line PD-1 inhibitor plus chemotherapy.

## Important note

These scripts are a reproducible analysis template. Before public release, run every script on the locked analysis dataset and verify that all estimates, p values, sample sizes, and figures match the final manuscript.

## Suggested run order

1. `00_packages_and_setup.R`
2. `01_cutoff_maxstat_and_ROC.R`
3. `02_Figure1_PDL1.R` through `06_Figure5_joint_RIME_PDL1.R`
4. `07_Supplementary_dynamic.R` through `10_Supplementary_tables.R`

## Expected analysis dataset

The input file is not included because it contains individual-level clinical data. The scripts expect a de-identified CSV file with at least: `id`, `arm`, `rime`, `pdl1`, `OS_time`, `OS_event`, `PFS_time`, `PFS_event`, and `ORR`. The input `rime` value is the final baseline RIME expression value used in the manuscript. In the locked source workbook, this is the natural-log-transformed averaged baseline measurement: `rime = ln(poly_RNA_ave)` (stored as `logpoly_RNA_ave`).

`arm` should be coded as `Control` or `Experimental`; survival times are in months; `OS_event` and `PFS_event` use 1 for an event and 0 for censoring; and `ORR` uses 1 for objective response and 0 otherwise. The exploratory grouping cutoff is derived using maximally selected rank statistics based on PFS time and event status; the reported analysis cutoff is 2.23. Time-dependent ROC analysis at 6 months is used for descriptive assessment in Figure S1.

Dynamic analyses additionally require `data/longitudinal_rime.csv`, with one row per patient and the columns `id`, `rime_baseline`, `rime_post`, `PFS_time`, `PFS_event`, `OS_time`, and `OS_event`. Dynamic change is defined as `(ONGOING - Baseline) / Baseline`, with a prespecified exploratory cutoff of -16.2%.

For adjusted analyses, include the prespecified clinical covariates: age, sex, ECOG, disease status, lymph-node status, metastatic burden, prior chemoradiotherapy, and PD-L1 status.

## Data availability

Individual-level clinical data are not publicly available because of patient privacy and institutional restrictions. De-identified data may be available from the Lead Contact upon reasonable request and subject to institutional approval.

## Software

The R version and package versions used for this release are recorded in `sessionInfo.txt`.

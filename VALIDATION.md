# Completed local validation — 5 October 2026

The production secondary-data analysis was executed on the actual downloaded source workbooks. No synthetic inputs were used to produce the reported results.

| Check | Result |
|---|---|
| Source input identity | Three publisher XLSX files verified against SHA-256 values in `config.json` |
| Source record selection | Genuine POAG Cross-ancestry records; blank separator rows excluded |
| Source QC and posterior bounds | QC flags limited to 0/1; posterior values within 0–1 |
| Baseline count cross-check | 76 loci and 228 gene nominations agree with Supplementary Data 13 in Hamel et al. |
| Distinct method-agreement definition | The source gene-level intersection has 23 genes across 18 loci; this pilot's stricter same-QTL/same-tissue requirement gives 12 genes across 11 loci |
| Atlas identity | Unique text gene labels retained; 26 non-text labels in each workbook audited and excluded |
| Atlas values | Finite, nonnegative published averages; expected cluster labels present |
| Missing mappings | 162 nominated genes match the atlas; 66 remain explicitly unmatched |
| Exploratory tests | 10,000 seeded draws; nonzero empirical P values; one BH family of nine classes |
| CCA comparison | Seven shared classes in both versions; same donors; unavailable classes not imputed as zero |
| Runtime verification | Analysis completed, figures rendered and all four visually inspected |
| Unit tests | Nine passed |
| Cache behavior | Validated analysis reused when figure code was revised |
| GitHub CI | Workflow supplied; remote execution not yet completed |

Software used: Python 3.12.14; NumPy 2.3.5; pandas 2.2.3; SciPy 1.17.0; matplotlib 3.10.8; openpyxl 3.1.5.

Automated checks do not establish causal validity or independent disease replication. Biological conclusions remain restricted to secondary evidence and healthy-retina expression context.

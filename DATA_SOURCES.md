# Data sources and attribution

| Source | Files / sheets used | Role |
|---|---|---|
| Hamel AR et al. (2024), Nature Communications 15:396, [doi:10.1038/s41467-023-44380-y](https://doi.org/10.1038/s41467-023-44380-y) | Supplementary workbook, Data 7 and 10 | Published cross-ancestry POAG eCAVIAR and enloc e/sQTL result tables. Their authors fitted these models. |
| Lukowski SW et al. (2019), The EMBO Journal 38:e100811, [doi:10.15252/embj.2018100811](https://doi.org/10.15252/embj.2018100811) | Dataset EV1 and EV2 | Published cluster-average expression, before and after CCA, from 20,009 retinal cells from three healthy donors. Both files describe the same biological samples. |
| Gharahkhani P et al. (2021), Nature Communications 12:1258, [doi:10.1038/s41467-020-20851-4](https://doi.org/10.1038/s41467-020-20851-4) | Underlying GWAS cited through Hamel's tables | Original POAG GWAS: 34,179 cases and 349,321 controls across ancestries, 127 loci. No individual-level data or genome-wide summary statistics are analyzed in this pilot. |

Exact source download URLs and SHA-256 values are in `config.json`. Original workbooks are downloaded unchanged and cached locally. `data/raw/` is excluded from Git so that the GitHub project distributes code and derived results rather than rehosting third-party workbooks. A separate private offline data cache may accompany the deliverable.

The original publications and input workbooks retain their original authorship and applicable publisher terms; the code's MIT license does not relicense these sources. The Hamel publication is CC BY 4.0. Consult each original publication for its current terms before redistributing source files.

The derived tables are explicitly labelled secondary analyses. No cell-level identifiers, protected GTEx records, clinical records, patient identities, credentials, or individual-level genotypes are included.

The input Excel atlas contains non-text gene labels, including dates produced by Excel gene-name conversion. These are excluded without ambiguous name reconstruction. Every exclusion records the source workbook, Excel row and input label.

The 2024 source workbook contains thresholded colocalization result tables, not a complete record of every tested gene-tissue pair. An absent entry is treated as unavailable/unsupported in the published table, never as a fitted zero posterior.

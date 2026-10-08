# Count-level single-cell extension of the glaucoma pilot

Author portfolio: Mila Desi Anasanti. Analysis run: 6 October 2026. Repository: https://github.com/MDAnusamandiriGroup/glaucoma-cell-context

## What was actually executed

This module reanalyses the published filtered retinal UMI count matrix of Lukowski et al. (2019), rather than relying only on the earlier pilot's published cell-class averages. The earlier 228-gene colocalisation evidence table is frozen and used solely to define the six retina-supported targets. No GWAS association, colocalisation or Mendelian-randomisation model is refitted here. Source files were checked against the published MD5 checksums. All 20,009 cell barcodes, count totals and detected-gene counts match the published metadata exactly.

The input contains 20,009 cells, 21,989 genes, three healthy donors and five libraries. Primary QC retained 19,865 cells (99.3%); stricter QC retained 14,836. Total counts were scaled to 10,000 and log1p transformed, with integer UMIs retained separately. Library-aware selection produced 2,000 highly variable genes; 30 principal components, a 15-neighbour graph, Leiden clustering and UMAP were computed. Leiden resolutions 0.5 and 1.0 produced 14 and 15 clusters respectively. No integrated or batch-corrected embedding was fitted.

## Main results

All 6 nominated targets map by exact symbol to the count matrix; 163/228 genes from the complete earlier evidence table map exactly. Published original cell labels remain the primary labels for expression summaries. Independently fitted clusters were audited against canonical marker panels and published-label composition: 10/14 clusters at resolution 0.5 had the same provisional marker-panel label as their dominant published label. This audit is descriptive and does not replace expert annotation.

With the predefined minimum of 20 cells per donor and cell type, the cell types eligible in all three donors are: Rod, Cone, Bipolar, RGC, Muller glia. Rare cell types with inadequate coverage are reported as unavailable for the primary three-donor comparison. Their measured expression is retained in the full donor table for descriptive inspection; missing groups are never treated as zero expression.

The primary ranking uses the equal-donor mean of log1p pseudobulk CPM. UMIs are summed within donor and published cell type, divided by all-gene UMIs in that group, and transformed with log1p. Each donor has equal weight. The highest context below means the highest score among eligible cell types; it is not a causal target-cell assignment.

| Gene | Highest context (eligible types) | LODO folds agreeing | Strict-QC comparison |
|---|---|---:|---|
| COL8A2 | RGC | 2/3 | Not evaluable |
| DGKG | Muller glia | 3/3 | Same |
| NPC2 | Muller glia | 3/3 | Same |
| PSMC3 | RGC | 1/3 | Same |
| RP11-466F5.8 | Rod | 2/3 | Same |
| SLC2A12 | Muller glia | 3/3 | Same |

The highest context remains the same in all three leave-one-donor-out (LODO) reaggregations for 3/6 genes: DGKG, NPC2, SLC2A12. LODO removes one donor from the expression summary and reaggregates the two retained donor pseudobulks. It does not rerun QC, feature selection, clustering or UMAP. The set of eligible cell types remains fixed across folds.

Several target genes have sparse counts across the nine assigned cell types: COL8A2: 6 UMIs, DGKG: 31 UMIs, SLC2A12: 19 UMIs. These totals exclude Unassigned cells. COL8A2 has no detected UMIs among the eligible comparison types after strict QC, so its strict-QC ranking is not evaluable. LODO rank agreement for a low-count gene does not imply a precisely estimated expression level or strong biological evidence.

Strict-QC comparisons use only cell types eligible in all three donors under both QC rules: Rod, Cone, Bipolar, RGC, Muller glia. The highest context agrees for 5/6 genes: DGKG, NPC2, PSMC3, RP11-466F5.8, SLC2A12. The separately reported minimum-10-cell sensitivity includes: Rod, Cone, Bipolar, Horizontal, RGC, Muller glia, Microglia. These sensitivity outputs are retained even where they differ from the primary result; thresholds were not chosen to maximise agreement.

## Interpretation and limitations

This is a reproducible technical pilot using real single-cell counts and donor-level descriptive sensitivity analyses. It demonstrates count-matrix ingestion, QC, normalisation, feature selection, dimensionality reduction, graph clustering, marker auditing and provenance-aware reporting. It is not a new glaucoma case-control single-cell study, an ageing association study, or a completed spatial-transcriptomics analysis. The donors are healthy, donor number is three, age is confounded with donor, and cell composition is unequal. Cells are not independent biological replicates. No differential-expression p-values, disease/age estimates, spatial mapping, or new doublet classifier are produced.

The published source matrix was already filtered by its authors. This workflow starts from their filtered UMI matrix, not FASTQ files or unfiltered droplets. Published original labels are retained; provisional labels based on marker panels are only an audit aid. Library-aware feature selection does not remove donor or technical batch effects. Candidate expression and LODO stability do not establish causality, independent replication, therapeutic efficacy, or publication readiness for a Q1 journal. Changes from the earlier atlas-average ranking can reflect donor balancing, different estimands and exclusion of underrepresented classes.

## Reproduce and resume

Install Python 3.12 and the packages in `single_cell/requirements.txt`, then run `python single_cell/run_single_cell.py`. After source files have been cached, `python single_cell/run_single_cell.py --offline` validates and reuses completed stages. Stage fingerprints include inputs, relevant code, parameters and software versions; output SHA256 values are validated before reuse. Partial runs are preserved in separate attempts, existing results are never overwritten, and a process lock prevents concurrent writers. Frozen computed matrices are intentionally excluded from ordinary GitHub upload packages.

## Sources

Lukowski SW et al. A single-cell transcriptome atlas of the adult human retina. The EMBO Journal (2019). https://doi.org/10.15252/embj.2018100811

Original filtered UMI count matrix and metadata: https://zenodo.org/records/5515631

Scanpy official preprocessing and clustering tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html

The earlier candidate table and its original scientific sources remain documented in the parent repository.

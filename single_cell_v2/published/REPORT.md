# Glaucoma cell-context V2: executed Python pilot

V2 expands the frozen published-genetics pilot to all 228 nominations, technical sampling robustness, independent-study RNA context, predefined metabolic programs and regulatory triangulation. It uses 19,865 V1 whole-cell profiles from three healthy donors and 51,645 nuclear profiles from four different published donors. The two studies are analysed separately.

## Executed results

- 163/228 candidates mapped exactly in V1; 156/228 in the original Wang count distribution.
- 94/163 V1 candidates retained the same unique top type in all three leave-one-donor-out folds; 150 retained their top type under strict QC.
- 2,400 technical draws: 300 per scenario across eight frozen scenarios. Sampling stays within donor/type; comparison types are fixed within each scenario.
- 78/139 jointly evaluable candidates had the same highest-expression type in the fixed five-type comparison. With seven types, 70/140 were concordant.
- 16 genes passed the descriptive context screening rule. This is not significance or a causal probability.
- Six predefined programs were scored with expression-matched controls. Only 3/6 had the same relative top type across studies when using identical shared panel/control pairs.
- 25 index variants overlapped reproducible Wang ATAC peaks; 50 were represented in the published glaucoma annotation table; 14 candidates had a source-reported gene link. Computed cell-type peak intersections reproduced all 50 published annotations exactly.
- Four index variants overlapped the recent Yuan 2026 source-listed retina/macula cCRE table. All four overlaps survived plausible one-base endpoint conventions. Microglia were excluded from this supplied atlas table because of low cell numbers.

## Six prespecified cases

| Gene | Five types: V1 / Wang | Seven types: V1 / Wang | V1 assigned-cell UMIs |
|---|---|---|---:|
| COL8A2 | RGC / Muller glia | Microglia / Microglia | 6 |
| DGKG | Muller glia / RGC | Muller glia / RGC | 31 |
| NPC2 | Muller glia / Muller glia | Microglia / Microglia | 3642 |
| PSMC3 | RGC / RGC | Horizontal / Horizontal | 3662 |
| RP11-466F5.8 | Rod / nan | Rod / nan | 103 |
| SLC2A12 | Muller glia / Muller glia | Muller glia / Muller glia | 19 |

NPC2 is technically stable in both comparison sets and concordant across studies. Its five-type maximum is Müller glia; its seven-type maximum is microglia. The same index variant, rs754458, overlaps a microglia ATAC peak. However, the published retinal eQTL table names NPC2 and LTBP2, the HiChIP table links YLPM1, and BPNet does not call this variant high-effect. This is a locus with reproducible RNA context and unresolved target-gene mechanism. It is not a new discovery of NPC2, a microglia-specific eQTL, or a proven causal NPC2 effect.

DGKG changes context across studies. SLC2A12 agrees across RNA studies but has low V1 counts and less than 90% agreement under equal-cell sampling. PSMC3 agrees across studies but is unstable under V1 equal-cell sampling. COL8A2 is especially sparse. RP11-466F5.8 is not exactly symbol-mapped in the external distribution; no alias was inferred.

## Methods and limits

Raw UMIs are summed within donor and published cell type, divided by all-gene UMIs in that group, converted to log1p CPM and averaged with equal donor weights. Primary types are Rod, Cone, Bipolar, RGC and Müller glia (at least 20 cells per donor/type). The extended sensitivity adds Horizontal and Microglia (at least 10 cells per V1 donor/type). Donor omissions reaggregate saved counts; they do not refit clustering.

Technical draws use cell sampling without replacement or binomial thinning of selected-gene and remaining UMIs together. They describe uncertainty conditional on the observed samples, not biological population confidence intervals. The 90% screening cutoff is descriptive.

Six MSigDB Hallmark 2025.1.Hs programs were selected before module execution. Control genes came from 20 expression bins, with up to 20 controls per panel gene and all six panels excluded from controls. Cross-study scoring uses identical intersected panel/control pairs and weights. The scores describe RNA programs, not metabolic flux, mitochondrial function or disease-specific activation. Whole-cell versus nuclear capture, donor characteristics and laboratory effects can contribute to differences.

Regulatory evidence uses exact index rsIDs and GRCh38 coordinates. Nearest-gene marker peaks, target predictions, promoter co-accessibility, HiChIP and retinal eQTL are separate evidence categories. Published BPNet scores and negative calls are imported; no model was trained. Hamel and Wang reuse overlapping retinal eQTL resources, so this is not independent genetic replication. No new GWAS, fine-mapping, colocalization or Mendelian randomization has been fitted.

Yuan's separate peak/gene matrices were not paired by row order without explicit identifiers. A validated enhancer-gene edge list is required for that additional ABC analysis. New colocalization/MR would require complete harmonized GWAS/QTL summary statistics; thresholded published nomination tables are insufficient.

## Reproducibility and disclosure

Python stages are checkpointed with source, software and code hashes. V1 models/checkpoints are read-only inputs. Failed attempts remain diagnostic and are excluded from handover. Code implementation and reporting were AI-assisted; executed outputs and numerical checks are supplied. This is a computational pilot for review, not a manuscript accepted for publication or an experimentally validated mechanism.

## Sources

- Lukowski et al. 2019, DOI https://doi.org/10.15252/embj.2018100811; counts https://zenodo.org/records/5515631
- Hamel et al. 2024, DOI https://doi.org/10.1038/s41467-023-44380-y
- Wang et al. 2022, DOI https://doi.org/10.1016/j.xgen.2022.100164; GSE196235 and CellxGene collection https://cellxgene.cziscience.com/collections/348da6dc-5bf6-435d-adc5-37747b9ae38a
- Yuan et al. 2026, DOI https://doi.org/10.1126/sciadv.adv9162
- MSigDB Hallmark 2025.1.Hs: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
- MitoCarta3.0: https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
- GENCODE v32 GRCh38: https://www.gencodegenes.org/human/release_32.html

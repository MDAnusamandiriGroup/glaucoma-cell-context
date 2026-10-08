# Sources and provenance

Retrieved 6 October 2026. Full URLs, sizes, local SHA-256 hashes and publisher checksums are in [data/source_registry.json](data/source_registry.json). Supporting method references and the inspected supplementary-methods file hash are in [data/method_sources.json](data/method_sources.json).

## EyeGEx retinal eQTL

Ratnapriya et al. 2019, *Retinal transcriptome and eQTL analyses identify genes associated with age-related macular degeneration*.

- Original publication: https://doi.org/10.1038/s41588-019-0351-9
- Author replacement archive: https://doi.org/10.5281/zenodo.18331341 (CC BY 4.0)
- Full nominal chromosome 14 data, 113,356,961 compressed bytes: https://zenodo.org/records/18331341/files/new_nominal.eQTL.chr14.with_thresholds.txt.gz?download=1
- Full gene permutation statistics: https://zenodo.org/records/18331341/files/new_perms_full.txt?download=1
- Author column descriptions: https://zenodo.org/records/18331341/files/README_Retina_eQTL_files_021125.txt?download=1
- Inspected supplementary methods: https://pmc-oa-opendata.s3.amazonaws.com/PMC6441365.1/NIHMS1518564-supplement-1.docx

The archive description mentions significant pairs and SE, but the actual nominal file contains full gene rows, including large P values, and no native SE column. Actual columns and count checks control this analysis. Coordinates in the explicit chromosome/position columns are GRCh38; embedded variant-ID positions can be GRCh37. `V12` in the nominal file is P and `V13` is beta. In the permutation file, `V12` is nominal dof1 and `V13` is optimized dof2. These schemas are different.

## Primary open-angle glaucoma GWAS

Gharahkhani et al. 2021, *Genome-wide meta-analysis identifies 127 open-angle glaucoma loci with consistent effect across ancestries*: https://doi.org/10.1038/s41467-020-20851-4

Primary European study GCST90011766: 16,677 cases and 199,580 controls.

https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011766/harmonised/GCST90011766.h.tsv.gz

Combined-ancestry sensitivity GCST90011770: 34,179 cases and 349,321 controls; overlapping with the primary analysis.

https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/GCST90011770/harmonised/GCST90011770.h.tsv.gz

The harmonized files use GRCh38. Regional extraction used the remotely indexed files and locally verified index MD5 checksums. Full GWAS files were not downloaded or locally hash-verified; published full-file checksums are retained. The regional TSVs used in the analysis have local SHA-256 hashes.

## Method and previous nomination

- Giambartolomei et al. 2014 ABF colocalization method: https://doi.org/10.1371/journal.pgen.1004383
- Official coloc input and prior documentation: https://chr1swallace.github.io/coloc/articles/a02_data.html
- Official QTLtools cis nominal/permutation documentation: https://qtltools.github.io/qtltools/pages/QTLtools-cis.1.html
- Earlier glaucoma gene nomination using these source cohorts, Hamel et al. 2024: https://doi.org/10.1038/s41467-023-44380-y

The new targeted fit does not constitute independent genetic replication. The healthy single-cell/single-nucleus expression component of the broader pilot is separate from bulk-retina eQTL mapping.

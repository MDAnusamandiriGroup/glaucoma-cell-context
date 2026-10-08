# Validation of the count-level extension

Executed on 6 October 2026 with Python 3.12.14 and Scanpy 1.11.5. The complete recorded scientific environment is in `published/provenance/run_index.json`; package pins are in `requirements.txt`.

## Production data checks

- All three source files matched published MD5 values.
- Matrix dimensions and unique symbols/barcodes were checked.
- Count-matrix barcodes match both metadata files exactly.
- UMI totals and detected genes match the published per-cell QC metadata for every one of the 20,009 cells.
- Normalised expression, PCA and UMAP contain finite values; cell order is identical between cluster and expression outputs.
- Target symbol mapping, donor counts, missing groups, eligibility and low-count expression are reported explicitly.
- All seven stage manifests and 46 completed output files passed SHA256 validation.
- A full offline cache-hit run reused all seven stages, made no analytical recomputations and left every completed output modification time unchanged.
- The process lock was released after that run.
- The earlier pilot's frozen analysis manifest was unchanged: SHA256 `406d81fbac36cf7a67ec2aee546ced6ce2f511f4f4d75fe107ed0745c3225380`.

## Nine focused fixture tests passed locally

1. Compatible cache reuse does not recompute or rewrite outputs.
2. Corrupt checkpoints fail closed and are preserved.
3. Partial stages are preserved and a separate attempt is used.
4. Existing published outputs cannot be overwritten.
5. The gzip/TAR count-matrix container and missing leading header field are read correctly.
6. Duplicate gene symbols are rejected.
7. Metadata join uses exact barcodes, preserves requested order and rejects mismatches/duplicates.
8. Missing cell groups and absent genes stay missing rather than becoming zero; donor pseudobulk CPM and strict-QC eligibility are checked on a fixture.
9. Donors have equal weight; LODO uses a fixed comparison set; the minimum-10-cell sensitivity is distinct from the primary comparison.

Every new Python module passed syntax parsing. Production models were not rerun for fixture tests. Results were rendered from saved tables and embeddings. Figures and document layouts were visually inspected. These checks establish technical reproducibility of this run; they are not independent biological validation or evidence that all scientific assumptions are correct. No GitHub CI success is claimed.

The full validation JSON and local test output are retained under `published/provenance/`.

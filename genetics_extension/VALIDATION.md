# Completed verification

Eight unit tests passed on 6 October 2026 with Python 3.12.14, NumPy 2.3.5, pandas 2.2.3 and SciPy 1.17.0.

- Recovery of known regression SE from Student-t P values, including a small-tail P value.
- Agreement of H0–H4 with independent explicit enumeration of all distinct-variant pairs in a small fixture.
- Invariance of unsigned ABFs to effect sign and of posteriors to variant order.
- Finite, normalized outputs with extreme shared and distinct signals.
- Rejection of invalid SE, P and region inputs.
- Cache reuse without calling the producer.
- Preservation and rejection of tampered output.
- Preservation and rejection of an interrupted stage.

The separate `validate_results.py` passed on the real six comparisons / 216 sensitivity settings. It checks processed-input and completed-output hashes, probability sums, variant keys, df/SE labels, effect orientation, full nominal gene counts, conditional H4 sets and Python syntax. An independent long-double calculation of the actual six primary posteriors agrees to a maximum absolute error of **1.887379141862766e−15**.

A second entry-point call reported `REUSE` for preprocessing and fitting, with unchanged manifests. The log is [published/resume_verification.txt](published/resume_verification.txt).

The four-page PDF was rendered and visually checked. Numerical and cache verification does not establish biological validity or validate approximate reconstructed SE against unavailable native SE. The R `coloc` package was not executed and no empirical R parity check was performed. The interpretation requires explicit single-causal-variant and prior-scale assumptions.

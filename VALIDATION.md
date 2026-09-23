# Validation

Tested on 23 September 2026 with Python 3.12 and the pinned requirements.

- Full pipeline ran in a fresh output directory using only the 12 specified
  H5AD files and the included archived statistics.
- All 31 generated PNG files matched the existing manuscript exports
  byte-for-byte. PDF exports were also generated for all 31 figures.
- All 512 local composition rows, 1,104 ACED abundance rows, and 24 local
  composition p/q comparisons matched the original results numerically.
- Whole-tissue and patient-level test tables matched the original results.
- The 200 region-level window summaries, 16 cohort-level window tests,
  and all 136 CN OLS results matched the original tables numerically.
- Six unit tests passed for stored-neighborhood occurrence counting,
  zero denominators, missing grade labels, and patient mapping.

Repeat the image check after a full run:

```bash
python verify_figures.py
```

PNG checksums require identical rendering; the statistical matching above
was checked independently against the original numerical tables.

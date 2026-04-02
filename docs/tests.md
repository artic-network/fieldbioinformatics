---
title: tests
summary: Description of the available tests
authors:
  - Will Rowe
  - Nick Loman
  - Sam Wilkinson
date: 2024-11-11
---

# The test suite

All tests run as part of CI via GitHub Actions. To run them locally, clone the repository and activate the conda environment:

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
conda env create -f environment.yml
conda activate artic
```

## Unit tests

Unit tests cover the core Python modules. Tests are located in `tests/` and are discovered automatically via the `pyproject.toml` configuration. To run all unit tests:

```sh
pytest
```

To exclude slow tests (see below):

```sh
pytest -m "not slow"
```

The unit tests cover:

- `align_trim` — amplicon soft clipping and alignment filtering
- `artic_vcf_filter` — variant quality and allele frequency filtering
- `artic_mask` — coverage masking
- `artic_get_models` — model download logic
- `utils` — scheme fetching, model selection, primer direction, and related utilities
- `pipeline` — CLI argument parsing and pipeline exit code handling

## Pipeline tests

To test the full pipeline end-to-end, use the `test-runner.sh` script:

```sh
./test-runner.sh clair3
```

This downloads a small subset of an Ebola virus amplicon sequencing run (flongle) from [artic.s3.climb.ac.uk](http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz) and runs the complete pipeline against it.

## Variant validation tests

The validation tests download reference SARS-CoV-2 datasets, run the full workflow, and validate the reported variants and consensus sequences against known truth sets. These tests require internet access and are marked `slow` in the pytest configuration, so they are excluded from the default `pytest` run.

To run them explicitly:

```sh
pytest -m slow
```

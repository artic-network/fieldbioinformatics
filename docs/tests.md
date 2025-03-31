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

All of the `artic` tests are run as part of the as github actions. To run any of the tests yourself, it is assumed that you have downloaded the codebase and are using an appropiate environment:

```
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
conda env create -f environment.yml
conda activate artic
```

## Unit tests

We have begun writing unit tests for the ARTIC Python modules. These currently includes tests for the `align_trim` module, which performs amplicon soft clipping and alignment filtering, and `vcftagprimersites` which processes the ARTIC primer scheme. To run all available unit tests:

```
pytest -s artic/*_unit_test.py
```

## Pipeline tests

To test the core pipeline, you can use the `test-runner.sh` bash script:

```
./test-runner.sh clair3
```

The pipeline tests use a small subset of an Ebola virus amplicon sequencing run (flongle) which is downloaded from [here](http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz) when the test is called.

## Variant validation tests

Finally, we have also included some validation tests that will download several reference SARS-CoV-2 datasets, run the workflow, and then validate the reported variants and consensus sequences. To run all of the available validation datasets:

```
pytest -s artic/minion_validator.py
```
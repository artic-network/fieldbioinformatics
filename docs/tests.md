---
title: tests
summary: Description of the available tests
authors:
  - Will Rowe
  - Nick Loman
date: 2020-03-30
---

# The test suite

All of the `artic` tests are run as part of the [Travis Continuous Integration](https://travis-ci.org/github/artic-network/fieldbioinformatics) testing. To run any of the tests yourself, it is assumed that you have downloaded the codebase and are using an appropiate environment:

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

To test the core pipeline, you can use the `test-runner.sh` bash script. You can test both the medaka and nanopolish workflows:

```
./test-runner.sh medaka
./test-runner.sh nanopolish
```

The pipeline tests use a small subset of an Ebola virus amplicon sequencing run (flongle) which is downloaded from [here](http://artic.s3.climb.ac.uk/run-folders/EBOV_Amplicons_flongle.tar.gz) when the test is called.

## Variant validation tests

Finally, we have also included some validation tests that will download several reference SARS-CoV-2 datasets, run the nanopolish/medaka workflows and then validate the reported variants and consensus sequences. To run all of the available validation datasets:

```
pytest -s artic/minion_validator.py
```

Or you can specify which workflow and how many datasets to validate against:

```
pytest -s artic/minion_validator.py --workFlow medaka --numValidations 2
```

> use --workFlows to specify workflow (medaka|nanopolish)

> use --numValidations to specify how many datasets to download and validate against (specifing -1 or too high a number will just run all the datasets)

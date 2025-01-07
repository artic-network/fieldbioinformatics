<div align="center">
    <img src="docs/artic-logo.png?raw=true?" alt="artic-logo" width="250">
    <h1>ARTIC</h1>
    <h3>a bioinformatics pipeline for working with virus sequencing data sequenced with nanopore</h3>
    <hr>
    <a href="https://github.com/artic-network/fieldbioinformatics/actions/workflows/unittests.yml"><img src="https://github.com/artic-network/fieldbioinformatics/actions/workflows/unittests.yml/badge.svg" alt="Unit tests"></a>
    <a href='http://artic.readthedocs.io/en/latest/?badge=latest'><img src='https://readthedocs.org/projects/artic/badge/?version=latest' alt='Documentation Status'></a>
    <a href="https://bioconda.github.io/recipes/artic/README.html"><img src="https://anaconda.org/bioconda/artic/badges/downloads.svg" alt="bioconda"></a>
    <a href="https://github.com/artic-network/fieldbioinformatics/blob/master/LICENSE"><img src="https://img.shields.io/badge/license-MIT-orange.svg" alt="License"></a>
    <a href="https://github.com/artic-network/fieldbioinformatics/actions/workflows/docker-build-push.yml"><img src="https://github.com/artic-network/fieldbioinformatics/actions/workflows/docker-build-push.yml/badge.svg" alt="Docker Image"></a>
</div>

---

## Overview

`artic` is a pipeline and set of accompanying tools for working with viral nanopore sequencing data, generated from tiling amplicon schemes.

It is designed to help run the artic bioinformatics protocols; for example the [MPXV protocol](https://artic.network/mpxv/mpxv-bioinformatics-sop.html).

Features include:

- read filtering
- primer trimming
- amplicon coverage normalisation
- variant calling
- consensus building

## Installation

### Via conda

```sh
conda install -c bioconda -c conda-forge artic
```

Please make sure you are using either mamba or conda version >= 23.10.0 where libmamba solver was made the default conda solver.

### Via source

#### 1. downloading the source:

<!-- Download a [release](https://github.com/artic-network/fieldbioinformatics/releases) or use the latest master (which tracks the current release): -->

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
```

#### 2. installing dependencies:

The `artic pipeline` has several [software dependencies](https://github.com/artic-network/fieldbioinformatics/blob/master/environment.yml). You can solve these dependencies using the minimal conda environment we have provided, Please make sure you are using either mamba or conda version >= 23.10.0 where libmamba solver was made the default conda solver.
:

```sh
conda env create -f environment.yml
conda activate artic
```

#### 3. installing the pipeline:

```sh
pip install .
```

#### 4. testing the pipeline:

First check the pipeline can be called.

```
artic -v
```

You can try the pipeline tests.

```
./test-runner.sh clair3
```

For further tests, such as the variant validation tests, check [the documentation](http://artic.readthedocs.io/en/latest/tests?badge=latest).

## Documentation

Documentation for the `artic pipeline` is available via [read the docs](http://artic.readthedocs.io/en/latest/?badge=latest).

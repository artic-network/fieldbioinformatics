---
title: installation
summary: The installation guide.
authors:
  - Will Rowe
  - Nick Loman
  - Sam Wilkinson
date: 2020-03-30
---

# Installation

As of [release 1.4.0](https://github.com/artic-network/fieldbioinformatics/releases/tag/1.4.0), we provide a docker image [available here](https://quay.io/repository/artic/fieldbioinformatics) and a conda package. You may also wish to install the package from source after installing the dependencies via Conda yourself.

## Via conda

```sh
conda install -c bioconda artic
```

## Via source

### 1. installing dependencies

The `artic pipeline` has several [software dependencies](https://github.com/artic-network/fieldbioinformatics/blob/master/environment.yml). You can solve these dependencies using the minimal conda environment we have provided, we strongly recommend that you use either the mamba solver or conda version >= 23.10.0 since libmamba solver is now the default conda solver:

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
conda env create -f environment.yml
```

### 2. installing the pipeline

```sh
conda activate artic
pip install .
```

### 3. testing the pipeline

First check the pipeline can be called:

```
artic -v
```

To check that you have all the required dependencies, you can try the pipeline tests like so:

```
./test-runner.sh clair3
```

For further tests, such as the variant validation tests, see [here](http://artic.readthedocs.io/en/latest/tests?badge=latest).
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

The pipeline is available as a Docker image, a conda package, or can be installed from source after resolving dependencies via conda.

## Via Docker

A Docker image is available at [quay.io/artic/fieldbioinformatics](https://quay.io/repository/artic/fieldbioinformatics). All Clair3 models are pre-bundled — `artic_get_models` does not need to be run separately.

```sh
docker pull quay.io/artic/fieldbioinformatics

docker run --rm \
  -v $(pwd):/data \
  -w /data \
  quay.io/artic/fieldbioinformatics:latest \
  artic minion \
    --scheme-name artic-inrb-mpox \
    --scheme-version v1.0.0 \
    --scheme-length 2500 \
    --read-file my_sample.fastq \
    my_sample
```

## Via conda

```sh
conda create -n artic -c bioconda artic
conda activate artic
```

After installation, download the Clair3 models:

```sh
artic_get_models
```

## Via source

### 1. Install dependencies

The pipeline has several [software dependencies](https://github.com/artic-network/fieldbioinformatics/blob/master/environment.yml). Use the provided conda environment file to resolve them. We strongly recommend conda >= 23.10.0 (which uses the libmamba solver by default):

```sh
git clone https://github.com/artic-network/fieldbioinformatics
cd fieldbioinformatics
conda env create -f environment.yml
```

### 2. Install the pipeline

```sh
conda activate artic
pip install .
```

### 3. Download Clair3 models

```sh
artic_get_models
```

### 4. Verify the installation

Check the pipeline can be called:

```sh
artic -v
```

To verify all dependencies are present, run the pipeline integration tests:

```sh
./test-runner.sh clair3
```

For further tests see [Tests](./tests.md).

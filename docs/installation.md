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

Images are published to [quay.io/artic/fieldbioinformatics](https://quay.io/repository/artic/fieldbioinformatics) for both `amd64` and `aarch64` in two variants:

| Tag | Clair3 models | Use when |
| --- | --- | --- |
| `latest` / `vX.Y.Z` | Not included | You will mount a pre-downloaded model directory |
| `latest-models-included` / `vX.Y.Z-models-included` | Bundled | You want a fully self-contained image |

### With models pre-bundled

```sh
docker pull quay.io/artic/fieldbioinformatics:latest-models-included

docker run --rm \
  -v $(pwd):/data \
  -w /data \
  quay.io/artic/fieldbioinformatics:latest-models-included \
  artic minion \
    --scheme-name artic-inrb-mpox \
    --scheme-version v1.0.0 \
    --scheme-length 2500 \
    --read-file my_sample.fastq \
    my_sample
```

### Without models (bring your own)

Download models first (see [Clair3 Models](./clair3-models.md)), then mount the directory at run time:

```sh
docker pull quay.io/artic/fieldbioinformatics

docker run --rm \
  -v $(pwd):/data \
  -v /path/to/models:/models \
  -w /data \
  quay.io/artic/fieldbioinformatics:latest \
  artic minion \
    --model-dir /models \
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

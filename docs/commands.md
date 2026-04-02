---
title: commands
summary: The available artic pipeline commands.
authors:
  - Sam Wilkinson
  - Will Rowe
  - Nick Loman
date: 2024-11-11
---

# Commands

This page documents the available commands via the `artic` command line interface.

## guppyplex

### Overview

Aggregate pre-demultiplexed reads from a directory (or directories) produced by MinKNOW into a single FASTQ file, optionally filtering by length and quality.

### Input

- directory or directories containing FASTQ files to aggregate

### Output

- a single aggregated FASTQ file (plain or gzip-compressed — determined by the output filename)

### Usage example

```bash
artic guppyplex --directory ./fastq_pass --prefix my_sample

# Write gzip-compressed output directly:
artic guppyplex --directory ./fastq_pass --output my_sample.fastq.gz
```

| Argument name(s)       | Required | Default value | Description                                                                                   |
| :--------------------- | :------- | :------------ | :-------------------------------------------------------------------------------------------- |
| `--directory`          | Y        | NA            | Directory containing FASTQ files to aggregate                                                 |
| `--prefix`             | N        | NA            | Prefix for the auto-generated output filename (`<prefix>_<dir>.fastq`)                        |
| `--output`             | N        | NA            | Write aggregated reads to this path (alternative to `--prefix`). If the path ends in `.gz` the output is gzip-compressed |
| `--max-length`         | N        | NA            | Remove reads greater than this length (bp)                                                    |
| `--min-length`         | N        | NA            | Remove reads less than this length (bp)                                                       |
| `--quality`            | N        | 7             | Remove reads below this mean quality score                                                    |
| `--sample`             | N        | 1             | Sampling fraction for random subsampling (e.g. `0.5` keeps ~50 % of reads)                   |
| `--skip-quality-check` | N        | False         | Do not filter on quality score (speeds up processing)                                         |
| `--threads`            | N        | 1             | Number of worker processes for parallel file processing                                       |

---

## minion

### Overview

Run the full alignment / variant-calling / consensus pipeline on a sample. See the [Core Pipeline](./minion.md) page for a detailed description of each stage.

### Input

- a FASTQ read file and a primer scheme (either provided directly via `--bed`/`--ref`, or fetched automatically via `--scheme-name`/`--scheme-version`)

### Output

- primer-trimmed alignments, variant calls, and a consensus sequence

### Usage example

```bash
artic minion \
  --scheme-name artic-inrb-mpox \
  --scheme-version v1.0.0 \
  --scheme-length 2500 \
  --read-file my_sample.fastq \
  my_sample
```

| Argument name(s)           | Required | Default value             | Description                                                                                                                                                                                       |
| :------------------------- | :------- | :------------------------ | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`                   | Y        | NA                        | The name of the sample; used as the prefix for all output files                                                                                                                                   |
| `--read-file`              | Y        | NA                        | Path to the input FASTQ (or FASTA) file                                                                                                                                                           |
| `--scheme-name`            | N        | NA                        | Name of the scheme to fetch from the primerschemes repository (e.g. `artic-inrb-mpox`, `sars-cov-2`)                                                                                             |
| `--scheme-version`         | N        | NA                        | Version of the scheme (e.g. `v1.0.0`). For schemes with automatic reference selection the correct suffix will be appended automatically; supply the full suffixed version to skip auto-selection  |
| `--scheme-length`          | N        | NA                        | Amplicon length of the scheme (required when a scheme has multiple lengths)                                                                                                                       |
| `--scheme-directory`       | N        | `./primer-schemes`        | Local directory to cache downloaded schemes                                                                                                                                                       |
| `--bed`                    | N        | NA                        | Path to a local primer scheme BED file (alternative to `--scheme-name`)                                                                                                                           |
| `--ref`                    | N        | NA                        | Path to a local reference FASTA file (alternative to `--scheme-name`)                                                                                                                             |
| `--model`                  | N        | NA                        | Clair3 model to use. If not provided the pipeline reads the `basecall_model_version_id` tag from the first read header and selects a model automatically. See [Clair3 Models](./clair3-models.md) |
| `--model-dir`              | N        | `$CONDA_PREFIX/bin/models/` | Directory containing Clair3 model files                                                                                                                                                         |
| `--min-depth`              | N        | 20                        | Minimum read depth required for a position to be included in the consensus (positions below this threshold are masked to `N`)                                                                     |
| `--normalise`              | N        | 100                       | Downsample each amplicon to this mean depth before variant calling. Set to `0` to disable normalisation                                                                                           |
| `--min-mapq`               | N        | 20                        | Discard reads with a mapping quality below this value                                                                                                                                             |
| `--primer-match-threshold` | N        | 35                        | Maximum allowed distance (bp) between a read end and a primer site for amplicon assignment                                                                                                        |
| `--allow-mismatched-primers` | N      | False                     | Retain reads whose primer pairs do not match (e.g. amplicon read-through). By default such reads are removed                                                                                      |
| `--threads`                | N        | 8                         | Number of threads for steps that support parallelism                                                                                                                                              |
| `--no-indels`              | N        | False                     | Exclude all insertions and deletions from variant calling                                                                                                                                         |
| `--no-frameshifts`         | N        | False                     | Exclude frameshift indels (lengths not divisible by 3) from the consensus                                                                                                                         |
| `--align-consensus`        | N        | False                     | After consensus generation, produce a MAFFT alignment of the consensus against the reference (`$SAMPLE.aligned.fasta`)                                                                            |
| `--linearise-fasta`        | N        | False                     | Write consensus FASTA with the sequence on a single line (no line wrapping)                                                                                                                       |
| `--dry-run`                | N        | False                     | Print all commands to the log without executing them                                                                                                                                              |

---

## artic filter

### Overview

Filter a FASTQ file by read length, writing passing reads to stdout.

### Usage example

```bash
artic filter --min-length 400 --max-length 700 reads.fastq > filtered.fastq
```

| Argument name(s) | Required | Default value | Description                          |
| :--------------- | :------- | :------------ | :----------------------------------- |
| `filename`       | Y        | NA            | Input FASTQ file                     |
| `--min-length`   | N        | NA            | Discard reads shorter than this (bp) |
| `--max-length`   | N        | NA            | Discard reads longer than this (bp)  |

---

## artic_get_models

### Overview

Download the Clair3 models required by the pipeline. Under a standard conda installation models are stored in `$CONDA_PREFIX/bin/models/`. Models are distributed as PyTorch checkpoints (`pileup.pt` and `full_alignment.pt`) from the [HKU Clair3 model repository](https://www.bio8.cs.hku.hk/clair3/).

!!! note
    The Docker image bundles all models — `artic_get_models` only needs to be run for conda or source installs.

### Usage example

```bash
artic_get_models
# or, for a non-conda install:
artic_get_models --model-dir /path/to/models
```

| Argument name(s) | Required | Default value             | Description                                                         |
| :--------------- | :------- | :------------------------ | :------------------------------------------------------------------ |
| `--model-dir`    | N        | `$CONDA_PREFIX/bin/models/` | Directory to download models into                                 |

For the full list of available models and how they are selected automatically at runtime, see [Clair3 Models](./clair3-models.md).

---

## artic_get_scheme

### Overview

Download a primer scheme independently of running the full pipeline. Writes `primer.bed` and `reference.fasta` to the current directory (or to `--scheme-directory` if provided).

### Usage example

```bash
artic_get_scheme \
  --scheme-name artic-inrb-mpox \
  --scheme-version v1.0.0 \
  --scheme-length 2500
```

| Argument name(s)     | Required | Default value      | Description                                                              |
| :------------------- | :------- | :----------------- | :----------------------------------------------------------------------- |
| `--scheme-name`      | Y        | NA                 | Name of the scheme (e.g. `artic-inrb-mpox`, `sars-cov-2`)               |
| `--scheme-version`   | Y        | NA                 | Version string (e.g. `v1.0.0`)                                           |
| `--scheme-length`    | N        | NA                 | Amplicon length (required when a scheme has multiple lengths)            |
| `--scheme-directory` | N        | `./primer-schemes` | Directory to write the downloaded scheme files into                      |
| `--read-file`        | N        | NA                 | FASTQ file used for automatic reference selection (where supported)      |

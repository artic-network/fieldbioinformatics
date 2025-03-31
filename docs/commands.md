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

Aggregate pre-demultiplexed reads from MinKNOW/Guppy

### Input

- director(y/ies) to aggregate from

### Output

- directory of aggregated files

### Usage example

```bash
artic guppyplex --directory ./
```

| Argument name(s)     | Required | Default value | Description                                                       |
| :------------------- | :------- | :------------ | :---------------------------------------------------------------- |
| --directory          | Y        | NA            | The director[y/ies] to gather files from                          |
| prefix               | Y        | NA            | Prefix for guppyplex files                                        |
| --max-length         | N        | NA            | Remove reads greater than max-length                              |
| --min-length         | N        | NA            | Remove reads less than than min-length                            |
| --quality quality    | N        | 7             | Remove reads against this quality filter                          |
| --sample sample      | N        | 1             | Sampling frequency for random sample of sequence to reduce excess |
| --skip-quality-check | N        | NA            | Do not filter on quality score (speeds up)                        |

---

## minion

### Overview

Run the alignment/variant-call/consensus pipeline

### Input

- a primer scheme and a sample directory

### Output

- trimmed alignments, variants calls and consensus sequence

### Usage example

```bash
artic minion <scheme> <sample>
```

| Argument name(s)         | Required | Default value             | Description                                                                                                                                                                              |
| :----------------------- | :------- | :------------------------ | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| sample                   | Y        | NA                        | The name of the sample                                                                                                                                                                   |
| --read-file              | Y        | NA                        | Use alternative FASTA/FASTQ file to <sample>.fasta                                                                                                                                       |
| --model-dir              | Y        | $CONDA_PREFIX/bin/models/ | Path containing clair3 models, defaults to models packaged with conda installation                                                                                                       |
| --primer-match-threshold | Y        | 35                        | Allow fuzzy primer position matching within this threshold                                                                                                                               |
| --min-depth              | Y        | 20                        | Minimum coverage required for a position to be included in the consensus sequence                                                                                                        |
| --model                  | N        | NA                        | Clair3 model to use, if not provided the pipeline will attempt to determine the appropriate model based on the `basecall_model_version_id` field in the input FASTQ header (recommended) |
| --normalise              | N        | 100                       | Reduce reference coverage to this mean depth per amplicon, deactivate with `--normalise 0`                                                                                               |
| --threads                | N        | 8                         | Number of threads to utilise during steps which support multiprocessing                                                                                                                  |
| --scheme-directory       | N        | ./primer-schemes          | Download schemes from the primerschemes repository (https://github.com/quick-lab/primerschemes/) to this directory                                                                       |
| --scheme-name            | N        |                           | Name of scheme to fetch from the primerschemes repository                                                                                                                                |
| --scheme-length          | N        |                           | Length of scheme to fetch from the primerschemes repository                                                                                                                              |
| --scheme-version         | N        |                           | Version of the scheme to fetch from the primerschemes repository                                                                                                                         |
| --bed                    | N        |                           | Bed file path                                                                                                                                                                            |
| --ref                    | N        |                           | Reference fasta path                                                                                                                                                                     |
| --min-mapq               | Y        | 20                        | Remove reads which map to the reference with a lower mapping quality than this                                                                                                           |
| --no-indels              | N        | False                     | Ignore insertions and deletions during variant calling, maintains the co-ordinates of the ref                                                                                            |
| --no-frameshifts         | N        | False                     | Do not allow frameshift variants (indels of lengths which are non divisible be 3 ) to be added to the consensus                                                                          |
| --align-consensus        | N        | False                     | Use MAFFT to produce an alignment of the produced consensus sequence against the reference                                                                                               |
| --linearise-fasta        | N        | False                     | Output linearised (unwrapped) FASTA consensus files                                                                                                                                      |
| --dry-run                | N        | False                     | Perform a dry run of the minion pipeline, outputing commands to a log but not executing them                                                                                             |

---

## artic_get_models

### Overview

Get the pre-trained Clair3 models provided by ONT in the [Rerio repository](https://github.com/nanoporetech/rerio/tree/master/clair3_models) 

### Usage example

```bash
artic_get_models
```

| Argument name(s) | Required | Default value             | Description                                                                                  |
| :--------------- | :------- | :------------------------ | :------------------------------------------------------------------------------------------- |
| --model-dir      | Y        | $CONDA_PREFIX/bin/models/ | Path containing clair3 models, defaults to path of models packaged with Clair3 conda package |
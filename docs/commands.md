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

- director[y/ies] to aggregate from

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

| Argument name(s)     | Required | Default value  | Description                                                                                  |
| :------------------- | :------- | :------------- | :------------------------------------------------------------------------------------------- |
| scheme               | Y        | NA             | The name of the primer scheme                                                                |
| sample               | Y        | NA             | The name of the sample                                                                       |
| --model       | Y        | NA             | Medaka or Clair3 model to use                                                                       |
| --normalise          | N        | 100            | Normalise down to moderate coverage to save runtime                                          |
| --threads            | N        | 8              | Number of threads                                                                            |
| --scheme-directory   | N        | /artic/schemes | Default scheme directory                                                                     |
| --scheme-name     | N        |                | Name of scheme to fetch from the primerschemes repository     |
| --scheme-length     | N        |                | Length of scheme to fetch from the primerschemes repository     |
| --scheme-version     | N        |                | Version of the scheme to fetch from the primerschemes repository     |
| --bed                | N        |                | Bed file path            |
| --ref                | N        |                | Reference fasta path     |
| --read-file          | N        | NA             | Use alternative FASTA/FASTQ file to <sample>.fasta                                           |
| --min-mapq           | Y        | 20             | Remove reads which map to the reference with a lower mapping quality than this               |
| --no-indels          | N        | False          | Ignore insertions and deletions during variant calling, maintains the co-ordinates of the ref|
| --no-frameshifts     | N        | False          | Do not allow frameshift variants (indels of lengths which are non divisible be 3 ) to be added to the consensus |
| --use-muscle         | N        | False          | Use muscle to produce an alignment of the produced consensus sequence against the reference  |
| --dry-run            | N        | False          | Perform a dry run of the minion pipeline, outputing commands to a log but not executing them |


---

## rampart

### Overview

Interactive prompts to start RAMPART

### Input

### Output

### Usage example

```bash

```

---

## run

### Overview

Process an entire run folder interactively

### Input

### Output

### Usage example

```bash

```

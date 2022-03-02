---
title: commands
summary: The available artic pipeline commands.
authors:
  - Will Rowe
  - Nick Loman
date: 2020-03-30
---

# Commands

This page documents the available commands via the `artic` command line interface.

---

## basecaller

### Overview

Display basecallers in files

### Input

### Output

### Usage example

```bash
artic basecaller <directory>
```

---

## demultiplex

### Overview

Run demultiplex

### Input

- undemultiplexed FASTA file

### Output

- demultiplexed FASTA file(s)

### Usage example

```bash
artic demultiplex <fasta>
```

| Argument name(s)      | Required | Default value | Description                    |
| :-------------------- | :------- | :------------ | :----------------------------- |
| fasta                 | Y        | NA            | The undemultiplexed FASTA file |
| --threads             | N        | 8             | The number of threads          |
| --prefix              | N        | NA            | Prefix for demultiplexed files |
| --no-remove-directory | N        | NA            | Don't remove the directory     |

---

## export

### Overview

The export command is used to make a redistributable package of data for re-analysis. This includes the FASTQ file, the sequencing summary and the FAST5 file. The selection of reads to be used comes from a BAM file, and only aligned reads are used.

### Input

- a completed minion pipeline run

### Output

- a redistributable package of data

### Usage example

```bash
artic export <prefix> <bamfile> <sequencing_summary> <fast5_directory> <output_directory>
```

| Argument name(s)   | Required | Default value | Description                          |
| :----------------- | :------- | :------------ | :----------------------------------- |
| prefix             | Y        | NA            | The run prefix                       |
| bamfile            | Y        | NA            | The BAM file to export reads from    |
| sequencing_summary | Y        | NA            | Path to Guppy sequencing summary     |
| fast5_directory    | Y        | NA            | The path to directory of FAST5 files |
| output_directory   | Y        | NA            | The path to export the data to       |

---

## extract

### Overview

Create an empty poredb database

### Input

- na

### Output

- an initialised poredb database

### Usage example

```bash
artic extract <directory>
```

| Argument name(s) | Required | Default value                    | Description                |
| :--------------- | :------- | :------------------------------- | :------------------------- |
| directory        | Y        | NA                               | The name of the database   |
| --basecalller    | N        | ONT Albacore Sequencing Software | The name of the basecaller |

---

## filter

### Overview

Filter FASTQ files by length

### Input

- unfiltered reads

### Output

- filtered reads

### Usage example

```bash
artic filter --max-length 500 --min-length 50 <filename>
```

| Argument name(s) | Required | Default value | Description                            |
| :--------------- | :------- | :------------ | :------------------------------------- |
| filename         | Y        | NA            | The reads to filter                    |
| --max-length     | N        | NA            | Remove reads greater than max-length   |
| --min-length     | N        | NA            | Remove reads less than than min-length |

---

## gather

### Overview

Gather up demultiplexed files

### Input

- director[y/ies] to gather from

### Output

- directory of gathered files

### Usage example

```bash
artic gather --directory ./
```

| Argument name(s)   | Required | Default value | Description                               |
| :----------------- | :------- | :------------ | :---------------------------------------- |
| --directory        | Y        | NA            | The director[y/ies] to gather files from  |
| prefix             | Y        | NA            | Prefix for gathered files                 |
| --max-length       | N        | NA            | Remove reads greater than max-length      |
| --min-length       | N        | NA            | Remove reads less than than min-length    |
| --prompt-directory | N        | NA            | The run directory for interactive prompts |
| --fast5-directory  | N        | NA            | The directory with fast5 files            |
| --no-fast5s        | N        | NA            | Do not use fast5s and nanopolish          |

---

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

| Argument name(s)            | Required | Default value  | Description                                                                                    |
| :-------------------------- | :------- | :------------- | :--------------------------------------------------------------------------------------------- |
| scheme                      | Y        | NA             | The name of the primer scheme                                                                  |
| sample                      | Y        | NA             | The name of the sample                                                                         |
| --medaka                    | N        | False          | Use medaka instead of nanopolish for variants                                                  |
| --medaka-model              | *        | NA             | Medaka model to use (required if --medaka set)                                                 |
| --minimap2                  | N        | True           | Use minimap2                                                                                   |
| --bwa                       | N        | False          | Use bwa instead of minimap2                                                                    |
| --normalise                 | N        | 100            | Normalise down to moderate coverage to save runtime                                            |
| --threads                   | N        | 8              | Number of threads                                                                              |
| --scheme-directory          | N        | /artic/schemes | Default scheme directory                                                                       |
| --max-haplotypes            | N        | 1000000        | Max-haplotypes value for nanopolish                                                            |
| --read-file                 | N        | NA             | Use alternative FASTA/FASTQ file to <sample>.fasta                                             |
| --fast5-directory           | N        | NA             | FAST5 Directory                                                                                |
| --sequencing-summary        | N        | NA             | Path to Guppy sequencing summary                                                               |
| --skip-nanopolish           | N        | False          | Skip nanopolish                                                                                |
| --min-depth                 | N        | 20             | Minimum read depth per-base for consensus masking                                              |
| --min-qual                  | N        | 20             | Minimum variant quality for consensus masking                                                  |
| --nanopolish-qual-cov-ratio | N        | 3.0            | The quality/coverage cutoff for filtering nanopolish variants                                  |
| --strict                    | N        | False          | Run withn full-length amplicon alignments & strict filtering of variants against primer scheme |
| --dry-run                   | N        | False          | Perform a dry run of the minion pipeline, outputing commands to a log but not executing them   |

* `--medaka-model` is required if `--medaka` is set.

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

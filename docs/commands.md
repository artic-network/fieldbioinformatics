---
title: commands
summary: The available artic pipeline commands.
authors:
  - Sam Wilkinson
  - Will Rowe
  - Nick Loman
date: 2024-08-16
---

# Commands

This page documents the available commands via the `artic` command line interface.

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
| --clair3             | N        | False          | Use clair3 instead of medaka for variants  (experimental feature from v1.4.0)                |
| --model       | Y        | NA             | Medaka or Clair3 model to use                                                                       |
| --normalise          | N        | 100            | Normalise down to moderate coverage to save runtime                                          |
| --threads            | N        | 8              | Number of threads                                                                            |
| --scheme-directory   | N        | /artic/schemes | Default scheme directory                                                                     |
| --max-haplotypes     | N        | 1000000        | Max-haplotypes value for nanopolish                                                          |
| --read-file          | N        | NA             | Use alternative FASTA/FASTQ file to <sample>.fasta                                           |
| --no-longshot        | N        | False          | Use medaka variant instead of longshot (experimental feautre from v1.2.0)                    |
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

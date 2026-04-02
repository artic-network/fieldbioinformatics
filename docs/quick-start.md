---
title: quick-start
summary: A quick start guide for using the artic pipeline.
authors:
  - Will Rowe
  - Nick Loman
  - Sam Wilkinson
date: 2024-11-11
---

# Quick start

This guide walks through a complete end-to-end run of the artic pipeline, from raw reads to consensus sequence. The example uses MPXV data with the `artic-inrb-mpox` scheme, but the steps are identical for any supported virus.

## Prerequisites

- The pipeline is installed (see [Installation](./installation.md))
- For conda/source installs: Clair3 models have been downloaded (see below)
- Your reads are in a single FASTQ file, basecalled with Dorado hac or sup mode

!!! note
    The pipeline requires reads from a tiling amplicon protocol. It will not produce meaningful results on metagenomic or whole-genome shotgun data.

## 1. Download Clair3 models (conda/source installs only)

```bash
artic_get_models
```

This downloads all supported Clair3 models to `$CONDA_PREFIX/bin/models/`. Only needs to be run once per environment. Skip this step if you are using the Docker image — models are pre-bundled.

## 2. Aggregate reads with guppyplex

If your reads are spread across multiple FASTQ files in a MinKNOW output directory, aggregate them first:

```bash
artic guppyplex \
  --directory ./fastq_pass \
  --min-length 1000 \
  --max-length 2000 \
  --prefix my_sample
```

This produces `my_sample_fastq_pass.fastq`. Adjust `--min-length` and `--max-length` to match the expected amplicon size for your scheme. If you already have a single FASTQ file you can skip this step.

## 3. Run the pipeline

```bash
artic minion \
  --scheme-name artic-inrb-mpox \
  --scheme-version v1.0.0 \
  --scheme-length 2500 \
  --read-file my_sample_fastq_pass.fastq \
  my_sample
```

The pipeline will:

1. Download and cache the primer scheme (first run only)
2. Select the appropriate Clair3 model automatically from the read headers
3. Align reads, trim to primer boundaries, and assign read groups
4. Call variants per primer pool with Clair3
5. Merge, filter, and normalise variants
6. Mask low-coverage positions and generate a consensus sequence

For schemes that support multiple reference sequences (e.g. multi-clade schemes), the pipeline will automatically select the best-matching reference for your sample.

## 4. Key output files

| File | Description |
| ---- | ----------- |
| `my_sample.consensus.fasta` | Final consensus sequence — the primary output |
| `my_sample.normalised.vcf.gz` | Normalised PASS variants applied to the consensus |
| `my_sample.pass.vcf` | PASS variants before normalisation |
| `my_sample.fail.vcf` | Variants that did not pass quality filtering |
| `my_sample.primertrimmed.rg.sorted.bam` | Primer-trimmed alignment |
| `my_sample.sorted.bam` | Raw alignment |
| `my_sample.minion.log.txt` | Full command log with runtimes |

For a complete description of all output files, see [Core Pipeline](./minion.md).

## Common options

**Specify the Clair3 model explicitly** (if auto-detection fails or you want to override):
```bash
artic minion ... --model r1041_e82_400bps_hac_v520 my_sample
```

**Produce a MAFFT alignment of the consensus against the reference:**
```bash
artic minion ... --align-consensus my_sample
```

**Increase coverage depth requirement:**
```bash
artic minion ... --min-depth 40 my_sample
```

**Dry run** (print commands without executing):
```bash
artic minion ... --dry-run my_sample
```

## Troubleshooting

If the pipeline exits unexpectedly, see [Troubleshooting](./troubleshooting.md) for a description of exit codes and common error causes.

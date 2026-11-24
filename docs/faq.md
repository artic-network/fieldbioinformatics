---
title: faq
summary: The FAQ.
authors:
  - Will Rowe
  - Nick Loman
  - Sam Wilkinson
date: 2020-03-30
---

# FAQ

## How do I process MPXV data?

A set of resources for processing MPXV sequencing data may be found [here](https://artic.network/mpxv), including running this pipeline on the command line and via the artic MPXV Nextflow pipelines in epi2me.

## Where can I find the SOP for SARS-CoV-2?

The standard operating procedure for the ARTIC Network SARS-CoV-2 bioinformatics can be found [here](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

## The pipeline exited with an error code — what does it mean?

See [Troubleshooting](./troubleshooting.md) for a full description of exit codes.

## I am using Dorado-basecalled reads but the pipeline fails with a model error

The pipeline reads the `basecall_model_version_id` tag from the first read header to select a Clair3 model automatically. If the tag is present but no compatible Clair3 model exists, the pipeline will exit with code 6.

A common cause is data basecalled with Dorado **fast** mode: no versioned Clair3 fast models exist for R10.4.1 data. Re-basecall with `hac` or `sup`.
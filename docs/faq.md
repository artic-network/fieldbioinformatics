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

A set of resources for processing MPXV sequencing data may be found [here](https://artic.network/mpxv), this includes running this pipeline on the command line and the artic MPXV nextflow pipelines via epi2me.

## Where can I find the SOP for SARS-CoV-2?

The standard operating proceedure for the ARTIC Network SARS-SoV-2 bioinformatics can be found [here](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

## Lab-on-an-SSD

Please refer to [the ARTIC website](https://artic.network/lab-on-an-SSD) for more information about lab-on-SSD.

## Adding this repository as a submodule of another

Within the parent repo add the submodule:

```
git submodule add https://github.com/artic-network/fieldbioinformatics.git
```

Commit the change and push:

```
git commit -m "adding submodule"
git push origin master
```

To update all submodules:

```
git submodule update --remote
```

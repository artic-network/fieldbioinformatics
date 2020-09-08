---
title: faq
summary: The FAQ.
authors:
  - Will Rowe
  - Nick Loman
date: 2020-03-30
---

# FAQ

## Where can I find the SOP for SARS-CoV-2

The standard operating proceedure for the ARTIC Network SARS-SoV-2 bioinformatics can be found [here](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html).

## Should I use the nanopolish or medaka workflow

We currently recommend the nanopolish workflow (if you have signal data available) as we have spent more time validating and supporting this workflow. That being said, both tend to give consistent results with our test datasets so the choice is yours.

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

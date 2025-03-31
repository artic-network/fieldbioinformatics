---
title: Primer Schemes
summary: Outline of primer scheme handling.
authors:
  - Will Rowe
  - Nick Loman
  - Sam Wilkinson
date: 2025-03-31
---

# Primer Schemes

## About

Both `artic tools` and the `artic pipeline` are designed to work with tiling amplicon schemes for targeted enrichment of viral genomes. The protocol is described in the [Kent C et al. 2024](https://doi.org/10.1101/2024.12.20.629611) primalscheme paper and  [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066) Nature Protocols paper and the schemes are produced by [Primal Scheme](https://primalscheme.com/), details of available schemes .

Tiling amplicon schemes contain primers that are designed using a greedy algorithm against multiple reference genomes. Tiling refers to the fact that neighbouring amplicons overlap one another, ensuring that complete coverage of the target genome is acheived. To prevent short overlap products being produced between neighbouring amplicons, 2 primer pools are used to alternate primer pairs. The following figure from the [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066) paper shows how primers from 2 pools are used to tile across a reference genome:

![primal-scheme-fig](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnprot.2017.066/MediaObjects/41596_2017_Article_BFnprot2017066_Fig3_HTML.jpg?as=webp)

> [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066)
> (a) Schematic showing the regions amplified in pools 1 (upper track) and 2 (lower track), and the intended overlap between pools (as determined in Step 1)...

## Availability

Frequently used primer schemes for viral genome sequencing are found in the [Quick-lab primer scheme repo](https://github.com/quick-lab/primerschemes) and details of available schemes in [primalscheme labs](https://labs.primalscheme.com/).
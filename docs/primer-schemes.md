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

Both `artic tools` and the `artic pipeline` are designed to work with tiling amplicon schemes for targeted enrichment of viral genomes. The protocol is described in the [Kent C et al. 2024](https://doi.org/10.1101/2024.12.20.629611) primalscheme paper and [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066) Nature Protocols paper. Schemes are produced by [Primal Scheme](https://primalscheme.com/).

Tiling amplicon schemes contain primers designed against multiple reference genomes using a greedy algorithm. Neighbouring amplicons overlap, ensuring complete coverage of the target genome. Two primer pools are used to alternate primer pairs, preventing short overlap products between adjacent amplicons. The figure below illustrates how pools tile across a genome:

![primal-scheme-fig](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnprot.2017.066/MediaObjects/41596_2017_Article_BFnprot2017066_Fig3_HTML.jpg?as=webp)

> [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066)
> (a) Schematic showing the regions amplified in pools 1 (upper track) and 2 (lower track), and the intended overlap between pools...

## Availability

Frequently used primer schemes for viral genome sequencing are found in the [Quick-lab primer scheme repo](https://github.com/quick-lab/primerschemes). A browsable index is available at [primalscheme labs](https://labs.primalscheme.com/).

## Scheme fetching and caching

When `--scheme-name` and `--scheme-version` are provided, the pipeline fetches the scheme manifest from:

```
https://raw.githubusercontent.com/quick-lab/primerschemes/main/index.json
```

Downloaded scheme files are cached in `--scheme-directory` (default: `./primer-schemes`). If the network is unavailable but a previously cached `index.json` and `aliases.json` exist in that directory, the pipeline will use the local cache and print a warning.

## Scheme aliases

The pipeline also consults an `aliases.json` file from the same repository. This allows short or alternative names to resolve transparently to canonical scheme names. For example, a user may supply `--scheme-name sars-cov-2` and have it resolve to the canonical name automatically. No extra action is required — alias resolution happens silently during scheme lookup.

## Automatic reference selection

Some schemes support multiple closely related reference sequences (for example, schemes covering multiple clades of a virus). When the scheme manifest contains a `refselect` key, the pipeline automatically determines the best-matching reference for a given sample:

1. Up to 10,000 reads are subsampled from the input FASTQ.
2. These reads are aligned against a multi-reference FASTA using minimap2 with 4 threads.
3. The reference with the most aligned reads is selected.
4. The selected reference suffix is appended to `--scheme-version` (e.g. `v1.0.0` becomes `v1.0.0-clade_i`).

To bypass automatic selection, provide the full suffixed version directly:

```
--scheme-version v1.0.0-clade_i
```

## BED file format versions

The pipeline auto-detects the primer BED format version by inspecting the first non-comment line. Three versions are recognised:

| Version | Columns | Primer name format | Example |
| ------- | ------- | ------------------ | ------- |
| V1 | 6 | `SCHEME_NUM_LEFT` / `SCHEME_NUM_RIGHT` (with optional `_alt` suffix) | `nCoV-2019_1_LEFT` |
| V2 | 7 | Same name pattern as V1 | `nCoV-2019_1_LEFT` |
| V3 | 7 | `SCHEME_NUM_LEFT_POOL` / `SCHEME_NUM_RIGHT_POOL` | `nCoV-2019_1_LEFT_1` |

The format version determines how amplicons and primer pools are inferred during alignment post-processing.

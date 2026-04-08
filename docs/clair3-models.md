---
title: Clair3 Models
summary: Available Clair3 models and how they are selected.
authors:
  - Sam Wilkinson
date: 2024-11-11
---

# Clair3 Models

## Overview

The pipeline uses [Clair3](https://github.com/HKU-BAL/Clair3) for variant calling. Clair3 requires a model that matches the pore chemistry and basecalling preset used to produce the input reads. Models are downloaded from the [HKU Clair3 model repository](https://www.bio8.cs.hku.hk/clair3/) and stored as PyTorch checkpoints (`pileup.pt` and `full_alignment.pt`).

## Downloading models

```sh
artic_get_models --models r1041_e82_400bps_sup_v400 r1041_e82_400bps_sup_v420 ...
```

Will fetch the models provided as a space-delimited list to the `--models` parameter, if no `--models` parameter is provided all models will be fetched.

Models are stored in `$CONDA_PREFIX/bin/models/` by default. For a custom location:

```sh
artic_get_models --model-dir /path/to/models
```

The `*-models-included` Docker image variant bundles all models. The default Docker image does not include models — use `--model-dir` to mount a pre-downloaded model directory. `artic_get_models` is needed for conda/source installs and to prepare models for use with the default Docker image.

## Automatic model selection

When `--model` is not provided, the pipeline reads the `basecall_model_version_id` tag from the first read in the input FASTQ and selects a model automatically. The tag is written by Dorado into each read header, for example:

```
@read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.2.0 ...
```

The selection algorithm:

1. Splits the tag on `_` to extract pore type, kit ID, speed tier, and basecaller version.
2. Filters available models by pore type and speed tier.
3. Matches the basecaller version string (e.g. `v520` for Dorado 5.2.0).
4. If exactly one model matches at each step it is selected and a warning is printed. If no match is found at any step the pipeline exits with code 6.

Use `--model` to override auto-selection at any time.

## Guppy vs Dorado compatibility

Clair3 models whose names contain a `g`-prefixed version suffix (e.g. `_g632`, `_g5014`) were trained on data basecalled with **Guppy**. These models are **not compatible with Dorado-basecalled data** and will produce incorrect variant calls if used with it.

If you see exit code 6 for a `fast`-mode Dorado model: no versioned Clair3 fast model exists for R10.4.1 data. Re-basecall with `hac` or `sup`.

## Available models

### Dorado models

These models correspond to data basecalled with [Dorado](https://github.com/nanoporetech/dorado). The Dorado basecall model column shows the `basecall_model_version_id` tag value that will auto-select each model.

#### R10.4.1 — 400 bps

| Clair3 model                | Dorado basecall model                |
| --------------------------- | ------------------------------------ |
| `r1041_e82_400bps_hac_v520` | `dna_r10.4.1_e8.2_400bps_hac@v5.2.0` |
| `r1041_e82_400bps_sup_v520` | `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` |
| `r1041_e82_400bps_hac_v500` | `dna_r10.4.1_e8.2_400bps_hac@v5.0.0` |
| `r1041_e82_400bps_sup_v500` | `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` |
| `r1041_e82_400bps_hac_v430` | `dna_r10.4.1_e8.2_400bps_hac@v4.3.0` |
| `r1041_e82_400bps_sup_v430` | `dna_r10.4.1_e8.2_400bps_sup@v4.3.0` |
| `r1041_e82_400bps_hac_v420` | `dna_r10.4.1_e8.2_400bps_hac@v4.2.0` |
| `r1041_e82_400bps_sup_v420` | `dna_r10.4.1_e8.2_400bps_sup@v4.2.0` |
| `r1041_e82_400bps_hac_v410` | `dna_r10.4.1_e8.2_400bps_hac@v4.1.0` |
| `r1041_e82_400bps_sup_v410` | `dna_r10.4.1_e8.2_400bps_sup@v4.1.0` |
| `r1041_e82_400bps_hac_v400` | `dna_r10.4.1_e8.2_400bps_hac@v4.0.0` |
| `r1041_e82_400bps_sup_v400` | `dna_r10.4.1_e8.2_400bps_sup@v4.0.0` |

!!! warning
    No versioned Clair3 fast models exist for R10.4.1 data. Reads basecalled with `dna_r10.4.1_e8.2_400bps_fast@vX.Y.Z` cannot be automatically matched to a model. Re-basecall with `hac` or `sup`.

#### R10.4.1 — 260 bps

| Clair3 model                | Dorado basecall model                |
| --------------------------- | ------------------------------------ |
| `r1041_e82_260bps_hac_v410` | `dna_r10.4.1_e8.2_260bps_hac@v4.1.0` |
| `r1041_e82_260bps_sup_v410` | `dna_r10.4.1_e8.2_260bps_sup@v4.1.0` |
| `r1041_e82_260bps_hac_v400` | `dna_r10.4.1_e8.2_260bps_hac@v4.0.0` |
| `r1041_e82_260bps_sup_v400` | `dna_r10.4.1_e8.2_260bps_sup@v4.0.0` |

#### R9.4.1

| Clair3 model              | Dorado basecall model    |
| ------------------------- | ------------------------ |
| `r941_prom_hac_g360+g422` | `dna_r9.4.1_e8_hac@vX.X` |
| `r941_prom_sup_g5014`     | `dna_r9.4.1_e8_sup@vX.X` |

R9.4.1 has only two models (one per preset). Auto-selection returns the correct model for any R9.4.1 hac or sup tag regardless of basecaller version.

---

### Guppy models (legacy)

These models were trained on data basecalled with [Guppy](https://community.nanoporetech.com/downloads). They are **not compatible with Dorado-basecalled data**. Use them only if your reads were produced by Guppy.

#### R10.4.1 — 400 bps

| Clair3 model                 | Guppy version | Preset |
| ---------------------------- | ------------- | ------ |
| `r1041_e82_400bps_hac_g632`  | g632          | hac    |
| `r1041_e82_400bps_hac_g615`  | g615          | hac    |
| `r1041_e82_400bps_sup_g615`  | g615          | sup    |
| `r1041_e82_400bps_fast_g632` | g632          | fast   |
| `r1041_e82_400bps_fast_g615` | g615          | fast   |

#### R10.4.1 — 260 bps

| Clair3 model                 | Guppy version | Preset |
| ---------------------------- | ------------- | ------ |
| `r1041_e82_260bps_hac_g632`  | g632          | hac    |
| `r1041_e82_260bps_sup_g632`  | g632          | sup    |
| `r1041_e82_260bps_fast_g632` | g632          | fast   |

#### R10.4 (E8.1)

| Clair3 model         | Guppy version | Preset |
| -------------------- | ------------- | ------ |
| `r104_e81_hac_g5015` | g5015         | hac    |
| `r104_e81_sup_g5015` | g5015         | sup    |

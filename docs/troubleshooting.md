---
title: Troubleshooting
summary: Exit codes and common error causes.
authors:
  - Sam Wilkinson
date: 2024-11-11
---

# Troubleshooting

## Exit codes

| Exit code | Cause | Resolution |
| --------- | ----- | ---------- |
| **2** | No reads aligned to the reference | Check that `--scheme-name`/`--scheme-version` match your data, and that your reads are from a tiling amplicon protocol against the same target |
| **3** | Primer scheme BED file could not be parsed | The BED file is malformed or uses a format not recognised by primalbedtools. The error message will include the specific parse failure. Check the file is a valid primer scheme BED and matches the expected format (see [Primer Schemes](./primer-schemes.md)) |
| **4** | Input FASTQ is empty | The read file contains no reads. Check the file exists and is not empty. If using `artic guppyplex`, check that the source directory contained FASTQ files passing your length/quality filters |
| **6** | Clair3 model selection failed | The `basecall_model_version_id` tag was missing, unrecognised, or matched only Guppy-era models that are incompatible with Dorado data. Use `--model` to specify a model explicitly — see [Clair3 Models](./clair3-models.md) |
| **127** | Required tool not found in PATH | Ensure all dependencies are installed and the correct conda environment is activated |
| **137** | Out of memory (OOM) | Try reducing `--threads`, enabling or lowering `--normalise`, or running on a machine with more memory |
| **143** | Process killed by signal | Typically a user interrupt (`Ctrl+C`) or a system-level kill. Re-run when resources are available |

## Common issues

### "No reads aligned to the reference" (exit 2)

The most frequent causes are:

- Wrong scheme: the primer scheme does not match the virus or amplicon design used in the library prep.
- Wrong reference: the reference sequence is too divergent from the sequenced sample.
- Very low depth: fewer than `--min-depth` reads mapped to any amplicon. Check raw read counts and quality.

### Model selection fails (exit 6)

The pipeline reads `basecall_model_version_id` from the first read header. If the tag is absent (e.g. the FASTQ was produced by an older basecaller or reformatted), auto-selection is not possible. Provide the model name directly:

```sh
artic minion --model r1041_e82_400bps_hac_v520 ... my_sample
```

If your reads were basecalled with Dorado **fast** mode on R10.4.1 data, no compatible Clair3 model exists. Re-basecall with `hac` or `sup`. See [Clair3 Models](./clair3-models.md) for the full model list.

### Out of memory (exit 137)

Clair3 is memory-intensive, particularly on long references or at high depth. Options:

- Reduce `--threads` (less parallelism, less peak memory).
- Enable depth normalisation: `--normalise 100` (the default) caps coverage per amplicon before variant calling.
- Use a machine with more RAM.

### Some variants are absent from the merged VCF

Variants called within primer binding sites are separated into `$SAMPLE.primers.vcf` and excluded from `$SAMPLE.merged.vcf`. This is expected behaviour — these positions are unreliable due to primer softmasking. Check `$SAMPLE.primersitereport.txt` for details.

### Pipeline produces no consensus or an all-N consensus

This usually means coverage is below `--min-depth` across most of the genome. Check `$SAMPLE.amplicon_depths.tsv` to see per-amplicon depth. Common causes:

- Very low input read count.
- Amplicon dropout (one or more amplicons failed in the library prep).
- Normalisation set too low (`--normalise` value less than actual mean depth, discarding too many reads).

### Checking the full command log

Every command executed by the pipeline, along with its wall-clock runtime, is written to `$SAMPLE.minion.log.txt`. This is the first place to look when diagnosing an unexpected failure.

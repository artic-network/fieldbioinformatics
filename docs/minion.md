---
title: minion
summary: Outline of minion workflows.
authors:
  - Will Rowe
  - Nick Loman
  - Sam Wilkinson
date: 2024-11-11
---

# Core pipeline

## About

This page describes the core pipeline which is run via the `artic minion` command.

Output files are grouped below into **outputs** (files you should use downstream) and **intermediates** (files left on disk after the run which may be useful for troubleshooting but are not intended as primary outputs).

## Stages

### Input validation

The pipeline fetches primer scheme and reference files from the [Quick-lab primerschemes repository](https://github.com/quick-lab/primerschemes) when `--scheme-name`, `--scheme-version`, and (where required) `--scheme-length` are provided. Downloaded files are cached in the directory given by `--scheme-directory` (default: `./primer-schemes`).

For example, to fetch the [artic-inrb-mpox/2500/v1.0.0](https://github.com/quick-lab/primerschemes/tree/main/primerschemes/artic-inrb-mpox/2500/v1.0.0) scheme:

```
--scheme-name artic-inrb-mpox
--scheme-length 2500
--scheme-version v1.0.0
```

Alternatively, provide `--bed` and `--ref` to use a local primer BED file and reference FASTA directly.

For further detail on scheme fetching, aliases, and automatic reference selection, see [Primer Schemes](./primer-schemes.md).

### Reference alignment and post-processing

The pipeline aligns basecalled reads against the reference using [minimap2](https://github.com/lh3/minimap2) with the ONT preset (`map-ont`). Alignments are filtered to remove unmapped reads, then sorted and indexed with samtools.

The `align_trim` module then post-processes the alignments to:

- assign each read to a derived amplicon
- assign each read a **read group** based on its primer pool
- softmask alignments within their amplicon boundaries

Optionally it can also:

- remove primer sequence by further softmasking (`--primer-match-threshold` controls fuzzy matching)
- downsample reads per amplicon to `--normalise` depth (set to `0` to disable)
- remove reads with mismatched primer pairs, e.g. from amplicon read-through (use `--allow-mismatched-primers` to retain these)

Softmasking adjusts the CIGAR of each [alignment segment](https://samtools.github.io/hts-specs/SAMv1.pdf) so that soft clips replace any reference- or query-consuming operations outside primer boundaries. The leftmost mapping position is updated accordingly.

More information on how the primer scheme is used to infer amplicons can be found in [Primer Schemes](./primer-schemes.md).

#### Outputs

| File | Description |
| ---- | ----------- |
| `$SAMPLE.sorted.bam` / `.bai` | Raw alignment of reads to the reference |
| `$SAMPLE.primertrimmed.rg.sorted.bam` / `.bai` | Primer-trimmed, read-group-annotated alignment used for all downstream steps |

#### Intermediates

| File | Description |
| ---- | ----------- |
| `$SAMPLE.alignreport.tsv` | Per-read primer assignment report from align_trim |
| `$SAMPLE.amplicon_depths.tsv` | Per-amplicon read depth report |

### Variant calling

Variant calling is performed per read group using [Clair3](https://github.com/HKU-BAL/Clair3) (`run_clair3.sh`). Per-pool VCFs are merged into a single file using `artic_vcf_merge`.

During merging, any variants that fall within primer binding sites are separated into `$SAMPLE.primers.vcf` and excluded from the merged output — these positions are unreliable due to primer softmasking. The details are written to `$SAMPLE.primersitereport.txt`.

The merged variants are then filtered by `artic_vcf_filter` into three output files. Each failing variant is classified as either **mask** or **discard**:

- **Mask** — the position is genuinely ambiguous; it will be replaced with `N` in the consensus. This applies when the variant has low quality, mixed reads (low allele frequency above the ignore threshold), a frameshift, or insufficient depth.
- **Discard** — the variant call is likely an artefact but the position itself is fine; the reference base is kept. This applies when the allele frequency is below the ignore threshold (suggesting noise rather than a real mixed site) or when the absolute alt read count is too low despite adequate quality.

| Filter | Default threshold | Result on failure |
| ------ | ----------------- | ----------------- |
| Variant quality (QUAL) | ≥ 10 (`--min-variant-quality`) | mask |
| Allele frequency (AF) — lower bound | ≥ 0.1 (`--min-mask-allele-frequency`) | discard (below this, variant is noise) |
| Allele frequency (AF) — upper bound | ≥ 0.6 (`--min-allele-frequency`) | mask (between bounds, position is ambiguous) |
| Frameshift indel quality | ≥ 50 (`--min-frameshift-quality`), or always excluded with `--no-frameshifts` | mask |
| Read depth (DP) | ≥ `--min-depth` | mask |
| Alt allele read count (AD) | ≥ 5 (`--min-minor-allele-count`) | discard (too few reads despite good quality) |
| All indels | excluded when `--no-indels` is set | — |

Finally, passing variants are normalised against the pre-consensus using [bcftools norm](https://github.com/samtools/bcftools) to ensure REF alleles match the masked reference before consensus generation.

#### Outputs

| File | Description |
| ---- | ----------- |
| `$SAMPLE.normalised.vcf.gz` / `.tbi` | Normalised PASS variants — these are the variants applied to produce the consensus |
| `$SAMPLE.pass.vcf` | PASS variants before normalisation |
| `$SAMPLE.fail.vcf` | Mask variants — positions that will be replaced with `N` in the consensus |
| `$SAMPLE.ignore.vcf` | Discarded variants — rejected as likely artefacts; the reference base is kept at these positions |

#### Intermediates

| File | Description |
| ---- | ----------- |
| `$SAMPLE.merged.vcf` | Pre-filter merged variants from all read groups |
| `$SAMPLE.$POOL.vcf` | Raw per-pool Clair3 output (one file per primer pool) |
| `$SAMPLE_rg_$POOL/` | Full Clair3 output directory per pool |
| `$SAMPLE.primers.vcf` | Variants at primer binding sites, excluded from the main merge |
| `$SAMPLE.primersitereport.txt` | Report of primer-site variant handling |

### Consensus building

Each position in the reference is checked for read depth against the value of `--min-depth` (default: 20) using `artic_make_depth_mask`. Positions below this threshold are recorded in the coverage mask.

`artic_mask` then produces a pre-consensus sequence by applying low-coverage masking and the mask variants (`$SAMPLE.fail.vcf`) to the reference, replacing low-confidence positions with `N`. Discarded variants (`$SAMPLE.ignore.vcf`) are not applied — the reference base is used at those positions. `bcftools consensus` is then run against the pre-consensus using the normalised PASS variants to produce the final consensus sequence. The consensus header is annotated with the artic workflow identifier by `artic_fasta_header`.

If `--align-consensus` is provided, the consensus is aligned against the reference using:

```
mafft --6merpair --addfragments $SAMPLE.consensus.fasta $REF > $SAMPLE.aligned.fasta
```

#### Outputs

| File | Description |
| ---- | ----------- |
| `$SAMPLE.consensus.fasta` | Final consensus sequence |
| `$SAMPLE.aligned.fasta` | MAFFT alignment of consensus against reference (only produced with `--align-consensus`) |
| `$SAMPLE.primer.bed` | Copy of the primer scheme BED used for this run |
| `$SAMPLE.reference.fasta` | Copy of the reference FASTA used for this run |
| `$SAMPLE.minion.log.txt` | Log of all commands executed with wall-clock runtimes |

#### Intermediates

| File | Description |
| ---- | ----------- |
| `$SAMPLE.coverage_mask.txt` | BED-format record of positions masked below `--min-depth` |
| `$SAMPLE.preconsensus.fasta` | N-masked reference with FAIL variants applied; input to `bcftools consensus` |

!!! note
    `$SAMPLE.primer.bed` and `$SAMPLE.reference.fasta` record which scheme was used for this run. In future, when per-segment reference selection is supported, these will become primary outputs of greater significance.

## Summary of pipeline modules

| Module | Function |
| ------ | -------- |
| `align_trim` | Amplicon assignment, softmasking, read-group annotation, and normalisation |
| `artic_vcf_merge` | Combines per-pool VCFs; separates primer-site variants |
| `artic_vcf_filter` | Filters merged VCF into PASS and FAIL files |
| `artic_make_depth_mask` | Produces a coverage mask from the post-processed alignment |
| `artic_mask` | Combines reference, FAIL variants, and coverage mask into a pre-consensus |
| `artic_fasta_header` | Applies the artic workflow identifier to the consensus header |
| `artic_get_models` | Downloads Clair3 model files |
| `artic_get_scheme` | Downloads a primer scheme independently of the full pipeline |

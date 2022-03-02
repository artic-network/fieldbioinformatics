---
title: minion
summary: Outline of minion workflows.
authors:
  - Will Rowe
  - Nick Loman
date: 2020-09-01
---

# Core pipeline

## About

This page describes the core pipeline which is run via the `artic minion` command.

There are **2 workflows** baked into the core pipeline, one which uses signal data (via [nanopolish](https://github.com/jts/nanopolish)) and one that does not (via [medaka](https://github.com/nanoporetech/medaka)). As the workflows are identical in many ways, this page will describe the pipeline as whole and notify the reader when there is dfferent behaviour between the two workflows.
It should be noted here that by default the `nanopolish` workflow is selected; you need to specify `--medaka` (and `--medaka-model`) if you want the medaka workflow enabled.

> **NOTE**: It is very important that you select the appropriate value for `--medaka-model`.

At the end of each stage, we list here the "useful" stage output files which are kept. There will also be some additional files leftover at the end of the pipeline but these can be ignored (and are hopefully quite intuitively named).

## Stages

### Input validation

The supplied `--scheme-directory` will be checked for a **reference sequence** and **primer scheme**. If `--scheme-version` is not supplied, version 1 will be assumed.

As of version 1.2.0, the pipeline will try and download the reference and scheme from the [artic primer scheme](https://github.com/artic-network/primer-schemes) repository if it is not found/provided.

Once the scheme and reference have been found, the pipeline will validate the scheme and extract the primer pool information.

If running the **nanopolish workflow**, the pipeline will check the reference only contains one sequence and then run `nanopolish index` to map basecalled reads to the signal data.

#### stage output

- primer validation log (optional)
- nanopolish index (workflow dependant)

### Reference alignment and post-processing

The pipeline will then perform a reference alignment of the basecalled reads against the specified reference sequence. By default [minimap](https://github.com/lh3/minimap2) is used but [bwa](https://github.com/lh3/bwa) can be chosen as an alternative. Both aligners use their respective ONT presets. The alignments are filtered to keep only mapped reads, and then sorted and indexed.

We then use the `align_trim` module to post-process the aligments.

The purpose of alignment post-processing is:

- assign each read alignment to a derived amplicon
- using the derived amplicon, assign each read a **read group** based on the primer pool
- softmask read alignments within their derived amplicon

Also, there is the option to:

- remove primer sequence by further softmasking the read alignments
- normalise/reduce the number of read alignments to each amplicon
- remove reads with imperfect primer pairing, e.g. from amplicon read through

By softmasking, we refer to the process of adjusting the CIGAR of each [alignment segment](https://samtools.github.io/hts-specs/SAMv1.pdf) such that **soft clips** replace any reference or query consuming operations in regions of an alignment that fall outside of primer boundaries. The leftmost mapping position of alignments are also updated during softmasking.

More information on how the primer scheme is used to infer amplicons can be found [here](./primer-schemes.md#querying-schemes).

#### stage output

| file name                             | description                                                                          |
| ------------------------------------- | ------------------------------------------------------------------------------------ |
| `$SAMPLE.sorted.bam`                  | the raw alignment of sample reads to reference genome                                |
| `$SAMPLE.trimmed.rg.sorted.bam`       | the post-processed alignment                                                         |
| `$SAMPLE.primertrimmed.rg.sorted.bam` | the post-processed alignment with additional softmasking to exclude primer sequences |

### Variant calling

Once alignments have been softmasked, sorted and indexed (again), they are used for variant calling. This is where the two workflows actually differ.

For the **medaka workflow**, we use the following commands on the `$SAMPLE.primertrimmed.rg.sorted.bam` alignment:

- [medaka consensus](https://github.com/nanoporetech/medaka)
- medaka variant or snps (if pipeline has been told not to detect INDELS via `--no-indels`)
- medaka tools annotate

And for the **nanopolish workflow** we use the following command on the `$SAMPLE.trimmed.rg.sorted.bam` alignment:

- [nanopolish variants](https://github.com/jts/nanopolish)

For both workflows, the variant calling steps are run for each **read group** in turn. We then merge variants reported per read group into a single file using the `artic_vcf_merge` module.

> Note: we use the `$SAMPLE.trimmed.rg.sorted.bam` alignment for the **nanopolish workflow** as nanopolish requires some leading sequence to make variant calls; we use the primer sequence for this purpose in order to call variants at the start of amplicons.

Opionally, we can check the merged variant file against the primer scheme. This will allow us to detect variants called in primer and amplicon overlap regions.

Finally, we use the `artic_vcf_filter` module to filter the merged variant file through a set of workflow specific checks and assign all variants as either PASS or FAIL. The final PASS file is subsequently indexed ready for the next stage.

#### stage output

| file name                | description                                                     |
| ------------------------ | --------------------------------------------------------------- |
| `$SAMPLE.$READGROUP.vcf` | the raw variants detected (one file per primer pool)            |
| `$SAMPLE.merged.vcf.gz`  | the raw variants detected merged into one file                  |
| `$SAMPLE.vcfreport.txt`  | a report evaluating reported variants against the primer scheme |
| `$SAMPLE.fail.vcf`       | variants deemed too low quality                                 |
| `$SAMPLE.pass.vcf.gz`    | detected variants (indexed)                                     |

### Consensus building

Prior to building a consensus, we use the post-processed alignment from the previous step to check each position of the reference sequence for sample coverage. Any poition that is not covered by at least 20 reads from either read group are marked as low coverage. We use the `artic_make_depth_mask` module for this, which produces coverage information for each read group and also produces a coverage mask to tell us which coordinates in the reference sequence failed the coverage threshold.

Next, to build a consensus sequence for a sample, we require a pre-consensus sequence based on the input reference sequence. The preconsensus has low quality sites masked out with `N`'s using the coverage mask and the `$SAMPLE.fail.vcf` file. We then use `bcftools consensus` to combine the preconsensus with the `$SAMPLE.pass.vcf` variants to produce a consensus sequence for the sample. The consensus sequence has the artic workflow written to its header.

Finally, the consensus sequence is aligned against the reference sequence using `muscle`.

#### stage output

| file name                  | description                                                           |
| -------------------------- | --------------------------------------------------------------------- |
| `$SAMPLE.*_mqc.json`       | stats files which MultiQC can use to make a report                    |
| `$SAMPLE.consensus.fasta`  | the consensus sequence for the input sample                           |
| `$SAMPLE.muscle.out.fasta` | an alignment of the consensus sequence against the reference sequence |

## Summary of pipeline modules

| module                | function                                                                                             |
| --------------------- | ---------------------------------------------------------------------------------------------------- |
| align_trim            | alignment post processing (amplicon assignment, softmasking, normalisation)                          |
| artic_vcf_merge       | combines VCF files from multiple read groups                                                         |
| artic_vcf_filter      | filters a combined VCF into PASS and FAIL variant files                                              |
| artic_make_depth_mask | create a coverage mask from the post-processed alignment                                             |
| artic_mask            | combines the reference sequence, FAIL variants and coverage mask to produce a pre-consensus sequence |
| artic_fasta_header    | applies the artic workflow and identifier to the consensus sequence header                           |
| artic_get_stats       | collects stats from the pipeline output files and generate summary files for use by MultiQC          |

## Optional pipeline report

As of version 1.3.0, you can run MultiQC (which should be installed as part of the artic conda environment) on the pipeline output directory and this will produce a report containing amplicon coverage plots and variant call information. To generate a report from within your pipeline output directory:

```
multiqc .
```
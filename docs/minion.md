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

At the end of each stage, we list here the "useful" stage output files which are kept. There will also be some additional files leftover at the end of the pipeline but these can be ignored (and are hopefully quite intuitively named).

## Stages

### Input validation

As of version 1.2.0, the pipeline will try and download the reference and scheme from the [artic primer scheme](https://github.com/artic-network/primer-schemes) repository if it is not found/provided. It will be downloaded to the directory provided by `--scheme-directory` which defaults to a subfolder of the current working directory `primer-schemes/`

This is done using the `--scheme-name` `--scheme-length` and `--scheme-version` parameters, for example, to fetch the [artic-inrb-mpox/2500/v1.0.0](https://github.com/quick-lab/primerschemes/tree/main/primerschemes/artic-inrb-mpox/2500/v1.0.0) scheme you should provide the following arguments:
* `--scheme-name artic-inrb-mpox`
* `--scheme-length 2500`
* `--scheme-version v1.0.0`

Alternatively you may provide the `--bed` and `--ref` arguments to point directly towards the primer bed file and reference fasta you wish to use.

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

We use the following commands on the `$SAMPLE.primertrimmed.rg.sorted.bam` alignment:

- [run_clair3.sh](https://github.com/HKU-BAL/Clair3)

The variant calling steps are run for each **read group** in turn. We then merge variants reported per read group into a single file using the `artic_vcf_merge` module.

Opionally, we can check the merged variant file against the primer scheme. This will allow us to detect variants called in primer and amplicon overlap regions.

We then use the `artic_vcf_filter` module to filter the merged variant file through a set of workflow specific checks and assign all variants as either PASS or FAIL. The final PASS file is subsequently indexed ready for the next stage.

Finally the variants which have passed filtering are normalised against the depth masked pre-consensus using [bcftools norm](https://github.com/samtools/bcftools) to ensure that the REF column in the VCF matches the pre-consensus for the final output.

#### stage output

| file name                | description                                                     |
| ------------------------ | --------------------------------------------------------------- |
| `$SAMPLE.$READGROUP.vcf` | the raw variants detected (one file per primer pool)            |
| `$SAMPLE.merged.vcf`     | the raw variants detected merged into one file                  |
| `$SAMPLE.vcfreport.txt`  | a report evaluating reported variants against the primer scheme |
| `$SAMPLE.fail.vcf`       | variants deemed too low quality                                 |
| `$SAMPLE.pass.vcf.gz`    | detected variants (indexed)                                     |
| `$SAMPLE.normalised.vcf.gz` | normalised variants (indexed)                                     |

### Consensus building

Prior to building a consensus, we use the post-processed alignment from the previous step to check each position of the reference sequence for sample coverage. Any poition that is not covered by at least 20 reads from either read group are marked as low coverage. We use the `artic_make_depth_mask` module for this, which produces coverage information for each read group and also produces a coverage mask to tell us which coordinates in the reference sequence failed the coverage threshold.

Next, to build a consensus sequence for a sample, we require a pre-consensus sequence based on the input reference sequence. The preconsensus has low quality sites masked out with `N`'s using the coverage mask and the `$SAMPLE.fail.vcf` file. We then use `bcftools consensus` to combine the preconsensus with the `$SAMPLE.normalised.vcf.gz` variants to produce a consensus sequence for the sample. The consensus sequence has the artic workflow written to its header.

Finally, the consensus sequence is aligned against the reference sequence using `muscle` if the `--use-muscle` parameter is provided.

#### stage output

| file name                  | description                                                           |
| -------------------------- | --------------------------------------------------------------------- |
| `$SAMPLE.*_mqc.json`       | stats files which MultiQC can use to make a report                    |
| `$SAMPLE.consensus.fasta`  | the consensus sequence for the input sample                           |
| `$SAMPLE.muscle.out.fasta` | an alignment of the consensus sequence against the reference sequence |

## Summary of pipeline modules

| module                    | function                                                                                             |
| ------------------------- | ---------------------------------------------------------------------------------------------------- |
| align_trim                | alignment post processing (amplicon assignment, softmasking, normalisation)                          |
| artic_vcf_merge           | combines VCF files from multiple read groups                                                         |
| artic_vcf_filter          | filters a combined VCF into PASS and FAIL variant files                                              |
| artic_make_depth_mask     | create a coverage mask from the post-processed alignment                                             |
| artic_mask                | combines the reference sequence, FAIL variants and coverage mask to produce a pre-consensus sequence |
| artic_fasta_header        | applies the artic workflow and identifier to the consensus sequence header                           |
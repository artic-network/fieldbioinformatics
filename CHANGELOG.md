# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.12.0] - 2026-11-19

### Added

- `artic minion` now writes the amplicon QC reports from `align_trim` 2.0.0, alongside the existing coverage reports:
    - `$SAMPLE.summary.tsv`: run summary, with fragments kept and the number removed for each reason
    - `$SAMPLE.amplicons.tsv`: what happened to the fragments of each amplicon
    - `$SAMPLE.primer_pairs.tsv`: the primer pair combinations seen
    - `$SAMPLE.split_reads.tsv`: split (supplementary) alignments and chimeras
    - `$SAMPLE.overlaps.tsv`: concordance between overlapping amplicons

### Changed

- Requires `align_trim` >= 2.0.0. `artic minion` passes `--qc-prefix` in place of the removed `--genome-coverage-report`. Until 2.0.0 is on bioconda, `environment.yml` installs it from PyPI.
- The `mean_depth` column in `$SAMPLE.amplicon_depths.tsv` is renamed `mean_fragment_depth`.
- Reads lying entirely inside the overlap of two neighbouring amplicons are now always removed, even with `--allow-mismatched-primers`, as they cannot be assigned to either amplicon.
- The changelog uses the Keep a Changelog format, and missing releases have been filled in from the git history and GitHub release notes.

### Fixed

- Clair3 no longer receives paths that contain spaces. Because Clair3 splits some of its input paths on whitespace, Clair3 >= 2.0.5 failed when the sample prefix, reference or model path contained a space, and earlier versions quietly called no variants. The per-pool BAM, the reference and the model are now placed in a temporary directory with no spaces in its path.

## [1.11.2] - 2026-08-26

### Added

- `artic_get_models` falls back to the CLIMB S3 mirror when the HKU FTP server is unreachable or times out ([#199](https://github.com/artic-network/fieldbioinformatics/pull/199)).

## [1.11.1] - 2026-07-16

### Changed

- `artic guppyplex` is about 10x faster and uses less memory. Large single FASTQ files are split so they can be processed in parallel ([#196](https://github.com/artic-network/fieldbioinformatics/pull/196)).

## [1.11.0] - 2026-06-29

### Added

- Two per-base coverage depth reports show the effect of normalisation: `$SAMPLE.pre-normalisation.coverage.tsv` and `$SAMPLE.post-normalisation.coverage.tsv`. They are useful for comparing specimens against negative controls (NTCs): at very high depth, normalisation can make an NTC and a specimen look much more alike than they are ([#194](https://github.com/artic-network/fieldbioinformatics/pull/194)).

### Changed

- Requires `align_trim` >= 1.2.0.

## [1.10.3] - 2026-06-15

### Fixed

- Fixed an off-by-one error that masked the wrong positions in the pre-consensus when applying failed variants ([#193](https://github.com/artic-network/fieldbioinformatics/pull/193)).

## [1.10.2] - 2026-06-10

### Added

- Support for the Dorado v6.0.0 HAC Clair3 model (`r1041_e82_400bps_hac_v600`).
- The basecall model can now be detected from the newer Dorado FASTQ header format, which stores it in a SAM `RG:Z:` tag ([#191](https://github.com/artic-network/fieldbioinformatics/pull/191)).

## [1.10.1] - 2026-04-28

### Added

- Variant filtering thresholds can now be set from the command line. See the [commands documentation](https://github.com/artic-network/fieldbioinformatics/blob/master/docs/commands.md#minion).
- Variants with a very low allele frequency are now discarded rather than masked.
- Tests for `artic guppyplex`.

### Changed

- The minimum required Clair3 version has been raised.

### Fixed

- `artic guppyplex` no longer fails on hidden filesystem metadata files that share a name with a FASTQ file.
- Clair3 `RefCall` records are now discarded. Previously they could be wrongly N-masked in the consensus sequence ([#189](https://github.com/artic-network/fieldbioinformatics/pull/189)).

## [1.10.0] - 2026-04-10

### Added

- Support for Clair3 2.0.0 and its PyTorch-based models.
- `artic guppyplex` runs on multiple processes and can write gzipped FASTQ.
- Docker images for `aarch64`, with an optional variant that bundles the Clair3 models.
- Better handling of read files with no reads, for example a negative control.
- Much better test coverage.

### Changed

- `artic_make_depth_mask` is much faster, especially on data that hasn't been normalised and on large genomes. Thanks to @ChrisgKent ([#185](https://github.com/artic-network/fieldbioinformatics/pull/185)).
- Documentation updated.

### Removed

- The `cyvcf2` dependency, which also makes the Docker image smaller.

### Fixed

- Paths containing spaces (for example under WSL or EPI2ME) are now handled throughout the pipeline.
- Fixed an off-by-one error in variant handling ([#187](https://github.com/artic-network/fieldbioinformatics/pull/187)).

## [1.9.1] - 2026-03-19

### Fixed

- Clair3 is pinned to < 2.0.0 after Clair3 switched from TensorFlow to PyTorch ([#183](https://github.com/artic-network/fieldbioinformatics/issues/183)).

## [1.9.0] - 2026-02-24

### Changed

- Better error handling. The pipeline no longer exits with code 20 for every failure: it returns specific exit codes and prints descriptive messages for common failures such as out-of-memory, termination, missing tools, and no reads aligning. Pipeline commands now run through `subprocess` rather than `os.system` ([#181](https://github.com/artic-network/fieldbioinformatics/pull/181)).

## [1.8.5] - 2025-10-17

### Fixed

- Deletions are now correctly identified as in-frame or out-of-frame.

## [1.8.4] - 2025-09-11

### Fixed

- Primer scheme selection hotfix: `artic_get_scheme` now writes the selected scheme to `primer.bed` and `reference.fasta`.
- The Clair3 variant filter no longer fails on records that have no `DP` field.

## [1.8.3] - 2025-09-11

### Added

- `$SAMPLE.primer.bed` and `$SAMPLE.reference.fasta` outputs record the primer scheme and reference used. This is useful when the pipeline selects the reference itself, for example with the `artic-inrb-mpox` scheme.

## [1.8.2] - 2025-09-10

### Changed

- Variant calling and filtering were tuned using EBOV data from the DRC outbreak ([#177](https://github.com/artic-network/fieldbioinformatics/pull/177)):
    - Clair3 now runs with `--haploid_sensitive` and `--enable_variant_calling_at_sequence_head_and_tail`.
    - Added an allele frequency filter.
    - Frameshift variants need a higher variant quality to pass.
- Per-pool primer-trimmed BAMs are deleted once variant calling has finished.

## [1.8.1] - 2025-08-29

### Removed

- The untrimmed BAM output, which served no purpose and caused confusion.

### Fixed

- Fixed the Clair3 variant calling command ([#175](https://github.com/artic-network/fieldbioinformatics/pull/175)).

## [1.8.0] - 2025-08-26

### Changed

- `align_trim` is now a [standalone tool](https://github.com/artic-network/align_trim) and a new dependency. If you installed from source, add it to your existing conda environment.
- Primer scheme BED parsing now uses [primalbedtools](https://github.com/ChrisgKent/primalbedtools) ([#173](https://github.com/artic-network/fieldbioinformatics/pull/173)).

### Removed

- Support for [V1 primer BED files](https://github.com/artic-network/primal-page?tab=readme-ov-file#artic-primerbed-versions), which have been deprecated for several years.
- The `pandas` dependency in `artic_mask`.

## [1.7.5] - 2025-07-17

### Changed

- Clair3 model handling has been improved.
- Packaging moved from `setup.cfg` to `pyproject.toml` ([#172](https://github.com/artic-network/fieldbioinformatics/pull/172)).

## [1.7.4] - 2025-07-15

### Changed

- `align_trim` verbose messages are clearer and include the read name wherever possible.

### Fixed

- Fixed primer trimming in `align_trim` ([artic-network/amplicon-nf#2](https://github.com/artic-network/amplicon-nf/issues/2)).
- Fixed a `NoneType` error when handling segments ([#170](https://github.com/artic-network/fieldbioinformatics/pull/170)).

## [1.7.3] - 2025-06-12

### Fixed

- Removed redundant Clair3 calls, which greatly increased runtime ([#169](https://github.com/artic-network/fieldbioinformatics/pull/169)).

## [1.7.2] - 2025-06-11

### Changed

- New consensus FASTA header format: the record ID is now the sample name, and the description gives the reference sequence and the fieldbioinformatics version (`>SAMPLE CHROM_artic-network/fieldbioinformatics_VERSION`). The variant caller is no longer included ([#168](https://github.com/artic-network/fieldbioinformatics/pull/168)).

## [1.7.1] - 2025-06-06

### Fixed

- Fixed the handling of reads in the `unmatched` read group by `artic_make_depth_mask`.

### Removed

- `artic minion` no longer writes per-read-group depth files (`--store-rg-depths`).

## [1.7.0] - 2025-06-06

### Fixed

- Depth masking for segmented viruses now covers every reference. Previously only the first reference was considered, so positions in the other segments reverted to the reference sequence ([#166](https://github.com/artic-network/fieldbioinformatics/pull/166)).

## [1.6.5] - 2025-06-06

### Added

- `--allow-mismatched-primers` option for `artic minion`. When it is set, `align_trim` keeps reads that appear to have primers from different amplicons.

## [1.6.4] - 2025-05-22

### Added

- Support for the [Dorado v5.2 models](https://github.com/nanoporetech/dorado/releases/tag/v1.0.0).

### Fixed

- `artic_get_models` now copies the models bundled with the Clair3 conda package (`r941_prom_sup_g5014` and `r941_prom_hac_g360+g422`) to a user-defined `--model-dir`.
- Model tarballs whose top-level directory name doesn't match the model name are now handled ([#165](https://github.com/artic-network/fieldbioinformatics/pull/165)).

## [1.6.3] - 2025-05-08

### Added

- Paired-end read support in `align_trim` ([#163](https://github.com/artic-network/fieldbioinformatics/pull/163)).

## [1.6.2] - 2025-03-09

### Added

- `--linearise-fasta` option to write linearised (unwrapped) FASTA output ([#160](https://github.com/artic-network/fieldbioinformatics/issues/160), [#161](https://github.com/artic-network/fieldbioinformatics/pull/161)).

## [1.6.1] - 2025-01-23

### Changed

- `--min-depth` now sets the depth threshold for the coverage mask and is no longer passed to Clair3. Its help text has been corrected.
- `--normalise` values above 10000 are rejected with a suggestion to use `--normalise 0` to turn off normalisation.

### Fixed

- Large `--normalise` values gave incorrect depths when building the coverage mask ([#155](https://github.com/artic-network/fieldbioinformatics/issues/155), [#157](https://github.com/artic-network/fieldbioinformatics/pull/157)).

## [1.6.0] - 2025-01-06

### Added

- Support for segmented viruses and multi-target amplicon panels, meaning primer schemes with more than one reference sequence ([#146](https://github.com/artic-network/fieldbioinformatics/pull/146)).

## [1.5.8] - 2024-12-10

### Added

- `seqtk` dependency.

## [1.5.7] - 2024-12-04

### Added

- For primer schemes with several references, the most appropriate reference is selected automatically. Currently this applies only to `artic-inrb-mpox` ([#147](https://github.com/artic-network/fieldbioinformatics/pull/147)).

### Changed

- Clair3 version pinned.

## [1.5.6] - 2024-11-28

### Changed

- `align_trim` prints much less by default, so its stderr can be shown while the pipeline runs. The `$SAMPLE.alignreport.er` file is gone; the per-read alignment report holds the same information in more detail.
- `align_trim` no longer pipes straight into `samtools`, so its exceptions are shown clearly, for example on a malformed primer BED file ([#144](https://github.com/artic-network/fieldbioinformatics/issues/144), [#145](https://github.com/artic-network/fieldbioinformatics/pull/145)).

### Removed

- The leftover `pyfaidx` version pin.

## [1.5.5] - 2024-11-21

### Changed

- MAFFT replaces MUSCLE for consensus alignment, which also helps ARM compatibility ([#143](https://github.com/artic-network/fieldbioinformatics/pull/143)).

## [1.5.4] - 2024-11-21

### Fixed

- Removed the explicit platform tag from the Docker image build, which caused image pull problems in Nextflow.

## [1.5.3] - 2024-11-14

### Fixed

- `artic --version`.

## [1.5.2] - 2024-11-13

### Changed

- Smaller Docker images and updated entry points.

### Fixed

- Fixed `setup.cfg`, including its Python version requirement.

## [1.5.1] - 2024-11-13

Earlier versions of this changelog listed these changes under v1.4.5. The 1.3.x and 1.4.x development versions were released together as 1.5.1 ([#137](https://github.com/artic-network/fieldbioinformatics/pull/137)).

### Added

- Primer trimming and normalisation for rapid barcoded (fragmented) data.
- Clair3 is now the default variant caller. Only the r9.4.1 models are included by default. `artic_get_models` fetches the ONT r10.4.1 models listed in the [rerio](https://github.com/nanoporetech/rerio) repository.
- The pipeline tries to pick a suitable Clair3 model from the `basecall_model_version_id` field that ONT sequencers add to read headers.
- The primer scheme fetcher now uses the [Quick lab primerschemes repository](https://github.com/quick-lab/primerschemes). Schemes not in that repository can be given directly with `--bed` and `--ref`.
- `--min-mapq`, `--min-depth` and `--min-variant-quality` options.
- Per-amplicon mean depth report (`$SAMPLE.amplicon_depths.tsv`).
- Docker images are built automatically and pushed to [quay.io](https://quay.io/repository/artic/fieldbioinformatics) for use in Nextflow pipelines and elsewhere.
- CI reimplemented as a GitHub Action.

### Changed

- The minion `--model-path` option is renamed `--model-dir`.

### Removed

- Nanopolish, because of several compatibility issues.
- Medaka, which drops long indels in a way that can't be changed.
- Longshot, which also drops long variants and isn't needed with Clair3.
- Old functionality that is no longer relevant, such as basecalling and `gather`.

### Fixed

- Overlapping variants are now normalised against the pre-consensus with `bcftools norm` ([#125](https://github.com/artic-network/fieldbioinformatics/issues/125)).

## [1.2.4] - 2023-08-22

### Changed

- Updated the pinned dependency versions, including Medaka (>= 1.8.0).

### Fixed

- Fixed the right-primer check in `vcftagprimersites.py`. Thanks to @mattloose ([#123](https://github.com/artic-network/fieldbioinformatics/pull/123)).

## [1.2.3] - 2022-09-22

### Changed

- Updated the Medaka version pin, giving access to more models.
- Longshot 0.4.5 ([#117](https://github.com/artic-network/fieldbioinformatics/issues/117)).

### Fixed

- Fixed a dependency bug that installed two forks of `pyvcf` side by side ([#119](https://github.com/artic-network/fieldbioinformatics/issues/119)).

## [1.2.2] - 2022-02-23

### Fixed

- The primer scheme downloader is reimplemented in Python. This fixes a bug in the `artic-tools` downloader with V4.1 of the ARTIC primers, whose version number is a float rather than an integer.

## [1.2.1] - 2021-01-14

### Added

- With `--strict`, the pipeline generates stats from the `align_trim` and `vcf_check` reports as a JSON file for MultiQC. This replaces the dependency on an ARTIC MultiQC module.

### Changed

- Uses upstream MultiQC instead of the ARTIC fork, which fixes pip dependency installation problems.
- Documented the missing `minion` options.

### Fixed

- `--strict` failed because an unindexed VCF was passed to `artic-tools`.

## [1.2.0] - 2021-01-12

### Added

- `--no-longshot` option for `minion`, which skips Longshot in the Medaka workflow and uses `medaka annotate` instead.
- `--strict` option for `minion`. It checks the primer scheme against ARTIC format conventions and removes merged variants that appear only once in amplicon overlaps or that fall in primer binding sites.
- `--scheme-version` option for `minion`. If a scheme isn't found locally, it is downloaded from the ARTIC repositories.
- MultiQC (ARTIC fork) in the conda environment, for amplicon coverage plots and variant quality information computed from the `align_trim` report ([#51](https://github.com/artic-network/fieldbioinformatics/issues/51)).
- [artic-tools](https://github.com/will-rowe/artic-tools) dependency.

### Changed

- Updated Medaka to v1.0.3.
- Before Medaka VCFs are filtered, a prefilter drops variants supported by fewer than 2 reads or with quality below 20, so they no longer go into the fail VCF or the consensus mask.
- The Longshot and Medaka filters in `vcf_filter` are merged into one.
- `--medaka-model` is required when `--medaka` is set.
- The pipeline splits BAMs by read group less often, so fewer BAM copies are made.
- The default scheme location is now the current directory.
- `primer.bed` schemes are handled as well as `scheme.bed` schemes.
- `align_trim` enforces correct primer pairing in both the Medaka and Nanopolish workflows.
- Repository history was cleaned to remove a large fast5 file. Anyone opening a pull request needs a fresh clone or fork.

### Deprecated

- The `plot_amplicon_depth` script, replaced by MultiQC.

### Fixed

- The input read file is checked for existence when the pipeline starts.

## [1.1.3] - 2020-06-09

### Fixed

- Merging right-hand primers produced the intersection rather than the union of alternate and canonical primers. Reported by Shirlee Wohl.

## [1.1.2] - 2020-06-03

### Fixed

- Primer parsing for schemes with no non-alt partners. Thanks to @m-bull.
- The amplicon plotting script now sorts before histogram binning, fixing the monotonic sort error with unsorted primer schemes.
- Fixed a `tqdm` import error under conda.

## [1.1.1] - 2020-05-27

Combines the two 1.1.0 release candidates with a few patches.

### Added

- `minion` produces amplicon depth plots.
- `artic_make_depth_mask` can warn about regions with low read-group coverage.

### Changed

- Masking uses pool depths.
- Test data is downloaded during the tests rather than shipped with the repository.
- Dependency version bumps, including Nanopolish 0.13.2.

### Fixed

- `artic_make_depth_mask` ignored the coverage of the first read group ([#34](https://github.com/artic-network/fieldbioinformatics/issues/34)).

## [1.1.0-rc2] - 2020-04-15

### Added

- `--no-indels` option.
- `--no-frameshifts` option for `vcf_filter`.
- `--output` option for `guppyplex`.
- Each command's execution time is recorded in the minion log.

### Changed

- Updated Medaka to 0.12.1 and Nanopolish to 0.13.0.
- `artic_make_depth_mask` uses pysam pileups in place of `samtools depth`.
- Faster Medaka steps.

### Fixed

- Longshot no longer drops variants because of the strand bias filter in the Medaka pipeline.
- Out-of-frame indels are failed and masked correctly.
- A variant that passed in one pool and failed in another was written only to the pass VCF ([#21](https://github.com/artic-network/fieldbioinformatics/issues/21)).
- Fixed the quality filter logic.

## [1.1.0-rc] - 2020-03-30

### Added

- Read groups for primer pools, taken from the BED file (for example `nCoV2019_1` and `nCoV2019_2`). They can be viewed in the Tablet viewer. Variants are called separately for each group, so mutations that disagree between groups (for example from contamination) can be detected.
- Indels supported with `bcftools consensus`.
- Longshot added to the Medaka workflow so the VCF can be filtered on depth.
- `artic guppyplex` for demultiplexed Guppy output, replacing `artic gather`.
- Support for the amplicons v3 format (Will Rowe).
- Test suite and documentation (Will Rowe).

### Changed

- New workflow model: every command runs once per sample, and `nanopolish index` is scripted.
- Updated Nanopolish to 0.12.5 for more informative VCF output, with nanopore filtering updated to use it (for example strand bias).

### Removed

- Old deprecated code (Will Rowe).

### Fixed

- `align_trim` removes invalid CIGAR strings that could make Medaka fail.

## [1.0.0] - 2018-06-09

### Changed

- Ported to Python 3 and packaged as the `artic` pipeline with setuptools entry points.
- `gather` supports sequencing summaries and Guppy output.
- `align_trim` output is more informative.

### Removed

- The ZiBRA Illumina pipeline, old shell scripts, the ZiBRA project schemes, and the Docker implementation.

## [0.13] - 2017-12-13

### Changed

- ARTIC fork.

## [0.12]

### Changed

- Updated Nanopolish to v0.8.4 for Albacore 2.0+ support. Run `nanopolish -d /path/to/fast5 input.fasta` before `zibra.py minion`.

## [0.11-dev]

### Added

- Stringent pipeline.

## [0.11] - 2017-07-25

### Added

- `align_trim_fasta.py` and `reconstitute.py` as the basis for more principled demultiplexing: reads are aligned first, then trimmed to the known amplicon end points +/- 40 bases. This helps stop chimeric reads from being assigned the wrong barcode.

### Changed

- Back on Ryan Wick's Porechop repository now that `--untrimmed` is implemented.

## [0.10] - 2017-07-15

### Added

- Higher `max-haplotypes` values can be chosen, to support more divergent references.

### Changed

- Confident deletions are now N-masked rather than ignored.
- Updated Nanopolish to HEAD, with fixes for its latest VCF format changes.

## [0.9]

### Added

- New `zibra.py` command line interface to make extraction and demultiplexing easier with Albacore and multiple groups.

## [0.8]

### Changed

- Tracks the latest Porechop, which has better adapter trimming.

## [0.7]

### Added

- Illumina wrapper scripts.

## [0.6]

### Changed

- Updated poretools.

## [0.5]

### Added

- Lassa virus scheme.

### Changed

- Updated Nanopolish to HEAD (7de633d01cc35a58e5537af5dd1024ae0040d15c) to fix Albacore 1.1 support.

## [0.4]

### Added

- Initial support for Yellow Fever (500 and 1000 bp schemes).

### Changed

- Updated Nanopolish to HEAD (aaaf90beb4efe37243c651e62fcfd11374fe2176) for Albacore 1 support.
- Extraction uses poretools with named basecallers.

## [0.3]

### Added

- Porechop.

## [0.2]

### Changed

- Updated Nanopolish to 0.6, which adds default support for various nanopore models.
- Updated samtools to 1.4.

[Unreleased]: https://github.com/artic-network/fieldbioinformatics/compare/v1.12.0...HEAD
[1.12.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.11.2...v1.12.0
[1.11.2]: https://github.com/artic-network/fieldbioinformatics/compare/v1.11.1...v1.11.2
[1.11.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.11.0...v1.11.1
[1.11.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.10.3...v1.11.0
[1.10.3]: https://github.com/artic-network/fieldbioinformatics/compare/v1.10.2...v1.10.3
[1.10.2]: https://github.com/artic-network/fieldbioinformatics/compare/v1.10.1...v1.10.2
[1.10.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.10.0...v1.10.1
[1.10.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.9.1...v1.10.0
[1.9.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.9.0...v1.9.1
[1.9.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.8.5...v1.9.0
[1.8.5]: https://github.com/artic-network/fieldbioinformatics/compare/v1.8.4...v1.8.5
[1.8.4]: https://github.com/artic-network/fieldbioinformatics/compare/v1.8.3...v1.8.4
[1.8.3]: https://github.com/artic-network/fieldbioinformatics/compare/v1.8.2...v1.8.3
[1.8.2]: https://github.com/artic-network/fieldbioinformatics/compare/v1.8.1...v1.8.2
[1.8.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.8.0...v1.8.1
[1.8.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.7.5...v1.8.0
[1.7.5]: https://github.com/artic-network/fieldbioinformatics/compare/v1.7.4...v1.7.5
[1.7.4]: https://github.com/artic-network/fieldbioinformatics/compare/v1.7.3...v1.7.4
[1.7.3]: https://github.com/artic-network/fieldbioinformatics/compare/v1.7.2...v1.7.3
[1.7.2]: https://github.com/artic-network/fieldbioinformatics/compare/v1.7.1...v1.7.2
[1.7.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.7.0...v1.7.1
[1.7.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.6.5...v1.7.0
[1.6.5]: https://github.com/artic-network/fieldbioinformatics/compare/v1.6.4...v1.6.5
[1.6.4]: https://github.com/artic-network/fieldbioinformatics/compare/v1.6.3...v1.6.4
[1.6.3]: https://github.com/artic-network/fieldbioinformatics/compare/v1.6.2...v1.6.3
[1.6.2]: https://github.com/artic-network/fieldbioinformatics/compare/v1.6.1...v1.6.2
[1.6.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.6.0...v1.6.1
[1.6.0]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.8...v1.6.0
[1.5.8]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.7...v1.5.8
[1.5.7]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.6...v1.5.7
[1.5.6]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.5...v1.5.6
[1.5.5]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.4...v1.5.5
[1.5.4]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.3...v1.5.4
[1.5.3]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.2...v1.5.3
[1.5.2]: https://github.com/artic-network/fieldbioinformatics/compare/v1.5.1...v1.5.2
[1.5.1]: https://github.com/artic-network/fieldbioinformatics/compare/v1.2.4...v1.5.1
[1.2.4]: https://github.com/artic-network/fieldbioinformatics/compare/v1.2.3...v1.2.4
[1.2.3]: https://github.com/artic-network/fieldbioinformatics/compare/1.2.2...v1.2.3
[1.2.2]: https://github.com/artic-network/fieldbioinformatics/compare/1.2.1...1.2.2
[1.2.1]: https://github.com/artic-network/fieldbioinformatics/compare/1.2.0...1.2.1
[1.2.0]: https://github.com/artic-network/fieldbioinformatics/compare/1.1.3...1.2.0
[1.1.3]: https://github.com/artic-network/fieldbioinformatics/compare/1.1.2...1.1.3
[1.1.2]: https://github.com/artic-network/fieldbioinformatics/compare/1.1.1...1.1.2
[1.1.1]: https://github.com/artic-network/fieldbioinformatics/compare/1.1.0-rc2...1.1.1
[1.1.0-rc2]: https://github.com/artic-network/fieldbioinformatics/compare/1.1.0-rc...1.1.0-rc2
[1.1.0-rc]: https://github.com/artic-network/fieldbioinformatics/compare/v1.0.0...1.1.0-rc
[1.0.0]: https://github.com/artic-network/fieldbioinformatics/compare/v0.13...v1.0.0
[0.13]: https://github.com/artic-network/fieldbioinformatics/compare/v0.11...v0.13
[0.11]: https://github.com/artic-network/fieldbioinformatics/compare/v0.10...v0.11
[0.10]: https://github.com/artic-network/fieldbioinformatics/releases/tag/v0.10

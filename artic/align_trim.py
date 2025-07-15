#!/usr/bin/env python

from copy import copy
import csv
import pysam
import sys
import numpy as np
import random
import argparse
from collections import defaultdict
from artic.utils import read_bed_file

# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]


def find_primer(bed, pos, direction, chrom, threshold=35):
    """Given a reference position and a direction of travel, walk out and find the nearest primer site.

    Parameters
    ----------
    bed : list
        A list of dictionaries, where each dictionary contains a row of bedfile data
    pos : int
        The position in the reference sequence to start from
    direction : string
        The direction to search along the reference sequence

    Returns
    -------
    tuple[int, int, dict] | bool
        A tuple containing the distance to the primer, the relative position of the primer, and the primer site, or False if no primer found
    """
    from operator import itemgetter

    if direction == "+":
        primer_distances = [
            (abs(p["start"] - pos), p["start"] - pos, p)
            for p in bed
            if (p["direction"] == direction)
            and (pos >= (p["start"] - threshold))
            and chrom == p["chrom"]
        ]

    else:
        primer_distances = [
            (abs(p["end"] - pos), p["end"] - pos, p)
            for p in bed
            if (p["direction"] == direction)
            and (pos <= (p["end"] + threshold))
            and chrom == p["chrom"]
        ]

    if not primer_distances:
        return False

    closest = min(
        primer_distances,
        key=itemgetter(0),
    )

    return closest


def trim(segment, primer_pos, end, verbose=False):
    """Soft mask an alignment to fit within primer start/end sites.

    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    vernose : bool
        If True, will print soft masking info during trimming
    """
    if verbose:
        print(
            f"Trimming segment {segment.query_name} {'end' if end else 'start'} to primer position {primer_pos}",
            file=sys.stderr,
        )
    # get a copy of the cigar tuples to work with
    cigar = copy(segment.cigartuples)

    # get the segment position in the reference (depends on if start or end of the segment is being processed)
    if not end:
        pos = segment.pos
    else:
        pos = segment.reference_end

    # process the CIGAR to determine how much softmasking is required
    eaten = 0
    while 1:

        # chomp CIGAR operations from the start/end of the CIGAR
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)
            if verbose:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
            if verbose:
                print(
                    "Ran out of cigar during soft masking - completely masked read will be ignored",
                    file=sys.stderr,
                )
            break

        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if consumesReference[flag]:
            if not end:
                pos += length
            else:
                pos -= length

        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if consumesQuery[flag]:
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= primer_pos and flag == 0:
            break
        if end and pos <= primer_pos and flag == 0:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if verbose:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if verbose:
            print("Inserted a %s, %s" % (0, extra), file=sys.stderr)
        if end:
            cigar.append((0, extra))
        else:
            cigar.insert(0, (0, extra))
        eaten -= extra

    # softmask the left primer
    if not end:

        # update the position of the leftmost mappinng base
        segment.pos = pos - extra
        if verbose:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if verbose:
                print(
                    "softmask created a leading deletion in the CIGAR, shuffling the alignment",
                    file=sys.stderr,
                )
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                segment.pos += length

        # now add the leading softmask
        cigar.insert(0, (4, eaten))

    # softmask the right primer
    else:
        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        raise ("invalid cigar operation created - possibly due to INDEL in primer")
    segment.cigartuples = cigar
    return


def handle_segment(
    segment: pysam.AlignedSegment,
    bed: dict,
    args: argparse.Namespace,
    min_mapq: int,
    report_writer: csv.DictWriter = False,
):
    """Handle the alignment segment including

    Args:
        segment (pysam.AlignedSegment): The alignment segment to process
        bed (dict): The primer scheme
        reportfh (typing.IO): The report file handle
        args (argparse.Namespace): The command line arguments

    Returns:
        tuple [int, pysam.AlignedSegment] | bool: A tuple containing the amplicon number and the alignment segment, or False if the segment is to be skipped
    """

    # filter out unmapped and supplementary alignment segments
    if segment.is_unmapped:
        if args.verbose:
            print("%s skipped as unmapped" % (segment.query_name), file=sys.stderr)
        return False
    if segment.is_supplementary:
        if args.verbose:
            print("%s skipped as supplementary" % (segment.query_name), file=sys.stderr)
        return False
    if segment.mapping_quality < min_mapq:
        if args.verbose:
            print(
                "%s skipped as mapping quality below threshold" % (segment.query_name),
                file=sys.stderr,
            )
        return False

    # locate the nearest primers to this alignment segment
    # p1 = find_primer(bed, segment.reference_start, "+", segment.reference_name, args.primer_match_threshold)
    p1 = find_primer(
        bed=bed,
        pos=segment.reference_start,
        direction="+",
        chrom=segment.reference_name,
        threshold=args.primer_match_threshold,
    )
    # p2 = find_primer(bed, segment.reference_end, "-", segment.reference_name, args.primer_match_threshold)
    p2 = find_primer(
        bed=bed,
        pos=segment.reference_end,
        direction="-",
        chrom=segment.reference_name,
        threshold=args.primer_match_threshold,
    )

    if not p1 or not p2:
        if args.verbose:
            print(
                "%s skipped as no primer found for segment" % (segment.query_name),
                file=sys.stderr,
            )
        return False

    # check if primers are correctly paired and then assign read group
    # NOTE: removed this as a function as only called once
    # TODO: will try improving this / moving it to the primer scheme processing code
    correctly_paired = (
        p1[2]["Primer_ID"].split("_")[1] == p2[2]["Primer_ID"].split("_")[1]
    )

    if not args.no_read_groups:
        if correctly_paired:
            segment.set_tag("RG", p1[2]["PoolName"])
        else:
            segment.set_tag("RG", "unmatched")

    # get the amplicon number
    amplicon = p1[2]["Primer_ID"].split("_")[1]

    if args.report:
        # update the report with this alignment segment + primer details
        report = {
            "chrom": segment.reference_name,
            "QueryName": segment.query_name,
            "ReferenceStart": segment.reference_start,
            "ReferenceEnd": segment.reference_end,
            "PrimerPair": f"{p1[2]['Primer_ID']}_{p2[2]['Primer_ID']}",
            "Primer1": p1[2]["Primer_ID"],
            "Primer1Start": p1[2]["start"],
            "Primer2": p2[2]["Primer_ID"],
            "Primer2Start": p2[2]["start"],
            "IsSecondary": segment.is_secondary,
            "IsSupplementary": segment.is_supplementary,
            "Start": p1[2]["start"],
            "End": p2[2]["end"],
            "CorrectlyPaired": correctly_paired,
        }

        report_writer.writerow(report)

    if args.remove_incorrect_pairs and not correctly_paired:
        if args.verbose:
            print(
                "%s skipped as not correctly paired" % (segment.query_name),
                file=sys.stderr,
            )
        return False

    if args.verbose:
        # Dont screw with the order of the dict
        report_str = "\t".join(str(x) for x in report.values())
        print(report_str, file=sys.stderr)

    # get the primer positions
    if args.trim_primers:
        p1_position = p1[2]["end"]
        p2_position = p2[2]["start"]
    else:
        p1_position = p1[2]["start"]
        p2_position = p2[2]["end"]

    # softmask the alignment if left primer start/end inside alignment
    if args.verbose:
        print(
            f"Trimming segment (start) {segment.query_name} refstart: {segment.reference_start} primer_end: {p1_position}",
            file=sys.stderr,
        )

    if segment.reference_start < p1_position:
        try:
            trim(segment, p1_position, False, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment.reference_start, p1_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking left primer in {} (error: {}), skipping".format(
                    segment.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # softmask the alignment if right primer start/end inside alignment
    if args.verbose:
        print(
            f"Trimming segment (end) {segment.query_name} refstart: {segment.reference_end} primer_start: {p2_position}",
            file=sys.stderr,
        )

    if segment.reference_end > p2_position:
        try:
            trim(segment, p2_position, True, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment.reference_start, p2_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking right primer in {} (error: {}), skipping".format(
                    segment.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # check the the alignment still contains bases matching the reference
    if "M" not in segment.cigarstring:
        if args.verbose:
            print(
                "%s dropped as does not match reference post masking"
                % (segment.query_name),
                file=sys.stderr,
            )
        return False

    return (amplicon, segment)


def handle_paired_segment(
    segments: tuple[pysam.AlignedSegment, pysam.AlignedSegment],
    bed: dict,
    args: argparse.Namespace,
    min_mapq: int,
    report_writer: csv.DictWriter = False,
):
    """Handle the alignment segment including

    Args:
        segment (pysam.AlignedSegment): The alignment segment to process
        bed (dict): The primer scheme
        reportfh (typing.IO): The report file handle
        args (argparse.Namespace): The command line arguments

    Returns:
        tuple [int, pysam.AlignedSegment] | bool: A tuple containing the amplicon number and the alignment segment, or False if the segment is to be skipped
    """

    segment1, segment2 = segments

    if not segment1 or not segment2:
        if args.verbose:
            print(
                "Segment pair skipped as at least one segment in pair does not exist",
                file=sys.stderr,
            )
        return False

    # filter out unmapped and supplementary alignment segments
    if segment1.is_unmapped or segment2.is_unmapped:
        if args.verbose:
            print(
                "Segment pair: %s skipped as unmapped" % (segment1.query_name),
                file=sys.stderr,
            )
        return False

    if segment1.is_supplementary or segment2.is_supplementary:
        if args.verbose:
            print(
                "Segment pair: %s skipped as supplementary" % (segment1.query_name),
                file=sys.stderr,
            )
        return False

    if segment1.mapping_quality < min_mapq or segment2.mapping_quality < min_mapq:
        if args.verbose:
            print(
                "Segment pair: %s skipped as mapping quality below threshold"
                % (segment1.query_name),
                file=sys.stderr,
            )
        return False

    # locate the nearest primers to this alignment segment
    p1 = find_primer(
        bed=bed,
        pos=segment1.reference_start,
        direction="+",
        chrom=segment1.reference_name,
        threshold=args.primer_match_threshold,
    )
    p2 = find_primer(
        bed=bed,
        pos=segment2.reference_end,
        direction="-",
        chrom=segment2.reference_name,
        threshold=args.primer_match_threshold,
    )

    if not p1 or not p2:
        if args.verbose:
            print(
                "Paired segment: %s skipped as no primer found for segment"
                % (segment1.query_name),
                file=sys.stderr,
            )
        return False

    # check if primers are correctly paired and then assign read group
    # NOTE: removed this as a function as only called once
    # TODO: will try improving this / moving it to the primer scheme processing code
    correctly_paired = p1[2]["Primer_ID"].replace("_LEFT", "") == p2[2][
        "Primer_ID"
    ].replace("_RIGHT", "")

    if not args.no_read_groups:
        if correctly_paired:
            segment1.set_tag("RG", p1[2]["PoolName"])
            segment2.set_tag("RG", p1[2]["PoolName"])
        else:
            segment1.set_tag("RG", "unmatched")
            segment2.set_tag("RG", "unmatched")

    # get the amplicon number
    amplicon = p1[2]["Primer_ID"].split("_")[1]

    if args.report:
        # update the report with this alignment segment + primer details
        report = {
            "chrom": segment1.reference_name,
            "QueryName": segment1.query_name,
            "ReferenceStart": segment1.reference_start,
            "ReferenceEnd": segment2.reference_end,
            "PrimerPair": f"{p1[2]['Primer_ID']}_{p2[2]['Primer_ID']}",
            "Primer1": p1[2]["Primer_ID"],
            "Primer1Start": p1[2]["start"],
            "Primer2": p2[2]["Primer_ID"],
            "Primer2Start": p2[2]["start"],
            "IsSecondary": segment1.is_secondary,
            "IsSupplementary": segment1.is_supplementary,
            "Start": p1[2]["start"],
            "End": p2[2]["end"],
            "CorrectlyPaired": correctly_paired,
        }

        report_writer.writerow(report)

    if args.remove_incorrect_pairs and not correctly_paired:
        if args.verbose:
            print(
                "Paired segment: %s skipped as not correctly paired"
                % (segment1.query_name),
                file=sys.stderr,
            )
        return False

    # get the primer positions
    if args.trim_primers:
        p1_position = p1[2]["end"]
        p2_position = p2[2]["start"]
    else:
        p1_position = p1[2]["start"]
        p2_position = p2[2]["end"]

    # softmask the alignment if left primer start/end inside alignment
    if args.verbose:
        print(
            f"Trimming segment (start) {segment1.query_name} refstart: {segment1.reference_start} refend: {segment1.reference_end} primer_end: {p1_position}",
            file=sys.stderr,
        )

    if segment1.reference_start < p1_position:
        try:
            trim(segment1, p1_position, False, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment1.reference_start, p1_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking left primer in {} (error: {}), skipping".format(
                    segment1.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    elif segment1.reference_end > p2_position:
        try:
            trim(segment1, p2_position, True, args.verbose)
            if args.verbose:
                print(
                    "ref end %s >= primer_position %s"
                    % (segment1.reference_end, p2_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking right primer in {} (error: {}), skipping".format(
                    segment1.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # softmask the alignment if right primer start/end inside alignment
    if args.verbose:
        print(
            f"Trimming segment (end) {segment2.query_name} refstart: {segment1.reference_start} refend: {segment1.reference_end} primer_start: {p2_position}",
            file=sys.stderr,
        )
    if segment2.reference_end > p2_position:
        try:
            trim(segment2, p2_position, True, args.verbose)
            if args.verbose:
                print(
                    "ref start %s >= primer_position %s"
                    % (segment2.reference_start, p2_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking right primer in {} (error: {}), skipping".format(
                    segment1.query_name, e
                ),
                file=sys.stderr,
            )
            return False
    elif segment2.reference_start < p1_position:
        try:
            trim(segment2, p1_position, False, args.verbose)
            if args.verbose:
                print(
                    "ref end %s >= primer_position %s"
                    % (segment2.reference_end, p1_position),
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                "problem soft masking left primer in {} (error: {}), skipping".format(
                    segment1.query_name, e
                ),
                file=sys.stderr,
            )
            return False

    # check the the alignment still contains bases matching the reference
    if "M" not in segment1.cigarstring or "M" not in segment2.cigarstring:
        if args.verbose:
            print(
                "Paired segment: %s dropped as does not match reference post masking"
                % (segment1.query_name),
                file=sys.stderr,
            )
        return False

    return (amplicon, segments)


def generate_amplicons(bed: list):
    """Generate a dictionary of amplicons from a primer scheme list (generated by vcftagprimersites/read_bed_file)

    Args:
        bed (list): A list of dictionaries, where each dictionary contains a row of bedfile data (generated by vcftagprimersites/read_bed_file), assumes that all redundant primers have been expanded

    Raises:
        ValueError: Primer direction not recognised

    Returns:
        dict: A dictionary of amplicons, where each key is the amplicon number and the value is a dictionary containing the primer start, primer end, insert start, insert end, length and circularity
    """

    amplicons = {}
    for primer in bed:

        amplicon = primer["Primer_ID"].split("_")[1]

        amplicons.setdefault(primer["chrom"], {})
        amplicons[primer["chrom"]].setdefault(amplicon, {})

        if primer["direction"] == "+":
            amplicons[primer["chrom"]][amplicon]["p_start"] = primer["start"]
            amplicons[primer["chrom"]][amplicon]["start"] = primer["end"] + 1

        elif primer["direction"] == "-":
            amplicons[primer["chrom"]][amplicon]["p_end"] = primer["end"]
            amplicons[primer["chrom"]][amplicon]["end"] = primer["start"] - 1

        else:
            raise ValueError("Primer direction not recognised")
    for chrom, amplicons_dict in amplicons.items():
        for amplicon in amplicons_dict:
            if not all([x in amplicons_dict[amplicon] for x in ["p_start", "p_end"]]):
                raise ValueError(
                    f"Primer scheme for amplicon {amplicon} for reference {chrom} is incomplete"
                )

            # Check if primer runs accross reference start / end -> circular virus
            amplicons_dict[amplicon]["circular"] = (
                amplicons_dict[amplicon]["p_start"] > amplicons_dict[amplicon]["p_end"]
            )

            # Calculate amplicon length considering that the "length" may be negative if the genome is circular
            amplicons_dict[amplicon]["length"] = abs(
                amplicons_dict[amplicon]["p_end"] - amplicons_dict[amplicon]["p_start"]
            )

    return amplicons


def normalise(trimmed_segments: dict, normalise: int, bed: list, verbose: bool = False):
    """Normalise the depth of the trimmed segments to a given value. Perform per-amplicon normalisation using numpy vector maths to determine whether the segment in question would take the depth closer to the desired depth accross the amplicon.

    Args:
        trimmed_segments (dict): Dict containing amplicon number as key and list of pysam.AlignedSegment as value
        normalise (int): Desired normalised depth
        bed (list): Primer scheme list (generated by vcftagprimersites/read_bed_file)
        trim_primers (bool): Whether to trim primers from the reads

    Raises:
        ValueError: Amplicon assigned to segment not found in primer scheme file

    Returns:
       list : List of pysam.AlignedSegment to output
    """

    amplicons = generate_amplicons(bed)

    output_segments = []

    # mean_depths = {x: {} for x in amplicons}
    mean_depths = {}

    for chrom, amplicon_dict in trimmed_segments.items():
        for amplicon, segments in amplicon_dict.items():
            if amplicon not in amplicons[chrom]:
                raise ValueError(f"Segment {amplicon} not found in primer scheme file")

            desired_depth = np.full_like(
                (amplicons[chrom][amplicon]["length"],), normalise, dtype=int
            )

            amplicon_depth = np.zeros(
                (amplicons[chrom][amplicon]["length"],), dtype=int
            )

            if not segments:
                if verbose:
                    print(
                        f"No segments assigned to amplicon {amplicon}, skipping",
                        file=sys.stderr,
                    )
                continue

            random.shuffle(segments)

            distance = np.mean(np.abs(amplicon_depth - desired_depth))

            for segment in segments:
                test_depths = np.copy(amplicon_depth)

                relative_start = (
                    segment.reference_start - amplicons[chrom][amplicon]["p_start"]
                )

                if relative_start < 0:
                    relative_start = 0

                relative_end = (
                    segment.reference_end - amplicons[chrom][amplicon]["p_start"]
                )

                test_depths[relative_start:relative_end] += 1

                test_distance = np.mean(np.abs(test_depths - desired_depth))

                if test_distance < distance:
                    amplicon_depth = test_depths
                    distance = test_distance
                    output_segments.append(segment)

            mean_depths[(chrom, amplicon)] = np.mean(amplicon_depth)

    return output_segments, mean_depths


def normalise_paired(trimmed_segments: dict, normalise: int, bed: list):
    """Normalise the depth of the trimmed segments to a given value. Perform per-amplicon normalisation using numpy vector maths to determine whether the segment in question would take the depth closer to the desired depth accross the amplicon.

    Args:
        trimmed_segments (dict): Dict containing amplicon number as key and list of tuples liek: [pysam.AlignedSegment, pysam.AlignedSegment] as value
        normalise (int): Desired normalised depth
        bed (list): Primer scheme list (generated by vcftagprimersites/read_bed_file)
        trim_primers (bool): Whether to trim primers from the reads

    Raises:
        ValueError: Amplicon assigned to segment not found in primer scheme file

    Returns:
       list : List of pysam.AlignedSegment to output
    """

    amplicons = generate_amplicons(bed)

    output_segments = []

    mean_depths = {}

    for chrom, amplicon_dict in trimmed_segments.items():
        for amplicon, segments in amplicon_dict.items():
            if amplicon not in amplicons[chrom]:
                raise ValueError(f"Segment {amplicon} not found in primer scheme file")

            desired_depth = np.full_like(
                (amplicons[chrom][amplicon]["length"],), normalise, dtype=int
            )

            amplicon_depth = np.zeros(
                (amplicons[chrom][amplicon]["length"],), dtype=int
            )

            if not segments:
                print(
                    f"No segments assigned to amplicon {amplicon}, skipping",
                    file=sys.stderr,
                )
                continue

            random.shuffle(segments)

            distance = np.mean(np.abs(amplicon_depth - desired_depth))

            for paired_segments in segments:

                test_depths = np.copy(amplicon_depth)

                segment1, segment2 = paired_segments

                for segment in (segment1, segment2):

                    relative_start = (
                        segment.reference_start - amplicons[chrom][amplicon]["p_start"]
                    )

                    if relative_start < 0:
                        relative_start = 0

                    relative_end = (
                        segment.reference_end - amplicons[chrom][amplicon]["p_start"]
                    )

                    test_depths[relative_start:relative_end] += 1

                test_distance = np.mean(np.abs(test_depths - desired_depth))

                if test_distance < distance:
                    amplicon_depth = test_depths
                    distance = test_distance
                    output_segments.append(segment1)
                    output_segments.append(segment2)

            mean_depths[(chrom, amplicon)] = np.mean(amplicon_depth)

    return output_segments, mean_depths


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def go(args):
    """Filter and soft mask an alignment file so that the alignment boundaries match the primer start and end sites.

    Based on the most likely primer position, based on the alignment coordinates.
    """
    # prepare the report outfile
    if args.report:
        reportfh = open(args.report, "w")
        report_headers = [
            "chrom",
            "QueryName",
            "ReferenceStart",
            "ReferenceEnd",
            "PrimerPair",
            "Primer1",
            "Primer1Start",
            "Primer2",
            "Primer2Start",
            "IsSecondary",
            "IsSupplementary",
            "Start",
            "End",
            "CorrectlyPaired",
        ]
        report_writer = csv.DictWriter(reportfh, fieldnames=report_headers)
        report_writer.writeheader()

    # open the primer scheme and get the pools
    bed = read_bed_file(args.bedfile)
    pools = set([row["PoolName"] for row in bed])
    chroms = set([row["chrom"] for row in bed])
    pools.add("unmatched")

    # open the input SAM file and process read groups
    infile = pysam.AlignmentFile("-", "rb")
    bam_header = infile.header.copy().to_dict()
    if not args.no_read_groups:
        bam_header["RG"] = []
        for pool in pools:
            read_group = {}
            read_group["ID"] = pool
            bam_header["RG"].append(read_group)

    # prepare the alignment outfile
    outfile = pysam.AlignmentFile("-", "wh", header=bam_header)

    trimmed_segments = {x: {} for x in chroms}

    if args.paired:
        read_pairs = read_pair_generator(infile)

        for segments in read_pairs:
            if args.report:
                trimming_tuple = handle_paired_segment(
                    segments=segments,
                    bed=bed,
                    args=args,
                    report_writer=report_writer,
                    min_mapq=args.min_mapq,
                )
            else:
                trimming_tuple = handle_paired_segment(
                    segments=segments,
                    bed=bed,
                    args=args,
                    min_mapq=args.min_mapq,
                )
            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_pair = trimming_tuple
            trimmed_segments[trimmed_pair[0].reference_name].setdefault(amplicon, [])

            if trimmed_segments:
                trimmed_segments[trimmed_pair[0].reference_name][amplicon].append(
                    trimmed_pair
                )

        if args.normalise:
            output_segments, mean_amp_depths = normalise_paired(
                trimmed_segments, args.normalise, bed
            )

            # write mean amplicon depths to file
            if args.amp_depth_report:
                with open(args.amp_depth_report, "w") as amp_depth_report_fh:
                    amp_depth_report_fh.write("chrom\tamplicon\tmean_depth\n")
                    for (chrom, amplicon), depth in mean_amp_depths.items():
                        amp_depth_report_fh.write(f"{chrom}\t{amplicon}\t{depth}\n")

            for output_segment in output_segments:
                outfile.write(output_segment)
    else:
        # iterate over the alignment segments in the input SAM file
        for segment in infile:
            if args.report:
                trimming_tuple = handle_segment(
                    segment=segment,
                    bed=bed,
                    args=args,
                    report_writer=report_writer,
                    min_mapq=args.min_mapq,
                )
            else:
                trimming_tuple = handle_segment(
                    segment=segment,
                    bed=bed,
                    args=args,
                    min_mapq=args.min_mapq,
                )
            if not trimming_tuple:
                continue

            # unpack the trimming tuple since segment passed trimming
            amplicon, trimmed_segment = trimming_tuple
            trimmed_segments[trimmed_segment.reference_name].setdefault(amplicon, [])

            if trimmed_segment:
                trimmed_segments[trimmed_segment.reference_name][amplicon].append(
                    trimmed_segment
                )

        # normalise if requested
        if args.normalise:
            output_segments, mean_amp_depths = normalise(
                trimmed_segments, args.normalise, bed, args.verbose
            )

            # write mean amplicon depths to file
            if args.amp_depth_report:
                with open(args.amp_depth_report, "w") as amp_depth_report_fh:
                    amp_depth_report_fh.write("chrom\tamplicon\tmean_depth\n")
                    for (chrom, amplicon), depth in mean_amp_depths.items():
                        amp_depth_report_fh.write(f"{chrom}\t{amplicon}\t{depth}\n")

            for output_segment in output_segments:
                outfile.write(output_segment)

        else:
            for chrom, amplicon_dict in trimmed_segments.items():
                for amplicon, segments in amplicon_dict.items():
                    for segment in segments:
                        outfile.write(segment)

    # close up the file handles
    infile.close()
    outfile.close()
    if args.report:
        reportfh.close()


def main():
    parser = argparse.ArgumentParser(
        description="Trim alignments from an amplicon scheme."
    )
    parser.add_argument("bedfile", help="BED file containing the amplicon scheme")
    parser.add_argument(
        "--normalise", type=int, help="Subsample to n coverage per strand"
    )
    parser.add_argument(
        "--min-mapq", type=int, default=20, help="Minimum mapping quality to keep"
    )
    parser.add_argument(
        "--primer-match-threshold",
        type=int,
        default=35,
        help="Fuzzy match primer positions within this threshold",
    )
    parser.add_argument("--report", type=str, help="Output report to file")
    parser.add_argument(
        "--amp-depth-report", type=str, help="Output amplicon depth tsv to path"
    )
    parser.add_argument(
        "--trim-primers",
        action="store_true",
        help="Trims primers from reads",
    )
    parser.add_argument(
        "--paired",
        action="store_true",
        help="Process paired-end reads",
    )
    parser.add_argument(
        "--no-read-groups",
        dest="no_read_groups",
        action="store_true",
        help="Do not divide reads into groups in SAM output",
    )
    parser.add_argument("--verbose", action="store_true", help="Debug mode")
    parser.add_argument("--remove-incorrect-pairs", action="store_true")

    args = parser.parse_args()

    go(args)


if __name__ == "__main__":
    main()

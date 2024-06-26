#!/usr/bin/env python

from copy import copy
import csv
import pysam
import sys
import numpy as np
import random
import typing
import argparse
from .vcftagprimersites import read_bed_file

# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]

# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]


def find_primer(bed, pos, direction) -> tuple[int, int, dict] | bool:
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
            if (p["direction"] == direction and p["start"] < pos)
        ]

        if not primer_distances:
            return False

        closest = min(
            primer_distances,
            key=itemgetter(0),
        )
    else:
        primer_distances = [
            (abs(p["end"] - pos), p["end"] - pos, p)
            for p in bed
            if (p["direction"] == direction and p["end"] > pos)
        ]

        if not primer_distances:
            return False

        closest = min(
            primer_distances,
            key=itemgetter(0),
        )

    return closest


def trim(segment, primer_pos, end, debug):
    """Soft mask an alignment to fit within primer start/end sites.

    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    debug : bool
        If True, will print soft masking info during trimming
    """
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
            if debug:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
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
    if debug:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if debug:
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
        if debug:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if debug:
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
) -> tuple[int, pysam.AlignedSegment] | bool:
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
        print("%s skipped as unmapped" % (segment.query_name), file=sys.stderr)
        return False
    if segment.is_supplementary:
        print("%s skipped as supplementary" % (segment.query_name), file=sys.stderr)
        return False
    if segment.mapping_quality < min_mapq:
        print(
            "%s skipped as mapping quality below threshold" % (segment.query_name),
            file=sys.stderr,
        )
        return False

    # locate the nearest primers to this alignment segment
    p1 = find_primer(bed, segment.reference_start, "+")
    p2 = find_primer(bed, segment.reference_end, "-")

    if not p1 or not p2:
        print(
            "%s skipped as no primer found for segment" % (segment.query_name),
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
            segment.set_tag("RG", p1[2]["PoolName"])
        else:
            segment.set_tag("RG", "unmatched")

    # get the amplicon number
    amplicon = p1[2]["Primer_ID"].split("_")[1]

    if args.report:
        # update the report with this alignment segment + primer details
        report = {
            "QueryName": segment.query_name,
            "ReferenceStart": segment.reference_start,
            "ReferenceEnd": segment.reference_end,
            "PrimerPair": f"{p1[2]['Primer_ID']}_{p2[2]['Primer_ID']}",
            "Primer1": p1[2]["Primer_ID"],
            "Primer1Start": abs(p1[1]),
            "Primer2": p2[2]["Primer_ID"],
            "Primer2Start": abs(p2[1]),
            "IsSecondary": segment.is_secondary,
            "IsSupplementary": segment.is_supplementary,
            "Start": p1[2]["start"],
            "End": p2[2]["end"],
            "CorrectlyPaired": correctly_paired,
        }

        report_writer.writerow(report)

    if args.remove_incorrect_pairs and not correctly_paired:
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
        print(
            "%s dropped as does not match reference post masking"
            % (segment.query_name),
            file=sys.stderr,
        )
        return False

    return (amplicon, segment)


def generate_amplicons(bed: list) -> dict:
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

        amplicons.setdefault(amplicon, {})

        if primer["direction"] == "+":
            amplicons[amplicon]["p_start"] = primer["start"]
            amplicons[amplicon]["start"] = primer["end"] + 1

        elif primer["direction"] == "-":
            amplicons[amplicon]["p_end"] = primer["end"]
            amplicons[amplicon]["end"] = primer["start"] - 1

        else:
            raise ValueError("Primer direction not recognised")

    for amplicon in amplicons:
        if not all([x in amplicons[amplicon] for x in ["p_start", "p_end"]]):
            raise ValueError(f"Primer scheme for amplicon {amplicon} is incomplete")

        # Check if primer runs accross reference start / end -> circular virus
        amplicons[amplicon]["circular"] = (
            amplicons[amplicon]["p_start"] > amplicons[amplicon]["p_end"]
        )

        # Calculate amplicon length considering that the "length" may be negative if the genome is circular
        amplicons[amplicon]["length"] = abs(
            amplicons[amplicon]["p_end"] - amplicons[amplicon]["p_start"]
        )

    return amplicons


def normalise(
    trimmed_segments: dict, normalise: int, bed: list, trim_primers: bool
) -> list:
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

    for amplicon, segments in trimmed_segments.items():
        if amplicon not in amplicons:
            raise ValueError(f"Segment {amplicon} not found in primer scheme file")

        desired_depth = np.full_like(
            (amplicons[amplicon]["length"],), normalise, dtype=int
        )

        if trim_primers:
            # We don't want to cover the primers in the normalisation
            desired_depth[
                amplicons[amplicon]["p_start"] : amplicons[amplicon]["start"] - 1
            ] = 0
            desired_depth[amplicons[amplicon]["end"] : amplicons[amplicon]["p_end"]] = 0

        amplicon_depth = np.zeros((amplicons[amplicon]["length"],), dtype=int)

        if not segments:
            print(
                f"No segments assigned to amplicon {amplicon}, skipping",
                file=sys.stderr,
            )
            continue

        random.shuffle(segments)

        distance = np.mean(np.abs(amplicon_depth - desired_depth))

        for segment in segments:
            test_depths = np.copy(amplicon_depth)

            relative_start = segment.reference_start - amplicons[amplicon]["start"]

            relative_end = segment.reference_end - amplicons[amplicon]["start"]

            test_depths[relative_start:relative_end] += 1

            test_distance = np.mean(np.abs(test_depths - desired_depth))

            if test_distance < distance:
                amplicon_depth = test_depths
                distance = test_distance
                output_segments.append(segment)

    return output_segments


def go(args):
    """Filter and soft mask an alignment file so that the alignment boundaries match the primer start and end sites.

    Based on the most likely primer position, based on the alignment coordinates.
    """
    # prepare the report outfile
    if args.report:
        reportfh = open(args.report, "w")
        report_headers = [
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

    trimmed_segments = {}

    # iterate over the alignment segments in the input SAM file
    for segment in infile:

        trimming_tuple = handle_segment(
            segment=segment,
            bed=bed,
            args=args,
            report_writer=report_writer,
            min_mapq=args.min_mapq,
        )
        if not trimming_tuple:
            continue

        # unpack the trimming tuple since segment passed trimming
        amplicon, trimmed_segment = trimming_tuple
        trimmed_segments.setdefault(amplicon, [])

        if trimmed_segment:
            trimmed_segments[amplicon].append(trimmed_segment)

    # normalise if requested
    if args.normalise:
        output_segments = normalise(
            trimmed_segments, args.normalise, bed, args.trim_primers
        )

        for output_segment in output_segments:
            outfile.write(output_segment)
    else:
        for amplicon, segments in trimmed_segments.items():
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
    parser.add_argument("--report", type=str, help="Output report to file")
    parser.add_argument(
        "--trim-primers",
        action="store_true",
        help="Trims primers from reads",
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

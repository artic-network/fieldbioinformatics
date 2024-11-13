# Written by Nick Loman (@pathogenomenick)

import os
import pandas as pd
from Bio import SeqIO


def run(parser, args):
    if not os.path.exists(args.output_directory):
        os.mkdir(args.output_directory)

    query_set = set()

    fqfn = "%s/%s.fastq" % (
        args.output_directory,
        args.prefix,
    )
    cmd = "samtools view -h -F 4 %s | samtools fastq - -0 %s" % (args.bamfile, fqfn)
    print(cmd)
    os.system(cmd)

    fn = "%s/%s.reads.list" % (
        args.output_directory,
        args.prefix,
    )
    fh = open(fn, "w")
    for rec in SeqIO.parse(open(fqfn), "fastq"):
        fh.write("%s\n" % (rec.id))
    fh.close()

    cmd = (
        'fast5_subset -r -i "%s" -l "%s/%s.reads.list" -f "%s_" -n 4000 -s "%s/fast5"'
        % (
            args.fast5_directory,
            args.output_directory,
            args.prefix,
            args.prefix,
            args.output_directory,
        )
    )
    os.system(cmd)

    mappingfn = "%s/fast5/filename_mapping.txt" % (args.output_directory)
    mapping_file = pd.read_csv(
        mappingfn,
        names=("read_id", "filename"),
        sep="\t",
        index_col="read_id",
        header=None,
    )

    summary_file = pd.read_csv(args.sequencing_summary, sep="\t", low_memory=False)
    filtered_summary_file = summary_file[summary_file.read_id.isin(query_set)].copy()

    filtered_summary_file["filename"] = filtered_summary_file.apply(
        lambda row: mapping_file.loc[row["read_id"]]["filename"], axis=1
    )

    summaryfn = "%s/%s_sequencing_summary.txt" % (args.output_directory, args.prefix)
    filtered_summary_file.to_csv(summaryfn, sep="\t", index=False)

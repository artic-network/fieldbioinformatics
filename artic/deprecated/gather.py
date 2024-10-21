import sys
from Bio import SeqIO
import tempfile
import os
import glob
import fnmatch
import shutil
import pandas as pd
from collections import defaultdict
import fnmatch

from . import rampart

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

# fast fastq code by Heng Li
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def write_fastq(fh, name, rec, qual):
    fh.write("@%s\n%s\n+\n%s\n" % (name, rec, qual))

def run(parser, args):
    if not args.directory:
        if not os.path.exists(args.prompt_directory):
            print ("Please specify a directory with --directory, artic gather --help for details.")
            raise SystemExit
        directories = os.listdir(args.prompt_directory)
        directories = [args.prompt_directory+'/'+d for d in directories if os.path.isdir(args.prompt_directory+'/'+d)]
        args.directory = [rampart.chooser(directories)]

    if not args.fast5_directory and not args.no_fast5s:
        print("Must supply a directory to fast5 files with --fast5-directory")
        print("If you do not want use fast5s with nanopolish use --no-fast5s instead")
        raise SystemExit(1)

    if isinstance(args.directory, list) and len(args.directory) > 1 and not args.prefix:
        print("Must supply a prefix if gathering multiple directories!", file=sys.stderr)
        raise SystemExit(1)

    if args.prefix:
        prefix = args.prefix
    else:
        prefix = os.path.split(args.directory[0])[-1]

    all_fastq_outfn = "%s_all.fastq" % (prefix)
    all_fastq_outfh = open(all_fastq_outfn, "w")

    summary_files = []

    fastq = defaultdict(list)
    for directory in args.directory:
        d = directory

        for root, dirs, files in os.walk(d):
            paths = os.path.split(root)
            barcode_directory = paths[-1]

            fastq[barcode_directory].extend([root+'/'+f for f in files if f.endswith('.fastq')])
            summary_files.extend([root+'/'+f for f in files if fnmatch.fnmatch(f, '*cing_summary*txt')])

    for barcode_directory, fastq in list(fastq.items()):
        if len(fastq):
            fastq_outfn = "%s_%s.fastq" % (prefix, barcode_directory)
            outfh = open(fastq_outfn, "w")
            print("Processing %s files in %s" % (len(fastq), barcode_directory), file=sys.stderr)

            dups = set()
            uniq = 0
            total = 0    
            limit_reached = False

            for f in fastq:
                for name, rec, qual in readfq(open(f)):
                    seq_length = len(rec)

                    if args.max_length and seq_length > args.max_length:
                        continue
                    if args.min_length and seq_length < args.min_length:
                        continue

                    total += 1
                    if name not in dups:
                        write_fastq(outfh, name, rec, qual)
                        write_fastq(all_fastq_outfh, name, rec, qual)

                        dups.add(name)
                        uniq += 1

                    if args.limit and uniq >= args.limit:
                        limit_reached = True
                        break

                if limit_reached:
                    break

            outfh.close()

            print("%s\t%s\t%s" % (fastq_outfn, total, uniq))

    all_fastq_outfh.close()

    print("Found the following summary files:\n", file=sys.stderr)
    for summaryfn in summary_files:
        print ("  " + summaryfn, file=sys.stderr)

    dfs = []

    for summaryfn in summary_files:
        df = pd.read_csv(summaryfn, sep="\t")
        # support for local basecalling
        if 'filename_fast5' in df.columns:
            df['filename'] = df['filename_fast5']
        dfs.append(df)

    summary_outfn = ""
    if dfs:
        summary_outfn = "%s_sequencing_summary.txt" % (prefix)
        summaryfh = open(summary_outfn, "w")
        pd.concat(dfs, sort=False).to_csv(summaryfh, sep="\t", index=False)
        summaryfh.close()
    else:
        print("No sequencing summary files found. This may be because the run is ongoing. You can proceed but nanopolish index will be slow and may not be able to use all of your data.")

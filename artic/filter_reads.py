import sys
from Bio import SeqIO

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID


def run(parser, args):
    dups = set()
    uniq = 0
    total = 0
    for rec in SeqIO.parse(open(args.filename), "fastq"):
        if args.max_length and len(rec) > args.max_length:
            continue
        if args.min_length and len(rec) < args.min_length:
            continue

        total += 1
        if rec.id not in dups:
            SeqIO.write([rec], sys.stdout, "fastq")

            dups.add(rec.id)
            uniq += 1

    print("%s\t%s\t%s" % (total, uniq), file=sys.stderr)

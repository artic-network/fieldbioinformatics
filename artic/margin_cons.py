#!/usr/bin/env python
from Bio import SeqIO
import sys
import vcf
import subprocess
from collections import defaultdict
import os.path
import operator
from .vcftagprimersites import read_bed_file

#MASKED_POSITIONS = [2282]
MASKED_POSITIONS = []

reference = sys.argv[1]
vcffile = sys.argv[2]
bamfile = sys.argv[3]

DEPTH_THRESHOLD = 20

def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print((bamfile, sys.stderr))

    p = subprocess.Popen(['samtools', 'depth', bamfile],
                             stdout=subprocess.PIPE)
    out, err = p.communicate()
    depths = defaultdict(dict)
    for ln in out.decode('utf-8').split("\n"):
            if ln:
                    contig, pos, depth = ln.split("\t")
                    depths[contig][int(pos)] = int(depth)
    return depths

depths = collect_depths(bamfile)

def report(r, status, allele):
    idfile = os.path.basename(vcffile).split(".")[0]
    print(("%s\t%s\tstatus\t%s" % (idfile, r.POS, status), sys.stderr))
    print(("%s\t%s\tdepth\t%s" % (idfile, r.POS, r.INFO.get('TotalReads', ['n/a'])), sys.stderr))
    print(("%s\t%s\tbasecalledfrac\t%s" % (idfile, r.POS, r.INFO.get('BaseCalledFraction', ['n/a'])), sys.stderr))
    print(("%s\t%s\tsupportfrac\t%s" % (idfile, r.POS, r.INFO.get('SupportFraction', ['n/a'])), sys.stderr))
    print(("%s\t%s\tallele\t%s" % (idfile, r.POS, allele), sys.stderr))
    print(("%s\t%s\tref\t%s" % (idfile, r.POS, r.REF), sys.stderr))

def main():
    cons = ''

    seq = list(SeqIO.parse(open(sys.argv[1]), "fasta"))[0]
    cons = list(seq.seq)

    for n, c in enumerate(cons):
        try:
            depth = depths[seq.id][n+1]
        except:
            depth = 0

        if depth < DEPTH_THRESHOLD:
            cons[n] = 'N'

    for mask in MASKED_POSITIONS:
        cons[mask-1] = 'N'

    sett = set()
    vcf_reader = vcf.Reader(open(vcffile, 'r'))
    for record in vcf_reader:
        if record.ALT[0] != '.':
            # variant call

            if record.POS in MASKED_POSITIONS:
                report(record, "masked_manual", "n")
                continue

            if 'PRIMER' in record.INFO:
                report(record, "primer_binding_site", "n")
                cons[record.POS-1] = 'N'
                continue

            support = float(record.INFO['SupportFraction'])
            total_reads = int(record.INFO['TotalReads'])
            qual = record.QUAL

            REF = record.REF
            ALT = str(record.ALT[0])

            if len(ALT) > len(REF):
                print(("Skipping insertion at position: %s" % (record.POS), sys.stderr))
                continue

            if qual >= 200 and total_reads >= 20:
                if len(REF) > len(ALT):
                    print(("N-masking confident deletion at %s" % (record.POS), sys.stderr))
                    for n in range(len(REF)):
                        cons[record.POS-1+n] = 'N'
                    continue

                report(record, "variant", ALT)
                sett.add(record.POS)
                if len(REF) > len(ALT):
                    print(("deletion", sys.stderr))
                    continue

                if len(ALT) > len(REF):
                    print(("insertion", sys.stderr))
                    continue

                cons[record.POS-1] = str(ALT)
            elif len(REF) > len(ALT):
                continue
            else:
                report(record, "low_qual_variant", "n")
                cons[record.POS-1] = 'N'
                continue    

    #print >>sys.stderr, str(sett)

    print((">%s" % (sys.argv[3])))
    print(("".join(cons)))


if __name__ == "__main__":
    main()


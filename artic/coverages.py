#!/usr/bin/env python

import vcf
import sys
import subprocess
import csv
import os
import numpy
from collections import defaultdict
from operator import attrgetter

def collect_depths(bamfile):
    if not os.path.exists(bamfile):
        raise SystemExit("bamfile %s doesn't exist" % (bamfile,))

    print(bamfile, file=sys.stderr)

    p = subprocess.Popen(['samtools', 'depth', bamfile],
                             stdout=subprocess.PIPE)
    depths = numpy.zeros(18959, dtype=numpy.int)
    out, err = p.communicate()
    for ln in out.split("\n"):
            if ln:
                    contig, pos, depth = ln.split("\t")
                    depths[int(pos)] = int(depth)
    return depths

bamfn = "EM_079517_%s_marginalign.sorted.bam" % (sys.argv[4])
depths = collect_depths(bamfn)
covered = len([a for a in depths if a >= 25])
print("UPDATE runs SET mean_cov = %s, median_cov = %s, covered = %s WHERE batch = '%s';" % (numpy.mean(depths), numpy.median(depths), covered, sys.argv[2]))

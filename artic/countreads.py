#!/usr/bin/env python

import subprocess
import sys

def get_aligned(bamfile):
    cmd = "samtools view -F 4 %s | cut -f 1  | sort | wc -l" % (bamfile)
    p = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
    return int(out)

cmd = "UPDATE runs SET num_reads_align = %s WHERE batch = '%s';" % ( 
     get_aligned("EM_079517_%s_hq_marginalign.sorted.bam" % (sys.argv[3],)),
     sys.argv[2]
)
print(cmd)


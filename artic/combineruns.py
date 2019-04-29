
import csv
import sys
from collections import defaultdict
from .runs import get_runs

runs = get_runs()

ref = sys.argv[1]
primer_scheme = sys.argv[2]

for run, samples in runs.items():
    print("cd %s" % (run,))
    print("go_zika.sh %s ../newdata/%s/downloads/pass %s" % (ref, run, primer_scheme))
    for sample, barcodes in samples.items():
        print("variants.sh %s %s %s %s" % (ref, sample, primer_scheme, " ".join(barcodes)))
    print("cd ..")

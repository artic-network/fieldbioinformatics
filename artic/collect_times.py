#!/usr/bin/env python

import sys
import subprocess
import csv
from time import gmtime, localtime, strftime

run=sys.argv[1]

def collect_times(directory):
    p = subprocess.Popen(['poretools', 'times', directory], stdout=subprocess.PIPE)
    stamps = [row['unix_timestamp'] for row in csv.DictReader(p.stdout, dialect='excel-tab')]

    return min(stamps), max(stamps), len(stamps)


# mean depth
# seq time
#  poretools times

t1 = collect_times('data/%s/pass' % (run,))
t2 = collect_times('data/%s/fail' % (run,))

with open("times/%s.times.txt" % (run,), "w") as fh:
    start_time = float(min(t1[0], t2[0]))
    end_time = float(max(t1[1], t2[1]))
    print("%s\t%s\t%s\t%s\t%s" % (
        t1[2],
        t2[2],
        strftime('%F %T', localtime(start_time)),
        strftime('%F %T', localtime(end_time)),
        end_time - start_time), file=fh)



import os
import os.path
import sys
import shutil
from collections import defaultdict
import re
from . import runs

lookup = dict([(i['Flowcell'], i) for i in runs.load_runs(sys.argv[2])])

unique = set()
flowcells = defaultdict(int)
for root, dirs, files in os.walk(sys.argv[1], topdown=False):
    for name in files:
        if name not in unique:
            m = re.search('_(FN.*?)_', name)
            if m:
                flowcells[m.group(1)] += 1
            unique.add(name)

for k, v in flowcells.items():
    if k in lookup:
        print("%s %s => %s" % (lookup[k]['Library'], k, v))
    else:
        print("No such flowcell %s" % (k,), file=sys.stderr)

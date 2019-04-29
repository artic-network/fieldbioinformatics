#!/usr/bin/env python
import sys
from Bio import SeqIO

id = 1
for ln in open(sys.argv[1]):
    if '00000000-0000-0000-0000-000000000000' in ln:
        ln = ln.replace('00000000-0000-0000-0000-000000000000', '00000000-0000-0000-0000-%012d' % (id))
        id += 1
    print(ln, end=' ')

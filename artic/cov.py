#!/usr/bin/env python

import subprocess
import sys
from tabulate import tabulate
from pandas import DataFrame
import collections
from .runs import get_runs
from operator import attrgetter
from Bio import SeqIO

class OrderedDefaultdict(collections.OrderedDict):
    """ A defaultdict with OrderedDict as its base class. """

    def __init__(self, default_factory=None, *args, **kwargs):
        if not (default_factory is None
                or isinstance(default_factory, collections.Callable)):
            raise TypeError('first argument must be callable or None')
        super(OrderedDefaultdict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory  # called by __missing__()

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key,)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else tuple()
        return self.__class__, args, None, None, iter(self.items())

    def __repr__(self):  # optional
        return '%s(%r, %r)' % (self.__class__.__name__, self.default_factory,
                               list(self.items()))

def shell(cmd):
    p = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE)
    out, err = p.communicate()
    return out

class Stat:
    def __init__(self, sample, reflen):
        cmd = "samtools view %s.sorted.bam | cut -f 1 | sort | uniq | wc -l" % (sample)
        self.reads = int(shell(cmd))

        cmd = "samtools view -F 4 %s.sorted.bam | cut -f 1 | sort | uniq | wc -l" % (sample)
        self.mapped = int(shell(cmd))

        cmd = "samtools depth %s.trimmed.sorted.bam | awk '($3>0)' | wc -l" % (sample)
        self.basescovered = int(shell(cmd))

        cmd = "samtools depth %s.trimmed.sorted.bam | awk '($3>=25)' | wc -l" % (sample)
        self.basescovered25x = int(shell(cmd))

    def hash(self):
        return collections.OrderedDict([('reads', self.reads),
                           ('mapped', self.mapped),
                           ('basescovered', self.basescovered),
                           ('basescovered25x', self.basescovered25x),
                           ('perc', "%.02d" % (100*self.basescovered25x/float(reflen)))])


reflen = len(list(SeqIO.parse(open(sys.argv[1]), "fasta"))[0])

table = []
OrderedDefaultdict(list)

runs = get_runs()
directory = '.'
for directory, samples in runs.items():
    for sample in list(samples.keys()):
        s = Stat("%s/%s" % (directory, sample), reflen)
        a = OrderedDefaultdict()
        a['run'] = directory
        a['sample'] = sample
        for k,v in s.hash().items():
            a[k] = v
        table.append(a)

headers = table[0]
print("\t".join(list(headers.keys())))
for row in table:
    print("\t".join([str(s) for s in list(row.values())]))


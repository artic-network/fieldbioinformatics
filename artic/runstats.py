#!/usr/bin/env python

import subprocess
import sys
from tabulate import tabulate
from pandas import DataFrame
import collections
from .runs import get_runs
from operator import attrgetter
from copy import copy

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
    def __init__(self, dir):
        cmd = "listdir %s | wc -l" % (dir,)
        self.uncalled = int(shell(cmd))

        cmd = "listdir %s/uploaded | wc -l" % (dir,)
        self.uploaded = int(shell(cmd))

        cmd = "find %s/downloads/pass | wc -l" % (dir,)
        self.passreads = int(shell(cmd))

        cmd = "listdir %s/downloads/fail | wc -l" % (dir,)
        self.failreads = int(shell(cmd))

    def hash(self):
        return collections.OrderedDict([('uncalled', self.uncalled),
                           ('uploaded', self.uploaded),
                           ('pass', self.passreads),
                           ('fail', self.failreads)])

table = []
OrderedDefaultdict(list)

#
#for directory in sys.argv[1:]:
#    for barcode in ['NB%02d' % (i,) for i in xrange(1,13)]:

runs = get_runs()
for directory in list(runs.keys()):
    s = Stat('newdata/'+directory)
    a = OrderedDefaultdict()
    a['directory'] = directory
    for k,v in s.hash().items():
        a[k] = v
    table.append(a)

headers = table[0]
print("\t".join(list(headers.keys())))
for row in table:
    print("\t".join([str(s) for s in list(row.values())]))

#print tabulate(table, tablefmt='pipe', headers='keys')

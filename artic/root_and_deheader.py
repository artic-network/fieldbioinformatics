#!/usr/bin/env python

from ete3 import Tree
import sys

tree = sys.argv[1]
root = sys.argv[2]

t = Tree(tree)
t.set_outgroup(t & root)

for leaf in t.iter_leaves():
    cols = leaf.name.split("|")
    if cols[0] == 'EBOV':
        leaf.name = cols[1]
    elif cols[1] == 'SLE':
        leaf.name = cols[0]

print(t.write())



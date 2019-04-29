#!/usr/bin/env python

from ete3 import Tree
import sys

tree = sys.argv[1]
root = sys.argv[2]

t = Tree(tree)
t.set_outgroup(t & root)

print(t.write())



#!/usr/bin/env python

import sys
import os.path
import glob
import re
from Bio import SeqIO

prefix = sys.argv[1]
dist = sys.argv[2]

def get_clusters(prefix, dist):
    s = set()
    clusters = glob.glob("%s_dist_%s_aligned_short.fasta-cluster*" % (prefix, dist))
    for c in clusters:
        m = re.match(r'%s_dist_%s_aligned_short.fasta-cluster(\d+)' % (prefix, dist), c)
        if m:
            s.add(int(m.group(1)))
    return sorted(list(s))

clusters = get_clusters(prefix, dist)

report = """
# Real-time Ebola surveillance in Guinea 2015

## Report date: 16th January 2016

### International context

These sequences can be viewed in international context via:

<http://ebola.nextflu.org>

### Data availability

Data is available via:

<http://github.com/nickloman/ebov>

And the EVIDENT-EU website:

<http://evident-project.eu>

### Understanding this report

For a guide on how to interpret phylogenetic trees, please refer to this useful guide by Andrew Rambaut:

<http://epidemic.bio.ed.ac.uk/how_to_read_a_phylogeny>

For questions or clarifications about how the sequencing is performed, please refer to:

<http://virological.org/t/real-time-sequencing-and-data-release-for-ebolavirus-genomic-surveillance-in-guinea/164/8>

### Contact

Please direct enquiries to Sophie Durraffour <sophieduraffour@yahoo.fr>

# Location map

![Location Map](map.png)

# Phylogenetic reconstruction

Clusters are defined with a phylogenetic distance of %s.

![Figure 1: Tree showing all Guinea surveillance isolates labeled high quality. Columns indicate a) sample identifier, b) prefecture, c) subprefecture, d) sampling date, e) cluster identifier, f) subcluster identifier](%s-bigtree.png)

Cluster 1 relates to previously described lineage 'GN1' and Cluster 2 relates to the previously described lineage 'SL3'.

Alignments of differences within clusters demonstrate only positions that vary between isolates shown in the cluster. However, the whole genome sequence is used to build the trees.

""" % (dist, prefix)

print(report)

def dump_cluster(c):
    if os.path.exists("%s_dist_%s_aligned_short.fasta-cluster%d.png" % (prefix, dist, c)):
        print("""
## Subcluster %d
        """ % (c))

        if c == 0:
            print("Cluster 0 represents isolates that do not cluster with any other isolates within the distance cut-off, i.e. singleton sequences. The sequences presented are unrelated.")

        print("""
![Subcluster %d](%s_dist_%s_aligned_short.fasta-cluster%d.pdf.png)
""" % (c, prefix, dist, c))
    else:
        print("""
## Subcluster %s

(Tree not shown for clusters with <5 isolates)

Isolates:

        """ % (c))
        for rec in SeqIO.parse(open("%s_dist_%s_aligned_short.fasta-cluster%d" % (prefix, dist, c)), "fasta"):
            print("  - %s" % (rec.id))


for c in clusters:
    if c == 0: continue
    dump_cluster(c)

dump_cluster(0)

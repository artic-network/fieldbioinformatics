#!/usr/bin/env python

import sys

reload(sys)
sys.setdefaultencoding( 'ISO8859-1' )

from ete3 import Tree, NodeStyle, TreeStyle, CircleFace, TextFace, PhyloTree, faces
from rulerface import RulerFace
import random
from ete3.treeview.faces import SequencePlotFace

import sys
import csv
import argparse

class MyTextFace(TextFace):
    def __init__(self, *args, **kwargs):
        TextFace.__init__(self, *args, **kwargs)
        self.margin_right = 10
        self.ftype = 'Arial'

class MyBorderFace(MyTextFace):
    def __init__(self, *args, **kwargs):
        MyTextFace.__init__(self, *args, **kwargs)
        self.border.width = 1
        self.margin_left = 2
        self.margin_right = 2

def reset_tree(t):
    nstyle = NodeStyle()
    nstyle['shape'] = 'sphere'
    nstyle['size'] = 0
    nstyle['fgcolor'] = 'black'

    for n in t.traverse():
        n.set_style(nstyle)

def read_clusters(fn):
    clusterdict = {}
    with open(fn) as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel-tab')
        for row in reader:
            clusterdict[row['node']] = row['group']
    return clusterdict

def read_positions(fn):
    positions = []
    with open(fn) as csvfile:
        for ln in csvfile:
            cols = ln.split("\t")
            print(cols[0])
            positions.append(int(cols[0]))
    return positions

def read_samples(csvfn):
    sampledict = {}
    with open(csvfn) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sampleid = row['id']
            sampledict[sampleid] = row
    return sampledict

def read_runs(n):
    runs = []
    with open('../runs.txt') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel-tab')
        runs = list(reader)
    return set([l['Sample'] for l in runs[-n:]])

seqcolors = {'A': '#5050ff', 'C': '#e00000', 'T': '#e6e600', 'G': '#00c000', 'N': '#b7b7b7'}
blackcolors = {'A': '#000000', 'C' : '#000000', 'T': '#000000', 'G' : '#000000', 'N': '#000000'}
runs = set()

def layout(n):
    nstyle = NodeStyle()
    nstyle['shape'] = 'sphere'
    nstyle['size'] = 0
    nstyle['fgcolor'] = 'black'
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_width"] = 2
    n.set_style(nstyle)

    clustering = None
    subclustering = None
    if args.clusters:
        clustering = read_clusters(args.clusters)
    if args.subclusters:
        subclustering = read_clusters(args.subclusters)

    if n.is_leaf():
        s = samples[n.name]
        if args.circular:
            radius=12
        else:
            radius=4
        c = CircleFace(radius=radius, color=s['prefec__colour'], style='circle')
        n.add_face(c, column=0)

        offset = 1
        if not args.hide_annotations:
            offset += 5
            tf = MyTextFace(s['id'], fsize=8, tight_text=False)
            n.add_face(tf, column=1, position="aligned")
            n.add_face(MyTextFace(s['prefec'], fsize=8, tight_text=False), column=2, position="aligned")
            n.add_face(MyTextFace(s['subprefec'], fsize=8, tight_text=False), column=3, position="aligned")
            n.add_face(MyTextFace(s['village'], fsize=8, tight_text=False), column=4, position="aligned")
            n.add_face(MyTextFace(s['date'], fsize=8, tight_text=False), column=5, position="aligned")

        if clustering:
            n.add_face(MyBorderFace(clustering[n.name], fsize=20, tight_text=False), column=offset, position="aligned")
            offset += 1

        if subclustering:
            n.add_face(MyBorderFace(subclustering[n.name], fsize=20, tight_text=False), column=offset, position="aligned")
            offset += 1

        if hasattr(n, 'sequence'):
            sf = faces.SequenceFace(n.sequence, "aa", fsize=10, fg_colors=blackcolors, bg_colors=seqcolors, col_w=11)
            n.add_face(sf, column=offset, position="aligned")
            offset += 1

def main(args):
    if args.alignment:
        t = PhyloTree(args.tree, alignment=args.alignment, alg_format='fasta')
    else:
        t = PhyloTree(args.tree)

    if args.highlight_new:
        runs = read_runs(args.highlight_new)

    t.set_outgroup('EM_079422')
    t.ladderize()

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.layout_fn = layout

    thick_hz_line = NodeStyle()
    thick_hz_line["hz_line_width"] = 8
    t.set_style(thick_hz_line)
    #t.children[0].set_style(thick_hz_line)
    #t.children[1].set_style(thick_hz_line)

    thick_vt_line = NodeStyle()
    thick_vt_line["vt_line_width"] = 4
    t.set_style(thick_vt_line)

    # header
    if not args.hide_annotations:
        ts.aligned_header.add_face(MyTextFace('Sample identifier', fstyle='Bold', fsize=8, tight_text=False), column = 1)
        ts.aligned_header.add_face(MyTextFace('Prefecture', fstyle='Bold', fsize=8, tight_text=False), column = 2)
        ts.aligned_header.add_face(MyTextFace('Sous-prefecture', fstyle='Bold', fsize=8, tight_text=False), column = 3)
        ts.aligned_header.add_face(MyTextFace('Village', fstyle='Bold', fsize=8, tight_text=False), column = 4)
        ts.aligned_header.add_face(MyTextFace('Sample received', fstyle='Bold', fsize=8, tight_text=False), column = 5)

    if args.positions:
        positions = read_positions(args.positions)

        alg_header = RulerFace(positions,
                              col_width=11,
                              height=0, # set to 0 if dont want to use values
                              kind="stick",
                              hlines = [0],
                              hlines_col = ["white"], # trick to hide hz line
                              )

        ts.aligned_header.add_face(alg_header, 6)

    #legend
    if args.legend:
        legend = {}
        for s in list(samples.values()):
            legend[s['prefec']] = s['prefec__colour']
        for p in sorted(legend.keys()):
            ts.legend.add_face(CircleFace(4, legend[p]), column=0)
            ts.legend.add_face(MyTextFace(p, fsize=6, tight_text=False), column=1)    
        ts.legend_position=1

    if args.circular:
        ts.mode = "c"
        ts.arc_start = -180 # 0 degrees = 3 o'clock
        ts.arc_span = 180

#    t.show(tree_style=ts)
    t.render(args.output, tree_style=ts, w=1024)

parser = argparse.ArgumentParser(description='Render a tree.')
parser.add_argument('tree')
parser.add_argument('metadata')
parser.add_argument('output')
parser.add_argument('--alignment')
parser.add_argument('--circular', action='store_true')
parser.add_argument('--hide-annotations', action='store_true')
parser.add_argument('--highlight-new', type=int)
parser.add_argument('--clusters')
parser.add_argument('--subclusters')
parser.add_argument('--positions')
parser.add_argument('--legend', action='store_true')

args = parser.parse_args()
samples = read_samples(args.metadata)
main(args)


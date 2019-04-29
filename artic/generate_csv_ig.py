#!/usr/bin/env python

import sys
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
import re

def clean(s):
    s = s.decode("utf-8").replace("\u00E9", "e").encode("utf-8")
    return s

public = 1
if sys.argv[3] == 'private':
    public = 0

#codes = {'GUI': 'Guinea', 'SLe': 'Sierra Leone', 'LIB': 'Liberia', '?': 'Unknown'}

locations = {}
for line in open(sys.argv[1], 'r'):
    cols = line.strip().split(',')
    cols[0] = clean(cols[0])
    locations[cols[0]] = (float(cols[1]), float(cols[2]))


colours = ['#E31A1C', '#33A02C', '#1F78B4', '#B15928', '#6A3D9A', '#FF7F00', '#FB9A99', '#B2DF8A', '#A6CEE3', '#FFFF99', '#CAB2D6', '#FDBF6F']
#colours = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928']

lines = [record.description.strip().split('|') for record in SeqIO.parse(open(sys.argv[2], 'r'), 'fasta')]
counts = defaultdict(int)

newlist = []
for line in lines:
    if line[1] == 'SLE':
        m = re.search('([A-Z].+)([A-Z].+)', line[2])
        if m:
            line.append(m.group(1) + ' ' + m.group(2))
        else:
            line.append(line[2])
        newlist.append(['EBOV', line[0], 'Goodfellow', 'SLE', line[2] + '--', line[3], line[4]])
    else:
        print(line, file=sys.stderr)
        locstring = line[4].split('-')
        if locstring[0] in list(locations.keys()):
            if '_'.join([locstring[0], locstring[1]]) in list(locations.keys()):
                line.append('_'.join([locstring[0], locstring[1]]))
            line.append(locstring[0])
        else:
            raise ValueError('Valid prefecture not found %s' %line[4])
        newlist.append(line)

for line in newlist:
    counts[line[6]] += 1

colourDict= {}
for n, (key, value) in enumerate(sorted(list(counts.items()), key=itemgetter(1), reverse=True)):
    colourDict[key] = colours[n]

header = 'tree_id,id,__latitude,__longitude,prefec,prefec__shape,prefec__colour,date,__day,__month,__year,all'
if public == 0:
    header = header + ',subprefec,village'
sys.stdout.write(header + '\r\n')

for line in newlist:
    # latitude needs to be negative
    locstring = line[4].split('-')
    infofield = "|".join(line[:-1])
    datecols = line[5].split('-')
    if line[6] in locations:
        prefec = line[6]
    else:
        prefec == 'Unknown'
    shape = 'square'
    if line[6] == 'Boke_Kamsar':
        shape = 'triangle'
    string = ','.join((line[1], line[1], '%.7f' %(locations[prefec][0]),
            '%.7f' %(locations[prefec][1]), prefec, shape, colourDict[prefec], \
            line[5], datecols[2], datecols[1], datecols[0], infofield))
    if public == 0:
        string = string + ','.join(('', locstring[1], locstring[2]))
    sys.stdout.write(string + '\r\n') # this may change on other OS

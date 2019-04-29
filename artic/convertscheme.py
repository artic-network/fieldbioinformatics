import sys

for ln in open(sys.argv[1]):
    cols = ln.rstrip().split("\t")

    if 'LEFT' in cols[3]:
        direction = '+'
    else:
        direction = '-'

    a,pair,b = cols[3].split('_')

    print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (cols[0], cols[1], cols[2], cols[3], 0, direction, pair))



#!/usr/bin/env python3

"""
Plot the mean read depth per amplicon.

This has been written for use in the ARTIC pipeline so there are no file checks - it assumes the following:
 * the primer scheme is in ARTIC format
 * the input depth files are in the format: `chrom\treadgroup\tposition\tdepth
 * readgroup equates to primer pool
 * the primer pairs in the scheme are sorted by amplicon number (i.e. readgroups are interleaved)
 * depth values are provided for all positions (see output of make_depth_mask.py for expected format)

"""
from .vcftagprimersites import read_bed_file
import sys
import pandas as pd
import numpy as np
import argparse
import os


os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def go(args):

    # get the primer scheme
    primerScheme = read_bed_file(args.primerScheme)

    # number the amplicons in the scheme and link them to primer start site
    ampliconCounter = 1

    # store the amplicon number and starts by read group dict
    rgAmplicons = {}
    rgStarts = {}

    # process the primers by readgroup
    for primer in primerScheme:
        poolName = primer['PoolName']
        if poolName not in rgAmplicons:
            rgAmplicons[poolName] = []
            rgStarts[poolName] = []
        if primer['direction'] == '+':
            rgAmplicons[poolName].append(ampliconCounter)
            rgStarts[poolName].append(primer['start'])
            ampliconCounter += 1

    # for pandas cut func to create bins, we need to add an extra value to the starts (just use inf)
    for startList in rgStarts.values():
        startList.append(np.inf)

    # process the depth files
    dfs = {}
    for depthFile in args.depthFiles:

        # read in the depth file
        df = pd.read_csv(depthFile, sep='\t', header=None,
                            names=['refName', 'readGroup',
                                 'position', 'depth'],
                            dtype={'refName': str, 'readGroup': str,
                                 'position': int, 'depth': int},
                            usecols=(0, 1, 2, 3),)

        # check that there aren't too many positions in the depth data for plotting
        # assert len(df.index) < 30000, "error: too many data points to plot"

        # check all ref positions have a depth value
        startPos = df["position"][0]
        endPos = df["position"][df.index[-1]]
        assert len(df.index) == ((endPos - startPos) +
                                 1), "error: depth needs to be reported at all positions"

        # check the primer scheme contains the readgroup
        rgList = df.readGroup.unique()
        assert len(rgList) == 1, "error: depth file has %d readgroups, need 1 (%s)" % (
            len(rgList), depthFile)
        rg = rgList[0]
        assert rg in rgAmplicons, "error: readgroup not found in provided primer scheme (%s)" % (
            rg)

        # get the amplicon starts for this readgroup
        amplicons = sorted(rgAmplicons[rg])
        starts = sorted(rgStarts[rg])

        # bin read depths by amplicon for this readgroup
        df['amplicon'] = pd.cut(
            x=df['position'], bins=starts, labels=amplicons)

        # store the mean of each bin
        bins = (df.groupby(['amplicon'])[
                'depth'].mean()).rename(depthFile.name)

        # add to the pile
        assert rg not in dfs, "error: readgroup present in multiple files (%s)" % (
            rg)
        dfs[rg] = bins

    # combine the series data from each input file
    newDF = pd.concat(dfs, axis=1)
    newDF.sort_index(axis=0, inplace=True)
    newDF.reset_index(inplace=True)

    # melt the DF for seaborn
    newDF = newDF.melt("amplicon", var_name="read group",
                       value_name="mean amplicon read depth")
    newDF = newDF.dropna()

    # plot the bar
    g = sns.catplot(data=newDF,
                        x="amplicon",
                        y="mean amplicon read depth",
                        hue="read group",
                        height=4,
                        aspect=3,
                        kind="bar",
                        dodge=False,
                        legend=False)
    g.set(yscale="log")
    g.fig.suptitle(args.sampleID)
    plt.legend(loc='upper right')
    plt.xticks(rotation=45, size=6)
    plt.savefig(args.outFilePrefix + "-barplot.png")
    plt.close()

    # plot the box
    g = sns.catplot(data=newDF,
                    x="read group",
                    y="mean amplicon read depth",
                    kind="box")
    g.fig.suptitle(args.sampleID)
    plt.savefig(args.outFilePrefix + "-boxplot.png")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--primerScheme', required=True,
                        help='the ARTIC primer scheme')
    parser.add_argument('--sampleID', required=True,
                        help='the sample ID for the provided depth files')
    parser.add_argument('--outFilePrefix', default="./amplicon-depth",
                        help='the prefix to give the output plot file')
    parser.add_argument(
        "depthFiles", type=argparse.FileType('r'), nargs='+', help='the depth files produced by make_depth_mask.py')
    args = parser.parse_args()
    go(args)


if __name__ == "__main__":
    main()

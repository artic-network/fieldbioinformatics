from Bio import SeqIO
import sys
import numpy

for record in SeqIO.parse(sys.argv[1], "fastq"):
    print(numpy.mean(record.letter_annotations["phred_quality"]))



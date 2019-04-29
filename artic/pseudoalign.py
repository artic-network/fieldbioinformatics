
from Bio import SeqIO, Seq
import sys

for rec in SeqIO.parse(open(sys.argv[1]), "fasta"):
    seqstr = str(rec.seq)
    seqstr = '-------------------------------------------' + seqstr + '-------------------------------------------------------------------------------------------------------------'
    rec.seq = Seq.Seq(seqstr)
    SeqIO.write([rec], sys.stdout, "fasta")


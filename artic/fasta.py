from poretools.Fast5File import Fast5FileSet
import sys

def extract_fast5(path, basecaller, flowcell_id):
    reads = set()

    for fast5 in Fast5FileSet(path, 0, basecaller):
        read_flowcell_id= fast5.get_flowcell_id()
        if flowcell_id != read_flowcell_id:
            print(("Skipping read from flowcell: %s" % (read_flowcell_id), sys.stderr))
            continue

        read_id = fast5.get_read_id()
        if read_id in reads:
            print(("Skipping duplicate read: %s" % (read_id), sys.stderr))
            continue

        reads.add(read_id)

        fas = fast5.get_fastas('fwd')
        for read in fas:
            print(read)
        fast5.close()

extract_fast5(sys.argv[1], 'ONT Albacore Sequencing Software=1.0.4', sys.argv[2])

# zibra.py 
#  run
#   --flowcell
#   --type 1d / 2d
#   --check-sample-name
#   --check-flowcell-name
#   --min-support-value
#   --min-depth
#   --min-log-likelihood
#   --normalised-depth
#   --use-indels
#   --trim-reads
#   <scheme> <sample> <directory>
#  list-schemes


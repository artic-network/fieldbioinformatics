from poretools.Fast5File import Fast5FileSet
import sys

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

def run(parser, args):
	flowcells = set()
	reads = set()
	i = 0
	basecaller_version = None

	for fast5 in Fast5FileSet(args.directory, None, args.basecaller):
#		if not basecaller_version:
#			basecaller_version = fast5.get_basecaller_version()
#		elif fast5.get_basecaller_version() != basecaller_version:
#			print >>sys.stderr, "ABORTED: More than one basecaller version found: %s, %s" % (basecaller_version, fast5.get_basecaller_version())
#			raise SystemExit
				

		if not fast5.is_open:
			print("Skipping read: %s" % (fast5.filename), file=sys.stderr)
			continue

		read_flowcell_id = fast5.get_flowcell_id()
		flowcells.add(read_flowcell_id)
		if len(flowcells) != 1:
			print("ABORTED: More than one flowcell found in dataset: %s" % (flowcells,), file=sys.stderr)
			raise SystemExit(1)

		#if flowcell_id != read_flowcell_id:
		#	print >>sys.stderr, "Skipping read from flowcell: %s" % (read_flowcell_id)
		#	continue

		read_id = fast5.get_read_id()
		if read_id in reads:
			print("Skipping duplicate read: %s" % (read_id), file=sys.stderr)
			continue

		reads.add(read_id)

		fas = fast5.get_fastas('fwd')
		for read in fas:
			if read:
				print(read)
		fast5.close()

		i += 1

		if i % 1000 == 0:
			print("Extracted %s reads" % (i,), file=sys.stderr)

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


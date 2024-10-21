from poretools.Fast5File import Fast5FileSet
import sys
import re
from collections import defaultdict

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

matcher = re.compile('Basecall_1D_(\d+)')

def get_basecaller_version(g):
	try:
		return g.attrs['chimaera version']
	except:
		pass

	try:
		return g.attrs['version']
	except:
		return None

def get_basecallers(fast5):
	basecallers = []
	analyses = fast5.hdf5file.get('Analyses')
	if analyses:
		for k, g in analyses.items():
			m = matcher.match(k)
			if m:
				basecaller_name = g.attrs['name']
				group = m.group(1)
				version = get_basecaller_version(g)

				basecallers.append([basecaller_name, group, version])
	return basecallers

def print_basecallers(basecallers):
	for key, item in basecallers.items():
		print("%s: %d reads" % (key, item))

def run(parser, args):
	flowcells = set()
	reads = set()
	basecallers = {}
	i = 0

	for fast5 in Fast5FileSet(args.directory):
		#, 0, basecaller):
		if not fast5.is_open:
			print("Skipping read: %s" % (fast5.filename), file=sys.stderr)
			continue

		bcs = get_basecallers(fast5)
		for bc in bcs:
			bcstring = "%s=%s" % (bc[0], bc[2])
			if bcstring not in basecallers:
				basecallers[bcstring] = 1
				print_basecallers(basecallers)
			else:
				basecallers[bcstring] += 1

		i += 1
		if i % 1000 == 0:
			print_basecallers(basecallers)

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


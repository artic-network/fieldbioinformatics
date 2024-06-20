import sys
import tempfile
import os
import shutil

from . import rampart

# extract with constraints:
#   -- only one group ever
#   -- only one flowcell ID ever
#   -- always unique read ID

def run(parser, args):
	tmpdir = tempfile.mkdtemp(dir='.')

	cmd = ("porechop --verbosity 2 --untrimmed -i \"%s\" -b %s --native_barcodes --discard_middle --require_two_barcodes --barcode_threshold 80 --threads %s --check_reads 10000 --barcode_diff 5 > %s.demultiplexreport.txt" % (args.fasta, tmpdir, args.threads, args.fasta))
	print(cmd, file=sys.stderr)
	os.system(cmd)

	a, b = os.path.split(args.fasta)
	prefix, ext = os.path.splitext(b)

	for fn in os.listdir(tmpdir):
		newfn = "%s-%s" % (prefix, os.path.basename(fn))
		shutil.move(tmpdir + '/' + fn, newfn)

		if newfn.endswith('.gz'):
			os.system("gunzip -f %s" % (newfn,))

		# build nanopolish index files for the demultiplexed reads if a readdb file exists for the input file
		master_readdb_fn = "%s.index.readdb" % (args.fasta)
		if os.path.exists(master_readdb_fn):
			# first we build the .fai files intentionally WITHOUT the fast5 directory argument
			# this will print an error that we ignore
			cmd = ("nanopolish index %s 2>/dev/null" % newfn)
			print(cmd, file=sys.stderr)
			os.system(cmd)

			new_readdb = "%s.index.readdb" % (newfn)
			os.symlink(master_readdb_fn, new_readdb)

	if not args.no_remove_directory:
		os.rmdir(tmpdir)


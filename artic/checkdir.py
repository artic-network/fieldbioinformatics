
import os.path
import sys

passexists = os.path.exists('data/%s/pass' % (sys.argv[1]))
failexists = os.path.exists('data/%s/fail' % (sys.argv[1]))

print("%s\t%s\t%s" % (sys.argv[1], passexists, failexists))

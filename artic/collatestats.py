import sys

headerprinted = False

for fn in sys.argv[1:]:
    fh = open(fn)
    headers = fh.readline()
    if not headerprinted:
        print("filename\t%s" % (headers), end=' ')
        headerprinted = True


    for ln in fh:
        print("%s\t%s" % (fn, ln), end=' ')

    


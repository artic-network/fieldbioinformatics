import os
import os.path
import sys
import shutil

input_dir = sys.argv[1]
albacore_dir = sys.argv[2]
process_dir = sys.argv[3]

basecalled_files = set()

for root, dirs, files in os.walk(albacore_dir, topdown=False):
    for name in files:
        basecalled_files.add(name)

for root, dirs, files in os.walk(input_dir, topdown=False):
        for name in files:
            if name in basecalled_files:
                albacore_root = root[len(input_dir):]
                # move it
                checkdir = process_dir + '/' + albacore_root
                if not os.path.exists(checkdir):
                    os.makedirs(checkdir)
                movefrom = input_dir + '/' + albacore_root + '/' + name
                moveto = process_dir + '/' + albacore_root + '/' + name
                print("Move %s to %s" % (movefrom, moveto))
                shutil.move(movefrom, moveto)





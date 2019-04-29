import os
import os.path
import sys
import shutil

input_dir = sys.argv[1]
output_dir = sys.argv[2]
matchpattern = sys.argv[3]

basecalled_files = set()

for root, dirs, files in os.walk(output_dir, topdown=False):
    for name in files:
        basecalled_files.add(name)

for root, dirs, files in os.walk(input_dir, topdown=False):
        for name in files:
            if name not in basecalled_files and matchpattern in name:
                albacore_root = root[len(input_dir):]
                # move it
                checkdir = output_dir + '/' + albacore_root
                if not os.path.exists(checkdir):
                    os.makedirs(checkdir)
                movefrom = input_dir + '/' + albacore_root + '/' + name
                moveto = output_dir + '/' + albacore_root + '/' + name
                print("Move %s to %s" % (movefrom, moveto))
                shutil.move(movefrom, moveto)





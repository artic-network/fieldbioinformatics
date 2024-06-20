# Written by Nick Loman (@pathogenomenick)

import os
import sys
from Bio import SeqIO
from clint.textui import colored, puts, indent


def chooser(directories):
    print("Found the following directories:")
    for n, d in enumerate(directories):
        print("  [%d]: %s" % (n + 1, d))
    print("Choose 0 to quit")

    while True:
        choice = input("=> ")
        try:
            choice_number = int(choice)
            if choice_number == 0:
                raise SystemExit(1)
            path = directories[choice_number - 1]
            break
        except Exception:
            print("Invalid choice, please select from the list above.")

    return path


def run(parser, args):
    directories = []

    for root, dirs, files in os.walk(args.run_directory, topdown=False):
        for d in dirs:
            if d == "fastq_pass":
                directories.append(root + "/" + d)

    basecalledPath = chooser(directories)

    directories = os.listdir(args.protocol_directory)
    directories = [
        d for d in directories if os.path.isdir(args.protocol_directory + "/" + d)
    ]

    protocolPath = chooser(directories)

    print(basecalledPath, protocolPath)

    skip_barcoding = False
    if os.path.exists("run_configuration.json"):
        reenter = input("Do you want to enter sample names again? (Y/N): ")
        if reenter.lower().startswith("n"):
            skip_barcoding = True

    if not skip_barcoding:
        barcodes = []
        print("Enter sample names:")
        for barcode in range(1, 25):
            print(" NB%02d: " % (barcode,), end="")
            barcodes.append(input("> "))

        fh = open("run_configuration.json", "w")
        fh.write(
            """{
	  \"samples\": [
		"""
        )
        first = True
        for n, b in enumerate(barcodes):
            if b:
                if not first:
                    fh.write(",\n")
                first = False
                fh.write(
                    """{
	      "name": "%s",
	      "description": "",
	      "barcodes": [ "NB%02d" ]
    }"""
                    % (barcodes[n].replace('"', '\\"'), n + 1)
                )
        fh.write("] }")
        fh.close()

    cmd = "rampart --basecalledPath %s --protocol %s/%s --clearAnnotated" % (
        basecalledPath,
        args.protocol_directory,
        protocolPath,
    )
    print(cmd)

    os.system(cmd)

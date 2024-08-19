from Bio import SeqIO


def fasta_header(fn, header):
    fh = open(fn)
    rec = list(SeqIO.parse(fh, "fasta"))
    if len(rec) > 1:
        print("Sorry, this script doesn't support multi-FASTA files!")
        raise SystemExit
    fh.close()

    fh = open(fn, "w")
    rec[0].id = header
    SeqIO.write([rec[0]], fh, "fasta")
    fh.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Trim alignments from an amplicon scheme."
    )
    parser.add_argument("filename")
    parser.add_argument("header")

    args = parser.parse_args()
    fasta_header(args.filename, args.header)


if __name__ == "__main__":
    main()

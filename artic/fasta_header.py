from Bio import SeqIO


def fasta_header(args):
    with open(args.filename) as fh:
        rec = list(SeqIO.parse(fh, "fasta"))

    fasta_format = "fasta-2line" if args.linearise_fasta else "fasta"

    with open(args.filename, "w") as fh:
        for record in rec:
            chrom = record.id
            record.id = f"{args.samplename}/{chrom}/ARTIC/{args.caller}"

            SeqIO.write(record, fh, fasta_format)


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Trim alignments from an amplicon scheme."
    )
    parser.add_argument("filename")
    parser.add_argument("samplename")
    parser.add_argument("caller")
    parser.add_argument("--linearise-fasta", action="store_true")

    args = parser.parse_args()
    fasta_header(args)


if __name__ == "__main__":
    main()

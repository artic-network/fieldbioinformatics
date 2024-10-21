import csv


def load_runs(fn):
    rows = []
    skipcomments = (row for row in open(fn) if not row.startswith("#"))
    for row in csv.DictReader(skipcomments, dialect="excel-tab"):
        if "Included" not in row:
            rows.append(row)
        elif int(row["Included"]):
            rows.append(row)

    return rows

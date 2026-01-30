#!/usr/bin/env python3
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--intron", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-p", "--prefix", required=True)
    args = ap.parse_args()

    counts = {}

    with open(args.intron) as f:
        for line in f:
            if line.strip():
                _, _, _, _, _, gene_cds = line.split()
                gene, exon = gene_cds.rsplit(".", 1)
                exon = int(exon)
                counts[gene] = max(counts.get(gene, 0), exon)

    out = open(args.output, "w")
    for g, c in counts.items():
        out.write(f"{g}\t{c}\t{args.prefix}\n")

    out.close()

if __name__ == "__main__":
    main()

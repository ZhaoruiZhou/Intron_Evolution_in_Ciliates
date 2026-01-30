#!/usr/bin/env python3
import argparse

def genome_size(path):
    total = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--genome", required=True)
    ap.add_argument("-i", "--intron", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    gsize = genome_size(args.genome)

    intron_bases = 0
    intron_count = 0

    with open(args.intron) as f:
        prev = {}
        for line in f:
            if line.strip():
                sc, _, s, e, strand, gene_cds = line.split()
                gene = gene_cds.rsplit(".", 1)[0]
                s = int(s)
                e = int(e)

                if gene in prev:
                    pe = prev[gene]
                    if s > pe + 1:
                        intron_bases += (s - pe - 1)
                        intron_count += 1

                prev[gene] = e

    out = open(args.output, "w")

    if gsize > 0:
        out.write(f"Intron bases / Genome size: {(intron_bases/gsize)*100:.4f}%\n")

    out.write(f"Intron count: {intron_count}\n")

    out.close()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--intron", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    stats = {0:0, 1:0, 2:0}
    prev = {}

    with open(args.intron) as f:
        for line in f:
            if line.strip():
                _, _, s, e, _, gene_cds = line.split()
                gene = gene_cds.rsplit(".", 1)[0]
                s = int(s)
                e = int(e)

                if gene in prev:
                    pe = prev[gene]
                    if s > pe + 1:
                        L = s - pe - 1
                        stats[L % 3] += 1

                prev[gene] = e

    total = sum(stats.values())

    out = open(args.output, "w")

    if total > 0:
        out.write(f"3n: {stats[0]/total*100:.2f}%\n")
        out.write(f"3n+1: {stats[1]/total*100:.2f}%\n")
        out.write(f"3n-1: {stats[2]/total*100:.2f}%\n")

    out.close()

if __name__ == "__main__":
    main()

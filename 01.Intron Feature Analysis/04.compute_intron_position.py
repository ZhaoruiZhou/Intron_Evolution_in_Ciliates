#!/usr/bin/env python3
import argparse

def parse_introns(path):
    cds = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                sc, _, s, e, strand, gene_cds = line.split()
                gene = gene_cds.rsplit(".", 1)[0]
                cds.setdefault(gene, []).append((int(s), int(e)))
    for g in cds:
        cds[g].sort()
    return cds

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--intron", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    cds = parse_introns(args.intron)
    out = open(args.output, "w")

    for gene, regions in cds.items():
        total_len = sum(e - s + 1 for s, e in regions)
        cum = 0
        positions = []

        for i in range(len(regions) - 1):
            s1, e1 = regions[i]
            s2, _ = regions[i + 1]
            cum += e1 - s1 + 1

            if s2 > e1 + 1:
                positions.append(cum / total_len)

        if positions:
            out.write(f"{gene}\t" + " ".join(f"{p:.3f}" for p in positions) + "\n")
        else:
            out.write(f"{gene}\tNone\n")

    out.close()

if __name__ == "__main__":
    main()

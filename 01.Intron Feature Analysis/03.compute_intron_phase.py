#!/usr/bin/env python3
import argparse

def parse_introns(path):
    cds = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                sc, _, s, e, strand, gene_cds = line.split()
                gene = gene_cds.rsplit(".", 1)[0]
                cds.setdefault(gene, []).append((int(s), int(e), strand))
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
    all_phases = []

    for gene, regions in cds.items():
        phases = []
        for i in range(len(regions) - 1):
            s1, e1, strand = regions[i]
            s2, e2, _ = regions[i + 1]

            if s2 > e1 + 1:
                phase = (e1 - s1 + 1) % 3
                phases.append(phase)
                all_phases.append(phase)

        if phases:
            out.write(f"{gene}\t" + " ".join(map(str, phases)) + "\n")
        else:
            out.write(f"{gene}\tNone\n")

    if all_phases:
        total = len(all_phases)
        for p in (0,1,2):
            pct = all_phases.count(p) / total * 100
            out.write(f"\nPhase {p}: {pct:.2f}%\n")

    out.close()

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse

def parse_fasta(path):
    genome = {}
    with open(path) as f:
        sid = None
        seq = []
        for line in f:
            if line.startswith(">"):
                if sid:
                    genome[sid] = "".join(seq)
                sid = line.strip().split()[0][1:]
                seq = []
            else:
                seq.append(line.strip())
        if sid:
            genome[sid] = "".join(seq)
    return genome

def parse_intron_file(path):
    cds = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                scaffold, _, start, end, strand, gene_cds = line.split()
                gene = gene_cds.rsplit(".", 1)[0]
                cds.setdefault(gene, []).append((scaffold, int(start), int(end), strand))
    for g in cds:
        cds[g].sort(key=lambda x: x[1])
    return cds

def revcomp(seq):
    comp = dict(A="T", T="A", C="G", G="C", N="N")
    return "".join(comp[b] for b in reversed(seq))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-g", "--genome", required=True)
    ap.add_argument("-i", "--intron", required=True)
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    genome = parse_fasta(args.genome)
    cds = parse_intron_file(args.intron)

    out = open(args.output, "w")

    for gene, regions in cds.items():
        for i in range(len(regions) - 1):
            sc, s1, e1, strand = regions[i]
            _, s2, _, _ = regions[i + 1]

            if s2 > e1 + 1:
                istart = e1 + 1
                iend = s2 - 1
                seq = genome[sc][istart - 1:iend]

                if strand == "-":
                    seq = revcomp(seq)

                out.write(f">{gene}_intron{i+1}\n{seq}\n")

    out.close()

if __name__ == "__main__":
    main()

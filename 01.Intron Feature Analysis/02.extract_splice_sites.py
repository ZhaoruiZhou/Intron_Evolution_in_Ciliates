#!/usr/bin/env python3
import argparse

def parse_fasta(path):
    genome = {}
    with open(path) as f:
        sid, seq = None, []
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

def parse_introns(path):
    cds = {}
    with open(path) as f:
        for line in f:
            if line.strip():
                sc, _, s, e, strand, gene_cds = line.split()
                gene = gene_cds.rsplit(".", 1)[0]
                cds.setdefault(gene, []).append((sc, int(s), int(e), strand))
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
    ap.add_argument("-o5", "--out5", required=True)
    ap.add_argument("-o3", "--out3", required=True)
    args = ap.parse_args()

    genome = parse_fasta(args.genome)
    cds = parse_introns(args.intron)

    out5 = open(args.out5, "w")
    out3 = open(args.out3, "w")

    for gene, regions in cds.items():
        for i in range(len(regions) - 1):
            sc, s1, e1, strand = regions[i]
            _, s2, _, _ = regions[i + 1]

            if s2 > e1 + 1:
                istart = e1 + 1
                iend = s2 - 1

                if strand == "+":
                    five = genome[sc][istart - 6:istart + 5]
                    three = genome[sc][iend - 5:iend + 6]
                else:
                    five = revcomp(genome[sc][iend - 5:iend + 6])
                    three = revcomp(genome[sc][istart - 6:istart + 5])

                out5.write(f">{gene}_intron{i+1}_5splice\n{five}\n")
                out3.write(f">{gene}_intron{i+1}_3splice\n{three}\n")

    out5.close()
    out3.close()

if __name__ == "__main__":
    main()

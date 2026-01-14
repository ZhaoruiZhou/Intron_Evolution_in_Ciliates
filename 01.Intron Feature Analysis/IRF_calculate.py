import pandas as pd
import pysam

def parse_cds(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None,
                     names=['chr', 'type', 'start', 'end', 'strand', 'gene_version'])
    df[['gene', 'rank']] = df['gene_version'].str.rsplit('.', n=1, expand=True)
    df['rank'] = df['rank'].astype(int)
    df = df.sort_values(by=['gene', 'rank'], ascending=[True, False])
    return df

def extract_introns(cds_df):
    introns = []
    for gene, group in cds_df.groupby('gene'):
        group = group.sort_values(by='start' if group.iloc[0]['strand'] == '+' else 'end')
        for i in range(len(group) - 1):
            exon1 = group.iloc[i]
            exon2 = group.iloc[i + 1]
            intron_start = exon1['end'] + 1
            intron_end = exon2['start'] - 1
            if intron_end > intron_start:
                introns.append({
                    'chr': exon1['chr'],
                    'start': intron_start,
                    'end': intron_end,
                    'strand': exon1['strand'],
                    'gene': exon1['gene'],
                    'intron_id': f"{exon1['gene']}_intron{i+1}",
                    'length': intron_end - intron_start + 1
                })
    return pd.DataFrame(introns)

def is_nir_read(read, intron_start, intron_end):
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op == 3:  # N = splice
            splice_start = ref_pos
            splice_end = ref_pos + length
            if splice_start <= intron_start and splice_end >= intron_end:
                return True
            ref_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X
            ref_pos += length
    return False

def is_ir_read(read, intron_start, intron_end, min_overlap=5):
    for op, _ in read.cigartuples:
        if op == 3:  # contains splice
            return False

    covered = 0
    for block_start, block_end in read.get_blocks():
        overlap_start = max(block_start, intron_start)
        overlap_end = min(block_end, intron_end)
        if overlap_end > overlap_start:
            covered += overlap_end - overlap_start

    return covered >= min_overlap

def count_ir_reads(bam, chrom, start, end, min_overlap=5, min_mapq=20):
    IR = 0
    NIR = 0

    try:
        for read in bam.fetch(chrom, max(0, start - 100), end + 100):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue

            if is_nir_read(read, start, end):
                NIR += 1
            elif is_ir_read(read, start, end, min_overlap=min_overlap):
                IR += 1

    except ValueError:
        pass

    return IR, NIR

def calculate_ir_ratio(bam_path, introns_df, min_overlap=5, min_mapq=20):
    bam = pysam.AlignmentFile(bam_path, 'rb')
    results = []

    for _, row in introns_df.iterrows():
        chrom = row['chr']
        start = int(row['start'])
        end = int(row['end'])
        intron_id = row['intron_id']
        gene = row['gene']
        length = row['length']

        ir, nir = count_ir_reads(bam, chrom, start, end, min_overlap, min_mapq)
        total = ir + nir
        ir_ratio = ir / total if total > 0 else 0

        results.append({
            'intron_id': intron_id,
            'gene': gene,
            'length': length,
            'IR_reads': ir,
            'NIR_reads': nir,
            'IR_ratio': round(ir_ratio, 4)
        })

    return pd.DataFrame(results)

if __name__ == '__main__':
    cds_file = 'test.txt'        # path to intron position file
    bam_file = 'RNA_sorted.bam'        # path to RNA-seq BAM file

    cds_df = parse_cds(cds_file)
    introns_df = extract_introns(cds_df)
    ir_df = calculate_ir_ratio(bam_file, introns_df)

    ir_df.to_csv('intron_retention_ratios.tsv', sep='\t', index=False)
    print(f"Processed {len(ir_df)} introns. Output: intron_retention_ratios.tsv")

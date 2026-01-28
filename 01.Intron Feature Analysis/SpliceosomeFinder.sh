#!/bin/bash
# Usage: ./SpliceosomeFinder.sh speciesA_proteins.fasta speciesA_genome.fasta

#######################################
### Initialize
#######################################
protein=$1
genome=$2
prefix=$(basename "${protein%.*}")  # Extract prefix from protein filename
human_ref="human_uniprot.fa"
pfam_db="/apps/users/zhouzhaorui/soft/Pfam/Pfam-A.hmm"
SEQKIT="/apps/users/zhouzhaorui/soft/conda1/envs/seqkit/bin/seqkit"

# Verify files exist
[ ! -f "$pfam_db" ] && echo "Error: Pfam database not found at $pfam_db" && exit 1
[ ! -f "$human_ref" ] && echo "Error: Human reference not found at $human_ref" && exit 1
[ ! -f "$protein" ] && echo "Error: SpeciesA protein file not found" && exit 1
[ ! -f "$genome" ] && echo "Error: SpeciesA genome file not found" && exit 1
[ ! -f "$SEQKIT" ] && echo "Error: seqkit not found at $SEQKIT" && exit 1

# Create prefixed directories
mkdir -p "${prefix}_tmp" "${prefix}_results"

#######################################
### Phase 1: RBH Extraction (EXACTLY as provided)
#######################################
makeblastdb -in $protein -input_type fasta -dbtype prot -parse_seqids -out "${prefix}_tmp/spdb"
blastp -query $human_ref -db "${prefix}_tmp/spdb" -evalue 1e-5 -outfmt 6 -num_threads 10 -out "${prefix}_tmp/blastp1.tab" -max_hsps 1 -num_alignments 1
makeblastdb -in $human_ref -input_type fasta -dbtype prot -parse_seqids -out "${prefix}_tmp/humandb"
blastp -query $protein -db "${prefix}_tmp/humandb" -evalue 1e-5 -outfmt 6 -num_threads 10 -out "${prefix}_tmp/blastp2.tab" -max_hsps 1 -num_alignments 1

# YOUR EXACT CODE BLOCK START
cut -f 1,2 "${prefix}_tmp/blastp1.tab" > "${prefix}_tmp/1.txt"
cut -f 2 "${prefix}_tmp/blastp2.tab" > "${prefix}_tmp/a.txt"
cut -f 1 "${prefix}_tmp/blastp2.tab" > "${prefix}_tmp/b.txt"
paste "${prefix}_tmp/a.txt" "${prefix}_tmp/b.txt" > "${prefix}_tmp/2.txt" &&
rm "${prefix}_tmp/a.txt" && rm "${prefix}_tmp/b.txt" &&
sort "${prefix}_tmp/1.txt" > "${prefix}_tmp/sort1.txt" && rm "${prefix}_tmp/1.txt"
sort "${prefix}_tmp/2.txt" > "${prefix}_tmp/sort2.txt" && rm "${prefix}_tmp/2.txt"
sed -i 's/\t//g' "${prefix}_tmp/sort1.txt"
sed -i 's/\t//g' "${prefix}_tmp/sort2.txt"

join "${prefix}_tmp/sort1.txt" "${prefix}_tmp/sort2.txt" > "${prefix}_results/rbh_pairs.txt"
sed -i 's/HUMAN/HUMAN\t/g' "${prefix}_results/rbh_pairs.txt"
# YOUR EXACT CODE BLOCK END

#######################################
### Phase 2: Extract Matched Sequences (Using specified seqkit)
#######################################
# Get matched human proteins
cut -f1 "${prefix}_results/rbh_pairs.txt" | sort -u > "${prefix}_tmp/matched_human_prots.txt"

# Get matched speciesA proteins
cut -f2 "${prefix}_results/rbh_pairs.txt" | sort -u > "${prefix}_tmp/matched_speciesA_prots.txt"

# Extract human sequences using specified seqkit
$SEQKIT grep -f "${prefix}_tmp/matched_human_prots.txt" $human_ref > "${prefix}_tmp/matched_human.fa"

# Extract speciesA sequences using specified seqkit
$SEQKIT grep -f "${prefix}_tmp/matched_speciesA_prots.txt" $protein > "${prefix}_tmp/matched_speciesA.fa"

#######################################
### Phase 3: Domain Analysis (ONLY matched sequences)
#######################################
echo "Running HMMER on matched sequences only..."
hmmsearch -E 0.01 --domtblout "${prefix}_results/human_pfam.domtbl" $pfam_db "${prefix}_tmp/matched_human.fa"  > /dev/null
awk '{print $1,$5}' "${prefix}_results/human_pfam.domtbl" | sort | grep -v '#' > "${prefix}_tmp/human_domains.txt"

hmmsearch -E 0.01 --domtblout "${prefix}_results/speciesA_pfam.domtbl" $pfam_db "${prefix}_tmp/matched_speciesA.fa"  > /dev/null
awk '{print $1,$5}' "${prefix}_results/speciesA_pfam.domtbl" | sort | grep -v '#' > "${prefix}_tmp/speciesA_domains.txt"

#######################################
### Phase 4: Create Domain Report
#######################################
echo -e "Human_Protein\tHuman_PFAMs\tSpeciesA_Protein\tSpeciesA_PFAMs" > "${prefix}_results/domain_report.tsv"

while read human_prot speciesA_prot; do
    # Get human PFAMs (comma-separated)
    human_pfams=$(grep "^$human_prot " "${prefix}_tmp/human_domains.txt" | cut -d' ' -f2 | sort | uniq | tr '\n' ',' | sed 's/,$//')
    [ -z "$human_pfams" ] && human_pfams="None"

    # Get speciesA PFAMs (comma-separated)
    speciesA_pfams=$(grep "^$speciesA_prot " "${prefix}_tmp/speciesA_domains.txt" | cut -d' ' -f2 | sort | uniq | tr '\n' ',' | sed 's/,$//')
    [ -z "$speciesA_pfams" ] && speciesA_pfams="None"

    echo -e "$human_prot\t$human_pfams\t$speciesA_prot\t$speciesA_pfams" >> "${prefix}_results/domain_report.tsv"
done < "${prefix}_results/rbh_pairs.txt"

#######################################
### Phase 5: Length Analysis (Using specified seqkit)
#######################################
# Create length comparison file
echo -e "Human_Protein\tSpeciesA_Protein\tHuman_Length\tSpeciesA_Length\tLength_Ratio" > "${prefix}_results/length_comparison.tsv"

while read human_prot speciesA_prot; do
    human_len=$($SEQKIT fx2tab -n -l $human_ref | awk -v prot="$human_prot" '$1==prot {print $2}')
    speciesA_len=$($SEQKIT fx2tab -n -l $protein | awk -v prot="$speciesA_prot" '$1==prot {print $2}')

    if [ -n "$human_len" ] && [ -n "$speciesA_len" ]; then
        ratio=$(awk -v sa_len="$speciesA_len" -v h_len="$human_len" 'BEGIN {printf "%.2f", sa_len/h_len}')
        echo -e "$human_prot\t$speciesA_prot\t$human_len\t$speciesA_len\t$ratio" >> "${prefix}_results/length_comparison.tsv"
    else
        echo -e "$human_prot\t$speciesA_prot\tNA\tNA\tNA" >> "${prefix}_results/length_comparison.tsv"
    fi
done < "${prefix}_results/rbh_pairs.txt"

#######################################
### Phase 6: Scoring (Updated logic)
#######################################
echo -e "Human_Protein\tScore\tEvidence" > "${prefix}_results/final_scores.tsv"

# Get all human proteins with BLASTP hits (from blastp1.tab)
cut -f1 "${prefix}_tmp/blastp1.tab" | sort -u > "${prefix}_tmp/all_blast_hits.txt"

# Process RBH pairs first
while read human_prot speciesA_prot; do
    # Get domains for both proteins
    human_pfams=$(grep "^$human_prot " "${prefix}_tmp/human_domains.txt" | cut -d' ' -f2 | sort | uniq)
    speciesA_pfams=$(grep "^$speciesA_prot " "${prefix}_tmp/speciesA_domains.txt" | cut -d' ' -f2 | sort | uniq)

    # Count PFAM domains
    human_pfam_count=$(echo "$human_pfams" | wc -w)
    speciesA_pfam_count=$(echo "$speciesA_pfams" | wc -w)

    # Check domain match based on new rules
    pfam_match=0

    if [ $human_pfam_count -gt 0 ]; then
        # Count common PFAMs
        common_pfams=$(comm -12 <(echo "$human_pfams") <(echo "$speciesA_pfams"))
        common_count=$(echo "$common_pfams" | wc -w)

        # Apply matching rules
        if [ $human_pfam_count -le 3 ] && [ $common_count -ge 1 ]; then
            pfam_match=1
        elif [ $human_pfam_count -gt 3 ] && [ $human_pfam_count -le 8 ] && [ $common_count -ge 2 ]; then
            pfam_match=1
        elif [ $human_pfam_count -gt 8 ] && [ $common_count -ge 3 ]; then
            pfam_match=1
        fi
    fi

    # Length check using the length comparison file (without bc)
    len_match=0
    ratio_line=$(grep "^$human_prot"$'\t'"$speciesA_prot" "${prefix}_results/length_comparison.tsv")
    if [ -n "$ratio_line" ]; then
        ratio=$(echo "$ratio_line" | awk '{print $5}')
        if [[ "$ratio" =~ ^[0-9.]+$ ]]; then
            # Using awk for floating point comparison
            len_match=$(awk -v r="$ratio" 'BEGIN {print (r >= 0.7 && r <= 1.3) ? 1 : 0}')
        fi
    fi

    # Scoring
    if [ $pfam_match -eq 1 ]; then
        if [ $len_match -eq 1 ]; then
            echo -e "$human_prot\t5\tRBH+Domain+Length" >> "${prefix}_results/final_scores.tsv"
        else
            echo -e "$human_prot\t4\tRBH+Domain" >> "${prefix}_results/final_scores.tsv"
        fi
    else
        echo -e "$human_prot\t2\tRBH_only" >> "${prefix}_results/final_scores.tsv"
    fi
done < "${prefix}_results/rbh_pairs.txt"

# Score non-RBH BLAST hits as 2 (BLAST-only)
cut -f1 "${prefix}_results/rbh_pairs.txt" | sort -u > "${prefix}_tmp/rbh_prots.txt"
comm -23 "${prefix}_tmp/all_blast_hits.txt" "${prefix}_tmp/rbh_prots.txt" > "${prefix}_tmp/non_rbh_hits.txt"

while read prot; do
    echo -e "$prot\t2\tBLAST_only" >> "${prefix}_results/final_scores.tsv"
done < "${prefix}_tmp/non_rbh_hits.txt"

#######################################
### Phase 7: tblastn for Unmatched (Only proteins without any BLASTP hits)
#######################################
echo "Processing unmatched proteins with tblastn..."
$SEQKIT seq -n $human_ref | sort > "${prefix}_tmp/all_human_prots.txt"
cut -f1 "${prefix}_tmp/blastp1.tab" | sort -u > "${prefix}_tmp/all_blast_hits.txt"
comm -23 "${prefix}_tmp/all_human_prots.txt" "${prefix}_tmp/all_blast_hits.txt" > "${prefix}_tmp/unmatched_prots.txt"

[ ! -s "${prefix}_tmp/unmatched_prots.txt" ] && echo "No unmatched proteins for tblastn" && exit 0

makeblastdb -in $genome -dbtype nucl -out "${prefix}_tmp/genome_db"

while read prot; do
    # Extract sequence using specified seqkit
    $SEQKIT grep -p $prot $human_ref > "${prefix}_tmp/query.fa"

    # Run tblastn
    tblastn -query "${prefix}_tmp/query.fa" -db "${prefix}_tmp/genome_db" -evalue 1e-5 \
             -outfmt '6 qseqid sseqid sseq' -num_threads 10 \
             -out "${prefix}_tmp/tblastn.out" -max_hsps 1 -num_alignments 1

    if [ -s "${prefix}_tmp/tblastn.out" ]; then
        # Extract protein sequence directly from output (column 3)
        translated_seq=$(cut -f3 "${prefix}_tmp/tblastn.out" | head -1)

        # Save translated sequence for HMMER
        echo -e ">${prot}_translated\n$translated_seq" > "${prefix}_tmp/translated.fa"

        # Domain check
        hmmsearch -E 0.01 --domtblout "${prefix}_tmp/trans_pfam.domtbl" $pfam_db "${prefix}_tmp/translated.fa"  > /dev/null
        trans_pfams=$(awk '{print $5}' "${prefix}_tmp/trans_pfam.domtbl" | sort | uniq | grep -v '#' | tr '\n' ',' | sed 's/,$//')
        human_pfams=$(grep "^$prot " "${prefix}_tmp/human_domains.txt" | cut -d' ' -f2 | sort | uniq | tr '\n' ',' | sed 's/,$//')

        # Check if at least one PFAM domain matches
        if [ -n "$human_pfams" ] && [ -n "$trans_pfams" ]; then
            common_pfams=$(comm -12 <(echo "$human_pfams" | tr ',' '\n') <(echo "$trans_pfams" | tr ',' '\n'))
            [ -n "$common_pfams" ] && pfam_match=1 || pfam_match=0
        else
            pfam_match=0
        fi

        if [ $pfam_match -eq 1 ]; then
            echo -e "$prot\t3\ttblastn+Domain" >> "${prefix}_results/final_scores.tsv"
        else
            echo -e "$prot\t1\ttblastn" >> "${prefix}_results/final_scores.tsv"
        fi
    else
        echo -e "$prot\t0\tNo_hit" >> "${prefix}_results/final_scores.tsv"
    fi
done < "${prefix}_tmp/unmatched_prots.txt"

#######################################
### Final Output
#######################################
sort -k2,2nr "${prefix}_results/final_scores.tsv" > "${prefix}_results/spliceosome_scores.ranked.tsv"

cut -f 2 "${prefix}_results/final_scores.tsv" > "${prefix}_tmp/one_zero.txt"
cut -f 1 "${prefix}_results/final_scores.tsv" > "${prefix}_tmp/tmp1.txt"
sed -i 's/sp//g' "${prefix}_tmp/tmp1.txt" &&
sed -i 's/tr//g' "${prefix}_tmp/tmp1.txt"
cat "${prefix}_tmp/tmp1.txt" | cut -c-6 > "${prefix}_tmp/tmp2.txt"
paste "${prefix}_tmp/tmp2.txt" "${prefix}_tmp/one_zero.txt" > "${prefix}_tmp/fig_data2.txt"
sort "${prefix}_tmp/fig_data2.txt" > "${prefix}_tmp/fig_data2_sort.txt"
join /apps/users/zhouzhaorui/work/spliceosome/try2/update/test/need/human_last_sort.txt "${prefix}_tmp/fig_data2_sort.txt" > "${prefix}_tmp/fig_data3.txt"
sort "${prefix}_tmp/fig_data3.txt" -k 3 > "${prefix}_tmp/fig_data3_sort.txt"
echo > "${prefix}_results/last_result.txt"
awk -F' ' -v prefix="$prefix" '{print $0" "prefix}' "${prefix}_tmp/fig_data3_sort.txt" >> "${prefix}_results/last_result.txt"
rm "${prefix}_tmp/tmp1.txt" "${prefix}_tmp/tmp2.txt" "${prefix}_tmp/one_zero.txt"
echo "Pipeline completed successfully for $prefix."
echo "Results:"
echo "1. Domain report: ${prefix}_results/domain_report.tsv"
echo "2. RBH pairs: ${prefix}_results/rbh_pairs.txt"
echo "3. Length comparison: ${prefix}_results/length_comparison.tsv"
echo "4. All scores: ${prefix}_results/final_scores.tsv"
echo "6. Ranked scores: ${prefix}_results/spliceosome_scores.ranked.tsv"

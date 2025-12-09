
## Hexenal Isomerases Lepidoptera

Scripts used in homolog identification and phylogenetic analysis of
hexenal isomerases in Lepidoptera

### 1. Identifying potential homologs with GMC_oxred_N and GMC_oxred_C motifs in 33 speices (1 dipteran, 4 trichopteran and 28 lepidopteran species)

    # prepare an HMM database for hmmscan
    cat GMC_oxred_C.hmm GMC_oxred_N.hmm > nc.hmm
    hmmpress nc.hmm

    # launch hmmscan
    bash sh_hmmscan.sh

    # select candidates with only 1 GMC_oxred_N motif and only 1 GMC_oxred_C motif
    cut -f5 hmmscan_dout | sort | uniq -c | sort -k1,1g | sed -r 's/^ +//g' | sed -r 's/ /\t/g' | grep -E '^2\b' | cut -f2 > candidate_list
    while read ID; do grep -E "$ID" hmmscan_dout >> candidate_list_hmm; done<<<$(cat candidate_list)

    # make sure GMC_oxred_N motif is upstream to GMC_oxred_C motif
    join -a 1 -a 2 -1 3 -2 3 -t $'\t' -e "NoHit" -o "1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5" \
    <(cat candidate_list_hmm | cut -f1,2,5,21,22 | grep 'GMC_oxred_N' | sort -t $'\t' -k3,3V) \
    <(cat candidate_list_hmm | cut -f1,2,5,21,22 | grep 'GMC_oxred_C' | sort -t $'\t' -k3,3V) | \
    grep -Ev 'NoHit' | awk 'BEGIN{FS=OFS="\t"}{if($5<$9) print $0}' | cut -f1,3,4,10 > asmid_aaid

    while IFS=$'\t' read ASM AA; do grep "$AA" ../411sp_aa/"$ASM".fa -A1 >> asmid_aaid.fa; done<<<$(cat asmid_aaid | cut -f1,2)

    # remove isoforms
    paste \
    <(cat asmid_aaid.fa | paste - - | cut -f1 | cut -d " " -f4 | sed -r 's/^gene://g') \
    asmid_aaid \
    <(cat asmid_aaid.fa | paste - - | awk 'BEGIN{FS=OFS="\t"}{print length($2), $0}') | awk 'BEGIN{FS=OFS="\t"}{if($3!~/PCG/) print $0}' | sort -t $'\t' -k1,1V -k6,6gr -k3,3V | sort -t $'\t' -k1,1V -u >> tmp_table_uniq_asmid_aaid

    # for *Chloridea virescens*, manually examine the isoforms
    paste \
    <(cat asmid_aaid.fa | paste - - | cut -f1 | cut -d " " -f4 | sed -r 's/^gene://g') \
    asmid_aaid \
    <(cat asmid_aaid.fa | paste - - | awk 'BEGIN{FS=OFS="\t"}{print length($2), $0}') | awk 'BEGIN{FS=OFS="\t"}{if($3~/PCG/) print $0}' | grep -Ev 'PCG78594.1|PCG78595.1' >> tmp_table_uniq_asmid_aaid

    # remove *Manduca sexta* ensembl entries
    grep -Ev "ENSMSXP.+JHU_Msex" tmp_table_uniq_asmid_aaid > tmp2_table_uniq_asmid_aaid

    # combine taxanomic information
    join -a2 -1 1 -2 2 -t $'\t' -e "NoHit" <(sort -t $'\t' -k1,1 asmid_ncbi_edirect_output) <(sort -t $'\t' -k2,2 tmp2_table_uniq_asmid_aaid) > table_uniq_asmid_aaid

    # cut to hmmscan env region
    while IFS=$'\t' read ID START END SEQ; do
      echo "$ID" | sed -r 's/^/>/g' >> cut_uniq_asmid_aaid.fa
      echo "$SEQ" | cut -c"$START"-"$END" >> cut_uniq_asmid_aaid.fa
    done<<<$(cat table_uniq_asmid_aaid | cut -f10,11,12,15)

    # add 3 fungal GMC as outgroup
    cat cut_uniq_asmid_aaid.fa /auto/brno2/home/bulah/nl/35sp/seq1204/outgroup_fungi/cut_outgroup.fa >> og_cut_uniq_asmid_aaid.fa

    # remove 14 fruit fly sequences
    cat og_cut_uniq_asmid_aaid.fa | paste - - | grep -Ev '>FBpp' | sed -r 's/\t/\n/g' > og_nofly_cut_uniq_asmid_aaid.fa

### 2. Identifying potential homologs with GMC_oxred_N and GMC_oxred_C motifs in 5 non-ditrysia speices

    # hmmscan search
    while read ID; do
      hmmscan -o out --domtblout tmp -E 1e-3 --domE 1e-3 --max --cpu "$PBS_NCPUS" nc.hmm ../getorf/"$ID"_aa.fa
      paste <(grep -Ev '^#' tmp | sed -r 's/ +/\t/g' | cut -f1-22) <(grep -Ev '^#' tmp | sed -r 's/ +/\t/g' | cut -f1-22 --complement | sed -r 's/\t/ /g') | sed -r "s/^/$ID\t/g" >> nd_hmmscan_dout
      rm -rf tmp out
    done<<<$(cut -f1 target | sed -r 's/000000//g')

    # select candidates with only 1 GMC_oxred_N motif and only 1 GMC_oxred_C motif
    cut -f5 nd_hmmscan_dout | sort | uniq -c | sort -k1,1g | sed -r 's/^ +//g' | sed -r 's/ /\t/g' | grep -E '^2\b' | cut -f2 > candidate_list
    while read ID; do grep -E "$ID" nd_hmmscan_dout >> candidate_list_hmm; done<<<$(cat candidate_list)

    # make sure GMC_oxred_N motif is upstream to GMC_oxred_C motif
    join -a 1 -a 2 -1 3 -2 3 -t $'\t' -e "NoHit" -o "1.1,1.2,1.3,1.4,1.5,2.1,2.2,2.3,2.4,2.5" \
    <(cat candidate_list_hmm | cut -f1,2,5,21,22 | grep 'GMC_oxred_N' | sort -t $'\t' -k3,3V) \
    <(cat candidate_list_hmm | cut -f1,2,5,21,22 | grep 'GMC_oxred_C' | sort -t $'\t' -k3,3V) | \
    grep -Ev 'NoHit' | awk 'BEGIN{FS=OFS="\t"}{if($5<$9) print $0}' | cut -f1,3,4,10 > asmid_aaid

    # cut to hmmscan env region
    while IFS=$'\t' read ASM AA START END; do
      grep "$AA" ../getorf/"$ASM"_aa.fa -A1 >> asmid_aaid.fa
      paste <(echo ">${AA}") <(grep "$AA" ../getorf/"$ASM"_aa.fa -A1 | tail -n 1 | cut -c"$START"-"$END") | sed -r 's/\t/\n/g' >> cut_asmid_aaid.fa
    done<<<$(cat asmid_aaid)

    # remove duplicates (isoforms)
    # 58 unique seq
    cat asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -u | wc -l
    # 19 duplicated seq
    cat asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -D | wc -l
    # 7 unique seq out of 19 duplicated seq
    cat asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -D | sort -t $'\t' -k2,2V -u | wc -l
    #
    cat <(cat asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -u) <(cat asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -D | sort -t $'\t' -k2,2V -u) | sed -r 's/\t/\n/g' > uniq_asmid_aaid.fa
    cat <(cat cut_asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -u) <(cat cut_asmid_aaid.fa | paste - - | sort -t $'\t' -k2,2V | uniq -f 1 -D | sort -t $'\t' -k2,2V -u) | sed -r 's/\t/\n/g' > uniq_cut_asmid_aaid.fa

### 3. Protein alignment, model selection, and tree building

    # concatenate all candidates
    cat og_cut_uniq_asmid_aaid_1185.fa non_ditrysia_uniq_cut_asmid_aaid.fa | paste - - | grep -Ev '^>GEOR01' | sed -r 's/\t/\n/g' > og_cut_uniq_asmid_aaid_1251.fa

    SEQ=og_cut_uniq_asmid_aaid_1251.fa

    einsi --thread 30 "$SEQ" > msa.fa

    clipkit msa.fa -s aa -l

    mv msa.fa.clipkit trimmed_msa.fa

    MSA=trimmed_msa.fa

    modeltest-ng -p 30 -d aa -i ${MSA} -h uigf -f ef -o ${MSA}.modeltestng

    MOD=$(tail ${MSA}.modeltestng.log | sed -r 's/^ +//g; s/ +/\t/g' | grep -E '^BIC' | cut -f2 | sed -r 's/ /x/g')

    iqtree2-mpi -T AUTO --threads-max ${PNS_NCPUS} -s ${MSA} -m ${MOD} -B 3000 -alrt 3000 --bnni --seed 12345 --bcor 0.99 -o AAF59929.2,XP_001727544.1,AAD01493.1

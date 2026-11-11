#!/bin/bash

#mkdir -p lists_newvep_2
cp /path/to/constraint/list/chen_constraint_z_Supplementary_Data_2.bed.gz .
gunzip chen_constraint_z_Supplementary_Data_2.bed.gz

for chr in {14..22}; do
    #input_file="/path/to/WGS/vep_annotations/ukb24310_c${chr}_ALL_indels_normalised.vep.gz"
    input_file="chr${chr}_merged_normalized_study_flag_corrected_processed.vep.gz"
    
    echo "Processing chromosome ${chr}..."

    module load dxpy
    dx download path/to/UKB/vep_annotations/${input_file}

## make different varlists for each of coding, proximal_utr, proximal_notutr, intronic, other ##
    zgrep -v "^##" "$input_file"|awk -v CHR="$chr" '
    BEGIN {
    FS  = "\t"   # read tabs
    OFS = "\t"   # write tabs as well (for reg)
    }
    {
        id = $1
        feat_type = $6
        conseq = $7
        mane = ($24 != "-" ? 1 : 0)
        biotype = $22
        chen = $89

        if (mane) {

            # ALL CODING -> coding_exons
            if (conseq ~ /missense_variant|start_lost|stop_lost|stop_gained|inframe_insertion|inframe_deletion|synonymous_variant|coding_sequence_variant|protein_altering_variant|splice_acceptor_variant|splice_donor_variant|frameshift_variant/) {
                fname = "lists_newvep_2/chr" CHR "_coding"
                print id >> fname
            }

            # UTR + UTR intronic -> proximal_utr
            if (conseq ~ /5_prime_UTR_variant|3_prime_UTR_variant/) {
                fname = "lists_newvep_2/chr" CHR "_proximal_utr"
                print id >> fname
            }

            # Proximal-intergenic -> proximal_notutr
            if (conseq ~ /upstream_gene_variant|downstream_gene_variant/) {
                fname = "lists_newvep_2/chr" CHR "_proximal_notutr"
                print id >> fname
            }

            # Intronic section of protein transcript -> introns (protein and RNA introns)
            if (conseq ~ /intron_variant/ && conseq !~ /coding_sequence_variant|splice_acceptor_variant|splice_donor_variant/) {
                fname = "lists_newvep_2/chr" CHR "_intronic"
                print id >> fname
            }

        }

        #Regulatory region -> regulatory
        # keep SNP id + CHEN col
        if (feat_type ~ /RegulatoryFeature/) {
            fname = "lists_newvep_2/chr" CHR "_regulatory"
            print id, chen >> fname 
        }

        # rna_exons
        if (biotype ~ /lncRNA|miRNA|rRNA|snoRNA|snRNA/ && conseq ~ /non_coding_transcript_exon_variant/) {
            fname = "lists_newvep_2/chr" CHR "_rna_exons"
            print id >> fname
        }

        # add RNA introns to _intronic
        if (biotype ~ /lncRNA|miRNA|rRNA|snoRNA|snRNA/ && conseq ~ /non_coding_transcript_variant/) {
            fname = "lists_newvep_2/chr" CHR "_intronic"
            print id >> fname
        }

        # Non-coding / intergenic categories (not subset to MANE Select) -> other
        if (conseq ~ /intron_variant/ && conseq !~ /splice_acceptor_variant|splice_donor_variant/) {
            fname = "lists_newvep_2/chr" CHR "_other"
            print id >> fname
        }

        if (conseq ~ /intergenic_variant|regulatory_region_variant|TF_binding_site_variant/) {
            fname = "lists_newvep_2/chr" CHR "_other"
            print id >> fname
        }

    }'

    echo "Extracted variants from VEP file"

## filter files ##

    # remove unconstrained regulatory from reg file using constraint col from VEP file, keep only SNP id col in reg file
    cat chen_constraint_z_Supplementary_Data_2.bed | awk -v c="chr${chr}" '$1==c && $4<4 {print ""c"_"$2"_"$3}' > chen_to_remove
    grep -vwf chen_to_remove lists_newvep_2/chr${chr}_regulatory | cut -f1 > reg_${chr}.tmp && mv reg_${chr}.tmp lists_newvep_2/chr${chr}_regulatory

    # remove duplicate lines for each varlist file
    for f in lists_newvep_2/chr${chr}_*; do
        sort -u "$f" > "$f".sorting.tmp && mv "$f".sorting.tmp "$f"
    done

    # remove duplicates as coding > rna_exons > proximal_utr > proximal_notutr > intronic > regulatory > other

    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_coding lists_newvep_2/chr${chr}_rna_exons > ${chr}_rna_exons.tmp && mv ${chr}_rna_exons.tmp lists_newvep_2/chr${chr}_rna_exons

    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_coding lists_newvep_2/chr${chr}_proximal_utr > ${chr}_proximal_utr.tmp && mv ${chr}_proximal_utr.tmp lists_newvep_2/chr${chr}_proximal_utr
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_rna_exons lists_newvep_2/chr${chr}_proximal_utr > ${chr}_proximal_utr.tmp && mv ${chr}_proximal_utr.tmp lists_newvep_2/chr${chr}_proximal_utr

    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_coding lists_newvep_2/chr${chr}_proximal_notutr > ${chr}_proximal_notutr.tmp && mv ${chr}_proximal_notutr.tmp lists_newvep_2/chr${chr}_proximal_notutr
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_rna_exons lists_newvep_2/chr${chr}_proximal_notutr > ${chr}_proximal_notutr.tmp && mv ${chr}_proximal_notutr.tmp lists_newvep_2/chr${chr}_proximal_notutr
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_utr lists_newvep_2/chr${chr}_proximal_notutr > ${chr}_proximal_notutr.tmp && mv ${chr}_proximal_notutr.tmp lists_newvep_2/chr${chr}_proximal_notutr

    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_coding lists_newvep_2/chr${chr}_intronic > ${chr}_intronic.tmp && mv ${chr}_intronic.tmp lists_newvep_2/chr${chr}_intronic
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_rna_exons lists_newvep_2/chr${chr}_intronic > ${chr}_intronic.tmp && mv ${chr}_intronic.tmp lists_newvep_2/chr${chr}_intronic
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_utr lists_newvep_2/chr${chr}_intronic > ${chr}_intronic.tmp && mv ${chr}_intronic.tmp lists_newvep_2/chr${chr}_intronic
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_notutr lists_newvep_2/chr${chr}_intronic > ${chr}_intronic.tmp && mv ${chr}_intronic.tmp lists_newvep_2/chr${chr}_intronic

    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_coding lists_newvep_2/chr${chr}_regulatory > ${chr}_regulatory.tmp && mv ${chr}_regulatory.tmp lists_newvep_2/chr${chr}_regulatory
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_rna_exons lists_newvep_2/chr${chr}_regulatory > ${chr}_regulatory.tmp && mv ${chr}_regulatory.tmp lists_newvep_2/chr${chr}_regulatory
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_utr lists_newvep_2/chr${chr}_regulatory > ${chr}_regulatory.tmp && mv ${chr}_regulatory.tmp lists_newvep_2/chr${chr}_regulatory
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_notutr lists_newvep_2/chr${chr}_regulatory > ${chr}_regulatory.tmp && mv ${chr}_regulatory.tmp lists_newvep_2/chr${chr}_regulatory
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_intronic lists_newvep_2/chr${chr}_regulatory > ${chr}_regulatory.tmp && mv ${chr}_regulatory.tmp lists_newvep_2/chr${chr}_regulatory

    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_coding lists_newvep_2/chr${chr}_other > ${chr}_other.tmp && mv ${chr}_other.tmp lists_newvep_2/chr${chr}_other
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_rna_exons lists_newvep_2/chr${chr}_other > ${chr}_other.tmp && mv ${chr}_other.tmp lists_newvep_2/chr${chr}_other
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_utr lists_newvep_2/chr${chr}_other > ${chr}_other.tmp && mv ${chr}_other.tmp lists_newvep_2/chr${chr}_other
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_proximal_notutr lists_newvep_2/chr${chr}_other > ${chr}_other.tmp && mv ${chr}_other.tmp lists_newvep_2/chr${chr}_other
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_intronic lists_newvep_2/chr${chr}_other > ${chr}_other.tmp && mv ${chr}_other.tmp lists_newvep_2/chr${chr}_other
    awk 'NR == FNR {a[$1]; next} !($1 in a) {print}' lists_newvep_2/chr${chr}_regulatory lists_newvep_2/chr${chr}_other > ${chr}_other.tmp && mv ${chr}_other.tmp lists_newvep_2/chr${chr}_other

    echo "Filtered lists"

    zstd -8 -T8 --rm lists_newvep_2/chr${chr}_*
    echo "Finished processing chromosome ${chr}"
done



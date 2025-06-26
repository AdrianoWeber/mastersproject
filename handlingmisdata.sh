#!/bin/bash

## This is a shell script completing missing data due to multi-samples vcf merging. It must have a input vcf, an output directory and a directory to bam files to query the coverage depth of every missing data with samtools depth.
## It returns a reannotated file named "cov[n_samples]samples_1KGP.vcf" with missing pos with coverage above 10 read (per base) assigned to reference. The script reannotate the first n samples.
##Beware that the script might be slow due to several reading an writing of files.

## This script is part of the master thesis of Adriano Weber at the University of Geneva (09/2024-09/2025).
# This script is licensed under the GNU General Public License v3.0 or later.
# See https://www.gnu.org/licenses/gpl-3.0.en.html for more details.

#Input arguments
INPUT_DIR="$1"
OUTPUT_DIR="$2"
BAM_DIR="$3"
N_COL="$4"

#USAGE
usage(){
  echo "This is a function to rewrite missing positions in regard of their coverage"
  echo "Usage: [VCF_to_rewrite] [output_directory] [directory_to_bam_files] [number_of_samples_to_rewrite]"
}

if [ "$#" -lt 4 ]; then
    echo "Error: Missing arguments."
    usage
    exit 1
fi

#Choosing an output file
vcf_out="$OUTPUT_DIR/cov"$N_COL"samples_1KGP.vcf"

#Finding individuals ID for the first n samples.
id_list=$(bcftools query -l "$INPUT_DIR" | head -n $N_COL)
echo "$id_list" > "$OUTPUT_DIR/samples-tmp.txt"
echo "Individuals to rewrite:$id_list"

# Estimer par individus puis par positions
while IFS= read -r line; do
    if [[ ! "$line" =~ ^#.* ]]; then
        modified_line="$line"
        for (( individu=1; individu<=N_COL; individu++ )); do
            id=$(echo "$id_list" | sed -n "${individu}p")
            if [[ $(awk -v col="$individu" '{print $(10 + $col)}' <<< "$line") =~ ^\./\..* ]]; then
                cov=$(samtools depth "$BAM_DIR/$id.deepvariant.haplotagged.bam" -r "$(awk '{print $1 ":" $2 "-" $2}' <<< "$line")" | awk -v col="$individu" '{print $3}')
                echo "The coverage depth is of: $cov."
                if [[ $cov -ge 10 ]]; then
                    modified_line=$(awk -v col="$individu" -v couv="$cov" 'BEGIN{FS=OFS="\t"} {$(10 + col)="0/0:.:"couv":."; print}' <<< "$modified_line")

                fi
            fi
        done
        echo "$modified_line" >> "$vcf_out"
    else
        echo "$line" >> "$vcf_out"
    fi
done < "$INPUT_DIR"
#Erasing temporary sample file
rm "$OUTPUT_DIR/samples-tmp.txt"
echo "Terminated. Output file to $vcf_out."


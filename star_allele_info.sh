#! /bin/bash

##The below script was written by Adriano Weber (2025)
##This is a script for extracting star alleles info from a vcf thanks to a .TSV file (previously imported from pharmGKB)
##The outputs are a vcf file with only star alleles region selected (.vcf) and a info table of star allele frequency in the selected population, MUS by default (.stat) 
##Usage: source star_allele_info.sh [tsv_file direction] [vcf_file direction] [gene name] [chromosome number]
##Beware when adatpting the scripts, some lines are made for a population named "MUS".

#Checking the presence of the software bcftools needed for this script
command -v bcftools >/dev/null 2>&1 || { echo "bcftools is required but not installed. Exiting."; exit 1; }

#Parameters
tsv_file=$1
vcf_file=$2
gene_name=$3
chromosome=$4
if [ -z "$tsv_file" ] || [ -z "$vcf_file" ] || [ -z "$gene_name" ] || [ -z "$chromosome" ]; then
    echo "Error: Missing required parameter(s)"
    echo "Usage: source star_allele_info.sh [tsv_file direction] [vcf_file direction] [gene name] [chromosome number]"
    exit 1
fi

#Variable declaration
pos_file="star_alleles${gene_name}.pos"
out_vcf="star_allele_${gene_name}.vcf"
tsv_info="info.tsv"
out_stat="MUS_${gene_name}_maf.stat"
pop_name="MUS"

#Making needed temporary files 
cat $tsv_file | grep -E "^${gene_name}\*[0-9]+\s" | cut -f 1,5,7,8 > $tsv_info #Print the main Star alleles number, start position, ref allele, alternative allele
cat $tsv_file | grep -E "^${gene_name}\*[0-9]+\s" | cut -f 5 | sed "s/^/chr${chromosome}\t/" | tail -n +2 > $pos_file # Making a position file for bcftools 

#Output vcf
bcftools view -R $pos_file -Ov -o $out_vcf $vcf_file

#Stat computation for MUSCAT (will have to change because of samples calling)
sample_list=$(bcftools view -h $out_vcf | tail -n 1 | awk '{for(i=1;i<=NF;i++) if($i ~ /^MUS/) print $i}' | tr "\n" "," | sed "s/.$//")
bcftools view -s $sample_list $out_vcf | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > "MUS${gene_name}_stat.1"
#Adding the star allele name
echo -e "Star_Allele\tCHROM\tPOS\tREF\tALT\tMAF"> $out_stat #Writing the first line of stat output
while IFS= read -r line; do
    pos=$(echo $line | cut -f 2 )
    star_all_name=$(grep $pos $tsv_info | cut -f 1 | head -n 1) #Here the head is to keep only the first correspondance with the position
    echo "$star_all_name\t$line" >> $out_stat
done < "MUS${gene_name}_stat.1"
#Remove the temporary files
rm "MUS${gene_name}_stat.1"
rm $tsv_info
rm $pos_file
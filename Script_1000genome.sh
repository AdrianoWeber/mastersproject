#!/bin/bash

###Check for depedencies
command -v bcftools >/dev/null 2>&1 || { echo "bcftools is required but not installed. Exiting."; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "wget is required but not installed. Exiting."; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "bgzip is required but not installed. Exiting."; exit 1; }

###USAGE
usage(){
  echo "This is a function to download data from the 1000Genomes database."
  echo "Usage: $0 -options <path_to_loci.zip>"
  echo "OPTIONS:"
  echo "-o: Where to write the output file(s)."
  echo "-l: Direction of the text file containing loci. The loci must be noted as *[chr]:[start_pos]-[end_pos]* with one locus per line"
  echo "-p: Comma seperated value of populations to extract. The populations accepted are: AFR: For african populations.\nAMR: For american populations.\nEAS: For east-asian populations.\nEUR: For european populations.\nSAS: For south-asian populations.\n"
  echo "-h: Help of the script."
}

###Initialisation of variables
OUTPUT_DIR=.
VCF_DIR="./1000genomes_vcf_files"
mkdir -p "$VCF_DIR"
base_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
panel_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

loci_file="0"
pop_names="0"
files_to_merge=()
#Population commented are excluded because they do not live in the continent to which they are assigned.
declare -A continent_populations
continent_populations["AFR"]="YRI,LWK,GWD,MSL,ESN" #ACB, ASW
continent_populations["AMR"]="PUR,CLM,PEL" #MXL
continent_populations["EAS"]="CHB,JPT,CHS,CDX,KHV"
continent_populations["EUR"]="TSI,FIN,GBR,IBS" #CEU
continent_populations["SAS"]="PJL,BEB" #GIH, ITU, STU

#Probablement ajouter un "files_to_delete" in fine

###Function to download vcf files
download_vcf(){
    chrom=$1
    if [[ "$chrom" =~ ^[0-9]+$ ]]; then
      vcf_file="ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    fi
    if [[ "$chrom" == 'X' ]]; then
      vcf_file="ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"
    fi
    if [[ "$chrom" == 'Y' ]]; then
      vcf_file="ALL.chr${chrom}.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
    fi
    if [[ "$chrom" == "M" ]]; then
      vcf_file="ALL.chr${chrom}T.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
    fi

    vcf_url="${base_url}/${vcf_file}"
    tbi_url="${vcf_url}.tbi"
    echo "Downloading VCF and index for chromosome ${chrom}..."
    echo "DIRECTION DE SORTIE: $VCF_DIR"
    wget -O "$VCF_DIR/datachr${chrom}_1Kgenome.vcf.gz" "$vcf_url" #Préciser VCF_DIR dans -O
    wget -O "$VCF_DIR/datachr${chrom}_1Kgenome.vcf.gz.tbi" "$tbi_url"
}

#Function to merge vcf files
merge_vcf(){
    vcf_list=("$@")
    merged_output="$OUTPUT_DIR/1000genomes_data_merged.vcf.gz"
    echo "Merging files..."
    echo "${vcf_list[@]}"
    bcftools view "${vcf_list[@]}" | bgzip -c > "$merged_output"
    echo "Merged VCF file created at $merged_output"
}

if [ "$#" -eq 0 ]; then
  echo "No parameter provided. Downloading all VCF from 1000 genomes...."
  chromosomes=( {1..22} X Y M ) #J'ai remplacé MT par M
   for chrom in "${chromosomes[@]}"
   do
        download_vcf "$chrom"
    done
  echo "DOWNLOADING COMPLETED!"
  exit 1
fi

###OPTIONS setting
while getopts ":ho:l:p:" opt; do 
    case "${opt}" in
          #help flag
      h)
        usage
        exit
        ;;
      o)
        OUTPUT_DIR="${OPTARG}"
        # Check if the direction exits and create it if not
        if [ ! -d "$OUTPUT_DIR" ]; then
            mkdir -p "$OUTPUT_DIR"  # Uses -p to create upper level direction
        fi
        ;;
      l)
        loci_file="${OPTARG}"
        ;;
      p)
        sample_list=""
        separator=""
        pop_names="${OPTARG}"
        pop_list=$(echo "$pop_names" | tr ',' ' ')
        wget -O "panel.txt" "$panel_url"
        for continent in $pop_list; do
          sample_list+="${separator}$(grep -E "$(echo "${continent_populations[${continent}]}" | tr ',' '|')" "panel.txt" | awk '{print $1}' | tr '\n' ',' | sed 's/,$//')"
          separator=","
        done
        echo "SAMPLE LIST: $sample_list"
        ;;
      #if another flag is used, error
      *)
        usage
        ;;
    esac
done

if [  "$loci_file" == "0" ]; then
  echo "No loci parameter provided. Downloading all VCF from 1000 genomes...."
  chromosomes=( {1..22} X Y MT )
   for chrom in "${chromosomes[@]}"
   do
        download_vcf "$chrom"
    done
  echo "DOWNLOADING COMPLETED!"
  exit 1
fi

###Loci selection
vcf_needed=$(awk -F '[:]' '{print $1}' "$loci_file" | sort -u)
echo "Downloading VCF needed for all loci provided..."
for chrom in $vcf_needed
do
    download_vcf "$chrom"
done
echo "Downloading completed!"
##Filtering
echo "Filtering VCF files to keep only the specified positions..."
while read -r line; do
    chrom=$(echo "$line" | awk -F '[:]' '{print $1}')
    positions=$(echo "$line" | awk -F '[:]' '{print $2}')
    echo "Processing chromosome: $chrom"
    echo "Positions: $positions"  
# Filtrer les fichiers VCF pour ne garder que les positions d'intérêt
    bcftools view --regions "${chrom}:${positions}" "$VCF_DIR/datachr${chrom}_1Kgenome.vcf.gz" \ | bgzip -c > "$OUTPUT_DIR/ALL.chr${chrom}.filtered.vcf.gz"
    echo "FILTERING POSITIONS..."
done < "$loci_file"
echo "Filtering completed! Filtered VCFs saved in $VCF_DIR"
files_to_merge+=("$OUTPUT_DIR/ALL.chr${chrom}.filtered.vcf.gz")

###Merging all data in 1 vcf
printf "Files to merge: %s\n" "${files_to_merge[@]}" 
echo "Writing a unique vcf..."
if  merge_vcf "${files_to_merge[@]}" ; then
    echo "Success!"
    else
        echo "ERROR during merging."
fi

###Populations selection
if [ "$pop_names" != "0" ]; then
    echo "FILTERING POPULATIONS..."
    tabix "$OUTPUT_DIR/1000genomes_data_merged.vcf.gz" 
    echo "SAMPLE LIST: $sample_list"
    bcftools view -s "$sample_list" "$OUTPUT_DIR/1000genomes_data_merged.vcf.gz" | bgzip -c > "$OUTPUT_DIR/1000genomes_data_merged_pop_filtered.vcf.gz"
    #Rajouter ici une mention de succès
    echo "Success! Program ended without problems."
fi

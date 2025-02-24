#!/bin/bash

###Check for depedencies
command -v bcftools >/dev/null 2>&1 || { echo "bcftools is required but not installed. Exiting."; exit 1; }
command -v wget >/dev/null 2>&1 || { echo "wget is required but not installed. Exiting."; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "bgzip is required but not installed. Exiting."; exit 1; }

###USAGE
usage(){
  echo "This is a function to download data from the 1000Genomes database from the 30xGrChr38."
  echo "This version download chromosome's data one at the time and delete it after usage. This is due to the size of data."
  echo "Usage: $0 [--options] <path_to_loci.zip>" #Not quite sure about <path_to_loci.zip>
  echo "OPTIONS:"
  echo "-o: Where to write the output file(s)."
  echo "-v: Where to download the 1KGP data."
  echo "-l: Direction of the text file containing loci. The loci must be noted as *[chr]:[start_pos]-[end_pos]* with one locus per line"
  echo -e "-p: Comma seperated value of populations to extract. The populations accepted are:\n\t- AFR: For african populations.\n\t- AMR: For american populations.\n\t- EAS: For east-asian populations.\n\t- EUR: For european populations.\n\t- SAS: For south-asian populations.\n"
  echo "-k: keep options. If provided then all data from 1KGP are kept."
  echo "-h: Display this usage manual for the script."
}

###Initialisation of variables
OUTPUT_DIR=.
VCF_DIR="./1KGP_vcf_files" # On pourrait le mettre dans la output dir plutot que wd.
#mkdir -p "$VCF_DIR"
# exec > >(tee -a $OUTPUT_DIR/fichier_log.log) 2>&1 #Test de si ça marche pour faire un fichier log, commented pour plus de simplicite pour l'instant.
base_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/"
index_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index" # Without the trios 

loci_file="0"
pop_names="0"
declare -a files_to_merge=()
keep=0
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
    #vcf_file="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
    if true; then
    if [[ "$chrom" =~ ^[0-9]+$ ]]; then
      vcf_file="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.shapeit2-duohmm-phased.vcf.gz"
    fi
    if [[ "$chrom" == 'X' ]]; then
      vcf_file="CCDG_14151_B01_GRM_WGS_2020-08-05_chr${chrom}.filtered.eagle2-phased.v2.vcf.gz"
    fi
    if [[ "$chrom" == 'Y' ]]; then
      echo "NO chromosome Y data so far."
    fi
    fi

    if [[ "$chrom" == "M" ]]; then
      echo "NO MT Data so far."
    fi

    vcf_url="${base_url}/${vcf_file}"
    tbi_url="${vcf_url}.tbi"
    if [ ! -e "$VCF_DIR/$vcf_file" ]; then #To avoid downloading multiple time
      echo "Downloading VCF and index for chromosome ${chrom}..."
      echo "DIRECTION DE SORTIE: $VCF_DIR"
      wget -O "$VCF_DIR/$vcf_file" "$vcf_url"
      sleep 7
      wget -O "$VCF_DIR/$vcf_file.tbi" "$tbi_url"
      sleep 5
      touch "$VCF_DIR/$vcf_file.tbi"
    fi
}

###Function to merge vcf files
merge_vcf(){
    declare -a vcf_list=()
    vcf_list=("$@")
    merged_output="$OUTPUT_DIR/1000genomes_data_merged.vcf.gz"
    echo "Merging files..."
    echo "${vcf_list[@]}"
    bcftools  concat -Da - Oz -o "$merged_output" "${vcf_list[@]}" #-D parameter is used to suppress doubles.
    tabix "$merged_output"
    echo "Merged VCF file created at $merged_output"
}

###OPTIONS setting
while getopts ":hv:o:l:p:k" opt; do 
    case "${opt}" in
          #help flag
      h)
        usage
        exit
        ;;

      v)
        VCF_DIR="${OPTARG}"
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
        wget -O "index.txt" "$index_url"
        echo "Writing a file without parental relations..."
        awk -F"\t" '!/^##/ {print $10 "\t" $11}' "index.txt" > samples.txt
        for continent in $pop_list; do
          if [[ ! "$continent" =~ ^(AFR|AMR|EAS|EUR|SAS)$ ]]; then
            echo "$continent is not a valid argument for -p option. See usage with ./Script_1000genome.sh -h."
            exit 1
          fi
          sample_list+="${separator}$(tail -n +2 "samples.txt" | grep -E "$(echo "${continent_populations[${continent}]}" | tr ',' '|')" | awk '{print $1}' | tr '\n' ',' | sed 's/,$//')"
          separator=","
        done
        echo "SAMPLE LIST: $sample_list"
        ;;
      #if another flag is used, error
      k)
        keep=1
        ;;
      :) 
        echo "Option -$OPTARG needs an argument."
        exit
        ;;
      *)
        echo "-$OPTARG is not a valid option."
        usage
        exit
        ;;
    esac
done

mkdir -p "$VCF_DIR"

###Manage if no  loci provided
if [  "$loci_file" = "0" ]; then
  echo "No loci parameter provided. Downloading all VCF from 1000 genomes...."
  chromosomes=( {1..22} X Y ) # M
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
  echo "DOWNLOADING DATA FOR CHROMOSOME $chrom ..."
    download_vcf "$chrom"
  echo "Downloading completed!"
  ##Filtering
  echo "Filtering VCF files to keep only the specified positions..."
  while read -r line; do
      if [ "$chrom" = "$(echo "$line" | awk -F '[:]' '{print $1}')" ];then
        positions=$(echo "$line" | awk -F '[:]' '{print $2}')
        echo "Processing chromosome: $chrom"
        echo "Filtering positions: $positions"  
        bcftools view --regions "chr${chrom}:${positions}" "$VCF_DIR/$vcf_file" | bgzip -c > "$OUTPUT_DIR/Chr${chrom}Pos${positions}.filtered.vcf.gz"
        tabix "$OUTPUT_DIR/Chr${chrom}Pos${positions}.filtered.vcf.gz"
        files_to_merge+=("$OUTPUT_DIR/Chr${chrom}Pos${positions}.filtered.vcf.gz")
        fi
  done < "$loci_file"

  #Option to remove data not enough space
  if [[ $keep -eq 0 ]]; then
    echo "REMOVING DATA FOR $chrom FOR MEMORY SANETY..."
    rm "$VCF_DIR/$vcf_file"
    rm "$VCF_DIR/$vcf_file"
  fi
done
echo "Filtering completed! Filtered VCFs saved in $VCF_DIR"

### Merging all data in 1 vcf
echo "Files to merge: ${files_to_merge[*]}"
echo "Writing a unique vcf..."
if  merge_vcf "${files_to_merge[@]}" ; then
    echo "Success! Merging completed."
    else
        echo "ERROR during merging."
        exit 1
fi

##(temporary commentary)##ICI ON INTEGRERA LA SELECTION DE POP A LA BOUCLE SI LES DATA SONT TROP GRANDES EGALEMENT
###Populations selection
if [ "$pop_names" != "0" ]; then
    echo "FILTERING POPULATIONS..."
    tabix "$OUTPUT_DIR/1000genomes_data_merged.vcf.gz" 
    if bcftools view --force-samples -s "$sample_list" "$OUTPUT_DIR/1000genomes_data_merged.vcf.gz" | bgzip -c > "$OUTPUT_DIR/1000genomes_data_merged_pop_filtered.vcf.gz"; then
      echo "Success! Program ended without major problems."
      else 
        echo "ERROR during population selection."
    fi
fi
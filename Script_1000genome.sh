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
continent_populations=(
    ["AFR"]="YRI,LWK,GWD,MSL,ESN" #ACB, ASW
    ["AMR"]="PUR,CLM,PEL" #MXL
    ["EAS"]="CHB,JPT,CHS,CDX,KHV"
    ["EUR"]="TSI,FIN,GBR,IBS" #CEU
    ["SAS"]="PJL,BEB" #GIH, ITU, STU
)
#Probablement ajouter un "files_to_delete" in fine



###Function to download vcf files
download_vcf(){
    chrom=$1
    vcf_file="ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    vcf_url="${base_url}/${vcf_file}"
    tbi_url="${vcf_url}.tbi"
    echo "Downloading VCF and index for chromosome ${chrom}..."
    wget -P "$VCF_DIR" "$vcf_url"
    wget -P "$VCF_DIR" "$tbi_url"
}

#Function to merge vcf files
merge_vcf(){
    vcf_list=("$@")
    if [ ${#vcf_list[@]} -lt 2 ]; then
    echo " ELEMENT VCF_LIST: ${vcf_list[*]}"
        echo "Less than 2 vcf files, skip merging."
        return
    fi
    merged_output="$OUTPUT_DIR/1000genomes_data_merged.vcf.gz"
    echo "Merging files..."
    bcftools "${vcf_list[@]}" | bgzip -c > "$merged_output"
    echo "Merged VCF file created at $merged_output"
}

if [ "$#" -eq 0 ]; then
  echo "No parameter provided. Downloading all VCF from 1000 genomes...."
  chromosomes=( {1..22} X Y MT )
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
            bcftools view --regions "${chrom}:${positions}" "$VCF_DIR/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" \ | bgzip -c > "$VCF_DIR/ALL.chr${chrom}.filtered.vcf.gz"
            echo "FILTERING..."
        done < "$loci_file"
        echo "Filtering completed! Filtered VCFs saved in $VCF_DIR"
        files_to_merge+=("$VCF_DIR/ALL.chr${chrom}.filtered.vcf.gz")
        ;;
      p)
        sample_list="" #Ligne également ajoutée par chatGPT
        pop_names="${OPTARG}"
        pop_list=$(echo "$pop_names" | tr ',' ' ')
        wget -O "panel.txt" "$panel_url"
        sample_list=$(grep -E "$(echo "${continent_populations[$pop_list]}" | tr ',' '|')" "panel.txt" | awk '{print $1}' | tr '\n' ',' | sed 's/,$//')
        echo "Sample list: $sample_list" #Je fais ça pour tester si j'ai la bonne liste 'échantillon (il me faut une liste séparée par des ",")
        #Ci-dessous de l'experimental
        #echo "$sample_list"
        #for pop in $pop_list; do
          #if [ -n "${continent_populations[$pop]}" ]; then
            #samples=$(grep -E "$(echo "${continent_populations[$pop]}" | tr ',' '|')" "panel.txt" | awk '{print $1}')
            #sample_list+="$samples\n"

            #echo "Sample list: $sample_list"
          #else
            #echo "Population $pop is not valid."
            #exit 1
          #fi
        #done
        #sample_list=$(echo -e "$sample_list" | tr '\n' ',' | sed 's/,$//')
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


###Populations selection
if [ "$pop_names" != "0" ]; then
    mkdir -p "./pop_subset"
    for file in "$VCF_DIR"/*.vcf.gz; do
        echo "Filtering..."
        bcftools view -s "$sample_list" "$file" | bgzip -c > "./pop_subset/$(basename "$file" .vcf.gz).filtered.vcf.gz"
        #Cette version ne fonctionne pas car elle cherche les noms des échantillons et je n'ai pas les noms correspondant à chaque indiv au sein des pops.
    done
    
#Adding files to "files_to_merge"
    for file in ./pop_subset/*.vcf.gz; do
        files_to_merge+=( "$file" )
    done
  else
    for file in "$VCF_DIR"/*.vcf.gz; do
      files_to_merge+=( "$file" )
    done
fi

# Création de l'output (merge des vcf en 1 seul)
printf "Files to merge: %s\n" "${files_to_merge[@]}" 
echo "Writing a unique vcf..."
if  merge_vcf "${files_to_merge[@]}" ; then
    echo "Success!"
    else
        echo "ERROR during merging."
fi
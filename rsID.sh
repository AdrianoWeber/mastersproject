#!/bin/bash

## This script rewrite a vcf by adding rsID from a database loaded locally.
## It is used as: rsID.sh [output_directory] [path_to_database] [-options...].
## See rsID.sh -h for more detail on usage.

## This script is part of the master thesis of Adriano Weber at the University of Geneva (09/2024-09/2025).
## This script is licensed under the GNU General Public License v3.0 or later.
## See https://www.gnu.org/licenses/gpl-3.0.en.html for more details.

#USAGE
usage(){
  echo "Usage: $0 -options <path_to_file.zip>"
  echo "OPTIONS:"
  echo "-o: Where to write the output file(s)."
  echo "-s: If provided then only the final annotated file will be kept"
  echo "-d: Where is the database since it's a local use."
  echo "-h: Help of the script."
}

if [ "$#" -eq 0 ]; then
  echo "No parameter provided."
  usage
  exit 1
fi

#Check if bcftools is installed
if ! command -v bcftools &> /dev/null; then
    echo "bcftools is not installed. Install it to proceed."
    exit 1
fi

#Initialisation of directory parameters
SILENT=0

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR=.
fi

if [ -z "$DB_DIR" ]; then
    DB_DIR=.
fi

#OPTIONS setting
while getopts ":ho:sd:" opt; do
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
      s)
        SILENT=1
        ;;
      d)
        DB_DIR="${OPTARG}"
        ;;
      *)
        usage
        ;;
    esac
done

shift $((OPTIND - 1)) ##Cette ligne est très sus

###Setting of variables
ZIP_FILE="$1" 
if [ ! -f "$ZIP_FILE" ]; then
    echo "Error with $ZIP_FILE no such file."
    usage
    exit 1
fi

BCF_LIST_FILE="bcf_files.txt"
tmp_files=()

###Opening of the zip file
echo "Unzipping input file..."
unzip "$ZIP_FILE" -d "$OUTPUT_DIR"

echo "Searching BCF files..."
find "$OUTPUT_DIR" -name "*.gz" > "$OUTPUT_DIR/$BCF_LIST_FILE" #This line read only the zipped vcf, it could be useful to extend to bcf files too
tmp_files+=( "$(find "$OUTPUT_DIR" -name "*.gz")" )
tmp_files+=( "$OUTPUT_DIR/$BCF_LIST_FILE" )

###Indexing BCF files
echo "Indexing BCF files..."
while IFS= read -r BCF_FILE; do
    echo "Indexation of $BCF_FILE..."
    if bcftools index "$BCF_FILE"; then
        tmp_files+=( "$BCF_FILE" "$BCF_FILE.csi" )
    else
        echo "Error during indexing $BCF_FILE"
        rm -f "${tmp_files[@]}"
        exit 1
    fi
done < "$OUTPUT_DIR/$BCF_LIST_FILE"

if [ ! -s "$OUTPUT_DIR/$BCF_LIST_FILE" ]; then
  echo "No BCF file was found."
  rm -f "${tmp_files[@]}"
  echo "All temporary files have been deleted."
  exit 1
fi

###Merging of BCF files
echo "Merging the BCF files using bcftools..."
if bcftools merge -o "$OUTPUT_DIR/merged_output.bcf" -O b -l "$OUTPUT_DIR/$BCF_LIST_FILE"; then
  tmp_files+=( "$OUTPUT_DIR/merged_output.bcf" )
  echo "Success ! The merged file was saved as 'merged_output.bcf'."
  tabix "$OUTPUT_DIR/merged_output.bcf"
  tmp_files+=("$OUTPUT_DIR/merged_output.bcf.csi")
else
  echo "Error during the merger of BCF files."
  rm -f "${tmp_files[@]}"
  echo "All temporary files have been deleted."
  exit 1
fi

###Rewriting of CHROM column
echo "Rewriting CHROM column to standarize it...."
if bcftools annotate --rename-chrs <(echo -e "chr1\t1\nchr2\t2\nchr3\t3\nchr4\t4\nchr5\t5\nchr6\t6\nchr7\t7\nchr8\t8\nchr9\t9\nchr10\t10\nchr11\t11\nchr12\t12\nchr13\t13\nchr14\t14\nchr15\t15\nchr16\t16\nchr17\t17\nchr18\t18\nchr19\t19\nchr20\t20\nchr21\t21\nchr22\t22\nchrX\tX\nchrY\tY\nchrM\tM") "$OUTPUT_DIR/merged_output.bcf" -o "$OUTPUT_DIR/rewrited_output.bcf" -O b; then
  tmp_files+=( "$OUTPUT_DIR/rewrited_output.bcf" )
  echo "Success! The rewrited files was saved as 'rewrited_output.bcf'"
  tabix "$OUTPUT_DIR/rewrited_output.bcf"
  tmp_files+=("$OUTPUT_DIR/rewrited_output.bcf.csi")
else
  echo "Error during the writing of BCF file."
  rm -f "${tmp_files[@]}"
  echo "All temporary files have been deleted."
  exit 1
fi


###Annotation of BCF
echo "Use of 'clinvar.vcf.gz' 2024-09-18 update from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
echo "Adding IDs with ClinVar...."

if bcftools annotate -a "$DB_DIR" -c ID -o "$OUTPUT_DIR/annotate_output.bcf" -O b "$OUTPUT_DIR/rewrited_output.bcf"; then
  echo "Success! The annotated file was saved as 'annotate_output.bcf'"
  tabix "$OUTPUT_DIR/annotate_output.bcf"
else
  echo "Error during the annotation of BCF file. Check access to the database."
  rm -f "${tmp_files[@]}"
  echo "All temporary files have been deleted."
  exit 1
fi

if [ "$SILENT" -eq 1 ]; then
  echo "Destruction of temporary files..."
  rm -f "${tmp_files[@]}"
  echo "All temporary files have been deleted."
fi
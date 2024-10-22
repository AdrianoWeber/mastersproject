#!/bin/bash
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <chemin_du_fichier_zip>"
  exit 1
fi

ZIP_FILE="$1"  #Ici cela doit être la direction du fichier zip
OUTPUT_DIR="${2:-.}" #Cela dit de prendre la deuxième option si fournie, sinon ce sera la direction actuelle
BCF_LIST_FILE="bcf_files.txt"

if ! command -v bcftools &> /dev/null
then
    echo "bcftools n'est pas installé. Veuillez l'installer avant de continuer."
    exit 1
fi

echo "Décompression du fichier zip..."
unzip "$ZIP_FILE" -d "$OUTPUT_DIR"

echo "Recherche des fichiers BCF..."
find "$OUTPUT_DIR" -name "*.gz" > "$BCF_LIST_FILE" #Ici modifier pour pouvoir utiliser des bcf

# Indexation des fichiers BCF
echo "Indexation des fichiers BCF..."
while IFS= read -r BCF_FILE; do
    echo "Indexation de $BCF_FILE..."
    bcftools index "$BCF_FILE"
    if [ $? -ne 0 ]; then
        echo "Erreur lors de l'indexation de $BCF_FILE"
        exit 1
    fi
done < "$OUTPUT_DIR/$BCF_LIST_FILE"

if [ ! -s "$BCF_LIST_FILE" ]; then
  echo "Aucun fichier BCF trouvé."
  exit 1
fi

echo "Fusion des fichiers BCF avec bcftools..."
bcftools merge -o merged_output.bcf -O b -l "$BCF_LIST_FILE"
if [ $? -eq 0 ]; then
  echo "Fusion terminée avec succès ! Le fichier final est 'merged_output.bcf'."
else
  echo "Erreur lors de la fusion des fichiers BCF."
  exit 1
fi
tabix merged_output.bcf

###Changement de la colonne CHROM pour correspondre à clinvar
echo "Réecriture de la colonne CHROM...."
bcftools annotate --rename-chrs <(echo -e "chr1\t1\nchr2\t2\nchr3\t3\nchr4\t4\nchr5\t5\nchr6\t6\nchr7\t7\nchr8\t8\nchr9\t9\nchr10\t10\nchr11\t11\nchr12\t12\nchr13\t13\nchr14\t14\nchr15\t15\nchr16\t16\nchr17\t17\nchr18\t18\nchr19\t19\nchr20\t20\nchr21\t21\nchr22\t22\nchrX\tX\nchrY\tY\nchrM\tM") merged_output.bcf -o rewrited_output.bcf -O b
if [ $? -eq 0 ]; then
  echo "Réécriture réussie! Le fichier réécrit a été enregistré sous 'rewrited_output.bcf'"
else
  echo "Erreur lors de la réécriture du fichier bcf."
  exit 1
fi
tabix rewrited_output.bcf

###Maintenant l'annotation du bcf avec des ID
echo "Usage de clinvar.vcf.gz update de 2024-09-18. https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
echo "Ajout des IDs selon ClinVar...."

bcftools annotate -a /mnt/c/Users/adria/Documents/Unige/Masters_project/Data_PacBio/clinvar.vcf.gz -c ID -o annotate_output.bcf -O b rewrited_output.bcf
if [ $? -eq 0 ]; then
  echo "Annotation réussie! Le fichier annoté a été enregistré sous 'annotate_output.bcf'"
else
  echo "Erreur lors de l'annotation du fichier bcf."
  exit 1
fi

tabix annotate_output.bcf
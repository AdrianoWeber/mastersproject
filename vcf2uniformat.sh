#! /bin/bash

##This script convert vcf file to UNIFORMAT style (format used by Gene[Rate])
## Usage: vcf2uniformat.sh [vcf_file] [output_file]

#Be aware that the software bcftools is needed for this script
command -v bcftools >/dev/null 2>&1 || { echo "bcftools is required but not installed. Exiting."; exit 1; }

#Loading vcf file thanks to options (first for input vcf file, second for output .unif file)
vcf_file=$1
out_filename=$2
#Getting samples list
IFS=$'\t' read -r -a samples <<< "$(bcftools view -h "$vcf_file" | tail -n 1 | cut -f10-)"
echo "Samples [1]: ${samples[*]}"

#Function to format the sample temporary file to an UNIFORMAT-like style. parameters ar file name and the length of haplotypes (comma separated)
format_sample(){
    sample_file=$1
    haplotypes=$2
    haplotypes=$(echo $haplotypes | tr "," "\t")
    sample_name=$(head -n 1 $sample_file | cut -d "]" -f2 | cut -d ":" -f1)
    start=2
    count=1
    out_line=$sample_name
    for line in $haplotypes; do
        hap1=$(tail -n "+$start" "$sample_file" | head -n "$line" | sed 's/\//\|/' | cut -d "|" -f1 | tr -d "\n" )
        hap2=$(tail -n "+$start" "$sample_file" | head -n "$line" | sed 's/\//\|/' | cut -d "|" -f2 | tr -d "\n" )
        out_line+="\tL$count*$hap1,L$count*$hap2"
        start=$((start+line+1))
        count=$((count+1))
    done
    echo "$out_line"
}

#Making tempo files with GT info
if [ -f "$out_filename.unif" ]; then
    read -p "Do you want to overwrite [y/n]? " answer
    if [[ "$answer" =~ ^[Yy]$ ]]; then
        echo "Overwriting..."
        rm $out_filename.unif
    else
        echo "Not overwriting. End of script"
        exit 1
    fi
fi
for samp in "${samples[@]}"; do
    echo "Sample [2]: $samp"
    bcftools query -f "[%GT]\n" -H -s $samp $vcf_file > "$samp.tmp"

    new_line=$(format_sample "$samp.tmp" "9,3,6")
    echo "Nouvelle ligne [3]: $new_line"
    echo -e $new_line >> $out_filename.tmp #the 9,3,6 part is the length of haplotypes for C19, C9 and C8. It must be changed if adapting the script.
    rm "$samp.tmp"
done

##Renaming the file with awk
awk -v info_out="$out_filename.info" '
BEGIN {
    label_count = 1
    FS = "\t"
    OFS = "\t"
}

{
    printf "%s", $1

    for (i=2; i<=NF; i++) {
        n = split($i, codes, ",")

        for (j=1; j<=n; j++) {
            if (match(codes[j], /^([^*]+\*)(.*)$/, arr)) {
                prefix = arr[1]
                code = arr[2]

                if (!(code in code_map)) {
                    label = "H" label_count++
                    code_map[code] = label
                    reverse_map[label] = code
                }

                codes[j] = prefix code_map[code]
            }
        }

        printf "%s%s", OFS, codes[1]
        for (j=2; j<=n; j++) {
            printf ",%s", codes[j]
        }
    }
    print ""
}

END {
    # Write Hn-to-code mapping
    for (label in reverse_map) {
        print label "\t" reverse_map[label] > info_out
    }
}
' $out_filename.tmp > $out_filename.unif

rm $out_filename.tmp
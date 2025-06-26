# mastersproject
This repository is part of my master's project at Univerity of Geneva.
It contains multiple scripts used during the project.
If you have any questions, please contact adriano.weber@etu.unige.ch

The main features are:

- [load_phased1KGP.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/load_phased1KGP.sh): A shell script to dowload high coverage data from the 1000 genomes project(1KGP) for the [phased VCF files](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). This script allows selection of loci and continent.
- [load_30x1KGP.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/load_30x1KGP.sh): A shell script to dowload high coverage data from the 1000 genomes project(1KGP) from [NYGC_GATK](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20190425_NYGC_GATK). This script allows selection of loci and continent.
- [load_release31KGP.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/load_release31KGP.sh): A shell script to dowload data from the 1000 genomes project(1KGP) for the [phase3 release of variant calls](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). This script allows selection of loci and continent.
- [vcf2uniformat.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/vcf2uniformat.sh): A shell script to transform a VCF file into `UNIFORMAT` style. This script is designed to get inpute files for the [Gene\[rate\] tool](https://hla-net.eu/tools/). **Beware if adapting the script, since it depends on the length of haplotypes.**
- [handlingmisdata.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/handlingmisdata.sh): A shell script completing missing data due to multi-samples vcf merging. It must have bam files to query the coverage depth of every missing data with <ins>samtools depth</ins>.
- [star_allele_info.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/star_allele_info.sh): A shell script that select star allele SNPs for a gene in a VCF file and compute the frequency of this allele in a population (MUS by default). **Beware if adapting the script, it is made for a population named "MUS".**
- [plot_plinkpca.R](https://github.com/AdrianoWeber/mastersproject/blob/main/plot_plinkpca.R): A R script to plot the output of a PLINK pca. 
- [rsID.sh](https://github.com/AdrianoWeber/mastersproject/blob/main/rsID.sh): A shell script to annotate VCF files with rsID (secondary so still not fully implemented).
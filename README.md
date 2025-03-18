# mastersproject
This branch is related to the data from 1000Genomes database for the 30xGrChr38 work.
This particuliary branch goes for [phased data](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/).
I am aware that multi branching isn't the cleanest way to keep multiple version of a script but due to time constraint I can't merge all versions in one.
The main features are:

- [Script_1000genomes.sh](https://github.com/AdrianoWeber/mastersproject/blob/1KG30x/Script_1000genome.sh): A shell script to dowload data from the 1000 genomes project(1KGP) from the 3 release. This script allows selection of loci and continent.
- [Script_handlingmisdata.sh](https://github.com/AdrianoWeber/mastersproject/blob/1KG30x/Script_handlingmisdata.sh): A shell script completing missing data due to multi-samples vcf merging. It must have bam files to query the coverage depth of every missing data with <ins>samtools depth</ins>.
- [star_allele_info.sh](https://github.com/AdrianoWeber/mastersproject/blob/1KG30x/star_allele_info.sh): A shell script that select star allele SNPs for a gene in a VCF file and compute the frequency of this allele in a population (MUS by default).

#!/bin/bash
#Script for running Local ancestry with RFMix2
#Juliana Acosta-Uribe 2023

#run nohup ./Local_ancestry.RFMix2.sh > Local_ancestry.RFMix2.log

# Set your variables:
SHAPEIT_HAPS='CLM.EUR.AFR.NAM' #Prefix for file.chr.haps.gz
QUERY_COHORT='TANGL_cohort.txt' 
REF_COHORT='AFR.EUR.NAT_cohorts.txt'
GENETIC_MAP='genetic_map_hg19.txt'
REF_COHORT_MAP=''  #file with the information of which sample belong to which ancestral group

# Needed software
#shapeit2: https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html
#vcftools: https://github.com/vcftools/vcftools
#RFMix2: https://github.com/slowkoni/rfmix

# Prep your files

## 1. Convert SHAPEIT2 output into vcf
for i in {1..22}; do
gunzip ${SHAPEIT_HAPS}.chr${i}.haps.gz
shapeit -convert --input-haps ${SHAPEIT_HAPS}.chr${i} --output-vcf ${SHAPEIT_HAPS}.chr${i}.vcf
bgzip ${SHAPEIT_HAPS}.chr${i}.haps; done

## 2. Extract the "query" and "reference" haplotypes which are to be analyzed. Files need to be indexed woth 'tabix'

for i in {1..22}; do
vcftools --vcf ${SHAPEIT_HAPS}.chr${i}.vcf --keep ${QUERY_COHORT} --recode --recode-INFO-all --stdout | bgzip > ${SHAPEIT_HAPS}.chr${i}_query.vcf.gz
tabix -p vcf ${SHAPEIT_HAPS}.chr${i}_query.vcf.gz 

vcftools --vcf ${SHAPEIT_HAPS}.chr${i}.vcf --keep ${REF_COHORT} --recode --recode-INFO-all --stdout | bgzip > ${SHAPEIT_HAPS}.chr${i}_reference.vcf.gz
tabix -p vcf ${SHAPEIT_HAPS}.chr${i}_reference.vcf.gz

# 3. Recode the genetic map according to RFMix2 recommendations:
## The genetic map file is tab delimited text containing at least 3 columns. The first 3 columns are intepreted as chromosome, physical position in bp, genetic position in cM.
## You can download a genetic map from https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/

for i in {1..22}; do
awk -v chr="${i}" '{ if( $1=chr) {print $1, $,2 $4}}' ${GENETIC_MAP} > genetic_map_chr${i}.RFMIX2.map; done

# 4. Run RFMix2
for i in 22 ; do
rfmix 
-f ${SHAPEIT_HAPS}.chr${i}_query.vcf.gz  
-r ${SHAPEIT_HAPS}.chr${i}_reference.vcf.gz 
-m ${REF_COHORT_MAP}
-g genetic_map_chr${i}.RFMIX2.map 
-o ${SHAPEIT_HAPS}.chr${i}_query.RFMIX2 
--chromosome=${i} 
-n 5 
-e 1 
--reanalyze-reference 
--random-seed=clock 
--n-threads=12; done 

# RFMix2 will generate the following output files:
# *output.msp.tsv*: The most likely assignment of subpopulations per CRF point. Produced by computing the viterbi algorithm on the CRF. The *.msp.tsv* file is condensed such that CRF windows are combined if all query samples are in the sample subpopulations for successive windows. 
# *output.fb.tsv*: The marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point. Produced by computing the forward-backward algorithm on the CRF
# *output.sis.tsv*
# *output.Q*: RFMix2 also reports global ancestry estimated based on local ancestries identified by their algorithm in their standard output file (∗results.rfmix.Q) corresponding to the .Q output files from ADMIXTURE that can be compared to ADMIXTURE generated global ancestry directly. 

#!/bin/bash
#Script for running Local ancestry with RFMix2
#Juliana Acosta-Uribe 2022

#run nohup ./Local_ancestry.RFMix2.sh > Local_ancestry.RFMix2.log

#Prep your files

#1. Convert SHAPEIT2 output into vcf
#for i in {1..22}; do
#gunzip ./SHAPEIT.output/CLM.EUR.AFR.NAM.chr${i}.haps.gz
#shapeit -convert --input-haps ./SHAPEIT.output/CLM.EUR.AFR.NAM.chr${i} --output-vcf ./SHAPEIT.vcf/CLM.EUR.AFR.NAM.chr${i}.vcf
#bgzip ./SHAPEIT.output/CLM.EUR.AFR.NAM.chr${i}.haps; done

#2. Extract the "query" and "reference" haplotypes which are to be analyzed. File needs to be indexed

#for i in {1..22}; do
#vcftools --vcf ./SHAPEIT.vcf/CLM.EUR.AFR.NAM.chr${i}.vcf --keep TANGL_cohort.txt --recode --recode-INFO-all --stdout | bgzip > CLM.EUR.AFR.NAM.chr${i}_Query.vcf.gz
#tabix -p vcf CLM.EUR.AFR.NAM.chr${i}_Query.vcf.gz 

#vcftools --vcf ./SHAPEIT.vcf/CLM.EUR.AFR.NAM.chr${i}.vcf --keep AFR.EUR.NAT_cohorts.txt --recode --recode-INFO-all --stdout | bgzip > CLM.EUR.AFR.NAM.chr${i}_Reference.vcf.gz
#tabix -p vcf CLM.EUR.AFR.NAM.chr${i}_Reference.vcf.gz; done

#3. Recode the genetic map according to RFMix2 recommendations:
#The genetic map file is tab delimited text containing at least 3 columns. The first 3 columns are intepreted as chromosome, physical position in bp, genetic position in cM.

#for i in {1..22}; do
#awk -v chr="${i}" '{print chr, $1, $3}' ~/1000GP/GeneticMaps/genetic_map_chr${i}_combined_b37.txt > genetic_map_chr${i}_combined_b37.RFMIX2.map
#sed -i "/position/d" genetic_map_chr${i}_combined_b37.RFMIX2.map; done

#4. Run RFMix2
for i in 22 ; do
rfmix \
-f CLM.EUR.AFR.NAM.chr${i}_Query.vcf.gz \
-r CLM.EUR.AFR.NAM.chr${i}_Reference.vcf.gz \
-m CLM.EUR.AFR.NAM.Reference.sample_map \
-g genetic_map_chr${i}_combined_b37.RFMIX2.map \
-o CLM.EUR.AFR.NAM.chr${i}_Query.RFMIX2 \
--chromosome=${i} \
-n 5 \
-e 1 \
--reanalyze-reference \
--random-seed=clock \
--n-threads=12; done 

#RFMix2 will generate the following output files:
# *output.msp.tsv*: The most likely assignment of subpopulations per CRF point. Produced by computing the viterbi algorithm on the CRF. The *.msp.tsv* file is condensed such that CRF windows are combined if all query samples are in the sample subpopulations for successive windows. 
# *output.fb.tsv*: The marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point. Produced by computing the forward-backward algorithm on the CRF
# *output.sis.tsv*
# *output.Q*:	 RFMix2 also reports global ancestry estimated based on local ancestries identified by their algorithm in their standard output file (∗results.rfmix.Q) corresponding to the .Q output files from ADMIXTURE that can be compared to ADMIXTURE generated global ancestry directly. 

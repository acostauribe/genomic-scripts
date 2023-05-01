#!/bin/bash
#Script for filtering variants annotated with ANNOVAR
#Juliana Acosta-Uribe 2022

# run as follows: 
#chmod u+x filter_variants.sh
#./filter_variants.sh | tee filter_variants.log

# 1. SET UP VARIABLES
## Set up the name of your vcf file. your 'ANNOVAR file' will be named based on the vcf.
vcf_file='redlat_exomes.neurodegeneration_genes148'
annovar_file=${vcf_file}.hg38_multianno
## Set up the path to the ANNOVAR database
annovar_database='/home/acostauribe/bin/annovar/humandb/'
## Set up reference genome
genome_version='hg38'

# 2. RUN ANNOVAR
echo $(date)
echo "Using ANNOVAR to annotate "${vcf_file}.vcf 

## A. If you want to remove monomorphic variants, run the next two lines
#mv ${vcf_file}.vcf > ${vcf_file}.original.vcf
#bcftools view --min-ac 1  ${vcf_file}.original.vcf >  ${vcf_file}.vcf

table_annovar.pl ${vcf_file}.vcf ${annovar_database} --buildver ${genome_version} --outfile ${vcf_file} --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

## B. Edit the ANNOVAR results header 
echo "Editing ANNOVAR header" 

## C. Get the sample Ids from the vcf
grep '#CHROM' ${annovar_file}.vcf > ${annovar_file}.chrom_line
head -1 ${annovar_file}.txt > ${annovar_file}.header
sed -i 's/ /\t/g' ${annovar_file}.header
cut -f -120 ${annovar_file}.header > ${annovar_file}.header_120
paste ${annovar_file}.header_120 ${annovar_file}.chrom_line > ${annovar_file}.header_new


## D. Edit values in ANNOVAR column IDs to make it more user friendly
awk -v OFS='\t' '{sub("AF","GnomAD_v3_all",$24);print}' ${annovar_file}.header_new > temp
mv temp ${annovar_file}.header_new 
sed -i 's/AF_raw/GnomAD_v3_raw/g' ${annovar_file}.header_new
sed -i 's/AF_male/GnomAD_v3_Male/g' ${annovar_file}.header_new
sed -i 's/AF_female/GnomAD_v3_Female/g' ${annovar_file}.header_new
sed -i 's/AF_afr/GnomAD_v3_African/g' ${annovar_file}.header_new
sed -i 's/AF_ami/GnomAD_v3_Amish/g' ${annovar_file}.header_new
sed -i 's/AF_amr/GnomAD_v3_American/g' ${annovar_file}.header_new
sed -i 's/AF_asj/GnomAD_v3_Ashkenazi/g' ${annovar_file}.header_new
sed -i 's/AF_eas/GnomAD_v3_EastAsian/g' ${annovar_file}.header_new
sed -i 's/AF_fin/GnomAD_v3_Finnish/g' ${annovar_file}.header_new
sed -i 's/AF_nfe/GnomAD_v3_European/g' ${annovar_file}.header_new
sed -i 's/AF_oth/GnomAD_v3_Other/g' ${annovar_file}.header_new
sed -i 's/AF_sas/GnomAD_v3_SouthAsian/g' ${annovar_file}.header_new
sed -i 's/ALL.sites.2015_08/1000GP_2015_All/g' ${annovar_file}.header_new
sed -i 's/AFR.sites.2015_08/1000GP_2015_African/g' ${annovar_file}.header_new
sed -i 's/EUR.sites.2015_08/1000GP_2015_European/g' ${annovar_file}.header_new
sed -i 's/EAS.sites.2015_08/1000GP_2015_EastAsian/g' ${annovar_file}.header_new
sed -i 's/SAS.sites.2015_08/1000GP_2015_SouthAsian/g' ${annovar_file}.header_new
sed -i 's/AMR.sites.2015_08/1000GP_2015_American/g' ${annovar_file}.header_new
sed -i 's/CLNALLELEID/ClinVar_AlleleID/g' ${annovar_file}.header_new
sed -i 's/CLNDN/ClinVar_Disease/g' ${annovar_file}.header_new
sed -i 's/CLNDISDB/ClinVar_DiseaseDatabase/g' ${annovar_file}.header_new
sed -i 's/CLNREVSTAT/ClinVar_ReviewStatus/g' ${annovar_file}.header_new
sed -i 's/CLNSIG/ClinVar_Significance/g' ${annovar_file}.header_new

## E. Paste edited header with annotations
sed '1d' ${annovar_file}.txt > ${annovar_file}.header_removed
cat ${annovar_file}.header_new ${annovar_file}.header_removed > ${annovar_file}.txt

## F. Count number of variants
echo "Number of variants in " ${annovar_file}.txt
wc -l ${annovar_file}.txt | awk '{print $1-1}' #the awk command substracts the line of the header

## G. Get protein altering variants
echo "Removing non protein coding variants" 
awk -v OFS='\t' 'NR==1; NR > 1{ if(($6 == "exonic") || ($6 == "exonic;splicing") || ($6 == "splicing")) { print } }' ${annovar_file}.txt > ${annovar_file}.coding.txt

echo "Number of variants in" ${annovar_file}.coding.txt
wc -l ${annovar_file}.coding.txt | awk '{print $1-1}'

## H. Exclude synonymous mutations
echo "Removing synonymous variants"
awk -v OFS='\t' 'NR==1; NR > 1{ if(($9 != "synonymous SNV") && ($6 != "unknown")) { print } }' ${annovar_file}.coding.txt > ${annovar_file}.coding.non-syn.txt

echo "Number of variants in" ${annovar_file}.coding.non-syn.txt
wc -l ${annovar_file}.coding.non-syn.txt | awk '{print $1-1}'

## I. Extract variants in the top 25 genes
#echo "Extracting variants in 25 genes for manual curation"
#awk -F "\t" 'NR==1; NR > 1{ if(($7=="CHCHD10") || ($7=="CHMP2B") || ($7=="CSF1R") || ($7=="FUS") || ($7=="GRN") || ($7=="HNRNPA1") || ($7=="HNRNPA2B1") || ($7=="LRRK2") || ($7=="OPTN") || ($7=="PRNP") || ($7=="SNCA") || ($7=="SQSTM1") || ($7=="TARDBP") || ($7=="UBQLN2") || ($7=="VCP") || ($7=="ABCA7") || ($7=="APOE") || ($7=="APP") || ($7=="MAPT") || ($7=="PSEN1") || ($7=="PSEN2") || ($7=="SORL1") || ($7=="TBK1") || ($7=="TREM2") || ($7=="APOE")) { print } }' ${annovar_file}.coding.non-syn.txt > ${annovar_file}.coding.non-syn.top25.txt

#echo "Number of variants in" ${annovar_file}.coding.non-syn.top25.txt
#wc -l ${annovar_file}.coding.non-syn.top25.txt | awk '{print $1-1}'

## Extract variants in the top 10 genes
echo "Extracting variants in 10 genes for manual curation"
awk -v OFS='\t''NR==1; NR > 1{ if(($7=="FUS") || ($7=="GRN") || ($7=="TARDBP") || ($7=="VCP") || ($7=="APOE") || ($7=="APP") || ($7=="MAPT") || ($7=="PSEN1") || ($7=="PSEN2") || ($7=="TBK1")) { print } }' ${annovar_file}.coding.non-syn.txt > ${annovar_file}.coding.non-syn.top10.txt

echo "Number of variants in" ${annovar_file}.coding.non-syn.top10.txt
wc -l ${annovar_file}.coding.non-syn.top10.txt | awk '{print $1-1}'

## Extract variants that have been reported as "Pathogenic" or "Likely pathogenic" in ClinVar
echo "Extracting variants that have been reported as pathogenic in ClinVar"
grep 'athogenic' ${annovar_file}.coding.non-syn.txt | grep -v 'Conflicting'> clinvar_pathogenic
cat ${annovar_file}.header_new clinvar_pathogenic > ${annovar_file}.coding.non-syn.ClinVar-Pathogenic.txt

echo "Number of variants in " ${annovar_file}.coding.non-syn.ClinVar-Pathogenic.txt
wc -l ${annovar_file}.coding.non-syn.ClinVar-Pathogenic.txt | awk '{print $1-1}'

## Extract variants with low allelic frequency and high CADD score
echo "Extracting variants for manual curation"
echo "Warning! Variants with missing MAF or CADD information will be retained for manual curation"
# MAF and CADD tresholds need to be defined inside the awk command
# If the value is a '.' or a blank, it will be printed out as well. when a '<=0' statement is present
# Columns are 16=MAF in ExAC_all, 25=MAF in 1000GP_all, 74=CADD_phred
awk -v OFS='\t' -v MAF='0.001' -v CADD='20' 'NR==1; NR > 1{ if(($16 <= MAF) && ($25 <= MAF) && ($74 >= CADD || $74<=0)) {print}}' \
${annovar_file}.coding.non-syn.txt > ${annovar_file}.coding.non-syn.manual-curation.txt

echo "Number of variants in " ${annovar_file}.coding.non-syn.manual-curation.txt
wc -l ${annovar_file}.coding.non-syn.manual-curation.txt | awk '{print $1-1}'

# Delete what we dont need
rm ${vcf_file}.avinput
rm ${annovar_file}.header*
rm ${annovar_file}.chrom_line
rm clinvar_pathogenic

# Create a directory
#echo "Creating a directory for annotations"
#mkdir ANNOVAR
#mv ${annovar_file}.* ./ANNOVAR
#mv filter_variants.* ./ANNOVAR
#mv ANNOVAR ANNOVAR_$(date "+%Y%m%d")

echo "Finished" $(date)
#!/bin/bash
#Script to extract and annotate a gene panel
#Juliana Acosta-Uribe 2023

# PART ONE: PREPARE YOUR FILES

## Set up the files that will be used
vcf_file='redlat_array' #name of vcf file. Must be gzipped and end in 'vcf.gz'
panel='neurodegeneration_panel_148.txt' #must be in bed format: 'chr', 'start bp', 'end bp'. Chromosome format must match vcf style
samples_passQC='PreimputationQC_pass_samples.one-column.txt' #one column with the IDs of samples that passed QC. Sample ID must match ID in vcf
annovar_file=${vcf_file}_panel.hg38_multianno
annovar_database_PATH='/home/acostauribe/bin/annovar/humandb/'
plink_file='redlat_array.qc' #give a prefix for the plink files that will be created (option A) or for the preexisting plink file (option B)
cohorts='sites.txt' # file where each line is the name of a file that contains the individuals from a site/group. Individual IDs need to be in plink FAM-ID IND-ID


echo "Starting script at" $(date)

# PART TWO: ANNOTATION OF THE VCF 

## 1. Extract target regions from vcf file

#vcftools --gzvcf ${vcf_file}.vcf.gz --bed ${panel} --keep  ${samples_passQC} --recode --recode-INFO-all --out ${vcf_file}_panel
#mv ${vcf_file}_panel.recode.vcf ${vcf_file}_panel.vcf 
#bcftools view --min-ac 1  ${vcf_file}_panel.vcf > tmp
#mv tmp ${vcf_file}_panel.vcf

### gzip and index
bgzip ${vcf_file}_panel.vcf
tabix -p vcf ${vcf_file}_panel.vcf.gz

## 2. Use ANNOVAR to annotate the vcf. [expects vcf aligned to hg38]

echo "Using ANNOVAR to annotate "${vcf_file}_panel.vcf 
table_annovar.pl ${vcf_file}_panel.vcf.gz ${annovar_database_PATH} --buildver hg38 --outfile ${vcf_file}_panel --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

### Edit the ANNOVAR results header 
echo "Editing ANNOVAR header" 

### Get the sample IDs from the vcf
grep '#CHROM' ${annovar_file}.vcf > ${annovar_file}.chrom_line
head -1 ${annovar_file}.txt > ${annovar_file}.header
sed -i 's/ /\t/g' ${annovar_file}.header
cut -f -120 ${annovar_file}.header > ${annovar_file}.header_120
paste ${annovar_file}.header_120 ${annovar_file}.chrom_line > ${annovar_file}.header_new

### Edit values in ANNOVAR column IDs to make it more user friendly
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

### Paste edited header with annotations
sed '1d' ${annovar_file}.txt > ${annovar_file}.header_removed
cat ${annovar_file}.header_new ${annovar_file}.header_removed > ${annovar_file}.txt

### Count number of columns
#head -1 ${annovar_file}.txt | awk -F'\t' '{print NF}' 
## Count number of variants
echo "Number of variants in " ${annovar_file}.txt
wc -l ${annovar_file}.txt | awk '{print $1-1}' #the awk command substracts the line of the header

## 3. Extract variants of interest

### Get Exonic and Splicing variants. 
echo "Retaining Exonic and Splicing variants" 
awk -v OFS='\t' 'NR==1; NR > 1{ if(($6 == "exonic") || ($6 == "exonic;splicing") || ($6 == "splicing")) { print } }' ${annovar_file}.txt > ${annovar_file}.exonic.txt
echo "Number of variants in" ${annovar_file}.exonic.txt
wc -l ${annovar_file}.exonic.txt | awk '{print $1-1}'

### Exclude synonymous mutations
echo "Removing synonymous variants"
awk -v OFS='\t' 'NR==1; NR > 1{ if(($9 != "synonymous SNV") && ($6 != "unknown")) { print } }' ${annovar_file}.exonic.txt > ${annovar_file}.exonic.non-syn.txt
echo "Number of variants in" ${annovar_file}.exonic.non-syn.txt
wc -l ${annovar_file}.exonic.non-syn.txt | awk '{print $1-1}'

### Extract variants in the top 25 genes of neurodegeneration
#echo "Extracting variants in 25 genes for manual curation"
#awk -v OFS='\t' 'NR==1; NR > 1{ if(($7=="CHCHD10") || ($7=="CHMP2B") || ($7=="CSF1R") || ($7=="FUS") || ($7=="GRN") || ($7=="HNRNPA1") || ($7=="HNRNPA2B1") || ($7=="LRRK2") || ($7=="OPTN") || ($7=="PRNP") || ($7=="SNCA") || ($7=="SQSTM1") || ($7=="TARDBP") || ($7=="UBQLN2") || ($7=="VCP") || ($7=="ABCA7") || ($7=="APOE") || ($7=="APP") || ($7=="MAPT") || ($7=="PSEN1") || ($7=="PSEN2") || ($7=="SORL1") || ($7=="TBK1") || ($7=="TREM2") || ($7=="APOE")) { print } }' ${annovar_file}.exonic.non-syn.txt > ${annovar_file}.exonic.non-syn.top25.txt
#echo "Number of variants in" ${annovar_file}.exonic.non-syn.top25.txt
#wc -l ${annovar_file}.exonic.non-syn.top25.txt | awk '{print $1-1}'

### Extract variants in the top 10 genes of neurodegeneration
#echo "Extracting variants in 10 genes for manual curation"
#awk -v OFS='\t' 'NR==1; NR > 1{ if(($7=="FUS") || ($7=="GRN") || ($7=="TARDBP") || ($7=="VCP") || ($7=="APOE") || ($7=="APP") || ($7=="MAPT") || ($7=="PSEN1") || ($7=="PSEN2") || ($7=="TBK1")) { print } }' ${annovar_file}.exonic.non-syn.txt > ${annovar_file}.exonic.non-syn.top10.txt
#echo "Number of variants in" ${annovar_file}.exonic.non-syn.top10.txt
#wc -l ${annovar_file}.exonic.non-syn.top10.txt | awk '{print $1-1}'

### Extract variants that have been reported as "Pathogenic" or "Likely pathogenic" in ClinVar
#echo "Extracting variants that have been reported as pathogenic in ClinVar"
#grep 'athogenic' ${annovar_file}.exonic.non-syn.txt | grep -v 'Conflicting'> clinvar_pathogenic
#cat ${annovar_file}.header_new clinvar_pathogenic > ${annovar_file}.exonic.non-syn.ClinVar-Pathogenic.txt
#echo "Number of variants in " ${annovar_file}.exonic.non-syn.ClinVar-Pathogenic.txt
#wc -l ${annovar_file}.exonic.non-syn.ClinVar-Pathogenic.txt | awk '{print $1-1}'

### Extract variants with low allelic frequency and high CADD score
#echo "Extracting variants for manual curation"
#echo "Warning! Variants with missing MAF or CADD information will be retained for manual curation"
### MAF and CADD tresholds need to be defined inside the awk command
### If the value is a '.' or a blank, it will be printed out as well. when a '<=0' statement is present
### Columns are 16=MAF in ExAC_all, 25=MAF in 1000GP_all, 74=CADD_phred
#awk -v OFS='\t' -v MAF='0.001' -v CADD='20' 'NR==1; NR > 1{ if(($16 <= MAF) && ($25 <= MAF) && ($74 >= CADD || $74<=0)) {print}}' \
#${annovar_file}.exonic.non-syn.txt > ${annovar_file}.exonic.non-syn.manual-curation.txt
#echo "Number of variants in " ${annovar_file}.exonic.non-syn.manual-curation.txt
#wc -l ${annovar_file}.exonic.non-syn.manual-curation.txt | awk '{print $1-1}'

## 4. Tidy up!

### Delete what we dont need
rm ${vcf_file}.avinput
rm ${annovar_file}.header*
rm ${annovar_file}.chrom_line
rm clinvar_pathogenic

### Store annonations in their own directory
#echo "Creating a directory for annotations"
#mkdir ANNOVAR
#mv ${annovar_file}.* ./ANNOVAR
#mv filter_variants.* ./ANNOVAR
#mv ANNOVAR ANNOVAR_$(date "+%Y%m%d")


# PART THREE: COUNT CASES AND CONTROLS WITH EACH VARIANT

## I will be counting cases and controls of ${annovar_file}.exonic.non-syn.txt

## 1. Import data into plink 
#====NEEDS TO BE FIXED

#plink --vcf ${vcf_file}.vcf --vcf-half-call m --allow-extra-chr --double-id --keep-allele-order --make-bed --out ${plink_file}
#.bim and .fam files were edited to for variant id and to have sex and phenotype
#plink --vcf ${plink_file} --extract range variants.txt --make-bed --out ${plink_file}_panel

###I had to make the variants.txt file manually
#===================


## Extract your cohorts according to sample name
#grep	'AF'	${plink_file}.fam	>	avila
#grep	'BE'	${plink_file}.fam	>	behrens
#grep	'BN'	${plink_file}.fam	>	bruno
#grep	'BS'	${plink_file}.fam	>	brusco
#grep	'CU'	${plink_file}.fam	>	custodio
#grep	'LO'	${plink_file}.fam	>	Lopera
#grep	'SL'	${plink_file}.fam	>	slachevski
#grep	'TA'	${plink_file}.fam	>	takada

## 2. Count SNPs in all the cohort

plink --bfile ${plink_file} --extract range variants.txt --model --allow-no-sex --out ${plink_file}_panel
#output ends in .model

### Retain 'GENO' counts
awk -v OFS='\t' 'NR==1; NR > 1{ if ($5 == "GENO") { print } }' ${plink_file}.all-variants.model > ${plink_file}.all-variants.geno
sed -i 's/UNAFF/All_unaffected/g' ${plink_file}.all-variants.geno
sed -i 's/AFF/All_affected/g' ${plink_file}.all-variants.geno
### Retain columns 2, 6 and 7
awk -v OFS='\t' '{print $2, $6, $7}' $plink_file.all-variants.geno > ${plink_file}.all-variants.cohort-comparison.txt 

## 3. Count SNPs in sub-cohorts
while read line; do 
    plink --bfile ${plink_file} \
    --extract vars.txt \
    --keep $line \
    --model \
    --allow-no-sex \
    --out ${plink_file}.$line \

    cat <(awk -v OFS='\t' 'BEGIN { print "AFF", "UNAFF"}') <(awk -v OFS='\t' '{ if ($5 == "GENO") {print $6,$7} }' ${plink_file}.$line.model) > ${plink_file}.$line.geno
    sed -i "s/UNAFF/$line.unaffected/g" ${plink_file}.$line.geno 
    sed -i "s/AFF/$line.affected/g" ${plink_file}.$line.geno
    #using double quotes allows to incorporate the variant
    #keep only the genotipic count values
    paste ${plink_file}.all-variants.cohort-comparison.txt ${plink_file}.$line.geno  > ${plink_file}.all-variants.cohort-comparison.temp 
    #add it to the counts from the entire cohort
    mv ${plink_file}.all-variants.cohort-comparison.temp ${plink_file}.all-variants.cohort-comparison.txt
    #delete temp file
    done < $cohorts


echo "Finished" $(date)
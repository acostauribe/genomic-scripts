#!/bin/bash
## Script to annotate a vcf using ANNOVAR
## Juliana Acosta-Uribe Sep. 2025

# This script is set up in three phases. You can turn "on" or "off" the different phases below.
# After setting up what you need, run this script as follows: 
## chmod u+x annotate_variants.sh
## ./annotate_variants.sh | tee annotate_variants.log\


# Provide a vcf:
vcf_file='LATAM_5000_ndeg.chr'
## Must be gzipped and end in 'vcf.gz'

## Chromosomes must be encoded as chr#. 
## If you need to change chromosome designation from # to chr# you need to create a file (e.g. chromosomes.txt) where every line has two columns: 'old name' 'new name' (e.g. 1 chr1)
#bcftools annotate --rename-chrs chromosomes.txt --output-type z file.vcf.gz > file.chr.vcf.gz


# ANNOVAR expects vcf aligned to hg38 and 'table_annovar.pl' must be in your path. 
## Provide a path to your ANNOVAR database
#annovar_database_PATH='/path/to/annovar/humandb/'
annovar_database_PATH='/home/acostauribe/bin/annovar/humandb/'





#================SCRIPT

## ANNOTATE VCF

echo "Starting script on $(date)"

## 1. Using ANNOVAR to annotate the vcf 
echo "Using ANNOVAR to annotate ${vcf_file}"
table_annovar.pl "${vcf_file}.vcf.gz" "${annovar_database_PATH}" --buildver hg38 --outfile "${vcf_file}" --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

## This will produce ${vcf_file}.hg38_multianno.vcf and ${vcf_file}.hg38_multianno.txt
annovar_file="${vcf_file}.hg38_multianno"

## 2. Edit the ANNOVAR results header 
echo "Editing ANNOVAR header" 

## Get the sample IDs from the vcf
grep '#CHROM' "${annovar_file}.vcf" > "${annovar_file}.chrom_line"
bgzip "${annovar_file}.vcf"

head -1 "${annovar_file}.txt" > "${annovar_file}.header"
sed -i 's/ /\t/g' "${annovar_file}.header"
cut -f -120 "${annovar_file}.header" > "${annovar_file}.header_120"
paste "${annovar_file}.header_120" "${annovar_file}.chrom_line" > "${annovar_file}.header_new"
rm "${annovar_file}.chrom_line"

## Edit strings in ANNOVAR column IDs to make it more user friendly
awk -v OFS='\t' '{sub("AF","GnomAD_v3_all",$24);print}' "${annovar_file}.header_new" > temp
mv temp "${annovar_file}.header_new"
## Using sed to make multiple replacements
sed -i -e 's/AF_raw/GnomAD_v3_raw/g' \
-e 's/AF_male/GnomAD_v3_Male/g' \
-e 's/AF_female/GnomAD_v3_Female/g' \
-e 's/AF_afr/GnomAD_v3_African/g' \
-e 's/AF_ami/GnomAD_v3_Amish/g' \
-e 's/AF_amr/GnomAD_v3_American/g' \
-e 's/AF_asj/GnomAD_v3_Ashkenazi/g' \
-e 's/AF_eas/GnomAD_v3_EastAsian/g' \
-e 's/AF_fin/GnomAD_v3_Finnish/g' \
-e 's/AF_nfe/GnomAD_v3_European/g' \
-e 's/AF_oth/GnomAD_v3_Other/g' \
-e 's/AF_sas/GnomAD_v3_SouthAsian/g' \
-e 's/ALL.sites.2015_08/1000GP_2015_All/g' \
-e 's/AFR.sites.2015_08/1000GP_2015_African/g' \
-e 's/EUR.sites.2015_08/1000GP_2015_European/g' \
-e 's/EAS.sites.2015_08/1000GP_2015_EastAsian/g' \
-e 's/SAS.sites.2015_08/1000GP_2015_SouthAsian/g' \
-e 's/AMR.sites.2015_08/1000GP_2015_American/g' \
-e 's/CLNALLELEID/ClinVar_AlleleID/g' \
-e 's/CLNDN/ClinVar_Disease/g' \
-e 's/CLNDISDB/ClinVar_DiseaseDatabase/g' \
-e 's/CLNREVSTAT/ClinVar_ReviewStatus/g' \
-e 's/CLNSIG/ClinVar_Significance/g' "${annovar_file}.header_new"

## Paste edited header with annotations
sed '1d' "${annovar_file}.txt" > "${annovar_file}.header_removed"
cat "${annovar_file}.header_new" "${annovar_file}.header_removed" > "${annovar_file}.txt"
rm "${annovar_file}.header"*

if [[ -f "${annovar_file}.txt" && -s "${annovar_file}.txt" ]]; then
  echo "ANNOVAR file was successfully generated: ${annovar_file}"
else
  echo "VCF annotation failed. Stopping script"
  return 1  # return with a status code of 1 to indicate an error.
fi

remove_genotypes="TRUE"
if [[ $remove_genotypes == "TRUE" ]]; then
  #echo "Removing genotypes from ANNOVAR file for easier filtering. Original file will be preserved as file.hg38_multianno.original.txt"
  annovar_file="${annovar_file%.hg38_multianno}"
  cut -f -128 "${annovar_file}.hg38_multianno.txt" > "${annovar_file}.no-geno.hg38_multianno.txt"
  mv "${annovar_file}.hg38_multianno.txt" "${annovar_file}.hg38_multianno.original.txt"
  annovar_file="${annovar_file}.no-geno.hg38_multianno"
fi

echo "Finished" $(date)


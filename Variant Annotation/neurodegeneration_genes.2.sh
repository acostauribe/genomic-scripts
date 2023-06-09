#!/bin/bash

## Set up the name of your vcf file
vcf_file='exome_neurodegeneration148'
plink_file=${vcf_file}.plink
sites='sites.txt' # file where each line is the name of a file that contains the individuals from a site/group. Individual IDs need to be in plink FAM-ID IND-ID
annovar_file=${vcf_file}.hg38_multianno
annovar_database='/home/acostauribe/bin/annovar/humandb/'

#Notice that this has been optimized for redLat chromosome format in the vcf "chrx". If the chromosome value in the vcf is just the number, line 64 needs to be edited.



#=============================

echo $(date)

#0. Prepare file
vcftools --gzvcf redlat_exomes.vcf.gz \
    --bed neurodegeneration_panel_148.txt \
    --remove exclude_qc_exome.txt \
    --mac 1 \
    --recode \
    --out ${vcf_file}

mv $vcf_file.recode.vcf ${vcf_file}.vcf

plink --vcf ${vcf_file}.vcf \
    --make-bed \
    --out ${plink_file}

plink --bfile ${plink_file} \
    --set-missing-var-ids @:\#[b38]\$1,\$2 \
    --allow-no-sex \
    --make-bed \
    --out ${plink_file}.var-id 

plink --bfile ${plink_file}.var-id \
    --update-ids exome_new_ids.txt \
    --allow-no-sex \
    --make-bed \
    --out ${plink_file}.id

plink --bfile ${plink_file}.id \
    --pheno exome_phenotypes.txt \
    --allow-no-sex \
    --make-bed \
    --out ${plink_file}.id.pheno


#1. Count SNPs in all the cohort

plink --bfile ${plink_file}.id.pheno \
    --mac 1 \
    --model \
    --allow-no-sex \
    --out ${plink_file}.all-variants 
   #output ends in .model

awk -v OFS='\t' 'NR==1; NR > 1{ if ($5 == "GENO") { print } }' ${plink_file}.all-variants.model > ${plink_file}.all-variants.geno
sed -i 's/UNAFF/All_unaffected/g' ${plink_file}.all-variants.geno
sed -i 's/AFF/All_affected/g' ${plink_file}.all-variants.geno

cat <(awk -v OFS='\t' 'BEGIN { print "ID","CHR", "BP", "A1", "A2"}') <(awk -v OFS='\t' '{print "chr"$1":"$4,$1,$4,$5,$6}' $plink_file.id.pheno.bim ) > ${plink_file}.id.pheno.bim.tmp
paste  <(awk -v OFS='\t' '{print $0}' $plink_file.id.pheno.bim.tmp) <(awk -v OFS='\t' '{print $2,$6,$7}' $plink_file.all-variants.geno)  > ${plink_file}.all-variants.cohort-comparison 
rm ${plink_file}.id.pheno.bim.tmp
sed -i 's/chr23/chrX/g' ${plink_file}.all-variants.cohort-comparison


#2. Count SNPs in sub-cohorts
while read line; do plink --bfile ${plink_file}.id.pheno \
    --model \
    --allow-no-sex \
    --keep $line \
    --out ${plink_file}.$line \
  
    cat <(awk -v OFS='\t' 'BEGIN { print "AFF", "UNAFF"}') <(awk -v OFS='\t' '{ if ($5 == "GENO") {print $6,$7} }' ${plink_file}.$line.model) > ${plink_file}.$line.geno
    sed -i "s/UNAFF/$line.unaffected/g" ${plink_file}.$line.geno 
    sed -i "s/AFF/$line.affected/g" ${plink_file}.$line.geno
    #using double quotes allows to incorporate the variant

    #keep only the genotipic count values
    paste ${plink_file}.all-variants.cohort-comparison ${plink_file}.$line.geno  > ${plink_file}.all-variants.cohort-comparison.temp 
    #add it to the counts from the entire cohort
    mv ${plink_file}.all-variants.cohort-comparison.temp ${plink_file}.all-variants.cohort-comparison
    #delete temp file
    done < $sites

#bcftools view --min-ac 1  redlat_array.neurodegeneration_genes148.vcf >  redlat_array.neurodegeneration_genes148.mac1.vcf

echo "Using ANNOVAR to annotate "${vcf_file}.vcf 

table_annovar.pl ${vcf_file}.vcf ${annovar_database} --buildver hg38 --outfile ${vcf_file} --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

## Edit the ANNOVAR results header 
echo "Editing ANNOVAR header" 

# Get the sample Ids from the vcf
grep '#CHROM' ${annovar_file}.vcf > ${annovar_file}.chrom_line
head -1 ${annovar_file}.txt > ${annovar_file}.header
sed -i 's/ /\t/g' ${annovar_file}.header
cut -f -120 ${annovar_file}.header > ${annovar_file}.header_120
paste ${annovar_file}.header_120 ${annovar_file}.chrom_line > ${annovar_file}.header_new


# Edit values in ANNOVAR column IDs to make it more user friendly
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

# Paste edited header with annotations
sed '1d' ${annovar_file}.txt > ${annovar_file}.header_removed
cat ${annovar_file}.header_new ${annovar_file}.header_removed > ${annovar_file}.txt

# Count number of columns
#head -1 ${annovar_file}.txt | awk -F'\t' '{print NF}' 

# Generate an ID per variant line
awk '{print $1":"$2}' ${annovar_file}.txt > ${annovar_file}_ids
# add it to the annovar file
paste ${annovar_file}.txt ${annovar_file}_ids > ${annovar_file}.txt
--> not working




## Count number of variants
echo "Number of variants in " ${annovar_file}.txt
wc -l ${annovar_file}.txt | awk '{print $1-1}' #the awk command substracts the line of the header

## Get protein altering variants
echo "Removing non protein coding variants" 
awk -v OFS='\t' 'NR==1; NR > 1{ if(($6 == "exonic") || ($6 == "exonic;splicing") || ($6 == "splicing")) { print } }' ${annovar_file}.txt > ${annovar_file}.coding.txt

echo "Number of variants in" ${annovar_file}.coding.txt
wc -l ${annovar_file}.coding.txt | awk '{print $1-1}'

## Exclude synonymous mutations
echo "Removing synonymous variants"
awk -v OFS='\t' 'NR==1; NR > 1{ if(($9 != "synonymous SNV") && ($6 != "unknown")) { print } }' ${annovar_file}.coding.txt > ${annovar_file}.coding.non-syn.txt

echo "Number of variants in" ${annovar_file}.coding.non-syn.txt
wc -l ${annovar_file}.coding.non-syn.txt | awk '{print $1-1}'

## Extract variants in the top 25 genes
#echo "Extracting variants in 25 genes for manual curation"
#awk -v OFS='\t' 'NR==1; NR > 1{ if(($7=="CHCHD10") || ($7=="CHMP2B") || ($7=="CSF1R") || ($7=="FUS") || ($7=="GRN") || ($7=="HNRNPA1") || ($7=="HNRNPA2B1") || ($7=="LRRK2") || ($7=="OPTN") || ($7=="PRNP") || ($7=="SNCA") || ($7=="SQSTM1") || ($7=="TARDBP") || ($7=="UBQLN2") || ($7=="VCP") || ($7=="ABCA7") || ($7=="APOE") || ($7=="APP") || ($7=="MAPT") || ($7=="PSEN1") || ($7=="PSEN2") || ($7=="SORL1") || ($7=="TBK1") || ($7=="TREM2") || ($7=="APOE")) { print } }' ${annovar_file}.coding.non-syn.txt > ${annovar_file}.coding.non-syn.top25.txt

#echo "Number of variants in" ${annovar_file}.coding.non-syn.top25.txt
#wc -l ${annovar_file}.coding.non-syn.top25.txt | awk '{print $1-1}'

## Extract variants in the top 10 genes
echo "Extracting variants in 10 genes for manual curation"
awk -v OFS='\t' 'NR==1; NR > 1{ if(($7=="FUS") || ($7=="GRN") || ($7=="TARDBP") || ($7=="VCP") || ($7=="APOE") || ($7=="APP") || ($7=="MAPT") || ($7=="PSEN1") || ($7=="PSEN2") || ($7=="TBK1")) { print } }' ${annovar_file}.coding.non-syn.txt > ${annovar_file}.coding.non-syn.top10.txt

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
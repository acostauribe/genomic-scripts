#!/bin/bash
## Script to extract and annotate a vcf aligned to hg38
## Juliana Acosta-Uribe December 2025

# Run this script as follows: 
## chmod u+x annotate_variants.sh
## ./annotate_variants.sh | tee annotate_variants.log\


## Provide a vcf
## Must be gzipped and end in 'vcf.gz'. Chromosomes must be encoded as chr#. 
vcf_file='test'

## Requires the following software: ANNOVAR and BCFtools (if you need to modify your vcf)
## Provide a path to your ANNOVAR 'table_annovar.pl' and humandb database. 
ANNOVAR="$HOME/bin/annovar/table_annovar.pl"
annovar_database_PATH="$HOME/bin/annovar/humandb/"


## Prepare your VCF 
# Provide path to bcftools
bcftools='/usr/bin/bcftools'

## If you need to change chromosome designation from # to chr# you need to create a file (e.g. chromosomes.txt) where every line has two columns: 'old name' 'new name' (e.g. 1 chr1)
${bcftools} annotate --rename-chrs ~/chromosomes.txt --output-type z $vcf_file.vcf.gz > $vcf_file.chr.vcf.gz
vcf_file=${vcf_file}'.chr'

## Align to reference and normalize. You need to provide a reference fasta file. Chromosome IDs mut match in the fasta and the vcf.
# fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
# The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html. It needs to be in the same directory than the fasta
fasta_file="$HOME/Utilities/hg38.chrm-string.fa.gz"
# index the VCF (requires bgzip-compressed VCF)
bcftools index -t "${vcf_file}.vcf.gz"
${bcftools} norm --check-ref x --fasta-ref ${fasta_file} --output-type z ${vcf_file}.vcf.gz > ${vcf_file}.ref.vcf.gz
# --check-ref warn (w), exclude (x), or set/fix (s)
vcf_file=${vcf_file}'.ref'

## Annotate with AC and AN if the file doesnt have it
${bcftools} +fill-tags ${vcf_file}.vcf.gz --output-type z --output ${vcf_file}.AC.vcf.gz -- -t AC,AN 
vcf_file=${vcf_file}'.AC'

#================SCRIPT

## ANNOTATE VCF

echo "Starting script on $(date)"

## Using ANNOVAR to annotate the vcf 
echo "Using ANNOVAR to annotate ${vcf_file}"

${ANNOVAR} "${vcf_file}.vcf.gz" "${annovar_database_PATH}" \
--buildver hg38 \
--outfile "${vcf_file}" \
--protocol refGene,ensGene,gnomad40_genome,allofus,avsnp151,clinvar_20250721,dbnsfp47a,dbscsnv11 \
--operation g,g,f,f,f,f,f,f \
--nastring . \
--vcfinput \
--remove 

## This will produce ${vcf_file}.hg38_multianno.vcf and ${vcf_file}.hg38_multianno.txt
annovar_file="${vcf_file}.hg38_multianno"

## Edit the ANNOVAR results header 
echo "Editing ANNOVAR header" 

# grab the VCF header line with the sample IDs
grep -m1 '^#CHROM' "${annovar_file}.vcf" > "${annovar_file}.chrom_line"

# compress vcf
bgzip -f "${annovar_file}.vcf"

# make annovar header tab-delimited (normalize any whitespace -> tabs)
head -1 "${annovar_file}.txt" \
| awk 'BEGIN{FS="[ \t]+"; OFS="\t"} {$1=$1}1' \
> "${annovar_file}.header"

# In my output I get three columns that I am still not sure what they are. My vcf starts again on column 'Otherinfo4'.
# find column index of Otherinfo4
col=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Otherinfo4"){print i; exit}}' "${annovar_file}.header")

# build new header: keep cols 1..col-1 from annovar header (annotations) + append VCF chrom line (sample IDs)
paste <(cut -f "1-$((col-1))" "${annovar_file}.header")  "${annovar_file}.chrom_line" > "${annovar_file}.header_new"

# ensure header_new is strictly tab-delimited
awk 'BEGIN{FS="[ \t]+"; OFS="\t"} {$1=$1}1' "${annovar_file}.header_new" > "${annovar_file}.header_new.tmp" 
mv "${annovar_file}.header_new.tmp" "${annovar_file}.header_new"

# cleanup
rm "${annovar_file}.chrom_line" "${annovar_file}.header"

sed -i -e 's/Func\.refGene/RefSeq_region/g' \
   -e 's/Gene\.refGene/RefSeq_gene/g' \
   -e 's/GeneDetail\.refGene/RefSeq_gene_detail/g' \
   -e 's/ExonicFunc\.refGene/RefSeq_exonic_effect/g' \
   -e 's/AAChange\.refGene/RefSeq_AA_change/g' \
   -e 's/Func\.ensGene/Ensembl_region/g' \
   -e 's/Gene\.ensGene/Ensembl_gene/g' \
   -e 's/GeneDetail\.ensGene/Ensembl_gene_detail/g' \
   -e 's/ExonicFunc\.ensGene/Ensembl_exonic_effect/g' \
   -e 's/AAChange\.ensGene/Ensembl_AA_change/g' \
   -e 's/gnomad40_genome_AF/gnomAD.v4_AF/g' \
   -e 's/gnomad40_genome_AF_raw/gnomAD.v4_AF_raw/g' \
   -e 's/gnomad40_genome_AF_XX/gnomAD.v4_AF_XX/g' \
   -e 's/gnomad40_genome_AF_XY/gnomAD.v4_AF_XY/g' \
   -e 's/gnomad40_genome_AF_grpmax/gnomAD.v4_max_AF/g' \
   -e 's/gnomad40_genome_faf95/gnomAD.v4_FAF95/g' \
   -e 's/gnomad40_genome_faf99/gnomAD.v4_FAF99/g' \
   -e 's/gnomad40_genome_fafmax_faf95_max/gnomAD.v4_max_FAF95/g' \
   -e 's/gnomad40_genome_fafmax_faf99_max/gnomAD.v4_max_FAF99/g' \
   -e 's/gnomad40_genome_AF_afr/gnomAD.v4_AF_AFR/g' \
   -e 's/gnomad40_genome_AF_ami/gnomAD.v4_AF_AMI/g' \
   -e 's/gnomad40_genome_AF_amr/gnomAD.v4_AF_AMR/g' \
   -e 's/gnomad40_genome_AF_asj/gnomAD.v4_AF_ASJ/g' \
   -e 's/gnomad40_genome_AF_eas/gnomAD.v4_AF_EAS/g' \
   -e 's/gnomad40_genome_AF_fin/gnomAD.v4_AF_FIN/g' \
   -e 's/gnomad40_genome_AF_mid/gnomAD.v4_AF_MID/g' \
   -e 's/gnomad40_genome_AF_nfe/gnomAD.v4_AF_NFE/g' \
   -e 's/gnomad40_genome_AF_remaining/gnomAD.v4_AF_Other/g' \
   -e 's/gnomad40_genome_AF_sas/gnomAD.v4_AF_SAS/g' \
   -e 's/gvs_all_af/AllOfUs_AF_all/g' \
   -e 's/gvs_max_af/AllOfUs_max_AF/g' \
   -e 's/gvs_afr_af/AllOfUs_AF_AFR/g' \
   -e 's/gvs_amr_af/AllOfUs_AF_AMR/g' \
   -e 's/gvs_eas_af/AllOfUs_AF_EAS/g' \
   -e 's/gvs_eur_af/AllOfUs_AF_EUR/g' \
   -e 's/gvs_mid_af/AllOfUs_AF_MID/g' \
   -e 's/gvs_oth_af/AllOfUs_AF_OTH/g' \
   -e 's/gvs_sas_af/AllOfUs_AF_SAS/g' \
   -e 's/avsnp151/dbSNP151_rsID/g' \
   -e 's/CLNALLELEID/ClinVar_allele_ID/g' \
   -e 's/CLNDN/ClinVar_disease_name/g' \
   -e 's/CLNDISDB/ClinVar_disease_IDs/g' \
   -e 's/CLNREVSTAT/ClinVar_review_status/g' \
   -e 's/CLNSIG/ClinVar_clinical_significance/g' \
   -e 's/ONCDN/Oncogenicity_disease_name/g' \
   -e 's/ONCDISDB/Oncogenicity_disease_IDs/g' \
   -e 's/ONCREVSTAT/Oncogenicity_review_status/g' \
   -e 's/ONC/Oncogenicity/g' \
   -e 's/SCIDN/Somatic_impact_disease_name/g' \
   -e 's/SCIDISDB/Somatic_impact_disease_IDs/g' \
   -e 's/SCIREVSTAT/Somatic_impact_review_status/g' \
   -e 's/SCI/Somatic_clinical_impact/g' \
   -e 's/Aloft_pred/ALoFT_pred/g' \
   "${annovar_file}.header_new"

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

remove_genotypes="FALSE"
if [[ $remove_genotypes == "TRUE" ]]; then
   echo "Removing genotypes from ANNOVAR file for easier filtering." 
   # Original file will be preserved as file.hg38_multianno.original.txt and file.no-geno.hg38_multianno.txt will be generated.
   
   # Get the index of the column that matches Otherinfo1
   Otherinfo1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Otherinfo1"){print i; exit}}' "${annovar_file}.txt")
  
   # If Genotypes are present, extract only the annotations
   if (( Otherinfo1 > 1 )); then
   annovar_file="${annovar_file%.hg38_multianno}"
   cut -f "1-$((Otherinfo1-1))" "${annovar_file}.hg38_multianno.txt" > "${annovar_file}.no-geno.hg38_multianno.txt"
   mv "${annovar_file}.hg38_multianno.txt" "${annovar_file}.hg38_multianno.original.txt"
   annovar_file="${annovar_file}.no-geno.hg38_multianno"
   else
   echo "Annovar file does not contain genotypes"
   fi

fi


echo "Finished" $(date)


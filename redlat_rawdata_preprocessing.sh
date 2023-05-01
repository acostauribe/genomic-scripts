#!/bin/bash
# This is a short script for the data preprocessing to harmonize the different types of data from the same cohort
# Juliana Acosta-Uribe 2022

RAW_dataset='ADFTD.vqsr.snp.indel'
seq='array'
#targets='xgen-exome-hyb-panel-v2-targets-hg38.bed' #make sure the vcf dile and the bed have the chromosomes using the same style
RAW_dataset_exclude_samples='excluded_samples_array.txt'
fasta_file='/home/acostauribe/public_html/Utilities/hg38.fa.gz'
new_ids='array_to_redlatIDs.txt'

# I. Extract targets and ReDLat samples from the original file
vcftools --gzvcf $RAW_dataset.vcf.gz 
--bed $targets 
--remove $RAW_dataset_exclude_samples 
--recode 
--recode-INFO-all 
--out $RAW_dataset.redlat

mv $RAW_dataset.redlat.recode.vcf $RAW_dataset.redlat.vcf
bgzip --threads 4 $RAW_dataset.redlat.vcf
# you can add --threads # to make it faster

# II. Rename samples according to ReDLat sequence IDs
# https://samtools.github.io/bcftools/bcftools.html#reheader

bcftools reheader --samples $new_ids $RAW_dataset.redlat.vcf.gz > redlat.$seq.vcf

# III. bgzip and index vcf
bgzip --threads 4 redlat_$seq.vcf
tabix -p vcf redlat_$seq.vcf #File should be indexed with Tabix

# IV Change chromosome designation from # to chr#
bcftools annotate --rename-chrs chromosomes.txt --output-type z redlat_$seq.vcf.gz > redlat_$seq.chr.vcf.gz 
tabix -p vcf redlat_$seq.chr.vcf.gz 

# IV. Check reference allele and normalize INDELs
bcftools norm --check-ref ws --fasta-ref $fasta_file --output-type z redlat.$seq.chr.vcf.gz > redlat_$seq.chr.temp.vcf
# --check-ref warn (w), exclude (x), or set/fix (s)
# --output-type compressed VCF (z) will bgzip the output
# fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
# The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html

mv redlat_$seq.temp.vcf.gz redlat_$seq.vcf.gz
tabix -p vcf redlat_$seq.vcf 

sftp://kosik.cnsi.ucsb.edu/home/acostauribe/public_html/Utilities/hg38.fa.gz.fai

bcftools view --output-type z --output joint_psomagen.redlat.vcf --threads 8 more each

 --gzvcf ./Exome_Seq/Joint_call_fixed/redlat_exomes.vcf.gz --out genome_vs_array --gzdiff ./ --diff-indv-discordance
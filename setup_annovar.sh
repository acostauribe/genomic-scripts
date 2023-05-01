#!/bin/bash
# This is a script to download, set up and run ANNOVAR to annotate vcf files
# Juliana Acosta-Uribe 2022

## 1. Download annovar- fill registration and it will give you the link.
# https://annovar.openbioinformatics.org/en/latest/user-guide/download/ 
# https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
# https://www.nature.com/articles/nprot.2015.105
# https://statgen.us/files/tutorials/FunctionalAnnotation_Annovar_final.pdf

## 2. Uncompress ANNOVAR
tar -xvfz annovar.latest.tar.gz
# Make sure that all the perl scripts are in the $PATH

## 3. Download databases for annotation:

# A. Gene based annotation
# https://annovar.openbioinformatics.org/en/latest/user-guide/gene/
# FASTA sequences for all annotated transcripts in RefSeq Gene [Updated 20200817 at UCSC]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
# FASTA sequences for all annotated transcripts in Gencode v41 Basic collection [Updated 20220712 at UCSC]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene humandb/

# B. Region based annotation
# https://annovar.openbioinformatics.org/en/latest/user-guide/region/
perl annotate_variation.pl -buildver hg38 -downdb cytoBand humandb/
# Transcription factor binding site
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar tfbsConsSites humandb/
# Identify variants disrupting microRNAs and snoRNAs
perl annotate_variation.pl -build hg38 -downdb wgRna humandb/

# C. Filter based annotation
# i. Population frequency:
# ExAC. 65000 exome allele allele frequency data for ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian)). version 0.3. Left normalization done. [Updated 20151129]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 humandb/ 
# gnomAD exome collection (v3), with "AF AF_popmax AF_male AF_female AF_raw AF_afr AF_sas AF_amr AF_eas AF_nfe AF_fin AF_asj AF_oth non_topmed_AF_popmax non_neuro_AF_popmax non_cancer_AF_popmax controls_AF_popmax" header [Updated 20191104]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad30_genome humandb/
# ii. Variant ID.
# dbSNP with allelic splitting and left-normalization [dbSNP150 Update 20170929]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 humandb/ 
# Clinvar [Updated 20220320]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20220320 humandb/
# iii. Deleteriousness scores
# whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 4 [Updated 20210710]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp42c humandb/
# dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest [Updated 20151218]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 humandb/
# Mendelian Clinically Applicable Pathogenicity (M-CAP) Score scores for non-synonymous variants. 
# http://bejerano.stanford.edu/mcap/index.html [updated 20161205]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar mcap13 humandb/
# REVEL scores for non-synonymous variants 
# https://sites.google.com/site/revelgenomics/about?pli=1 [updated 20161205]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar revel humandb/
# InterVar: clinical interpretation of missense variants [Updated 20180325]
perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118 humandb/


## 4. Annotate your file

# The table_annovar.pl argument is a command used to annotate your input and generate a tab-delimited output that contains representative columns for each of the annotations.
# Here is an example of a script to run the annotation

## Set up the name of your vcf file
vcf_file='joint_psomagen.dementiagenes'
annovar_file=${vcf_file}.hg38_multianno
annovar_database='/home/acostauribe/bin/annovar/humandb/'

# Run Annovar
table_annovar.pl ${vcf_file}.vcf ${annovar_database} --buildver hg38 --outfile ${vcf_file} --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

# The —protocol' argument should be followed by the exact names of the annotation sources
# Follow the '—operation' argument by the annotation type: 'g' for gene-based annotation, 'r' for region-based annotation and 'f' for filter-based annotation. 

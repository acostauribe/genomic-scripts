#!/bin/bash
# This is a script to download, set up and run ANNOVAR to annotate vcf files
# Juliana Acosta-Uribe. Updated Dec 3 2025

### 1. Download annovar- fill registration and it will give you the link.
# https://annovar.openbioinformatics.org/en/latest/user-guide/download/ 
# https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
# https://www.nature.com/articles/nprot.2015.105
# https://statgen.us/files/tutorials/FunctionalAnnotation_Annovar_final.pdf

### 2. Uncompress ANNOVAR
tar -xvfz annovar.latest.tar.gz
# Make sure that the /annovar directory with the scripts  is in the $PATH
# provide the path to the $humandb directory in the annovar directory
humandb='$HOME/bin/annovar/humandb'

### 3. Download databases for annotation:

## A. Gene based annotation
## https://annovar.openbioinformatics.org/en/latest/user-guide/gene/
# FASTA sequences for all annotated transcripts in RefSeq Gene [Updated 20200817 at UCSC]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene $humandb
# FASTA sequences for all annotated transcripts in Gencode v41 Basic collection [Updated 20240513 at UCSC]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene $humandb

## B. Region based annotation
## https://annovar.openbioinformatics.org/en/latest/user-guide/region/
# Chromosome cytology bands
annotate_variation.pl -buildver hg38 -downdb cytoBand $humandb
# Transcription factor binding site
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar tfbsConsSites $humandb
# Identify variants disrupting microRNAs and snoRNAs
annotate_variation.pl -build hg38 -downdb wgRna $humandb

## C. Filter based annotation
# i. Population frequency:
# ExAC. 65000 exome allele allele frequency data for ALL, AFR (African), AMR (Admixed American), EAS (East Asian), FIN (Finnish), NFE (Non-finnish European), OTH (other), SAS (South Asian)). version 0.3. Left normalization done. [Updated 20151129]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 $humandb 
# gnomAD exome collection (v4), with "AF AF_popmax AF_male AF_female AF_raw AF_afr AF_sas AF_amr AF_eas AF_nfe AF_fin AF_asj AF_oth non_topmed_AF_popmax non_neuro_AF_popmax non_cancer_AF_popmax controls_AF_popmax" header [Updated 20231127]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad40_genome $humandb
# All of Us whole-genome data from first ~250k srWGS
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar allofus $humandb

# ii. Variant ID.
# dbSNP with allelic splitting [dbSNP151 Update 20240525]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp151 $humandb 
# Clinvar [Updated 20250721]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20250721 $humandb
# Protein domain variants for dbnsfp47a
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47a_interpro $humandb

# iii. Deleteriousness scores
# whole-exome SIFT, PolyPhen2 HDIV, PolyPhen2 HVAR, LRT, MutationTaster, MutationAssessor, FATHMM, MetaSVM, MetaLR, VEST, CADD, GERP++, PhyloP and SiPhy scores from dbNSFP version 4.7 [Updated 20240525]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp47a $humandb
# dbscSNV version 1.1 for splice site prediction by AdaBoost and Random Forest [Updated 20151218]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbscsnv11 $humandb
# Mendelian Clinically Applicable Pathogenicity (M-CAP) Score scores for non-synonymous variants. 
# http://bejerano.stanford.edu/mcap/index.html [updated 20161205]
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar mcap13 $humandb


## 4. Annotate your file

# The table_annovar.pl argument is a command used to annotate your input and generate a tab-delimited output that contains representative columns for each of the annotations.
# Here is an example of a script to run the annotation

## Set up the name of your vcf file
vcf_file='joint_psomagen.dementiagenes'
annovar_file=${vcf_file}.hg38_multianno
annovar_database='/home/acostauribe/bin/annovar/$humandb'

# Run Annovar
table_annovar.pl ${vcf_file}.vcf ${annovar_database} --buildver hg38 --outfile ${vcf_file} --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

# The —protocol' argument should be followed by the exact names of the annotation sources
# Follow the '—operation' argument by the annotation type: 'g' for gene-based annotation, 'r' for region-based annotation and 'f' for filter-based annotation. 

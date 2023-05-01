#!/bin/bash
#Code for fast1 to vcf processing

#for file.fastq.gz
fastq='1' 

# STEP 1. Running FastQC as part of a pipeline
# download from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

#mkdir QC
#mkdir input 
mkdir output

fastqc ${fastq}_1.fastq.gz ${fastq}_2.fastq.gz -o output

# STEP 2. EDICION DE SECUENCIAS DE MALA CALIDAD
#download and install 

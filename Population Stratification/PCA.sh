#!/bin/bash 
#Script to perform a PCA  with a 1000gp subset
#Juliana Acosta Uribe 2022

#usage chmod u+x PCA.sh 
#./PCA.sh

file='1000G.toy'

#1000G.toy is a subset of the 1000Genomes project dataset
#remember that all individuals need to be unrelated. if needed do:
#king -b $file.bed --unrelated
#that will create a list of unrelated individuals that you can extract from $file using plink --keep

#1. Retain variants with MAF > 10%
#plink --bfile $file --maf 0.1 --make-bed --out $file.maf

#2. Calculate LD
#--indep-pairwise <window size>['kb'] <step size (variant ct)> <r^2 threshold>
#plink --bfile $file.maf --indep-pairwise 50 10 0.2 

#3. Retain variants not in LD
#plink --file $file.maf --extract plink.prune.in --make-bed --out $file.maf.ld

#4. Perform a PCA
#a dataset for a PCA should be > 100.000 variants
plink --file $file.maf.ld --pca tabs --out $file.maf.ld.pca

#Files to plot:
#scatter plot $file.maf.ld.pca.eigenvec
#screeplot $file.maf.ld.pca.eigenval

#5. Run ADMIXTURE
for k in {1..10}; do 
./admixture $file.maf.ld --cv ${k} | tee log${k}.out; done

#Files to plot:
# Stacked barplot $file.maf.ld.K.Q plot with Pong
# take cverror from logK.out and plot
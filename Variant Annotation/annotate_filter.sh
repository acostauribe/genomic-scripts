#!/bin/bash
## Script to extract and annotate a vcf 
## Juliana Acosta-Uribe April 2024

# This script is set up in three phases. You can turn "on" or "off" the different phases below.
# After setting up what you need, run this script as follows: 
## chmod u+x annotate_variants.sh
## ./annotate_filter.sh | tee aannotate_filter.log \

## PART ONE: ANNOTATE VCF
### Extracts areas and samples of interest, checking alignment and normalization, annotate vcf with ANNOVAR
### outputs are a annotated vcf "*.hg38_multianno.vcf.gz" and a tab delimited table with the variant information and genotypes "*.hg38_multianno.txt"
### Set to "TRUE" is you want to run this step
annotate_vcf=TRUE

## PART TWO: COUNT CASES AND CONTROLS CARRYING EACH VARIANT
### This step uses plink to add a count of the cases and controls of a given variant.
### Appends columns to the "*.hg38_multianno.txt" file
### Set to "TRUE" is you want to run this step
count_carriers=TRUE

## PART THREE: FILTER ANNOTATED VARIANTS
### Filter annotated variants in "*.hg38_multianno.txt" based on prefered settings. 
### Set to "TRUE" is you want to run this step
filter_annotated=TRUE


## In the following lines you need to specify the sub-analyses you want, and provide the name and path for needed files

# PART ONE: ANNOTATE VCF

# Provide a vcf:
vcf_file='joint_redlat_1-22-24.ndg-panel.redlat.id' 
## Must be gzipped and end in 'vcf.gz'

## Chromosomes must be encoded as chr#. 
## If you need to change chromosome designation from # to chr# you need to create a file (e.g. chromosomes.txt) where every line has two columns: 'old name' 'new name' (e.g. 1 chr1)
## bcftools annotate --rename-chrs chromosomes.txt --output-type z file.vcf.gz > file.chr.vcf.gz

# Filter your vcf:
## VCF files can get very big, you can extract and annotate a set of regions of interest.
## If you wantto extract a subset of regions, set extract_panel to "TRUE" and provide a file as described
extract_panel=FALSE
## Provide a file with your regions of interest (i.e. genes) 
## File must be in bed format: 'chr', 'start bp', 'end bp'. Chromosome format must match vcf style (eg. "chr1" vs. "1")
panel='' 

## If you want to extract a subset of samples from the vcf, set extract_samples to "TRUE" and provide a file as described
extract_samples=FALSE
## Provide a file with your samples of interest
## File must have one column with the IDs of samples that you want to retain. Sample ID must match an ID in the vcf
samples='' 

# ANNOVAR expects vcf aligned to hg38 and 'table_annovar.pl' must be in your path. 
## Provide a path to your ANNOVAR database
#annovar_database_PATH='/path/to/annovar/humandb/'
annovar_database_PATH='/home/acostauribe/bin/annovar/humandb/'

## Before annotation, the vcf will be checked for proper alignmnent and normalization
## Provide a path to the hg38 reference genome (fa.gz) and the respective index 
#fasta_file='/path/to/hg38.fa.gz'
fasta_file='/home/acostauribe/public_html/Utilities/hg38.fa.gz'
## fasta file was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
## The index file fasta.fai was created using http://www.htslib.org/doc/samtools-faidx.html

## This step requires the following software: ANNOVAR, Vcftools, BCFtools


# PART TWO: COUNT CASES AND CONTROLS WITH EACH VARIANT

## Provide a file with the case/control designation for your samples.
## File must be two columns. Fisrt column matching the sample ID on the VCF, second column the phenotype of the sample.
## Case/control phenotypes are expected to be encoded as 1=unaffected (control), 2=affected (case)
## -9 or 0 are accepted as missing value
phenotypes='phenotypes.txt'
## Variant counts are given in Homozygous alternate | Heterozygous | Homozygous reference format

## You can provide a file with a set of variants of interest to use for the case-control count (interest_variants)
interest_variants=''
## interest_variants must be a single column file where each variant is represented by their identifier in the vcf (e.g. rs number)
## Part one, Step 4 assigns an ID '%CHROM\_%POS\_%REF\_%FIRST_ALT'  to all variants that didnt have an rs number in the vcf
## If you do not provide this file, then the computer will perform the case control count in all variants. 

## you can also provide a 
cohorts='groups' # file where each line is the name of a file that contains the individuals from a cohort/group of interest. 

## This step requires to have Plink 1.9, Python3 and pandas installed
## You can install pandas as follows
#pip install pandas


# PART THREE: FILTER ANNOTATIONS

## If you want to filter your variants after annotating them, you can set up the following filterings:

## You can limit these searches to a subset of genes you are interested in.
extract_genes='TRUE'
## Provide a text file named gene_list with one gene per line.
gene_list="genes_interest.txt"
## In addition to the filtering set before, this will generate a subset with all the annotations in the gene_list

## Extract variants that have been reported as "Pathogenic" or "Likely pathogenic" in ClinVar.
## Generates annotated_file.ClinVar-Pathogenic.txt
extract_pathogenic='TRUE'
## You can extract reported pathogenic variants in all the data, or in a subset of genes (see line 82)

## Get Exonic and Splicing variants. 
## Generates annotated_file.exonic.txt
extract_coding='TRUE'

## Exclude synonymous mutations. 
## Generates annotated_file.non-syn.txt
remove_synonymous='TRUE'

## Extract variants with low allelic frequency and high CADD score for manual curation
## Generates annotated_file.manual-curation.txt 
manual_curation='TRUE'
## Define the CADD and MAF thresholds
MAF_filter='0.01'
CADD_filter='20'
## Variants with missing MAF or CADD information will be retained for manual curation"
## If the value is a '.' or a blank, it will be printed out as well.
## Ensure MAF_filter and CADD_filter are properly defined as numeric values

## The order of these filterings is extract_genes > extract_pathogenic ; extract_genes > extract_coding > remove_synonymous > manual_curation
## You can remove a step, but the order will be kept



#================SCRIPT

## PART ONE: ANNOTATE VCF

echo "Starting script on $(date)"

if [[ $annotate_vcf == "TRUE" ]]; then

    ## 0. Search for vcf
    if [ -f "${vcf_file}.vcf" && -s "${vcf_file}.vcf"]; then
        echo "gzipping ${vcf_file}.vcf"
        bgzip "${vcf_file}.vcf"
        echo "Starting analysis of ${vcf_file}.vcf.gz"
    elif [ -f "${vcf_file}.vcf.gz" && -s "${vcf_file}.vcf.gz"]; then
        echo "Starting analysis of ${vcf_file}.vcf.gz"
    else
        echo "Please provide a valid vcf file"
        return 1
    fi

    ## 1. Extract target regions and samples from vcf file
    if [[ $extract_panel == 'TRUE' || $extract_samples == 'TRUE' ]]; then
        if [[ $extract_panel == 'TRUE' && $extract_samples == 'TRUE' ]]; then
            echo "Extracting regions and samples of interest"
            vcftools --gzvcf "${vcf_file}.vcf.gz" --bed "$panel" --keep "$samples" --recode --recode-INFO-all --out "${vcf_file}_samples_panel"
            mv "${vcf_file}_samples_panel.recode.vcf" "${vcf_file}_samples_panel.vcf" 
            vcf_file="${vcf_file}_samples_panel" 

        elif [[ $extract_panel == 'TRUE' ]]; then
            echo "Extracting regions of interest"
            vcftools --gzvcf "${vcf_file}.vcf.gz" --bed "$panel" --recode --recode-INFO-all --out "${vcf_file}_panel"
            mv "${vcf_file}_panel.recode.vcf" "${vcf_file}_panel.vcf" 
            vcf_file="${vcf_file}_panel"

        elif [[ $extract_samples == 'TRUE' ]]; then
            echo "Extracting samples of interest"
            vcftools --gzvcf "${vcf_file}.vcf.gz" --keep "$samples" --recode --recode-INFO-all --out "${vcf_file}_samples"
            mv "${vcf_file}_samples.recode.vcf" "${vcf_file}_samples.vcf" 
            vcf_file="${vcf_file}_samples"
        fi
        bgzip "${vcf_file}.vcf"
        # Compressing the file at the end, instead of using '--stdout | bgzip' in vcftools, allows generating a log file useful for debugging.
    fi

    ## 2. Removing monomorphic variants
    echo "Removing monomorphic variants"
    bcftools view --min-ac 1 --output-type z "${vcf_file}.vcf.gz" > "${vcf_file}.tmp.vcf.gz"
    ## gzip and indexing
    tabix -p vcf "${vcf_file}.tmp.vcf.gz"

    ## 3. Check alignment and INDEL normalization
    echo "Checking alignment and INDEL normalization"
    bcftools norm --check-ref ws --fasta-ref "$fasta_file" --output-type z "${vcf_file}.tmp.vcf.gz" > "${vcf_file}_ref.vcf.gz"
    # --check-ref warn (w), exclude (x), or set/fix (s)
    # --output-type compressed VCF (z) will bgzip the output
    vcf_file="${vcf_file}_ref"
    rm "${vcf_file}.tmp."*

    ## 4. Assign variant identifiers.
    echo "Assigning variant identifiers"
    ## Making sure every variant has an rs ID. This is very useful for vcftools --> plink variant extraction
    ## Missing IDs will be set to chromosome_position_reference_first-alternate
    bcftools annotate --set-id '+%CHROM\_%POS\_%REF\_%FIRST_ALT' --output-type z "${vcf_file}.vcf.gz" > "${vcf_file}_rs.vcf.gz"
    vcf_file="${vcf_file}_rs"

    ## 5. Using ANNOVAR to annotate the vcf 
    echo "Using ANNOVAR to annotate ${vcf_file}"
    table_annovar.pl "${vcf_file}.vcf.gz" "${annovar_database_PATH}" --buildver hg38 --outfile "${vcf_file}" --protocol refGene,ensGene,exac03,gnomad30_genome,ALL.sites.2015_08,AFR.sites.2015_08,EUR.sites.2015_08,EAS.sites.2015_08,SAS.sites.2015_08,AMR.sites.2015_08,avsnp150,clinvar_20220320,dbnsfp33a,dbscsnv11,revel --operation g,g,f,f,f,f,f,f,f,f,f,f,f,f,f --nastring . --vcfinput --remove 

    ## This will produce ${vcf_file}.hg38_multianno.vcf and ${vcf_file}.hg38_multianno.txt
    annovar_file="${vcf_file}.hg38_multianno"

    ## 6. Edit the ANNOVAR results header 
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
fi


# PART TWO: COUNT CASES AND CONTROLS CARRYING EACH VARIANT

if [[ $count_carriers == "TRUE" ]]; then

    ## 1. Checking for ANNOVAR output files

    if [[ $annotate_vcf != "TRUE" ]]; then
        echo "Checking for ANNOVAR output files"
        annovar_file=($(find . -maxdepth 1 -type f -name "*.hg38_multianno.txt" ))
        # Check the number of files found
        if [[ ${#annovar_file[@]} -eq 1 ]]; then
            # There is one file found, assign to annovar_file
            annovar_file="${annovar_file[0]}"
            if [[ -s "$annovar_file" ]]; then
                echo "ANNOVAR file is properly formatted"
                annovar_file="${annovar_file%.txt}"
            else
                echo "$annovar_file is not readable. Stopping script"
                return 1  # return with a status code of 1 to indicate an error.
            fi
        elif [[ ${#annovar_file[@]} -gt 1 ]]; then
            # More than one file found, raise a warning
            echo "Warning: More than one ANNOVAR file found. Keep ony one file ending in '.hg38_multianno.txt' in the working directory"
            echo "Stopping script"
            return 1  # return with a status code of 1 to indicate an error.
        else
            # If no file was found, raise a warning
            echo "Provide a valid ANNOVAR file for variant filtering. Script requires one file ending in '.hg38_multianno.txt' in the working directory"
            echo "Stopping script"
            return 1  # return with a status code of 1 to indicate an error.
        fi
        annovar_vcf=($(find . -maxdepth 1 -type f -name "*.hg38_multianno.vcf.gz" ))
        # Check the number of files found
        if [[ ${#annovar_vcf[@]} -eq 1 ]]; then
            # There is one file found, assign to annovar_vcf
            annovar_vcf="${annovar_vcf[0]}"
            if [[ -s "$annovar_vcf" ]]; then
                echo "ANNOVAR annotated vcf found $annovar_vcf"
                annovar_vcf="${annovar_vcf%.vcf.gz}"
            else
                echo "$annovar_vcf is not readable. Stopping script"
                return 1  # return with a status code of 1 to indicate an error.
            fi
        elif [[ ${#annovar_vcf[@]} -gt 1 ]]; then
            # More than one file found, raise a warning
            echo "Warning: More than one ANNOVAR annotated vcf found. Keep ony one file ending in '.hg38_multianno.vcf.gz' in the working directory"
            echo "Stopping script"
            return 1  # return with a status code of 1 to indicate an error.
        else
            # If no file was found, raise a warning
            echo "Provide a valid ANNOVAR annotated vcf for variant counting. Script requires one file ending in '.hg38_multianno.vcf.gz' in the working directory"
            echo "Stopping script"
            return 1  # return with a status code of 1 to indicate an error.
        fi
    fi
    
    ## 2. Import data into plink 
    echo "Importing vcf into plink"
    plink_file="${annovar_vcf}.plink"

    plink --vcf "${annovar_vcf}.vcf.gz" \
        --vcf-half-call m \
        --allow-extra-chr \
        --double-id \
        --keep-allele-order \
        --make-bed \
        --out "$plink_file"

    ## 3. Generate a plink phenotype file
    if [[ -f "$phenotypes" && -s "$phenotypes" ]]; then
        awk '{print $1, $1, $2}' "$phenotypes" > "${phenotypes}_plink"
    else 
        echo "Please provide the phenotypes of the samples in your cohort"
    fi

    ## 4. Count SNPs in all the cohort. Output ends in ".model"
    echo "Using plink to count variants in cases and controls"
    if [[ -f "$interest_variants " && -s "$interest_variants" ]]; then
        plink --bfile "${plink_file}" \
        --extract "$interest_variants" \
        --pheno "${phenotypes}_plink" \
        --model \
        --allow-no-sex \
        --out "${plink_file}.all"
    else
        plink --bfile "${plink_file}" \
        --pheno "${phenotypes}_plink" \
        --model \
        --allow-no-sex \
        --out "${plink_file}.all"
    fi
    ## Retain 'GENO' counts
    awk -v OFS='\t' 'NR==1; NR > 1{ if ($5 == "GENO") { print } }' "${plink_file}.all.model" > "${plink_file}.all.geno"
    sed -i 's/UNAFF/All_unaffected/g' "${plink_file}.all.geno"
    sed -i 's/AFF/All_affected/g' "${plink_file}.all.geno"
    sed -i "s/SNP/ID/g" "${plink_file}.all.geno"
    ## Retain columns 2, 6 and 7
    awk -v OFS='\t' '{print $2, $6, $7}' "${plink_file}.all.geno" > "${plink_file}.cohort-comparison.txt"

    ## 5. Count SNPs in sub-cohorts
    if [[ -f "$cohorts" && -s "$cohorts" ]]; then
        echo "Using Plink to count variant carriers in cohorts"

        while read -r line; do
            awk '{ print $1, $1}' "$line" > "${line}_plink"
        done < "$cohorts"

        if [[ -f "$interest_variants" && -s "$interest_variants" ]]; then
            while read -r line; do 
                plink --bfile "$plink_file" \
                --extract "$interest_variants" \
                --keep "${line}_plink" \
                --pheno "${phenotypes}_plink" \
                --model \
                --allow-no-sex \
                --out "${plink_file}.$line" 
                
                cat <(awk -v OFS='\t' 'BEGIN { print "AFF", "UNAFF"}') <(awk -v OFS='\t' '{ if ($5 == "GENO") {print $6,$7} }' "${plink_file}.$line.model") > "${plink_file}.$line.geno"
                sed -i "s/UNAFF/$line.unaffected/g" "${plink_file}.$line.geno" 
                sed -i "s/AFF/$line.affected/g" "${plink_file}.$line.geno"
                #keep only the genotipic count values
                paste "${plink_file}.cohort-comparison.txt" "${plink_file}.$line.geno"  > "${plink_file}.cohort-comparison.temp"
                #add it to the counts from the entire cohort
                mv "${plink_file}.cohort-comparison.temp" "${plink_file}.cohort-comparison.txt"
            done < "$cohorts"
        else
            while read -r line; do 
                plink --bfile "$plink_file" \
                --keep "${line}_plink" \
                --pheno "${phenotypes}_plink" \
                --model \
                --allow-no-sex \
                --out "${plink_file}.$line" 
                
                cat <(awk -v OFS='\t' 'BEGIN { print "AFF", "UNAFF"}') <(awk -v OFS='\t' '{ if ($5 == "GENO") {print $6,$7} }' "${plink_file}.$line.model") > "${plink_file}.$line.geno"
                sed -i "s/UNAFF/$line.unaffected/g" "${plink_file}.$line.geno" 
                sed -i "s/AFF/$line.affected/g" "${plink_file}.$line.geno"
                #keep only the genotipic count values
                paste "${plink_file}.cohort-comparison.txt" "${plink_file}.$line.geno"  > "${plink_file}.cohort-comparison.temp"
                #add it to the counts from the entire cohort
                mv "${plink_file}.cohort-comparison.temp" "${plink_file}.cohort-comparison.txt"
            done < "$cohorts"
        fi
    fi
    ## 6. check output
    if [[ -s "${plink_file}.cohort-comparison.txt" ]]; then
        echo "Variant counts written to ${plink_file}.cohort-comparison.txt"
    else
        echo "${plink_file}.cohort-comparison.txt file was not properly generated. Stopping script"
        return 1  # return with a status code of 1 to indicate an error.
    fi

    ## 7. Use python to merge variant annotation () and variant count
    ## Notice that Pyhton sub-script cannot be indented. Leave as is
    ## Execute Python code using heredoc

python3 - <<EOF
import sys
import pandas as pd

# Read TSV files into pandas DataFrames, assuming they have headers
file_a = pd.read_table("${annovar_file}.txt")
file_b = pd.read_table("${plink_file}.cohort-comparison.txt")


# Merge the dataframes on 'ID'
# Right merge allows to keep all variants in vcf even if they are not imported into plink (This is common for multialleic vars).
merged_data = pd.merge(file_a, file_b, on='ID', how='left')

# Save the merged DataFrame to a TSV file with tab as the separator
merged_data.to_csv('annotated_variants_cohort_counts.tsv', sep='\t', index=False)
EOF

   if [[ -s "annotated_variants_cohort_counts.tsv" ]]; then
        mv "annotated_variants_cohort_counts.tsv" "${annovar_file}.counts.txt"
        annovar_file="${annovar_file}.counts"
        echo "Merging variant annotations and carrier counts is complete. Output saved to ${annovar_file}.counts"
    else
        echo "${annovar_file}.counts file was not properly generated. Stopping script"
        return 1  # return with a status code of 1 to indicate an error.
    fi
fi

## PART THREE: FILTER ANNOTATED VARIANTS

if [[ $filter_annotated == 'TRUE' ]]; then

    if [[ $annotate_vcf != "TRUE" && $count_carriers != "TRUE" ]]; then
        annovar_file=($(find . -maxdepth 1 -type f -name "*.hg38_multianno.txt" ))
        # Check the number of files found
        if [[ ${#annovar_file[@]} -eq 1 ]]; then
            # There is one file found, assign to annovar_file
            annovar_file="${annovar_file[0]}"
            if [[ -s "$annovar_file" ]]; then
                echo "ANNOVAR file is properly formatted"
                annovar_file="${annovar_file%.txt}"
            else
                echo "$annovar_file is not readable. Stopping script"
                return 1  # return with a status code of 1 to indicate an error.
            fi
        elif [[ ${#annovar_file[@]} -gt 1 ]]; then
            # More than one file found, raise a warning
            echo "Warning: More than one ANNOVAR file found. Keep ony one file ending in '.hg38_multianno.txt' in the working directory"
            echo "Stopping script"
            return 1  # return with a status code of 1 to indicate an error.
        else
            # If no file was found, raise a warning
            echo "Provide a valid ANNOVAR file for variant filtering. Script requires one file ending in '.hg38_multianno.txt' in the working directory"
            echo "Stopping script"
            return 1  # return with a status code of 1 to indicate an error.
        fi
    fi
    
    echo "Filtering annotated variants in ${$annovar_file}"

    ## 1. Count number of variants
    echo "Number of variants in ${annovar_file}.txt:"
    wc -l "${annovar_file}.txt" | awk '{print $1-1}' # Subtract the header line to get the number of variants

    ## 2. Extract variants in a subset of genes 
    if [[ $extract_genes == 'TRUE' ]]; then
        echo "Extracting variants in ${gene_list}"
        awk -v OFS='\t' 'NR==FNR {genes[$1]; next} NR > 1 {if ($7 in genes) print}' "$gene_list" "${annovar_file}.txt" > "${annovar_file}.gene_list.txt"
        annovar_file="${annovar_file}.gene_list"
        echo "Number of variants in ${annovar_file}.txt:"
        wc -l "${annovar_file}.txt" | awk '{print $1-1}'
    fi

    ## 3. Extract variants that have been reported as "Pathogenic" or "Likely pathogenic" in ClinVar
    if [[ $extract_pathogenic == 'TRUE' ]]; then
        echo "Extracting variants that have been reported as pathogenic in ClinVar"
        echo "Warning: Pathogenic variants are reported as of March 20 2022"
        awk -v OFS='\t' 'NR > 1 && /athogenic/ && !/Conflicting/' "${annovar_file}.txt" > "${annovar_file}.ClinVar-Pathogenic.txt"
        echo "Number of variants in ${annovar_file}.ClinVar-Pathogenic.txt:"
        wc -l "${annovar_file}.ClinVar-Pathogenic.txt" | awk '{print $1-1}'
    fi

    ## 4. Get Exonic and Splicing variants.
    if [[ $extract_coding == 'TRUE' ]]; then
        echo "Retaining Exonic and Splicing variants"
        awk -v OFS='\t' 'NR==1; NR > 1{ if(($6 == "exonic") || ($6 == "exonic;splicing") || ($6 == "splicing")) { print } }' "${annovar_file}.txt" > "${annovar_file}.exonic.txt"
        annovar_file="${annovar_file}.exonic"
        echo "Number of variants in ${annovar_file}.txt:"
        wc -l "${annovar_file}.txt" | awk '{print $1-1}'
    fi

    ## 5. Exclude synonymous mutations
    if [[ $remove_synonymous == 'TRUE' ]]; then
        echo "Removing synonymous variants"
        awk -v OFS='\t' 'NR==1; NR > 1{ if(($9 != "synonymous SNV") && ($6 != "unknown")) { print } }' "${annovar_file}.txt" > "${annovar_file}.non-syn.txt"
        annovar_file="${annovar_file}.non-syn"
        echo "Number of variants in ${annovar_file}.txt:"
        wc -l "${annovar_file}.txt" | awk '{print $1-1}'
    fi

    ## 6. Extract variants with low allelic frequency and high CADD score
    if [[ $manual_curation == 'TRUE' ]]; then
        echo "Extracting coding variants for manual curation"
        echo "Warning! Variants with missing MAF or CADD information will be retained for manual curation"
        # Ensure MAF_filter and CADD_filter are properly defined as numeric values
        awk -v OFS='\t' -v MAF="$MAF_filter" -v CADD="$CADD_filter" 'NR==1; NR > 1{ if((($16 == "." || $16 <= MAF) && ($25 == "." || $25 <= MAF)) && ($74 == "." || $74 >= CADD)) {print} }' \
        "${annovar_file}.txt" > "${annovar_file}.manual-curation.txt"
        echo "Number of variants in ${annovar_file}.manual-curation.txt:"
        wc -l "${annovar_file}.manual-curation.txt" | awk '{print $1-1}'
    fi
fi

## Tidy up!
rm 
*_plink

## Store annotations in their own directory
echo "Creating a directory for annotations"
mkdir ANNOVAR
mv ${annovar_file}.* ./ANNOVAR
#mv filter_variants.* ./ANNOVAR
mv ANNOVAR ANNOVAR_$(date "+%Y%m%d")

echo "Finished" $(date)

#!/bin/bash
set -euo pipefail

## Script to extract, normalize, and annotate a VCF
## Juliana Acosta-Uribe, May 2026
##
## Usage:
##   chmod u+x annotate_variants.sh
##   ./annotate_variants.sh | tee annotate_variants.log

# ==============================================================================
# CONFIGURATION
# ==============================================================================

vcf_file='LATAM5k_joint_call_11-14-25_dp10_gq20.ndeg'
bcftools='bcftools'
ANNOVAR=/home/acostauribe/bin/annovar/table_annovar.pl
annovar_database_PATH='/home/acostauribe/bin/annovar/humandb/'

# hg38 FASTA from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/
# Index (.fai) must be in the same directory (created with: samtools faidx)
fasta_file="$HOME/Utilities/hg38.chrm-string.fa"

# Set to TRUE to strip genotype columns from the final ANNOVAR .txt output
remove_genotypes="FALSE"

# ==============================================================================
# PHASE 1: RENAME CHROMOSOMES (# -> chr#)
# ==============================================================================

echo ">>> [$(date)] Starting script"
echo ">>> Renaming chromosomes..."

${bcftools} annotate \
    --rename-chrs ~/chromosomes.txt \
    --output-type z \
    -o "${vcf_file}.chr.vcf.gz" \
    "${vcf_file}.vcf.gz"

vcf_file="${vcf_file}.chr"

# ==============================================================================
# PHASE 2: INDEX, NORMALIZE, AND RE-INDEX
# ==============================================================================

echo ">>> Indexing..."
${bcftools} index -t -f "${vcf_file}.vcf.gz"

echo ">>> Normalizing and aligning to reference..."
# --check-ref w: warn only on mismatches (first-pass checkpoint)
# NOTE: --check-ref s does NOT fix strand issues — do not use it for that purpose
${bcftools} norm \
    --check-ref w \
    -f "${fasta_file}" \
    -Oz \
    -o "${vcf_file}.ref.vcf.gz" \
    "${vcf_file}.vcf.gz" \
    2> ref_mismatches.txt

mismatch_count=$(grep -c "REF_MISMATCH" ref_mismatches.txt || true)
echo ">>> REF mismatches found: ${mismatch_count}"

if [ "${mismatch_count}" -gt 0 ]; then
    echo ">>> Mismatches detected. Checking if source is PLINK (PR flag)..."
    pr_flag=$(${bcftools} view -h "${vcf_file}.vcf.gz" | grep -c 'ID=PR' || true)

    if [ "${pr_flag}" -gt 0 ]; then
        echo ">>> PLINK-derived VCF detected (provisional REF alleles). Fixing REF with --check-ref s..."
        ${bcftools} norm \
            --check-ref s \
            -f "${fasta_file}" \
            -Oz \
            -o "${vcf_file}.ref.vcf.gz" \
            "${vcf_file}.vcf.gz" \
            2> ref_fix.txt

        fixed=$(grep -c "Fixing" ref_fix.txt || true)
        remaining=$(grep -c "REF_MISMATCH" ref_fix.txt || true)
        echo ">>> REF alleles fixed: ${fixed}"
        echo ">>> Remaining mismatches after fix: ${remaining}"

        if [ "${remaining}" -gt 0 ]; then
            echo "WARNING: ${remaining} mismatches could not be fixed. Review ref_fix.txt before proceeding."
        fi
    else
        echo "WARNING: REF mismatches found but source is not PLINK. Review ref_mismatches.txt before proceeding."
        echo "Exiting. Re-run after inspecting mismatches."
        exit 1
    fi
fi

vcf_file="${vcf_file}.ref"

echo ">>> Indexing normalized VCF..."
${bcftools} index -t -f "${vcf_file}.vcf.gz"

# ==============================================================================
# PHASE 3: FILL AC/AN TAGS
# ==============================================================================

echo ">>> Filling AC and AN tags..."
${bcftools} +fill-tags \
    "${vcf_file}.vcf.gz" \
    --output-type z \
    --output "${vcf_file}.AC.vcf.gz" \
    -- -t AC,AN

vcf_file="${vcf_file}.AC"
echo ">>> Filled tags. Current file: ${vcf_file}.vcf.gz"

# ==============================================================================
# PHASE 4: ANNOTATE WITH ANNOVAR
# ==============================================================================

echo ">>> Annotating with ANNOVAR: ${vcf_file}.vcf.gz"
${ANNOVAR} "${vcf_file}.vcf.gz" "${annovar_database_PATH}" \
    --buildver hg38 \
    --outfile "${vcf_file}" \
    --protocol refGene,ensGene,gnomad40_genome,allofus,avsnp151,clinvar_20250721,dbnsfp47a,dbnsfp47a_interpro,dbscsnv11,revel \
    --operation  g,g,f,f,f,f,f,f,f,f \
    --nastring . \
    --vcfinput \
    --remove

# Outputs: ${vcf_file}.hg38_multianno.vcf and ${vcf_file}.hg38_multianno.txt
annovar_file="${vcf_file}.hg38_multianno"

# ==============================================================================
# PHASE 5: REFORMAT ANNOVAR HEADER
# ==============================================================================

echo ">>> Reformatting ANNOVAR header..."

# Compress VCF output
bgzip -f "${annovar_file}.vcf"

# Extract the VCF #CHROM line (contains sample IDs)
zgrep -m1 '^#CHROM' "${annovar_file}.vcf.gz" > "${annovar_file}.chrom_line"

# Normalize ANNOVAR header to strict tab-delimited
head -1 "${annovar_file}.txt" \
    | awk 'BEGIN{FS="[ \t]+"; OFS="\t"} {$1=$1}1' \
    > "${annovar_file}.header"

# ANNOVAR txt output includes mystery columns before the VCF resumes at Otherinfo4.
# Build a clean header: ANNOVAR annotations (cols 1..Otherinfo4-1) + VCF #CHROM line
col=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Otherinfo4"){print i; exit}}' "${annovar_file}.header")

paste \
    <(cut -f "1-$((col-1))" "${annovar_file}.header") \
    "${annovar_file}.chrom_line" \
    > "${annovar_file}.header_new"

# Ensure strictly tab-delimited
awk 'BEGIN{FS="[ \t]+"; OFS="\t"} {$1=$1}1' \
    "${annovar_file}.header_new" \
    > "${annovar_file}.header_new.tmp"
mv "${annovar_file}.header_new.tmp" "${annovar_file}.header_new"

# Cleanup intermediates
rm "${annovar_file}.chrom_line" "${annovar_file}.header"

# Rename columns to human-readable labels
sed -i \
    -e 's/Func\.refGene/RefSeq_region/g' \
    -e 's/Gene\.refGene/RefSeq_gene/g' \
    -e 's/GeneDetail\.refGene/RefSeq_gene_detail/g' \
    -e 's/ExonicFunc\.refGene/RefSeq_exonic_effect/g' \
    -e 's/AAChange\.refGene/RefSeq_AA_change/g' \
    -e 's/Func\.ensGene/Ensembl_region/g' \
    -e 's/Gene\.ensGene/Ensembl_gene/g' \
    -e 's/GeneDetail\.ensGene/Ensembl_gene_detail/g' \
    -e 's/ExonicFunc\.ensGene/Ensembl_exonic_effect/g' \
    -e 's/AAChange\.ensGene/Ensembl_AA_change/g' \
    -e 's/gnomad40_genome_AF\b/gnomAD.v4_AF/g' \
    -e 's/gnomad40_genome_AF_raw/gnomAD.v4_AF_raw/g' \
    -e 's/gnomad40_genome_AF_XX/gnomAD.v4_AF_XX/g' \
    -e 's/gnomad40_genome_AF_XY/gnomAD.v4_AF_XY/g' \
    -e 's/gnomad40_genome_AF_grpmax/gnomAD.v4_max_AF/g' \
    -e 's/gnomad40_genome_faf95\b/gnomAD.v4_FAF95/g' \
    -e 's/gnomad40_genome_faf99\b/gnomAD.v4_FAF99/g' \
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
    -e 's/CLNDN\b/ClinVar_disease_name/g' \
    -e 's/CLNDISDB/ClinVar_disease_IDs/g' \
    -e 's/CLNREVSTAT/ClinVar_review_status/g' \
    -e 's/CLNSIG/ClinVar_clinical_significance/g' \
    -e 's/ONCDN\b/Oncogenicity_disease_name/g' \
    -e 's/ONCDISDB/Oncogenicity_disease_IDs/g' \
    -e 's/ONCREVSTAT/Oncogenicity_review_status/g' \
    -e 's/ONC\b/Oncogenicity/g' \
    -e 's/SCIDN\b/Somatic_impact_disease_name/g' \
    -e 's/SCIDISDB/Somatic_impact_disease_IDs/g' \
    -e 's/SCIREVSTAT/Somatic_impact_review_status/g' \
    -e 's/SCI\b/Somatic_clinical_impact/g' \
    -e 's/Aloft_pred/ALoFT_pred/g' \
    "${annovar_file}.header_new"

# Replace original header with the cleaned, renamed one
sed '1d' "${annovar_file}.txt" > "${annovar_file}.header_removed"
cat "${annovar_file}.header_new" "${annovar_file}.header_removed" > "${annovar_file}.txt"
rm "${annovar_file}.header_new" "${annovar_file}.header_removed"

# Validate output
if [[ ! -f "${annovar_file}.txt" || ! -s "${annovar_file}.txt" ]]; then
    echo "ERROR: ANNOVAR annotation failed — output file missing or empty. Stopping."
    exit 1
fi
echo ">>> ANNOVAR annotation complete: ${annovar_file}.txt"

# ==============================================================================
# PHASE 6: OPTIONALLY REMOVE GENOTYPE COLUMNS
# ==============================================================================

if [[ "${remove_genotypes}" == "TRUE" ]]; then
    echo ">>> Removing genotype columns from ANNOVAR txt output..."

    Otherinfo1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Otherinfo1"){print i; exit}}' "${annovar_file}.txt")

    if (( Otherinfo1 > 1 )); then
        base="${annovar_file%.hg38_multianno}"
        cut -f "1-$((Otherinfo1-1))" "${base}.hg38_multianno.txt" \
            > "${base}.no-geno.hg38_multianno.txt"
        mv "${base}.hg38_multianno.txt" "${base}.hg38_multianno.original.txt"
        annovar_file="${base}.no-geno.hg38_multianno"
        echo ">>> Genotypes removed. Output: ${annovar_file}.txt"
    else
        echo ">>> No genotype columns found in ANNOVAR output — skipping."
    fi
fi

# ==============================================================================
echo ">>> [$(date)] Done. Final file: ${annovar_file}.txt"

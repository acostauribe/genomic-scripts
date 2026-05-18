## Takes the output from annotate_variants.sh and extracts a set of genes/samples

import subprocess
import sys

# ── Dependency check ──────────────────────────────────────────────────────────
try:
    import pandas as pd
except ImportError:
    print("pandas not found, installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
    import pandas as pd

# ── User-defined inputs ───────────────────────────────────────────────────────
annotated_vcf = "LATAM5k_joint_call_11-14-25_dp10_gq20.ndeg.chr.ref.hg38_multianno.txt"
output_file   = "LATAM5k_joint_call_11-14-25_dp10_gq20.ndeg.chr.ref.hg38_multianno.PSEN1.txt"   # name of the output file

# genome_list: sample IDs to retain (columns after INFO).
# Options:
#   - Path to a .txt file (one value per line): genome_list = "genome_list.txt"
#   - Inline list:                              genome_list = ["sample_A", "sample_B"]
#   - No filter (keep all):                     genome_list = None
genome_list = None

# gene_list: genes to retain (rows where RefSeq_gene or Ensembl_gene match).
# Options:
#   - Path to a .txt file (one value per line): gene_list = "gene_list.txt"
#   - Inline list:                              gene_list = ["BRCA1", "BRCA2"]
#   - No filter (keep all):                     gene_list = None
gene_list = ["PSEN1"]

# region_list: functional regions to retain (rows where RefSeq_region or Ensembl_region match).
# Options:
#   - Path to a .txt file (one value per line): region_list = "region_list.txt"
#   - Inline list:                              region_list = ["UTR5", "intronic", "exonic", "UTR3", "splicing"]
#   - No filter (keep all):                     region_list = None
region_list = ["UTR5", "exonic", "UTR3", "splicing"]

# ── User-defined options ──────────────────────────────────────────────────────
remove_polymorphic = True   # set to False to skip this filter

# ── Resolve lists (file path, inline list, or None) ───────────────────────────
def resolve_list(value):
    if value is None:
        return None
    if isinstance(value, list):
        return [v.strip() for v in value if str(v).strip()]
    with open(value) as f:
        return [line.strip() for line in f if line.strip()]

genome_list = resolve_list(genome_list)
gene_list   = resolve_list(gene_list)
region_list = resolve_list(region_list)

# ── Load ──────────────────────────────────────────────────────────────────────
df = pd.read_csv(annotated_vcf, sep="\t", low_memory=False)

# ── Column filtering ──────────────────────────────────────────────────────────
info_idx = df.columns.get_loc("INFO")
cols_before_info = df.columns[:info_idx + 2].tolist()  # +2 to include FORMAT after INFO

if genome_list is not None:
    cols_after_info = [c for c in df.columns[info_idx + 2:] if c in genome_list]
else:
    cols_after_info = df.columns[info_idx + 2:].tolist()

df = df[cols_before_info + cols_after_info]

# ── Row filtering: genes ──────────────────────────────────────────────────────
if gene_list is not None:
    mask = df["RefSeq_gene"].isin(gene_list) | df["Ensembl_gene"].isin(gene_list)
    df = df[mask]

# ── Row filtering: regions ────────────────────────────────────────────────────
if region_list is not None:
    mask = df["RefSeq_region"].isin(region_list) | df["Ensembl_region"].isin(region_list)
    df = df[mask]

# ── Remove non-polymorphic rows (optional) ────────────────────────────────────
# A row is kept if at least one sample column (after FORMAT) contains '0/1' or '1/1'
if remove_polymorphic:
    format_idx = df.columns.get_loc("FORMAT")
    sample_cols = df.columns[format_idx + 1:]
    gt_pattern = r"(^|:)(0/1|1/1)(:|$)"
    has_variant = df[sample_cols].apply(
        lambda col: col.astype(str).str.contains(gt_pattern, regex=True)
    ).any(axis=1)
    df = df[has_variant]

# ── Output ────────────────────────────────────────────────────────────────────
df.to_csv(output_file, sep="\t", index=False)

print(f"Done. {len(df)} rows retained.")
print(f"Output written to: {output_file}")

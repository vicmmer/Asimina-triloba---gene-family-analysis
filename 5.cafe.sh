#!/usr/bin/env bash
set -euo pipefail
#conda activate cafe5
# ---------- settings you can tweak ----------
INPUT="cafe_gene_families.tsv"
TREE="SpeciesTree_ultrametric.txt"
TOP_N=40              # drop top-N max–min outliers (try 20–50)
K=2                   # gamma categories (start with 2)
ITER=1000             # optimizer iterations
OUTDIR="results_k2_poisson"
# -------------------------------------------

echo "[1/4] Rank families by (max-min) and select top $TOP_N to drop..."
awk 'BEGIN{FS=OFS="\t"} NR==1{next}
{
  max=-1; min=1e18;
  for(i=3;i<=NF;i++){ if($i>max) max=$i; if($i<min) min=$i }
  print $2, max-min   # $2 = Family ID
}' "$INPUT" | sort -k2,2nr > diff_by_family.tsv

head -n "$TOP_N" diff_by_family.tsv | cut -f1 > exclude_ids.txt

echo "[2/4] Create filtered input..."
awk 'BEGIN{FS=OFS="\t"} NR==FNR{bad[$1]=1; next}
     NR==1 || !($2 in bad)' exclude_ids.txt "$INPUT" > cafe_input_filtered.tsv

echo "[3/4] Run CAFE5 (k=$K, Poisson root, I=$ITER)..."
rm -rf "$OUTDIR"
cafe5 -i cafe_input_filtered.tsv -t "$TREE" -k "$K" -p -I "$ITER" -o "$OUTDIR"

echo "[4/4] Done. Report file: $OUTDIR/Gamma_report.cafe"
~

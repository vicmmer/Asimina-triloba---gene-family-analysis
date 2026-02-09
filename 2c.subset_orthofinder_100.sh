#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash 2c.subset_orthofinder_100.sh /path/to/OrthoFinder/Results_XXXX 100 123
#
# Args:
#   1) OrthoFinder Results directory (e.g., OrthoFinder/Results_Feb09/ )
#   2) Number of orthogroups to keep (default 100)
#   3) Seed for random sampling (default 123). Use "NONE" to take first N deterministically.

RESULTS_DIR="${1:?Provide OrthoFinder Results dir}"
N_KEEP="${2:-100}"
SEED="${3:-123}"

# Normalize path (remove trailing slash)
RESULTS_DIR="${RESULTS_DIR%/}"

ORTHOGROUPS_DIR="$RESULTS_DIR/Orthogroups"
OGSEQ_DIR="$RESULTS_DIR/Orthogroup_Sequences"

if [[ ! -d "$ORTHOGROUPS_DIR" ]]; then
  echo "ERROR: Can't find Orthogroups folder at: $ORTHOGROUPS_DIR" >&2
  exit 1
fi
if [[ ! -d "$OGSEQ_DIR" ]]; then
  echo "ERROR: Can't find Orthogroup_Sequences folder at: $OGSEQ_DIR" >&2
  exit 1
fi

OUT_DIR="${RESULTS_DIR}_SUBSET${N_KEEP}"
mkdir -p "$OUT_DIR/Orthogroups" "$OUT_DIR/Orthogroup_Sequences"

# 1) Get list of orthogroup IDs from Orthogroups.tsv (expects OG0000001-style first column)
OG_TSV="$ORTHOGROUPS_DIR/Orthogroups.tsv"
if [[ ! -f "$OG_TSV" ]]; then
  echo "ERROR: Can't find $OG_TSV" >&2
  exit 1
fi

tmp_all="$(mktemp)"
tmp_pick="$(mktemp)"
trap 'rm -f "$tmp_all" "$tmp_pick"' EXIT

# grab OG IDs from first column, skipping header if present
awk 'NR==1 && $1 ~ /^Orthogroup/ {next} {print $1}' "$OG_TSV" | sed '/^$/d' > "$tmp_all"

TOTAL=$(wc -l < "$tmp_all" | tr -d ' ')
if (( TOTAL < N_KEEP )); then
  echo "ERROR: Only $TOTAL orthogroups found, cannot keep $N_KEEP" >&2
  exit 1
fi

# 2) Select N_KEEP orthogroups
if [[ "$SEED" == "NONE" ]]; then
  # deterministic: first N after sorting
  sort "$tmp_all" | head -n "$N_KEEP" > "$tmp_pick"
else
  # reproducible random sample using python seed
  python3 - <<'PY' "$tmp_all" "$tmp_pick" "$N_KEEP" "$SEED"
import sys, random
inp, outp, n, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
random.seed(seed)
with open(inp) as f:
    ogs = [ln.strip() for ln in f if ln.strip()]
random.shuffle(ogs)
with open(outp, "w") as f:
    f.write("\n".join(ogs[:n]) + "\n")
PY
fi

echo "Selected $(wc -l < "$tmp_pick" | tr -d ' ') orthogroups → $OUT_DIR"
cp "$tmp_pick" "$OUT_DIR/Orthogroups/subset_orthogroups.list"

# 3) Filter Orthogroups.tsv (keep header if present)
# Make a fast lookup table and filter rows by OG id.
awk 'NR==1 && $1 ~ /^Orthogroup/ {print; next} FNR==NR {keep[$1]=1; next} ($1 in keep) {print}' \
  "$tmp_pick" "$OG_TSV" > "$OUT_DIR/Orthogroups/Orthogroups.tsv"

# 4) Filter GeneCount file if it exists
GC_TSV="$ORTHOGROUPS_DIR/Orthogroups.GeneCount.tsv"
if [[ -f "$GC_TSV" ]]; then
  awk 'NR==1 {print; next} FNR==NR {keep[$1]=1; next} ($1 in keep) {print}' \
    "$tmp_pick" "$GC_TSV" > "$OUT_DIR/Orthogroups/Orthogroups.GeneCount.tsv"
fi

# 5) Copy only the orthogroup fasta files
# OrthoFinder usually names these like Orthogroup_Sequences/OG0000001.fa
while read -r og; do
  # handle .fa or .fasta (and gz variants) if present
  found=0
  for ext in fa fasta fa.gz fasta.gz; do
    f="$OGSEQ_DIR/${og}.${ext}"
    if [[ -f "$f" ]]; then
      cp "$f" "$OUT_DIR/Orthogroup_Sequences/"
      found=1
      break
    fi
  done
  if [[ $found -eq 0 ]]; then
    echo "WARNING: No fasta found for $og in $OGSEQ_DIR (skipping)" >&2
  fi
done < "$tmp_pick"

echo "Done. Subset results at: $OUT_DIR"
echo "Next: run your 2b/3/4/5 steps using $OUT_DIR instead of $RESULTS_DIR"

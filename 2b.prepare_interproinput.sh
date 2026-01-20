# Purpose:
# Prepare OrthoFinder orthogroup protein FASTA files for InterProScan by removing stop codons (*) from sequence lines while preserving FASTA headers.

#!/usr/bin/env bash
set -e

clean_dir="interproscan_input"
mkdir -p "$clean_dir"

# Auto-detect newest OrthoFinder Orthogroup_Sequences folder
orig_dir=$(ls -td orthofinder/Results_*/Orthogroup_Sequences | head -n 1)

echo " Cleaning OrthoFinder FASTAs for InterProScan "
echo "Input:  $orig_dir"
echo "Output: $clean_dir"
echo

for file in "$orig_dir"/*.fa; do
  out="$clean_dir/$(basename "$file")"
  echo "Cleaning $(basename "$file")"

  awk '
    /^>/ { print; next }
    { gsub(/\*/, "", $0); print }
  ' "$file" > "$out"
done

echo
echo "=== Done! Cleaned FASTAs are in: $clean_dir ==="

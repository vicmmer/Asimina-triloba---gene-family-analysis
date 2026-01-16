#!/usr/bin/env bash
set -euo pipefail

# ====== INPUT / OUTPUT DIRS (EDIT IF NEEDED) ======
orig_dir="orthofinder/Results_Dec03/Orthogroup_Sequences"
clean_dir="interproscan_input"

mkdir -p "$clean_dir"

echo "Cleaning FASTA files from: $orig_dir"
echo "Writing cleaned FASTAs to: $clean_dir"
echo

shopt -s nullglob
found_any=false

for file in "$orig_dir"/*.fa; do
    found_any=true
    base=$(basename "$file")
    out="${clean_dir}/${base}"

    echo "Cleaning $file -> $out"

    awk '
        # Header lines: keep as-is
        /^>/ { print; next }
        # Sequence lines: remove asterisks only
        {
            gsub(/\*/, "", $0)
            print
        }
    ' "$file" > "$out"
done

if [ "$found_any" = false ]; then
    echo "No .fa files found in $orig_dir"
    exit 1
fi

echo
echo "Done! Cleaned FASTA files are in: $clean_dir"

# Purpose: Run InterProScan on a single cleaned protein FASTA (from interproscan_input/) and write a TSV result to interproscan_output using the same basename.

#!/bin/bash
set -e

INPUT_DIR="interproscan_input"
OUTPUT_DIR="interproscan_output"
IPS_PATH="/home/vmartinez/my_interproscan/interproscan-5.75-106.0/interproscan.sh"

JOBS=40   # number of parallel jobs
CPU=1     # CPUs per InterProScan run

mkdir -p "$OUTPUT_DIR"

echo "=== Running InterProScan with GNU parallel ==="
echo "Input: $INPUT_DIR"
echo "Jobs:  $JOBS"
echo

ls "$INPUT_DIR"/*.fa | parallel -j "$JOBS" --eta --verbose '
  file={}
  base=$(basename "$file" .fa)
  out="'"$OUTPUT_DIR"'/${base}.tsv"

  '"$IPS_PATH"' \
    -i "$file" \
    -f tsv \
    -appl Pfam,PANTHER \
    --iprlookup \
    --goterms \
    -cpu '"$CPU"' \
    -o "$out"
'

echo " InterProScan finished. Outputs in: $OUTPUT_DIR "






#!/bin/bash
set -e

# Set lineage and mode
LINEAGE="embryophyta_odb10"
MODE="protein"
THREADS=30

# Parent results folder
OUT_PARENT="busco_results"
mkdir -p "$OUT_PARENT"

# Loop over each FASTA file
# (handles *.fa and *.fasta; skips if none)
for fasta in *.fa *.fasta; do
    [ -e "$fasta" ] || continue   # skip pattern when it matches nothing

    # Get file name only
    fname=$(basename "$fasta")

    # Strip extension (everything after the first dot)
    # e.g. Annona_cherimola.fa            -> Annona_cherimola
    #      Annona_muricata_annotated.fa   -> Annona_muricata_annotated
    BASENAME="${fname%%.*}"

    echo "Running BUSCO on $fasta (basename: $BASENAME)..."

    busco \
        -i "$fasta" \
        -l "$LINEAGE" \
        -m "$MODE" \
        -c "$THREADS" \
        -o "$BASENAME" \
        --out_path "$OUT_PARENT" \
        -f
done

echo "All BUSCO analyses complete."

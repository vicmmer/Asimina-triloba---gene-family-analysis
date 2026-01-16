# to run outside of script: conda activate orthofinder_env

#!/bin/bash

INPUT_DIR="protein_sequences"
OUTPUT_DIR="orthofinder"

THREADS=30

orthofinder -f "$INPUT_DIR" \
            -o "$OUTPUT_DIR" \
            -t "$THREADS" \
            -a 10 \
            -S diamond

echo "OrthoFinder finished! Check inside: $OUTPUT_DIR/OrthoFinder/"

#!/bin/bash

BUSCO_DIR="busco_results"
OUTFILE="busco_summary_table.tsv"

# header
echo -e "Species\tC_pct\tS_pct\tD_pct\tF_pct\tM_pct\tn\tC_count\tS_count\tD_count\tF_count\tM_count\tTotal_count" > "$OUTFILE"

for dir in "$BUSCO_DIR"/*; do
    if [ -d "$dir" ]; then
        species=$(basename "$dir")

        summary=$(ls "$dir"/short_summary*.txt 2>/dev/null)
        if [ -z "$summary" ]; then
            echo "Warning: no summary for $species"
            continue
        fi

        # ----- percent line -----
        line=$(grep "C:" "$summary" | head -n1 | sed 's/^[[:space:]]*//')
        if [ -z "$line" ]; then
            echo "Warning: could not parse percent line for $species"
            continue
        fi

        # extract percentages and n
        C_pct=$(echo "$line" | sed -n 's/.*C:\([0-9.]\+\)%.*/\1/p')
        S_pct=$(echo "$line" | sed -n 's/.*S:\([0-9.]\+\)%.*/\1/p')
        D_pct=$(echo "$line" | sed -n 's/.*D:\([0-9.]\+\)%.*/\1/p')
        F_pct=$(echo "$line" | sed -n 's/.*F:\([0-9.]\+\)%.*/\1/p')
        M_pct=$(echo "$line" | sed -n 's/.*M:\([0-9.]\+\)%.*/\1/p')
        n=$(echo  "$line" | sed -n 's/.*n:\([0-9]\+\).*/\1/p')

        # ----- counts -----
        C_count=$(grep "Complete BUSCOs (C)" "$summary" | awk '{print $1}')
        S_count=$(grep "Complete and single-copy BUSCOs (S)" "$summary" | awk '{print $1}')
        D_count=$(grep "Complete and duplicated BUSCOs (D)" "$summary" | awk '{print $1}')
        F_count=$(grep "Fragmented BUSCOs (F)" "$summary" | awk '{print $1}')
        M_count=$(grep "Missing BUSCOs (M)" "$summary" | awk '{print $1}')
        Total_count=$(grep "Total BUSCO groups searched" "$summary" | awk '{print $1}')

        echo -e "${species}\t${C_pct}\t${S_pct}\t${D_pct}\t${F_pct}\t${M_pct}\t${n}\t${C_count}\t${S_count}\t${D_count}\t${F_count}\t${M_count}\t${Total_count}" >> "$OUTFILE"
    fi
done

echo "BUSCO summary table written to $OUTFILE"

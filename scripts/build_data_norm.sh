#!/usr/bin/env bash
set -euo pipefail

RAW_DIR="data"
NORM_DIR="data_norm"

if [[ ! -d "$RAW_DIR" ]]; then
    echo "ERROR: raw data directory missing: $RAW_DIR"
    exit 1
fi

rm -rf "$NORM_DIR"
mkdir -p "$NORM_DIR"

for f in "$RAW_DIR"/*.fastq.gz; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")

    # Strip extensions
    stem="${base%.fastq.gz}"

    # Identify sample + read
    if [[ "$stem" =~ (.+)_R1_001$ ]]; then
        sample="${BASH_REMATCH[1]}"
        read="R1"
    elif [[ "$stem" =~ (.+)_R2_001$ ]]; then
        sample="${BASH_REMATCH[1]}"
        read="R2"
    elif [[ "$stem" =~ (.+)_1$ ]]; then
        sample="${BASH_REMATCH[1]}"
        read="R1"
    elif [[ "$stem" =~ (.+)_2$ ]]; then
        sample="${BASH_REMATCH[1]}"
        read="R2"
    elif [[ "$stem" =~ ^[^_]+$ ]]; then
        sample="$stem"
        read="SE"
    else
        echo "Skipping unrecognized file: $base"
        continue
    fi

    ln -sf "../$f" "${NORM_DIR}/${sample}_${read}.fastq.gz"
done

echo "Normalized files:"
ls -1 "$NORM_DIR" | head

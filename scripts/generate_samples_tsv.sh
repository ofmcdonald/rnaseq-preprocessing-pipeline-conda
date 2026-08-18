#!/usr/bin/env bash
set -euo pipefail

OUTPUT="samples.tsv"
DATA_DIR="data_norm"

echo "sample" > "$OUTPUT"

R1_COUNT=0
R2_COUNT=0
SE_COUNT=0

# Check that the data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
    echo "ERROR: normalized data directory missing: $DATA_DIR"
    exit 1
fi

# Count and validate R1 files
for f in "$DATA_DIR"/*_R1.fastq.gz; do
    [[ -e "$f" ]] || continue

    base=$(basename "$f")
    sample="${base%_R1.fastq.gz}"

    R1_COUNT=$((R1_COUNT + 1))

    if [[ ! -f "$DATA_DIR/${sample}_R2.fastq.gz" ]]; then
        echo "ERROR: missing R2 for paired-end sample $sample"
        exit 1
    fi
done

# Count R2 files
for f in "$DATA_DIR"/*_R2.fastq.gz; do
    [[ -e "$f" ]] || continue
    R2_COUNT=$((R2_COUNT + 1))
done

# Count single-end files
for f in "$DATA_DIR"/*_SE.fastq.gz; do
    [[ -e "$f" ]] || continue
    SE_COUNT=$((SE_COUNT + 1))
done

echo "Found $R1_COUNT R1 samples"
echo "Found $R2_COUNT R2 samples"
echo "Found $SE_COUNT single-end samples"

# Generate sample list
for f in "$DATA_DIR"/*.fastq.gz; do
    [[ -e "$f" ]] || continue

    base=$(basename "$f")

    if [[ "$base" == *_R1.fastq.gz ]]; then
        sample="${base%_R1.fastq.gz}"
        echo "$sample" >> "$OUTPUT"

    elif [[ "$base" == *_SE.fastq.gz ]]; then
        sample="${base%_SE.fastq.gz}"
        echo "$sample" >> "$OUTPUT"

    elif [[ "$base" == *_R2.fastq.gz ]]; then
        # R2 is already represented by the corresponding R1
        continue
    fi
done

# Remove duplicate samples and sort
{
    head -n 1 "$OUTPUT"
    tail -n +2 "$OUTPUT" | sort -u
} > "${OUTPUT}.tmp"

mv "${OUTPUT}.tmp" "$OUTPUT"

echo "Generated $OUTPUT"
column -t "$OUTPUT"

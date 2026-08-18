#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE_BIN="snakemake"

# Resolve repository root relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "$WORKDIR"

"${SNAKEMAKE_BIN}" \
  --profile profiles/slurm \
  --executor slurm \
  -j 100 \
  "$@"

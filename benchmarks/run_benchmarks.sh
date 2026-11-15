#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/home/jayaram/Projects/biocpy/genomicranges/.venv/bin/python}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

echo "1. Generating shared random genomic intervals (seed=42)..."
"$PYTHON_BIN" "$SCRIPT_DIR/generate_intervals.py"

echo ""
echo "2. Running BiocPy GenomicRanges benchmark (1 and 4 threads)..."
"$PYTHON_BIN" "$SCRIPT_DIR/benchmark_python.py"

echo ""
echo "3. Running R Bioconductor GenomicRanges benchmark..."
"$RSCRIPT_BIN" "$SCRIPT_DIR/benchmark_r.R"


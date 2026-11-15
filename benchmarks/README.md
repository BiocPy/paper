# GenomicRanges Overlap and Proximity Benchmarks

This directory contains reproducible scripts for generating the `GenomicRanges` runtime comparison (Table 2 in the manuscript) between R/Bioconductor (`GenomicRanges`) and Python (`genomicranges` / `iranges` backed by `nclist-cpp`).

## Design

To guarantee an exact comparison across languages and thread counts:
1. [`generate_intervals.py`](./generate_intervals.py) samples random genomic intervals (`seqnames`, `start`, `width`, `strand`) across the 24 human chromosomes (`chr1`–`chr22`, `chrX`, `chrY`) using a fixed seed (`SEED = 42`) and writes the subject and query tables to `data/*.csv`.
2. [`benchmark_python.py`](./benchmark_python.py) loads those CSV files into BiocPy `GenomicRanges` objects and times `find_overlaps()` and `nearest()` across `num_threads=1` and `num_threads=4` (mean of 3 runs), verifying that all thread counts return identical results.
3. [`benchmark_r.R`](./benchmark_r.R) loads the exact same CSV files into Bioconductor `GRanges` objects and times `findOverlaps()` and `nearest()` (mean of 3 runs), verifying identical overlap hit counts with Python.

## Running the Benchmarks

```bash
bash run_benchmarks.sh
```

Output summary CSVs are written to `results_python.csv` and `results_r.csv`.


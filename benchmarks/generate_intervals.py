#!/usr/bin/env python3
"""Generate reproducible random genomic intervals for R and Python benchmarks.

Saves paired subject and query interval tables as CSV files so that both
R (Bioconductor GenomicRanges) and Python (BiocPy genomicranges) evaluate
findOverlaps/find_overlaps and nearest on the exact same genomic coordinates.
"""

from pathlib import Path
import numpy as np
import pandas as pd

SEED = 42
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
STRANDS = ["+", "-", "*"]
BENCHMARK_PAIRS = [
    (100_000, 100_000),
    (1_000_000, 100_000),
    (1_000_000, 1_000_000),
]


def sample_intervals(rng: np.random.Generator, n: int) -> pd.DataFrame:
    """Sample n random genomic intervals across 24 human chromosomes."""
    seqnames = rng.choice(CHROMS, size=n)
    start = rng.integers(1, 100_000_000, size=n, dtype=np.int32)
    width = rng.integers(100, 2_000, size=n, dtype=np.int32)
    strand = rng.choice(STRANDS, size=n)
    return pd.DataFrame(
        {
            "seqnames": seqnames,
            "start": start,
            "width": width,
            "strand": strand,
        }
    )


def main() -> None:
    data_dir = Path(__file__).resolve().parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    for n_subj, n_query in BENCHMARK_PAIRS:
        tag = f"subj_{n_subj}_query_{n_query}"
        subj_path = data_dir / f"{tag}_subj.csv"
        query_path = data_dir / f"{tag}_query.csv"

        subj_df = sample_intervals(rng, n_subj)
        query_df = sample_intervals(rng, n_query)

        subj_df.to_csv(subj_path, index=False)
        query_df.to_csv(query_path, index=False)
        print(f"Generated {subj_path.name} ({n_subj} rows) and {query_path.name} ({n_query} rows)")


if __name__ == "__main__":
    main()


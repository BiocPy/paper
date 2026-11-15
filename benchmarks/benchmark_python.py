#!/usr/bin/env python3
"""Benchmark BiocPy GenomicRanges find_overlaps and nearest across threads."""

from pathlib import Path
import time
import numpy as np
import pandas as pd
from iranges import IRanges
from genomicranges import GenomicRanges

from generate_intervals import BENCHMARK_PAIRS, main as generate_data

N_RUNS = 3
THREADS = [1, 4]


def load_granges(csv_path: Path) -> GenomicRanges:
    """Load a CSV interval file into a BiocPy GenomicRanges object."""
    df = pd.read_csv(csv_path)
    return GenomicRanges(
        seqnames=df["seqnames"].tolist(),
        ranges=IRanges(
            df["start"].to_numpy(dtype=np.int32),
            df["width"].to_numpy(dtype=np.int32),
        ),
        strand=df["strand"].tolist(),
    )


def main() -> None:
    data_dir = Path(__file__).resolve().parent / "data"
    if not data_dir.exists() or not any(data_dir.glob("*.csv")):
        print("Interval CSVs not found; generating reproducible intervals first...")
        generate_data()

    print("=== BiocPy GenomicRanges Benchmark ===")
    results = []

    for n_subj, n_query in BENCHMARK_PAIRS:
        tag = f"subj_{n_subj}_query_{n_query}"
        subj_gr = load_granges(data_dir / f"{tag}_subj.csv")
        query_gr = load_granges(data_dir / f"{tag}_query.csv")

        row = {"n_subj": n_subj, "n_query": n_query}
        ref_hits_len = None
        ref_nearest = None

        for t in THREADS:
            t0 = time.perf_counter()
            for _ in range(N_RUNS):
                hits = subj_gr.find_overlaps(query_gr, num_threads=t)
            row[f"find_overlaps_{t}T"] = (time.perf_counter() - t0) / N_RUNS

            if ref_hits_len is None:
                ref_hits_len = len(hits)
                row["num_hits"] = ref_hits_len
            else:
                assert len(hits) == ref_hits_len, (
                    f"Hit count mismatch between threads: {len(hits)} vs {ref_hits_len}"
                )

            t0 = time.perf_counter()
            for _ in range(N_RUNS):
                near_idx = subj_gr.nearest(query_gr, num_threads=t)
            row[f"nearest_{t}T"] = (time.perf_counter() - t0) / N_RUNS

            if ref_nearest is None:
                ref_nearest = near_idx
            else:
                assert np.array_equal(near_idx, ref_nearest)

        results.append(row)
        print(
            f"Subject={n_subj:,} x Query={n_query:,} (hits={row['num_hits']:,}): "
            f"find_overlaps 1T={row['find_overlaps_1T']:.3f}s, 4T={row['find_overlaps_4T']:.3f}s | "
            f"nearest 1T={row['nearest_1T']:.3f}s, 4T={row['nearest_4T']:.3f}s"
        )

    out_csv = Path(__file__).resolve().parent / "results_python.csv"
    pd.DataFrame(results).to_csv(out_csv, index=False)
    print(f"Saved Python benchmark results to {out_csv}")


if __name__ == "__main__":
    main()

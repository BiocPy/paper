#!/usr/bin/env Rscript
# Benchmark Bioconductor GenomicRanges findOverlaps and nearest on the exact
# same interval datasets used by benchmark_python.py.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
data_dir <- file.path(script_dir, "data")

benchmark_pairs <- list(
  c(100000L, 100000L),
  c(1000000L, 100000L),
  c(1000000L, 1000000L)
)
n_runs <- 3L

load_granges <- function(csv_path) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  GRanges(
    seqnames = df$seqnames,
    ranges = IRanges(start = df$start, width = df$width),
    strand = df$strand
  )
}

cat(sprintf("=== R Bioconductor GenomicRanges (v%s) Benchmark ===\n",
            as.character(packageVersion("GenomicRanges"))))

results <- data.frame(
  n_subj = integer(),
  n_query = integer(),
  num_hits = integer(),
  findOverlaps_R = numeric(),
  nearest_R = numeric(),
  stringsAsFactors = FALSE
)

for (pair in benchmark_pairs) {
  n_subj <- pair[1]
  n_query <- pair[2]
  tag <- sprintf("subj_%d_query_%d", n_subj, n_query)

  subj_gr <- load_granges(file.path(data_dir, paste0(tag, "_subj.csv")))
  query_gr <- load_granges(file.path(data_dir, paste0(tag, "_query.csv")))

  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(n_runs)) {
    hits <- findOverlaps(query_gr, subj_gr)
  }
  t_fo <- (proc.time()[["elapsed"]] - t0) / n_runs

  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(n_runs)) {
    near_idx <- nearest(query_gr, subj_gr)
  }
  t_near <- (proc.time()[["elapsed"]] - t0) / n_runs

  cat(sprintf(
    "Subject=%s x Query=%s (hits=%s): findOverlaps=%.3fs | nearest=%.3fs\n",
    format(n_subj, big.mark = ",", scientific = FALSE),
    format(n_query, big.mark = ",", scientific = FALSE),
    format(length(hits), big.mark = ",", scientific = FALSE),
    t_fo,
    t_near
  ))

  results <- rbind(
    results,
    data.frame(
      n_subj = n_subj,
      n_query = n_query,
      num_hits = length(hits),
      findOverlaps_R = t_fo,
      nearest_R = t_near
    )
  )
}

out_csv <- file.path(script_dir, "results_r.csv")
write.csv(results, out_csv, row.names = FALSE)
cat(sprintf("Saved R benchmark results to %s\n", out_csv))


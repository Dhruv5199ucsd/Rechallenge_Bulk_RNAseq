#!/usr/bin/env Rscript

# Hashtag demultiplexing wrapper for the Kallies Sci Immuno 2025 Scanpy workflow.
# The analysis logic lives in src/hashtag_demultiplex.py so it can use the same
# AnnData/Scanpy environment as the surrounding notebooks.

repo_dir <- "/home/dk5299/Projects_31926/RNA-seq/Kallies_Sci_Immuno_2025"
setwd(repo_dir)
Sys.setenv(MPLCONFIGDIR = "/tmp/matplotlib-kallies-sci-immuno-2025")

message("Running hashtag demultiplexing from: ", repo_dir)

project_python <- file.path(repo_dir, ".venv", "bin", "python")
if (file.exists(project_python)) {
  status <- system2(
    project_python,
    args = c("-m", "src.hashtag_demultiplex"),
    stdout = "",
    stderr = ""
  )
} else {
  status <- system2(
    "uv",
    args = c("run", "python", "-m", "src.hashtag_demultiplex"),
    stdout = "",
    stderr = ""
  )
}

if (!identical(status, 0L)) {
  stop("Hashtag demultiplexing failed with exit status: ", status, call. = FALSE)
}

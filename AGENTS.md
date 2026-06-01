# Project Instructions

## Project Layout

This repository contains the code-facing side of the Kallies Sci Immuno 2025 single-cell RNA-seq analysis. Keep large data and derived `.h5ad` objects on SATA storage through `src/paths.py`.

- Repository: `/home/dk5299/Projects_31926/RNA-seq/Kallies_Sci_Immuno_2025`
- Data storage: `/mnt/sata4/Dhruv_2025/Projects_31926/RNA_Seq/Kallies_Sci_Immuno_2025`
- Raw 10x files: `DATA_DIR`
- Intermediate AnnData objects: `INTERMEDIATE_DIR`
- Tables and final outputs: `RESULTS_DIR`

Do not reorganize the project unless explicitly requested. Follow the existing `src/paths.py` constants instead of hard-coding alternate locations.

## Notebook Order

Run the analysis in numbered order, saving one object per step:

1. `notebooks/00_Load_Data.ipynb`
2. `notebooks/01_QC_Filtering.ipynb`
3. `notebooks/01b_hashtag_demultiplexing.R`
4. `notebooks/02_Normalization.ipynb`
5. `notebooks/03_clustering.ipynb`
6. `notebooks/04_Cell_Annotation.ipynb`
7. `notebooks/05_Memory_vs_Rechallenge.ipynb`
8. `notebooks/06_Marker_Gene_Dotplots_LogNormalized.ipynb`
9. `notebooks/07_Bulk_RNAseq.ipynb`

Notebook `07_Bulk_RNAseq.ipynb` is a companion bulk RNA-seq visualization workflow using the Von Hoesslin Sci Immuno 2022 project outputs. It is not part of the single-cell preprocessing chain.

## Hashtag Mapping

The local hashtag reference is `DATA_DIR / "GSE287643_hashtag_reference.csv"`. The raw feature matrix includes three `Antibody Capture` features: `HT1`, `HT2`, and `HT3`.

Use the Man 2025 processing workflow as the interpretation reference:

- `HT1` maps to memory small-intestinal IEL: `Mem_si_IEL`
- `HT2` maps to memory small-intestinal lamina propria: `Mem_si_LP`
- cells not assigned to `HT1` or `HT2` are treated as rechallenged small-intestinal IEL: `Rechallenged_si_IEL`

For this local WT-only matrix, use these analysis labels:

- `HT1` -> `sIEL`
- `HT2` -> `sLPL`
- `Negative` or `HT3` -> `rechallenge`

The demultiplexing step should write a metadata CSV and subset `.h5ad` objects so that downstream notebooks can load explicit sample groups rather than inferring them again.

## Implementation Conventions

- Use `uv run python ...` for Scanpy/AnnData work.
- Keep notebooks importable from the repo root by appending the repo path or using `src.paths`.
- Save large objects only under `INTERMEDIATE_DIR`.
- Save small summary tables under `RESULTS_DIR`.
- Preserve raw counts in a layer before normalization.
- Do not delete raw data or overwrite earlier workflow objects without a clear reason.

## Bulk RNA-seq Companion Workflow

Use `notebooks/07_Bulk_RNAseq.ipynb` to compare Kallies candidate genes against the Von Hoesslin bulk RNA-seq time course.

- External source project: `/home/dk5299/Projects_31926/RNA-seq/Von_Hoesslin_Sci_Immuno_2022`
- External DESeq2 tables: `results/tables/02_deseq`
- External count matrix: `count matrix/GSE208097_count_matrix_compiled`
- Local outputs: `RESULTS_DIR / "07_bulk_rnaseq"`
- Local web viewer: `web/bulk_rnaseq_viewer.html`

The bulk notebook estimates DESeq2-style median-ratio size factors across all Von Hoesslin samples, validates the normalized counts against the existing Von Hoesslin R-normalized table, and exports per-sample candidate-gene normalized counts for the web viewer.

When updating candidate genes, rerun at least:

1. the `candidate_genes = [...]` cell
2. the candidate DESeq2 stats export cell
3. the DESeq2-normalized expression export cell

The web viewer needs these files from `RESULTS_DIR / "07_bulk_rnaseq"`:

- required: `07_bulk_rnaseq_candidate_gene_deseq_normalized_counts_by_sample_d0_d2_d8.tsv`
- optional padj stats: `07_bulk_rnaseq_candidate_gene_d2_d8_deseq.tsv`

Use the web viewer gene selector for large gene lists instead of rendering all genes at once.

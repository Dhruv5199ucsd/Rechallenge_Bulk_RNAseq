# Kallies Sci Immuno 2025 scRNA-seq Analysis

This repository contains the code-facing workflow for the Kallies Sci Immuno
2025 single-cell RNA-seq analysis. Large AnnData objects and generated outputs
are stored on SATA storage through `src/paths.py`.

## Key Paths

- Repository: `/home/dk5299/Projects_31926/RNA-seq/Kallies_Sci_Immuno_2025`
- Storage root: `/mnt/sata4/Dhruv_2025/Projects_31926/RNA_Seq/Kallies_Sci_Immuno_2025`
- Raw/input data: `DATA_DIR`
- Intermediate `.h5ad` objects: `INTERMEDIATE_DIR`
- Tables and figures: `RESULTS_DIR`

Use `uv run python ...` for Python/Scanpy work.

## Active Notebook Workflow

Run the main workflow in numbered order:

1. `notebooks/00_Load_Data.ipynb`
2. `notebooks/01_QC_Filtering.ipynb`
3. `notebooks/01b_hashtag_demultiplexing.R`
4. `notebooks/02_Normalization.ipynb`
5. `notebooks/03_clustering.ipynb`
6. `notebooks/04_Cell_Annotation.ipynb`
7. `notebooks/05_Memory_vs_Rechallenge.ipynb`
8. `notebooks/06_Marker_Gene_Dotplots_LogNormalized.ipynb`
9. `notebooks/07_Bulk_RNAseq.ipynb`

Older exploratory notebooks were moved to:

```text
notebooks/archive_exploratory/
```

Older exploratory result folders/files were moved to:

```text
RESULTS_DIR / "archive_exploratory"
```

## Hashtag Mapping

The local WT-only matrix uses:

- `HT1` -> `sIEL`
- `HT2` -> `sLPL`
- `Negative` or `HT3` -> `rechallenge`

Downstream condition labels:

- `sIEL` -> `IEL_memory`
- `sLPL` -> `LPL_memory`
- `rechallenge` -> `rechallenge`

## Normalization Convention

`notebooks/02_Normalization.ipynb` was updated to fit the pipeline.

It now:

- reads the demultiplexed singlet and subset objects from notebook `01b`
- preserves raw counts in `layers["counts"]`
- normalizes and log-transforms expression
- stores full-gene log-normalized expression in `.raw`
- runs HVG selection, scaling, and PCA
- writes the normalized objects expected by notebook `03`

Important outputs:

- `INTERMEDIATE_DIR / "02_normalized_hvg_pca.h5ad"`
- `INTERMEDIATE_DIR / "02_normalized_hvg_pca_sIEL.h5ad"`
- `INTERMEDIATE_DIR / "02_normalized_hvg_pca_sLPL.h5ad"`
- `INTERMEDIATE_DIR / "02_normalized_hvg_pca_rechallenge.h5ad"`
- `RESULTS_DIR / "02_normalization_summary.csv"`

Notebook `05_Memory_vs_Rechallenge.ipynb` was also updated so DE rebuilds
log-normalized expression from `layers["counts"]` instead of assuming `.raw`
contains raw counts.

## Marker Gene Dot Plot Workflow

Use this for comparing marker genes between conditions on log-normalized values:

```text
notebooks/06_Marker_Gene_Dotplots_LogNormalized.ipynb
```

The notebook explicitly rebuilds expression from raw counts:

```python
expr_adata.X = expr_adata.layers["counts"].copy()
sc.pp.normalize_total(expr_adata, target_sum=1e4)
sc.pp.log1p(expr_adata)
```

It does not use scaled `.X` values.

### Marker Gene Sheet

Edit this CSV to change the plotted genes:

```text
data/marker_gene_panels.csv
```

Required columns:

```text
panel,gene,display_order
```

Current default marker sheet:

```text
egress,Ccr7,1
egress,S1pr1,2
egress,Klf2,3
retention,Cd69,4
retention,Itgae,5
retention,Ccr9,6
retention,Cxcr6,7
```

### Current Marker Dot Plot Outputs

All current marker dot plot outputs are written here:

```text
RESULTS_DIR / "06_marker_gene_dotplots_lognormalized"
```

Main figures:

- `marker_gene_dotplot_by_condition_lognormalized.png`
- `marker_gene_dotplot_by_condition_lognormalized.pdf`
- `marker_gene_dotplot_by_cluster_lognormalized.png`
- `marker_gene_dotplot_by_cluster_lognormalized.pdf`
- `marker_gene_dotplot_by_condition_cluster_lognormalized.png`
- `marker_gene_dotplot_by_condition_cluster_lognormalized.pdf`
- `marker_gene_umap_lognormalized_expression.png`
- `marker_gene_umap_lognormalized_expression.pdf`

Main tables:

- `marker_gene_lognorm_summary_by_condition.tsv`
- `marker_gene_lognorm_summary_by_cluster.tsv`
- `marker_gene_lognorm_summary_by_condition_cluster.tsv`
- `marker_gene_mean_lognorm_condition_deltas.tsv`

The dot plots encode:

- color = mean log1p-normalized expression
- dot size = percent of cells detected

## Notes on Interpretation

Use the marker dot plots and summary tables to compare expression between:

- `IEL_memory`
- `LPL_memory`
- `rechallenge`

Because these plots use log-normalized counts, they are appropriate for marker
visualization and condition summaries. They are not based on scaled expression
and should not be interpreted as PCA-scaled values.

## Bulk RNA-seq Companion Workflow

Use this notebook to compare candidate genes against the Von Hoesslin Sci Immuno
2022 bulk RNA-seq time course:

```text
notebooks/07_Bulk_RNAseq.ipynb
```

Inputs come from:

```text
/home/dk5299/Projects_31926/RNA-seq/Von_Hoesslin_Sci_Immuno_2022
```

The notebook uses:

- `results/tables/02_deseq/*_d2_ova_vs_d0_nr_results.tsv`
- `results/tables/02_deseq/*_d8_ova_vs_d0_nr_results.tsv`
- `count matrix/GSE208097_count_matrix_compiled`
- `results/tables/00_loading/00_sample_metadata.tsv`

Expression values are DESeq2-style normalized counts. The notebook estimates
median-ratio size factors across all Von Hoesslin samples, then validates the
result against the existing Von Hoesslin R output:

```text
results/tables/03_downstream/05_gene_barplots/05_selected_gene_normalized_expression.tsv
```

Current bulk outputs are written to:

```text
RESULTS_DIR / "07_bulk_rnaseq"
```

Key files:

- `07_bulk_rnaseq_candidate_gene_deseq_normalized_counts_by_sample_d0_d2_d8.tsv`
- `07_bulk_rnaseq_candidate_gene_d2_d8_deseq.tsv`
- `07_bulk_rnaseq_candidate_gene_padj_summary.tsv`
- `07_bulk_rnaseq_deseq_normalization_validation_against_existing_table.tsv`
- `07_bulk_rnaseq_DP_candidate_gene_timecourse_deseq_normalized_counts_d0_d2_d8.png`
- `07_bulk_rnaseq_DN_candidate_gene_timecourse_deseq_normalized_counts_d0_d2_d8.png`
- `07_bulk_rnaseq_DP_candidate_gene_scaled_dotplot_d0_d2_d8.png`
- `07_bulk_rnaseq_DN_candidate_gene_scaled_dotplot_d0_d2_d8.png`

### Updating Bulk Candidate Genes

Edit the `candidate_genes = [...]` cell in `notebooks/07_Bulk_RNAseq.ipynb`.
For a large list, such as 100 genes, rerun only the downstream export cells if
you plan to use the web viewer:

1. candidate DESeq2 stats export
2. DESeq2-normalized expression export
3. optional summary and notebook plot cells

The web viewer can then plot subsets without rerunning the notebook.

### Web Viewer

Open:

```text
web/bulk_rnaseq_viewer.html
```

Load these files from `RESULTS_DIR / "07_bulk_rnaseq"`:

- required expression file: `07_bulk_rnaseq_candidate_gene_deseq_normalized_counts_by_sample_d0_d2_d8.tsv`
- optional stats file: `07_bulk_rnaseq_candidate_gene_d2_d8_deseq.tsv`

The viewer supports:

- DP or DN selection
- time-course bar plots
- scaled dot plots
- one-gene selection
- custom gene lists through the `Genes to plot` text box
- padj stars when the optional stats file is loaded

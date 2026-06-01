# Kallies / Von Hoesslin RNA-seq Analysis and Website Pipeline

This document explains the full analysis pipeline behind the hosted RNA-seq gene
explorer. It covers both the Kallies single-cell RNA-seq analysis and the Von
Hoesslin bulk RNA-seq comparison, including why each major step was done.

## Overall Goal

The goal is to compare gene expression during memory and rechallenge/reinfection
states using two complementary datasets:

1. Kallies Sci Immuno 2025 single-cell RNA-seq
2. Von Hoesslin Sci Immuno 2022 bulk RNA-seq

The final website lets you enter any gene and view:

- bulk RNA-seq expression across `d0 NR`, `d2 OVA`, and `d8 OVA`
- single-cell pseudobulk expression across memory and rechallenge groups

The website is dynamic. It does not store plot images for every gene. It stores
static normalized expression tables, and JavaScript draws SVG plots in the
browser when a gene is selected.

## Project Layout

Repository:

```text
/home/dk5299/Projects_31926/RNA-seq/Kallies_Sci_Immuno_2025
```

Large data and generated analysis files:

```text
/mnt/sata4/Dhruv_2025/Projects_31926/RNA_Seq/Kallies_Sci_Immuno_2025
```

Path constants are defined in:

```text
src/paths.py
```

Important path constants:

```text
DATA_DIR
INTERMEDIATE_DIR
RESULTS_DIR
```

## Single-cell RNA-seq Pipeline

The single-cell pipeline follows the numbered notebooks in this repository.

```text
notebooks/00_Load_Data.ipynb
notebooks/01_QC_Filtering.ipynb
notebooks/01b_hashtag_demultiplexing.R
notebooks/02_Normalization.ipynb
notebooks/03_clustering.ipynb
notebooks/04_Cell_Annotation.ipynb
notebooks/05_Memory_vs_Rechallenge.ipynb
notebooks/06_Marker_Gene_Dotplots_LogNormalized.ipynb
```

### 1. Load Raw Data

Notebook:

```text
notebooks/00_Load_Data.ipynb
```

Input:

```text
DATA_DIR / matrix.mtx
DATA_DIR / features.tsv
DATA_DIR / barcodes.tsv
```

Output:

```text
INTERMEDIATE_DIR / 00_loaded_raw.h5ad
```

Purpose:

- read the raw 10x-style count matrix
- store cells, genes, and raw molecular counts in an AnnData object
- create a stable starting point for all downstream single-cell analysis

Why this step matters:

- raw counts are the only valid starting point for QC, normalization, and
  pseudobulk
- saving the loaded object avoids repeatedly parsing large raw files

### 2. Quality Control and Filtering

Notebook:

```text
notebooks/01_QC_Filtering.ipynb
```

Output:

```text
INTERMEDIATE_DIR / 01_filtered.h5ad
```

Typical QC metrics include:

- number of detected genes per cell
- total counts per cell
- mitochondrial count percentage
- genes detected across cells

Purpose:

- remove poor-quality cells
- remove empty or low-complexity droplets
- remove cells with unusually high mitochondrial content
- remove genes with insufficient detection

Why this step matters:

- low-quality cells can create artificial clusters
- dying or stressed cells can dominate downstream variation
- filtering improves clustering, marker detection, and condition comparisons

### 3. Hashtag Demultiplexing

Notebook/script:

```text
notebooks/01b_hashtag_demultiplexing.R
```

Outputs:

```text
INTERMEDIATE_DIR / 01b_demultiplexed_all.h5ad
INTERMEDIATE_DIR / 01b_demultiplexed_singlets.h5ad
INTERMEDIATE_DIR / 01b_sIEL.h5ad
INTERMEDIATE_DIR / 01b_sLPL.h5ad
INTERMEDIATE_DIR / 01b_rechallenge.h5ad
RESULTS_DIR / 01b_hashtag_demultiplexing_metadata.csv
RESULTS_DIR / 01b_hashtag_demultiplexing_summary.csv
```

Hashtag mapping:

```text
HT1 -> sIEL
HT2 -> sLPL
Negative or HT3 -> rechallenge
```

Downstream condition labels:

```text
sIEL -> IEL_memory
sLPL -> LPL_memory
rechallenge -> rechallenge
```

Purpose:

- assign cells to biological groups using hashtag counts
- separate memory intestinal epithelial lymphocyte, memory lamina propria, and
  rechallenge groups
- remove or flag ambiguous cells where appropriate

Why this step matters:

- condition labels are required for memory vs rechallenge comparisons
- incorrect demultiplexing would mix biological groups and obscure true signal
- saving explicit labels prevents downstream notebooks from re-inferring groups

### 4. Single-cell Normalization

Notebook:

```text
notebooks/02_Normalization.ipynb
```

Outputs:

```text
INTERMEDIATE_DIR / 02_normalized_hvg_pca.h5ad
INTERMEDIATE_DIR / 02_normalized_hvg_pca_sIEL.h5ad
INTERMEDIATE_DIR / 02_normalized_hvg_pca_sLPL.h5ad
INTERMEDIATE_DIR / 02_normalized_hvg_pca_rechallenge.h5ad
RESULTS_DIR / 02_normalization_summary.csv
```

Key convention:

```text
layers["counts"] = raw counts
```

Normalization steps:

1. preserve raw counts in `layers["counts"]`
2. normalize total counts per cell
3. log-transform normalized expression
4. identify highly variable genes
5. scale data for PCA
6. compute PCA

Purpose:

- correct for differences in sequencing depth between cells
- make cells more comparable
- reduce the influence of very highly expressed genes
- select informative genes for dimensionality reduction

Why normalize single-cell data:

- cells have different library sizes, meaning one cell may have many more
  captured molecules than another
- without normalization, clustering can reflect sequencing depth instead of
  biology
- log transformation compresses large expression differences and makes marker
  visualization more interpretable

Why preserve raw counts:

- raw counts are needed for pseudobulk
- raw counts are needed if normalization needs to be rebuilt later
- scaled PCA values should not be used as expression values for gene plots

### 5. Clustering and UMAP

Notebook:

```text
notebooks/03_clustering.ipynb
```

Outputs include:

```text
INTERMEDIATE_DIR / 03_clustered_umap.h5ad
INTERMEDIATE_DIR / 03_clustered_umap_sIEL.h5ad
INTERMEDIATE_DIR / 03_clustered_umap_sLPL.h5ad
INTERMEDIATE_DIR / 03_clustered_umap_rechallenge.h5ad
RESULTS_DIR / 03_clustered_umap_cluster_summary.csv
```

Purpose:

- group transcriptionally similar cells
- visualize major cell states using UMAP
- create cluster labels for marker and annotation analysis

Why this step matters:

- clustering helps identify whether the dataset contains distinct cell states
- UMAP gives an interpretable view of condition and cluster structure
- cluster labels can be used to summarize marker expression and cell abundance

### 6. Cell Annotation

Notebook:

```text
notebooks/04_Cell_Annotation.ipynb
```

Outputs:

```text
INTERMEDIATE_DIR / 04_cell_annotation.h5ad
INTERMEDIATE_DIR / 04_cd8_subset.h5ad
RESULTS_DIR / 04_cluster_markers_leiden_0_5.xlsx
RESULTS_DIR / 04_cd8_cluster_scores.csv
```

Purpose:

- identify cell populations using marker genes
- annotate clusters manually
- create the CD8 T-cell subset used for memory and rechallenge analysis

Why this step matters:

- downstream biological interpretation depends on analyzing the relevant cell
  type
- the website pseudobulk uses the CD8 subset, not all cells
- restricting to the relevant population reduces confounding from unrelated
  cell types

### 7. Memory vs Rechallenge Analysis

Notebook:

```text
notebooks/05_Memory_vs_Rechallenge.ipynb
```

Outputs:

```text
INTERMEDIATE_DIR / 05_memory_vs_rechallenge.h5ad
RESULTS_DIR / 05_DE_rechallenge_vs_IEL_memory_all_CD8.xlsx
RESULTS_DIR / 05_DE_rechallenge_vs_IEL_memory_top50_expression_zscores.csv
RESULTS_DIR / 05_condition_cluster_counts.csv
```

Purpose:

- compare CD8 cells from memory and rechallenge conditions
- identify genes and patterns associated with rechallenge
- summarize condition composition across clusters

Why this step matters:

- this is the main biological contrast in the single-cell dataset
- it links cell-state changes to rechallenge/reinfection biology

### 8. Marker Gene Dot Plots

Notebook:

```text
notebooks/06_Marker_Gene_Dotplots_LogNormalized.ipynb
```

Marker gene input:

```text
data/marker_gene_panels.csv
```

Outputs:

```text
RESULTS_DIR / 06_marker_gene_dotplots_lognormalized
```

Purpose:

- visualize selected marker genes by condition, cluster, or condition-cluster
- summarize mean expression and percent detected

Why this step matters:

- dot plots are useful for marker interpretation
- they show both expression intensity and prevalence across cells
- they use log-normalized expression rebuilt from raw counts, not scaled PCA
  values

## Single-cell Pseudobulk for Website

Script:

```text
scripts/build_website_data.py
```

Input:

```text
INTERMEDIATE_DIR / 04_cd8_subset.h5ad
```

The script uses:

```text
layers["counts"]
```

Groups:

```text
IEL_memory
LPL_memory
rechallenge
```

Output:

```text
data/sc_pseudobulk_all_gene_counts.tsv
docs/data/sc_pseudobulk_all_gene_counts.tsv
```

### What Pseudobulk Means

Pseudobulk means raw single-cell counts are summed across cells within each
condition.

For each gene:

```text
pseudobulk raw count = sum of raw counts from all cells in the condition
```

Then the condition-level pseudobulk count matrix is median-ratio normalized.

### Why Pseudobulk Was Used

Single-cell plots based on average log-normalized values can be useful, but they
are not the same as bulk-style expression.

Pseudobulk is used because:

- it starts from raw counts
- it reduces cell-level noise
- it makes the single-cell data more comparable to bulk RNA-seq
- it avoids treating thousands of cells as independent biological replicates

### Why Normalize Pseudobulk Counts

The condition groups have different numbers of cells and different total counts.

Raw summed counts would be larger simply because:

- a condition has more cells
- a condition has higher total captured molecules

Median-ratio normalization adjusts for these differences, so expression is more
comparable across pseudobulk conditions.

### Pseudobulk Interpretation Caveat

The current single-cell metadata does not expose true biological replicate IDs.

Therefore, the website shows condition-level pseudobulk expression. It is useful
for visualization, but it should not be interpreted as replicate-level
pseudobulk differential expression.

There are no single-cell pseudobulk SEM error bars because each condition is
represented as one summed profile.

## Bulk RNA-seq Pipeline

Bulk RNA-seq is used as an external comparison dataset from the Von Hoesslin Sci
Immuno 2022 project.

Source project:

```text
/home/dk5299/Projects_31926/RNA-seq/Von_Hoesslin_Sci_Immuno_2022
```

Raw count matrix:

```text
count matrix/GSE208097_count_matrix_compiled
```

Sample metadata:

```text
count matrix/GSE208097_sample_key.tsv
```

Original DESeq2 outputs:

```text
results/tables/02_deseq
```

Companion notebook in this repository:

```text
notebooks/07_Bulk_RNAseq.ipynb
```

### Bulk Samples Used by Website

The website uses OVA rechallenge and baseline samples:

```text
DP_nr_d0_*
DP_ova_d2_*
DP_ova_d8_*
DN_nr_d0_*
DN_ova_d2_*
DN_ova_d8_*
```

These are displayed as:

```text
d0 NR
d2 OVA
d8 OVA
```

for:

```text
DP
DN
```

### Bulk Normalization

The website builder estimates DESeq2-style median-ratio size factors across the
selected bulk samples.

Output:

```text
data/bulk_all_gene_summary.tsv
docs/data/bulk_all_gene_summary.tsv
```

Important columns:

```text
gene_id
gene_name
population
time_condition
mean_normalized_count
sem_normalized_count
n_samples
plot_expression
```

### Why Bulk Counts Were Normalized

Bulk RNA-seq samples have different sequencing depths and RNA composition.

Raw counts are not directly comparable between samples because a higher count
can reflect:

- more sequencing reads
- larger library size
- composition differences
- true biological expression

Median-ratio normalization adjusts sample-level size factors so expression can
be compared across samples more fairly.

### Why Summarize by Mean and SEM

Bulk samples have biological replicates.

For each gene, population, and time condition, the website reports:

- mean normalized count
- SEM across samples
- number of samples

This gives both the central expression trend and a simple measure of replicate
variation.

## Website Data Builder

Script:

```text
scripts/build_website_data.py
```

Run from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python scripts/build_website_data.py
```

The thread variables are used to keep computation conservative and avoid using
many CPU cores.

The script builds:

```text
docs/data/bulk_all_gene_summary.tsv
docs/data/sc_pseudobulk_all_gene_counts.tsv
```

The live GitHub Pages root also contains copies at:

```text
data/bulk_all_gene_summary.tsv
data/sc_pseudobulk_all_gene_counts.tsv
```

## Website Implementation

Main viewer:

```text
index.html
```

Backup/copy:

```text
docs/index.html
```

The browser loads:

```text
data/bulk_all_gene_summary.tsv
data/sc_pseudobulk_all_gene_counts.tsv
```

Then JavaScript:

1. parses the TSV files
2. builds the searchable gene list
3. filters rows for the selected gene or gene list
4. draws SVG bar plots dynamically
5. optionally expands the current plot in a larger modal

No server backend is used.

No pre-rendered gene plot images are stored.

## What the Plots Show

### Bulk RNA-seq Plot

Shows:

```text
d0 NR
d2 OVA
d8 OVA
```

for either:

```text
DP
DN
```

Bars:

```text
mean normalized count
```

Error bars:

```text
SEM across bulk biological samples
```

### Single-cell Pseudobulk Plot

Shows:

```text
IEL_memory
LPL_memory
rechallenge
```

Bars:

```text
condition-level pseudobulk normalized count
```

No error bars are shown because the current pseudobulk table has one summed
profile per condition.

## GitHub Pages Hosting

Current URL:

```text
https://dhruv5199ucsd.github.io/Rechallenge_Bulk_RNAseq/
```

Recommended GitHub Pages settings:

```text
Source: Deploy from a branch
Branch: main
Folder: /root
```

The root `index.html` is the actual interactive viewer.

If the site appears stale after a push, use a cache-busting URL:

```text
https://dhruv5199ucsd.github.io/Rechallenge_Bulk_RNAseq/?v=2
```

## How to Update the Website

If the analysis data change:

1. rerun the appropriate single-cell or bulk notebooks
2. rerun:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python scripts/build_website_data.py
```

3. copy updated files from `docs/data/` to `data/` if needed
4. commit and push

Example:

```bash
git add index.html docs/index.html data docs/data scripts/build_website_data.py WEBSITE_PIPELINE.md
git commit -m "Update RNA-seq gene explorer"
git push origin main
```

## Key Interpretation Points

- Bulk RNA-seq values are normalized counts summarized across biological
  replicates.
- Single-cell website values are pseudobulk normalized counts from summed raw
  single-cell counts.
- Single-cell pseudobulk is condition-level because replicate IDs are not
  available in the current metadata.
- Normalization is required because raw counts are affected by sequencing depth,
  cell number, and library composition.
- The website is for interactive visualization, not for running new differential
  expression tests in the browser.

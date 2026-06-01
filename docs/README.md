# RNA-seq Gene Explorer

This directory is the GitHub Pages site.

## Regenerate Website Data

Run from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python scripts/build_website_data.py
```

The site loads:

- `data/bulk_all_gene_summary.tsv`
- `data/sc_pseudobulk_all_gene_counts.tsv`

The single-cell table is pseudobulked from raw counts in `layers["counts"]`,
summed by condition, then median-ratio normalized.

## GitHub Pages

1. Push this repository to GitHub.
2. In GitHub, open `Settings` -> `Pages`.
3. Under `Build and deployment`, choose `Deploy from a branch`.
4. Choose the branch you pushed, then choose `/docs` as the folder.
5. Save. GitHub will publish `docs/index.html`.

The page is fully static and does not need a server backend.

#!/usr/bin/env python3
"""Build static website data for bulk RNA-seq and scRNA-seq pseudobulk plots."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.paths import INTERMEDIATE_DIR, REPO_DIR


VON_HOESSLIN_DIR = Path("/home/dk5299/Projects_31926/RNA-seq/Von_Hoesslin_Sci_Immuno_2022")
COUNT_MATRIX = VON_HOESSLIN_DIR / "count matrix" / "GSE208097_count_matrix_compiled"
SAMPLE_KEY = VON_HOESSLIN_DIR / "count matrix" / "GSE208097_sample_key.tsv"
SC_OBJECT = INTERMEDIATE_DIR / "04_cd8_subset.h5ad"

BULK_KEEP_GROUPS = {
    ("nr", "d0"): "d0 NR",
    ("ova", "d2"): "d2 OVA",
    ("ova", "d8"): "d8 OVA",
}

SC_CONDITION_LABELS = {
    "sIEL": "IEL_memory",
    "sLPL": "LPL_memory",
    "rechallenge": "rechallenge",
}


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    """Estimate DESeq2-style median-ratio size factors from genes with all counts > 0."""
    positive = counts.gt(0).all(axis=1)
    if not positive.any():
        positive = counts.sum(axis=1).gt(0)
    valid = counts.loc[positive].astype(float)
    log_geo_means = np.log(valid.replace(0, np.nan)).mean(axis=1, skipna=True)
    geo_means = np.exp(log_geo_means)
    ratios = valid.div(geo_means, axis=0).replace([np.inf, -np.inf], np.nan)
    size_factors = ratios.median(axis=0, skipna=True)
    size_factors = size_factors.replace(0, np.nan)
    if size_factors.isna().any():
        fallback = counts.sum(axis=0) / np.median(counts.sum(axis=0))
        size_factors = size_factors.fillna(fallback)
    return size_factors


def sem(values: pd.Series) -> float:
    if len(values) < 2:
        return 0.0
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def parse_bulk_sample(sample_id: str) -> dict[str, str] | None:
    parts = sample_id.split("_")
    if len(parts) != 4:
        return None
    population, group_code, timepoint, replicate = parts
    condition = BULK_KEEP_GROUPS.get((group_code, timepoint))
    if condition is None:
        return None
    return {
        "sample_id": sample_id,
        "population": population,
        "group_code": group_code,
        "timepoint": timepoint,
        "replicate": replicate,
        "time_condition": condition,
    }


def build_bulk(output_dir: Path) -> None:
    counts = pd.read_csv(COUNT_MATRIX, sep="\t")
    gene_info = counts[["gene_id", "gene_name"]].copy()
    sample_cols = [col for col in counts.columns if col not in {"gene_id", "gene_name"}]
    sample_meta = pd.DataFrame(
        row for col in sample_cols if (row := parse_bulk_sample(col)) is not None
    )
    selected_cols = sample_meta["sample_id"].tolist()
    selected_counts = counts[selected_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    size_factors = median_ratio_size_factors(selected_counts)
    normalized = selected_counts.div(size_factors, axis=1)
    normalized.insert(0, "gene_name", gene_info["gene_name"])
    normalized.insert(0, "gene_id", gene_info["gene_id"])

    long = normalized.melt(
        id_vars=["gene_id", "gene_name"],
        var_name="sample_id",
        value_name="normalized_count",
    ).merge(sample_meta, on="sample_id", how="left")
    long["plot_expression"] = long["normalized_count"] + 1

    summary = (
        long.groupby(["gene_id", "gene_name", "population", "time_condition"], observed=True)
        .agg(
            mean_normalized_count=("normalized_count", "mean"),
            sem_normalized_count=("normalized_count", sem),
            n_samples=("sample_id", "nunique"),
        )
        .reset_index()
    )
    summary["plot_expression"] = summary["mean_normalized_count"] + 1

    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "bulk_all_gene_summary.tsv", sep="\t", index=False)
    size_factors.rename("size_factor").to_frame().to_csv(
        output_dir / "bulk_all_gene_size_factors.tsv", sep="\t"
    )


def condition_series(adata: ad.AnnData) -> pd.Series:
    if "condition" in adata.obs:
        return adata.obs["condition"].astype(str)
    if "sample_group" not in adata.obs:
        raise KeyError("Expected either obs['condition'] or obs['sample_group'] in AnnData.")
    return adata.obs["sample_group"].astype(str).map(SC_CONDITION_LABELS).fillna(
        adata.obs["sample_group"].astype(str)
    )


def pseudobulk_counts(adata: ad.AnnData, conditions: pd.Series) -> pd.DataFrame:
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    matrix = matrix.tocsr()

    rows = []
    for condition in ["IEL_memory", "LPL_memory", "rechallenge"]:
        mask = conditions.eq(condition).to_numpy()
        if not mask.any():
            continue
        summed = np.asarray(matrix[mask].sum(axis=0)).ravel()
        rows.append(pd.Series(summed, name=condition))
    if not rows:
        raise ValueError("No cells matched the expected pseudobulk conditions.")
    counts = pd.DataFrame(rows).T
    counts.insert(0, "gene_name", adata.var_names.to_list())
    if "gene_ids" in adata.var:
        counts.insert(0, "gene_id", adata.var["gene_ids"].astype(str).to_list())
    else:
        counts.insert(0, "gene_id", adata.var_names.to_list())
    return counts


def build_single_cell_pseudobulk(output_dir: Path) -> None:
    adata = ad.read_h5ad(SC_OBJECT)
    conditions = condition_series(adata)
    cell_counts = conditions.value_counts().rename_axis("condition").reset_index(name="n_cells")
    raw_counts = pseudobulk_counts(adata, conditions)

    condition_cols = [col for col in raw_counts.columns if col not in {"gene_id", "gene_name"}]
    size_factors = median_ratio_size_factors(raw_counts[condition_cols])
    normalized = raw_counts[condition_cols].div(size_factors, axis=1)
    normalized.insert(0, "gene_name", raw_counts["gene_name"])
    normalized.insert(0, "gene_id", raw_counts["gene_id"])

    long = normalized.melt(
        id_vars=["gene_id", "gene_name"],
        var_name="condition",
        value_name="normalized_count",
    )
    raw_long = raw_counts.melt(
        id_vars=["gene_id", "gene_name"],
        var_name="condition",
        value_name="raw_count",
    )
    long = long.merge(raw_long, on=["gene_id", "gene_name", "condition"], how="left")
    long = long.merge(cell_counts, on="condition", how="left")
    long["plot_expression"] = long["normalized_count"] + 1

    output_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(output_dir / "sc_pseudobulk_all_gene_counts.tsv", sep="\t", index=False)
    cell_counts.to_csv(output_dir / "sc_pseudobulk_cell_counts.tsv", sep="\t", index=False)
    size_factors.rename("size_factor").to_frame().to_csv(
        output_dir / "sc_pseudobulk_size_factors.tsv", sep="\t"
    )


def write_manifest(output_dir: Path) -> None:
    manifest = pd.DataFrame(
        [
            {
                "dataset": "bulk",
                "file": "bulk_all_gene_summary.tsv",
                "description": "All-gene Von Hoesslin bulk RNA-seq DESeq2-style normalized summary for d0 NR, d2 OVA, and d8 OVA.",
            },
            {
                "dataset": "single_cell_pseudobulk",
                "file": "sc_pseudobulk_all_gene_counts.tsv",
                "description": "All-gene Kallies CD8 scRNA-seq pseudobulk counts summed by condition and median-ratio normalized.",
            },
        ]
    )
    manifest.to_csv(output_dir / "manifest.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=REPO_DIR / "docs" / "data")
    parser.add_argument("--bulk-only", action="store_true")
    parser.add_argument("--single-cell-only", action="store_true")
    args = parser.parse_args()

    if not args.single_cell_only:
        build_bulk(args.output_dir)
    if not args.bulk_only:
        build_single_cell_pseudobulk(args.output_dir)
    write_manifest(args.output_dir)


if __name__ == "__main__":
    main()

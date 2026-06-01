#!/usr/bin/env python
"""Score gene modules on normalized Kallies Sci Immuno 2025 AnnData objects."""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import scanpy as sc

REPO_DIR = Path(__file__).resolve().parents[1]
if str(REPO_DIR) not in sys.path:
    sys.path.append(str(REPO_DIR))

from src.paths import INTERMEDIATE_DIR, RESULTS_DIR  # noqa: E402


DEFAULT_MODULES: dict[str, list[str]] = {
    "t_cell_core": ["Cd3d", "Cd3e", "Trac", "Lck", "Lat"],
    "cd8_core": ["Cd8a", "Cd8b1"],
    "cytotoxic": ["Nkg7", "Gzmb", "Gzma", "Prf1", "Ccl5", "Ctsw"],
    "effector_inflammatory": ["Ifng", "Tnf", "Il2ra", "Cxcr3", "Irf1"],
    "memory_like": ["Tcf7", "Lef1", "Il7r", "Sell", "Ccr7", "Bcl2"],
    "tissue_residency": ["Itgae", "Cd69", "Cxcr6", "Zfp683", "Rgs1"],
    "activation_exhaustion": ["Pdcd1", "Lag3", "Tigit", "Havcr2", "Ctla4", "Tox"],
    "proliferation": ["Mki67", "Top2a", "Stmn1", "Tubb5", "Hmgb2"],
    "interferon_response": ["Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1", "Bst2"],
    "trafficking": ["S1pr1", "Klf2", "Cxcr3", "Ccr5", "Ccr9", "Itgb7"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score gene modules using Scanpy score_genes. By default the script "
            "reads 02_normalized_hvg_pca.h5ad, rebuilds log-normalized expression "
            "from layers['counts'], and writes score tables to RESULTS_DIR."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=INTERMEDIATE_DIR / "02_normalized_hvg_pca.h5ad",
        help="Input .h5ad file. Defaults to the full normalized object.",
    )
    parser.add_argument(
        "--modules",
        type=Path,
        help=(
            "Optional CSV with columns 'module' and 'gene'. If omitted, built-in "
            "CD8/T-cell modules are used."
        ),
    )
    parser.add_argument(
        "--output-prefix",
        default="06_module_scores",
        help="Prefix for output files in RESULTS_DIR.",
    )
    parser.add_argument(
        "--groupby",
        nargs="*",
        default=["sample_group", "leiden_0_5", "cell_type_manual"],
        help="obs columns to summarize when present.",
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Target sum for count normalization before log1p.",
    )
    parser.add_argument(
        "--ctrl-size",
        type=int,
        default=50,
        help="Control gene set size passed to sc.tl.score_genes.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random seed passed to sc.tl.score_genes.",
    )
    parser.add_argument(
        "--write-h5ad",
        action="store_true",
        help="Also save a scored AnnData object under INTERMEDIATE_DIR.",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        help="Optional output .h5ad path. Used only with --write-h5ad.",
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z]+", "_", value.strip()).strip("_").lower()
    if not cleaned:
        raise ValueError(f"Invalid empty module name from {value!r}")
    return cleaned


def load_modules(path: Path | None) -> dict[str, list[str]]:
    if path is None:
        return DEFAULT_MODULES

    modules_df = pd.read_csv(path)
    required = {"module", "gene"}
    missing = required.difference(modules_df.columns)
    if missing:
        raise ValueError(f"{path} is missing required column(s): {sorted(missing)}")

    modules: dict[str, list[str]] = {}
    for module, genes_df in modules_df.dropna(subset=["module", "gene"]).groupby("module"):
        genes = [str(gene).strip() for gene in genes_df["gene"] if str(gene).strip()]
        modules[str(module)] = list(dict.fromkeys(genes))
    return modules


def normalize_from_counts(adata: sc.AnnData, target_sum: float) -> sc.AnnData:
    scored = adata.copy()
    if "counts" in scored.layers:
        scored.X = scored.layers["counts"].copy()
    elif "raw_counts" in scored.layers:
        scored.X = scored.layers["raw_counts"].copy()
    else:
        print(
            "WARNING: no counts/raw_counts layer found; scoring from existing .X. "
            "This is only appropriate if .X is log-normalized expression.",
            file=sys.stderr,
        )

    scored.uns.pop("log1p", None)
    sc.pp.normalize_total(scored, target_sum=target_sum)
    sc.pp.log1p(scored)
    return scored


def present_genes(genes: Iterable[str], var_names: pd.Index) -> list[str]:
    return [gene for gene in genes if gene in var_names]


def score_modules(
    adata: sc.AnnData,
    modules: dict[str, list[str]],
    ctrl_size: int,
    random_state: int,
) -> tuple[list[str], pd.DataFrame]:
    score_columns: list[str] = []
    rows: list[dict[str, object]] = []

    for module, genes in modules.items():
        column = f"module_score_{safe_name(module)}"
        genes_present = present_genes(genes, adata.var_names)
        genes_missing = sorted(set(genes).difference(genes_present))

        rows.append(
            {
                "module": module,
                "score_column": column,
                "n_genes_requested": len(genes),
                "n_genes_present": len(genes_present),
                "n_genes_missing": len(genes_missing),
                "genes_present": ";".join(genes_present),
                "genes_missing": ";".join(genes_missing),
            }
        )

        if len(genes_present) < 2:
            print(
                f"Skipping {module}: fewer than 2 genes present "
                f"({len(genes_present)}/{len(genes)}).",
                file=sys.stderr,
            )
            continue

        sc.tl.score_genes(
            adata,
            gene_list=genes_present,
            score_name=column,
            ctrl_size=ctrl_size,
            random_state=random_state,
            use_raw=False,
        )
        score_columns.append(column)

    return score_columns, pd.DataFrame(rows)


def write_outputs(
    adata: sc.AnnData,
    score_columns: list[str],
    module_summary: pd.DataFrame,
    output_prefix: str,
    groupby_columns: list[str],
    write_h5ad: bool,
    output_h5ad: Path | None,
) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    module_file = RESULTS_DIR / f"{output_prefix}_module_gene_coverage.csv"
    module_summary.to_csv(module_file, index=False)
    print("Saved:", module_file)

    obs_columns = [column for column in ["sample_group", "hash_id", "demux_call"] if column in adata.obs]
    cell_scores = adata.obs[obs_columns + score_columns].copy()
    cell_file = RESULTS_DIR / f"{output_prefix}_cell_scores.csv"
    cell_scores.to_csv(cell_file, index=True, index_label="barcode")
    print("Saved:", cell_file)

    for groupby in groupby_columns:
        if groupby not in adata.obs:
            continue
        summary = (
            adata.obs.groupby(groupby, observed=True)[score_columns]
            .agg(["mean", "median", "std"])
            .reset_index()
        )
        summary.columns = [
            "_".join(str(part) for part in col if part) if isinstance(col, tuple) else col
            for col in summary.columns
        ]
        summary["n_cells"] = summary[groupby].map(adata.obs[groupby].value_counts())
        summary_file = RESULTS_DIR / f"{output_prefix}_by_{safe_name(groupby)}.csv"
        summary.to_csv(summary_file, index=False)
        print("Saved:", summary_file)

    if write_h5ad:
        if output_h5ad is None:
            output_h5ad = INTERMEDIATE_DIR / f"{output_prefix}.h5ad"
        output_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_h5ad)
        print("Saved:", output_h5ad)


def main() -> None:
    args = parse_args()
    modules = load_modules(args.modules)
    print(f"Loaded {len(modules)} module(s).")

    adata = sc.read_h5ad(args.input)
    print("Input:", args.input)
    print(adata)

    scored = normalize_from_counts(adata, target_sum=args.target_sum)
    score_columns, module_summary = score_modules(
        scored,
        modules=modules,
        ctrl_size=args.ctrl_size,
        random_state=args.random_state,
    )

    if not score_columns:
        raise RuntimeError("No module scores were computed; check module gene symbols.")

    write_outputs(
        scored,
        score_columns=score_columns,
        module_summary=module_summary,
        output_prefix=args.output_prefix,
        groupby_columns=args.groupby,
        write_h5ad=args.write_h5ad,
        output_h5ad=args.output_h5ad,
    )


if __name__ == "__main__":
    main()

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path("/home/dk5299/Projects_31926/RNA-seq/Kallies_Sci_Immuno_2025")
FUNG_DIR = Path("/home/dk5299/Projects_31926/RNA-seq/Fung_et_al/GSE185541")
OUT_DIRS = [REPO_DIR / "docs" / "data", REPO_DIR / "data"]

COUNTS_PATH = FUNG_DIR / "data" / "GSE185541_merged_gene_counts.txt.gz"
COL_DATA_PATH = FUNG_DIR / "results" / "tables" / "00_loading" / "GSE185541_col_data_for_deseq2.tsv"
SIZE_FACTORS_PATH = FUNG_DIR / "results" / "tables" / "00_loading" / "GSE185541_size_factors.tsv"

OUT_FILE_NAME = "fung_gse185541_lp_trm_gene_summary.tsv"
DATASET_NAME = "fung_gse185541_lp_trm"
DATASET_DESCRIPTION = (
    "Fung et al. GSE185541 bulk RNA-seq DESeq2-normalized mean expression and SEM "
    "for CD103+Tom+ and CD103-Tom- LP Trm cells at baseline and after infection."
)

SUBSET_LABELS = {
    "cd103_pos": "CD103+Tom+ LP Trm",
    "cd103_neg": "CD103-Tom- LP Trm",
}
INFECTION_LABELS = {
    "baseline": "Baseline",
    "infected": "Infected",
}


def sem(values: pd.Series) -> float:
    if len(values) < 2:
        return 0.0
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def build_summary() -> pd.DataFrame:
    counts = pd.read_csv(COUNTS_PATH, compression="gzip")
    col_data = pd.read_csv(COL_DATA_PATH, sep="\t")
    size_factors = pd.read_csv(SIZE_FACTORS_PATH, sep="\t").set_index("sample")["size_factor"]

    samples = col_data["sample"].tolist()
    missing_samples = [sample for sample in samples if sample not in counts.columns]
    if missing_samples:
        raise ValueError(f"Count matrix is missing samples: {missing_samples}")

    size_factors = size_factors.reindex(samples)
    if size_factors.isna().any():
        missing_factors = size_factors[size_factors.isna()].index.tolist()
        raise ValueError(f"Missing size factors for samples: {missing_factors}")

    sample_meta = col_data.loc[:, ["sample", "subset", "infection"]].copy()
    sample_meta["population"] = sample_meta["subset"].map(SUBSET_LABELS)
    sample_meta["time_condition"] = sample_meta["infection"].map(INFECTION_LABELS)

    if sample_meta["population"].isna().any() or sample_meta["time_condition"].isna().any():
        raise ValueError("Encountered unexpected Fung subset or infection labels.")

    normalized = (
        counts.loc[:, samples]
        .apply(pd.to_numeric, errors="coerce")
        .fillna(0.0)
        .div(size_factors, axis=1)
    )
    normalized.insert(0, "gene_name", counts["gene_name"].fillna(""))
    normalized.insert(0, "gene_id", counts["Geneid"])

    long = normalized.melt(
        id_vars=["gene_id", "gene_name"],
        var_name="sample",
        value_name="normalized_count",
    ).merge(sample_meta, on="sample", how="left")

    summary = (
        long.groupby(["gene_id", "gene_name", "population", "time_condition"], observed=True)
        .agg(
            mean_normalized_count=("normalized_count", "mean"),
            sem_normalized_count=("normalized_count", sem),
            n_samples=("sample", "nunique"),
        )
        .reset_index()
    )
    summary["plot_expression"] = summary["mean_normalized_count"] + 1

    population_order = list(SUBSET_LABELS.values())
    time_order = list(INFECTION_LABELS.values())
    summary["population"] = pd.Categorical(summary["population"], population_order, ordered=True)
    summary["time_condition"] = pd.Categorical(summary["time_condition"], time_order, ordered=True)
    summary = summary.sort_values(["gene_name", "gene_id", "population", "time_condition"]).reset_index(drop=True)
    summary["population"] = summary["population"].astype(str)
    summary["time_condition"] = summary["time_condition"].astype(str)
    return summary


def update_manifest(out_dir: Path) -> None:
    manifest_file = out_dir / "manifest.tsv"
    manifest = pd.read_csv(manifest_file, sep="\t")
    new_row = {
        "dataset": DATASET_NAME,
        "file": OUT_FILE_NAME,
        "description": DATASET_DESCRIPTION,
    }
    manifest = manifest.loc[manifest["dataset"].ne(DATASET_NAME)]
    manifest = pd.concat([manifest, pd.DataFrame([new_row])], ignore_index=True)
    manifest.to_csv(manifest_file, sep="\t", index=False)


def main() -> None:
    summary = build_summary()
    for out_dir in OUT_DIRS:
        out_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out_dir / OUT_FILE_NAME, sep="\t", index=False)
        update_manifest(out_dir)

    print(f"Saved {OUT_FILE_NAME} to {len(OUT_DIRS)} data directories")
    print(f"Rows: {len(summary)}")
    print(f"Genes: {summary[['gene_id', 'gene_name']].drop_duplicates().shape[0]}")
    print("Populations:", ", ".join(SUBSET_LABELS.values()))
    print("Conditions:", ", ".join(INFECTION_LABELS.values()))


if __name__ == "__main__":
    main()

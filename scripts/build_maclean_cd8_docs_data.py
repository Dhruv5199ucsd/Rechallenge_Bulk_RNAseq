from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


REPO_DIR = Path("/home/dk5299/Projects_31926/RNA-seq/Kallies_Sci_Immuno_2025")
MACLEAN_RESULTS = Path(
    "/mnt/sata4/Dhruv_2025/Projects_31926/RNA_Seq/MacLean_Immunity_2022_GSE194058/results/05_cd8_subset_pseudobulk"
)
OUT_DIR = REPO_DIR / "docs" / "data"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(
        MACLEAN_RESULTS / "05_cd8a_cd8b1_positive_pseudobulk_counts_by_sample.tsv",
        sep="\t",
        index_col=0,
    )
    metadata = pd.read_csv(MACLEAN_RESULTS / "05_cd8a_cd8b1_positive_metadata.csv", index_col=0)
    metadata = metadata.loc[metadata["condition"].isin(["resting", "rechallenged"])]
    n_cells = metadata.groupby("sample_id", observed=True).size()

    samples = [sample for sample in counts.columns if sample in set(metadata["sample_id"])]
    counts = counts[samples]
    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes.replace(0, np.nan), axis=1) * 1e6

    rows = []
    for condition in ["resting", "rechallenged"]:
        condition_samples = [sample for sample in samples if sample.startswith(f"{condition}_")]
        condition_cpm = cpm[condition_samples]
        condition_counts = counts[condition_samples]
        mean_cpm = condition_cpm.mean(axis=1)
        sem_cpm = condition_cpm.sem(axis=1)
        raw_sum = condition_counts.sum(axis=1)
        cell_sum = int(n_cells.reindex(condition_samples).fillna(0).sum())
        for gene in counts.index:
            rows.append(
                {
                    "gene_id": gene,
                    "gene_name": gene,
                    "condition": condition,
                    "normalized_count": mean_cpm.loc[gene],
                    "sem_normalized_count": sem_cpm.loc[gene],
                    "raw_count": raw_sum.loc[gene],
                    "n_cells": cell_sum,
                    "n_samples": len(condition_samples),
                    "plot_expression": mean_cpm.loc[gene] + 1,
                }
            )

    out_file = OUT_DIR / "maclean_cd8_resting_rechallenge_pseudobulk_gene_summary.tsv"
    pd.DataFrame(rows).to_csv(out_file, sep="\t", index=False)

    manifest_file = OUT_DIR / "manifest.tsv"
    manifest = pd.read_csv(manifest_file, sep="\t")
    new_row = {
        "dataset": "maclean_cd8_resting_rechallenge",
        "file": out_file.name,
        "description": "MacLean GSE194058 CD8a/Cd8b1-positive scRNA-seq pseudobulk mean CPM by resting and rechallenged condition; CLL excluded.",
    }
    manifest = manifest.loc[manifest["dataset"].ne(new_row["dataset"])]
    manifest = pd.concat([manifest, pd.DataFrame([new_row])], ignore_index=True)
    manifest.to_csv(manifest_file, sep="\t", index=False)

    print(f"Saved {out_file}")
    print(f"Rows: {len(rows)}")
    print(f"Samples: {', '.join(samples)}")


if __name__ == "__main__":
    main()

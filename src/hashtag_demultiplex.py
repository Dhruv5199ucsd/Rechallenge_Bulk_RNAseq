from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

from src.paths import DATA_DIR, INTERMEDIATE_DIR, RESULTS_DIR


HASH_TO_SAMPLE = {
    "HT1": "sIEL",
    "HT2": "sLPL",
    "HT3": "rechallenge",
    "Negative": "rechallenge",
}

HASH_TO_MAN_GROUP = {
    "HT1": "Mem_si_IEL",
    "HT2": "Mem_si_LP",
    "HT3": "Rechallenged_si_IEL",
    "Negative": "Rechallenged_si_IEL",
}


def _read_hashtag_reference() -> pd.DataFrame:
    reference_file = DATA_DIR / "GSE287643_hashtag_reference.csv"
    reference = pd.read_csv(reference_file)
    expected = {"id", "name", "sequence", "feature_type"}
    missing = expected - set(reference.columns)
    if missing:
        raise ValueError(f"Hashtag reference is missing columns: {sorted(missing)}")
    return reference


def _hto_counts(filtered: sc.AnnData, raw: sc.AnnData) -> pd.DataFrame:
    if "feature_types" not in raw.var:
        raise ValueError("Raw AnnData object lacks var['feature_types']; reload 10x with feature types.")

    hto_mask = raw.var["feature_types"].eq("Antibody Capture")
    if hto_mask.sum() == 0:
        raise ValueError("No Antibody Capture features were found in the raw AnnData object.")

    missing_cells = filtered.obs_names.difference(raw.obs_names)
    if len(missing_cells) > 0:
        raise ValueError(f"{len(missing_cells)} filtered cells are absent from the raw AnnData object.")

    return raw[filtered.obs_names, hto_mask].to_df()


def assign_hashes(
    counts: pd.DataFrame,
    min_count: float = 50,
    min_top_to_second_ratio: float = 2,
) -> pd.DataFrame:
    values = counts.to_numpy()
    order = np.argsort(values, axis=1)
    top_pos = order[:, -1]
    second_pos = order[:, -2] if values.shape[1] > 1 else order[:, -1]

    top_count = values[np.arange(values.shape[0]), top_pos]
    second_count = values[np.arange(values.shape[0]), second_pos]
    ratio = top_count / np.maximum(second_count, 1)
    top_hash = counts.columns[top_pos].astype(str)

    hash_id = pd.Series(top_hash, index=counts.index, dtype="object")
    hash_id.loc[top_count < min_count] = "Negative"
    hash_id.loc[(top_count >= min_count) & (ratio < min_top_to_second_ratio)] = "Doublet"

    metadata = pd.DataFrame(
        {
            "hash_id": hash_id,
            "hto_max": top_count,
            "hto_second": second_count,
            "hto_top_to_second_ratio": ratio,
            "hto_total": counts.sum(axis=1).to_numpy(),
        },
        index=counts.index,
    )
    metadata["sample_group"] = metadata["hash_id"].map(HASH_TO_SAMPLE).fillna("doublet")
    metadata["man_2025_group"] = metadata["hash_id"].map(HASH_TO_MAN_GROUP).fillna("Doublet")
    metadata["demux_call"] = np.where(metadata["hash_id"].eq("Doublet"), "Doublet", "Singlet")
    return metadata.join(counts.add_prefix("hto_"))


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    INTERMEDIATE_DIR.mkdir(parents=True, exist_ok=True)

    reference = _read_hashtag_reference()
    raw = sc.read_h5ad(INTERMEDIATE_DIR / "00_loaded_raw.h5ad")
    filtered = sc.read_h5ad(INTERMEDIATE_DIR / "01_filtered.h5ad")
    counts = _hto_counts(filtered, raw)

    missing_reference = sorted(set(counts.columns) - set(reference["id"]))
    if missing_reference:
        raise ValueError(f"Antibody Capture features missing from hashtag reference: {missing_reference}")

    demux = assign_hashes(counts)
    filtered.obs = filtered.obs.join(demux)

    metadata_file = RESULTS_DIR / "01b_hashtag_demultiplexing_metadata.csv"
    summary_file = RESULTS_DIR / "01b_hashtag_demultiplexing_summary.csv"
    all_file = INTERMEDIATE_DIR / "01b_demultiplexed_all.h5ad"
    singlet_file = INTERMEDIATE_DIR / "01b_demultiplexed_singlets.h5ad"

    demux.to_csv(metadata_file)
    summary = (
        filtered.obs.groupby(["hash_id", "sample_group", "man_2025_group", "demux_call"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["demux_call", "sample_group", "hash_id"])
    )
    summary.to_csv(summary_file, index=False)

    filtered.write_h5ad(all_file)
    singlets = filtered[~filtered.obs["hash_id"].eq("Doublet")].copy()
    singlets.write_h5ad(singlet_file)

    for sample_group in ["sIEL", "sLPL", "rechallenge"]:
        subset = singlets[singlets.obs["sample_group"].eq(sample_group)].copy()
        subset.write_h5ad(INTERMEDIATE_DIR / f"01b_{sample_group}.h5ad")

    print("Hashtag reference:")
    print(reference[["id", "name", "sequence", "feature_type"]].to_string(index=False))
    print("\nDemultiplexing summary:")
    print(summary.to_string(index=False))
    print(f"\nSaved metadata: {metadata_file}")
    print(f"Saved all cells: {all_file}")
    print(f"Saved singlets: {singlet_file}")


if __name__ == "__main__":
    main()

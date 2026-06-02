#!/usr/bin/env python3
"""Bulk RNA-seq PCA, DESeq2 volcano plots, and preranked GSEA."""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.append(str(Path(__file__).resolve().parents[1]))
from src.paths import RESULTS_DIR


VON_HOESSLIN_DIR = Path("/home/dk5299/Projects_31926/RNA-seq/Von_Hoesslin_Sci_Immuno_2022")
COUNT_MATRIX = VON_HOESSLIN_DIR / "count matrix" / "GSE208097_count_matrix_compiled"
DESEQ_DIR = VON_HOESSLIN_DIR / "results" / "tables" / "02_deseq"
OUT_DIR = RESULTS_DIR / "08_de_gsea" / "bulk"
GSEA_DIR = OUT_DIR / "gsea_prerank"
DEFAULT_GENE_SETS = [
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2023",
    "KEGG_2019_Mouse",
]
CONTRAST_LABELS = {
    "bulk_DP_d2_ova_vs_d0_nr": "DP d2",
    "bulk_DP_d8_ova_vs_d0_nr": "DP d8",
    "bulk_DN_d2_ova_vs_d0_nr": "DN d2",
    "bulk_DN_d8_ova_vs_d0_nr": "DN d8",
}

CONTRASTS = [
    ("DP", "d2", "02_DP_d2_ova_vs_d0_nr_results.tsv"),
    ("DP", "d8", "02_DP_d8_ova_vs_d0_nr_results.tsv"),
    ("DN", "d2", "02_DN_d2_ova_vs_d0_nr_results.tsv"),
    ("DN", "d8", "02_DN_d8_ova_vs_d0_nr_results.tsv"),
]


def parse_sample(sample_id: str) -> dict[str, str] | None:
    parts = sample_id.split("_")
    if len(parts) != 4:
        return None
    population, treatment, timepoint, replicate = parts
    return {
        "sample_id": sample_id,
        "population": population,
        "treatment": treatment,
        "timepoint": timepoint,
        "replicate": replicate,
        "condition": f"{timepoint}_{treatment}",
        "plot_group": f"{population} {timepoint.upper()} {treatment.upper()}",
    }


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    positive = counts.gt(0).all(axis=1)
    if not positive.any():
        positive = counts.sum(axis=1).gt(0)
    valid = counts.loc[positive].astype(float)
    log_geo_means = np.log(valid.replace(0, np.nan)).mean(axis=1, skipna=True)
    geo_means = np.exp(log_geo_means)
    ratios = valid.div(geo_means, axis=0).replace([np.inf, -np.inf], np.nan)
    size_factors = ratios.median(axis=0, skipna=True).replace(0, np.nan)
    if size_factors.isna().any():
        fallback = counts.sum(axis=0) / np.median(counts.sum(axis=0))
        size_factors = size_factors.fillna(fallback)
    return size_factors


def load_bulk_counts() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(COUNT_MATRIX, sep="\t")
    sample_cols = [c for c in counts.columns if c not in {"gene_id", "gene_name"}]
    sample_meta = pd.DataFrame([parse_sample(c) for c in sample_cols]).dropna()
    sample_cols = sample_meta["sample_id"].tolist()
    count_matrix = counts[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    count_matrix.index = counts["gene_name"]
    gene_info = counts[["gene_id", "gene_name"]].copy()
    return gene_info, count_matrix, sample_meta


def write_pca(count_matrix: pd.DataFrame, sample_meta: pd.DataFrame) -> None:
    keep = sample_meta[
        sample_meta["treatment"].isin(["nr", "ova"])
        & sample_meta["timepoint"].isin(["d0", "d2", "d8"])
    ].copy()
    keep = keep[
        ((keep["treatment"] == "nr") & (keep["timepoint"] == "d0"))
        | ((keep["treatment"] == "ova") & (keep["timepoint"].isin(["d2", "d8"])))
    ]
    selected = keep["sample_id"].tolist()
    filtered = count_matrix[selected]
    filtered = filtered.loc[filtered.sum(axis=1) >= 10]
    size_factors = median_ratio_size_factors(filtered)
    log_norm = np.log2(filtered.div(size_factors, axis=1) + 1)

    features = log_norm.var(axis=1).sort_values(ascending=False).head(5000).index
    x = log_norm.loc[features].T
    x_scaled = StandardScaler().fit_transform(x)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(x_scaled)
    pca_df = keep.set_index("sample_id").loc[x.index].copy()
    pca_df.index.name = "sample_id"
    pca_df = pca_df.reset_index()
    pca_df["PC1"] = coords[:, 0]
    pca_df["PC2"] = coords[:, 1]
    pca_df["PC1_percent"] = pca.explained_variance_ratio_[0] * 100
    pca_df["PC2_percent"] = pca.explained_variance_ratio_[1] * 100
    pca_df.to_csv(OUT_DIR / "bulk_pca_sample_coordinates.tsv", sep="\t", index=False)

    colors = {"DP": "#4b9ed8", "DN": "#df7b35"}
    markers = {"d0_nr": "o", "d2_ova": "s", "d8_ova": "^"}
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    for _, row in pca_df.iterrows():
        ax.scatter(
            row["PC1"],
            row["PC2"],
            s=72,
            color=colors.get(row["population"], "#777777"),
            marker=markers.get(row["condition"], "o"),
            edgecolor="black",
            linewidth=0.5,
        )
        ax.text(row["PC1"], row["PC2"], f" {row['sample_id']}", fontsize=7, va="center")
    ax.axhline(0, color="#d0d0d0", linewidth=0.7)
    ax.axvline(0, color="#d0d0d0", linewidth=0.7)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("Bulk RNA-seq PCA: d0 NR vs OVA rechallenge")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "bulk_pca_d0_nr_d2_d8_ova.png", dpi=300)
    fig.savefig(OUT_DIR / "bulk_pca_d0_nr_d2_d8_ova.pdf")
    plt.close(fig)


def read_deseq(population: str, timepoint: str, filename: str) -> pd.DataFrame:
    df = pd.read_csv(DESEQ_DIR / filename, sep="\t")
    df["population"] = population
    df["timepoint"] = timepoint
    df["contrast"] = f"{population}_{timepoint}_ova_vs_d0_nr"
    df["contrast_label"] = f"{population} {timepoint.upper()} OVA vs D0 NR"
    return df


def ranking_score(df: pd.DataFrame) -> pd.Series:
    p = pd.to_numeric(df["pvalue"], errors="coerce").clip(lower=np.finfo(float).tiny)
    lfc = pd.to_numeric(df["log2FoldChange"], errors="coerce")
    return np.sign(lfc) * -np.log10(p)


def write_ranked_list(df: pd.DataFrame, contrast: str) -> None:
    ranked = df[["gene_name", "rank_score"]].replace([np.inf, -np.inf], np.nan).dropna()
    ranked = ranked.sort_values("rank_score", ascending=False).drop_duplicates("gene_name")
    ranked.to_csv(OUT_DIR / f"{contrast}.rnk", sep="\t", index=False, header=False)


def write_volcano(df: pd.DataFrame, contrast: str, title: str) -> None:
    plot = df.copy()
    plot["log2FoldChange"] = pd.to_numeric(plot["log2FoldChange"], errors="coerce")
    plot["padj"] = pd.to_numeric(plot["padj"], errors="coerce")
    plot["pvalue"] = pd.to_numeric(plot["pvalue"], errors="coerce")
    plot = plot.dropna(subset=["log2FoldChange"])
    plot["neg_log10_padj"] = -np.log10(plot["padj"].fillna(1).clip(lower=np.finfo(float).tiny))
    plot["status"] = "not significant"
    plot.loc[(plot["padj"] < 0.05) & (plot["log2FoldChange"] >= 1), "status"] = "up"
    plot.loc[(plot["padj"] < 0.05) & (plot["log2FoldChange"] <= -1), "status"] = "down"

    colors = {"up": "#b2182b", "down": "#2166ac", "not significant": "#9aa0a6"}
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    for status, sub in plot.groupby("status"):
        ax.scatter(
            sub["log2FoldChange"],
            sub["neg_log10_padj"],
            s=7,
            alpha=0.55,
            color=colors[status],
            linewidths=0,
            label=f"{status} (n={len(sub)})",
        )
    ax.axvline(-1, color="#444444", linestyle="--", linewidth=0.8)
    ax.axvline(1, color="#444444", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(0.05), color="#444444", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)

    top = plot.sort_values(["padj", "pvalue"], na_position="last").head(12)
    for _, row in top.iterrows():
        if pd.notna(row["gene_name"]):
            ax.text(row["log2FoldChange"], row["neg_log10_padj"], f" {row['gene_name']}", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"{contrast}_volcano.png", dpi=300)
    fig.savefig(OUT_DIR / f"{contrast}_volcano.pdf")
    plt.close(fig)


def write_deseq_outputs() -> list[str]:
    all_tables = []
    contrasts = []
    for population, timepoint, filename in CONTRASTS:
        df = read_deseq(population, timepoint, filename)
        df["rank_score"] = ranking_score(df)
        contrast = f"bulk_{population}_{timepoint}_ova_vs_d0_nr"
        contrasts.append(contrast)
        df.to_csv(OUT_DIR / f"{contrast}_deseq_all_genes.tsv", sep="\t", index=False)
        write_ranked_list(df, contrast)
        write_volcano(df, contrast, df["contrast_label"].iloc[0])
        all_tables.append(df)

    combined = pd.concat(all_tables, ignore_index=True)
    combined.to_csv(OUT_DIR / "bulk_deseq_d0_vs_rechallenge_all_contrasts.tsv", sep="\t", index=False)

    summary = (
        combined.assign(
            significant=(
                (pd.to_numeric(combined["padj"], errors="coerce") < 0.05)
                & (pd.to_numeric(combined["log2FoldChange"], errors="coerce").abs() >= 1)
            )
        )
        .groupby("contrast", observed=True)
        .agg(
            n_genes=("gene_name", "count"),
            n_significant=("significant", "sum"),
            n_up=("log2FoldChange", lambda x: int(((combined.loc[x.index, "significant"]) & (x >= 1)).sum())),
            n_down=("log2FoldChange", lambda x: int(((combined.loc[x.index, "significant"]) & (x <= -1)).sum())),
        )
        .reset_index()
    )
    summary.to_csv(OUT_DIR / "bulk_deseq_d0_vs_rechallenge_summary.tsv", sep="\t", index=False)
    return contrasts


def gsea_gene_sets() -> list[str]:
    value = os.environ.get("GSEA_GENE_SETS", "").strip()
    if not value:
        return DEFAULT_GENE_SETS
    return [item.strip() for item in value.split(",") if item.strip()]


def safe_label(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in value).strip("_")


def normalize_gsea_result(res: pd.DataFrame, contrast: str, gene_set: str) -> pd.DataFrame:
    result = res.copy()
    if "Term" not in result.columns:
        result = result.reset_index().rename(columns={"index": "Term"})
    result.insert(0, "gene_set", gene_set)
    result.insert(0, "contrast", contrast)
    for col in ["ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]:
        if col in result.columns:
            result[col] = pd.to_numeric(result[col], errors="coerce")
    return result


def run_preranked_gsea(contrasts: list[str]) -> pd.DataFrame:
    try:
        import gseapy as gp
    except ImportError as exc:
        message = (
            "GSEA skipped because the Python package `gseapy` is not installed. "
            "Run `uv sync`, then rerun `uv run python scripts/bulk_pca_de_gsea.py`."
        )
        (OUT_DIR / "GSEA_SKIPPED.txt").write_text(message + "\n")
        print(message, file=sys.stderr)
        return pd.DataFrame()

    gene_sets = gsea_gene_sets()
    all_results = []
    for contrast in contrasts:
        rnk = OUT_DIR / f"{contrast}.rnk"
        for gene_set in gene_sets:
            label = safe_label(gene_set)
            outdir = GSEA_DIR / contrast / label
            outdir.mkdir(parents=True, exist_ok=True)
            try:
                prerank = gp.prerank(
                    rnk=str(rnk),
                    gene_sets=gene_set,
                    outdir=str(outdir),
                    min_size=15,
                    max_size=500,
                    permutation_num=1000,
                    seed=6,
                    threads=1,
                    verbose=False,
                    no_plot=True,
                )
            except Exception as exc:
                failure = outdir / "GSEA_FAILED.txt"
                failure.write_text(f"{type(exc).__name__}: {exc}\n")
                print(f"GSEA failed for {contrast} / {gene_set}: {exc}", file=sys.stderr)
                continue
            result = normalize_gsea_result(prerank.res2d, contrast, gene_set)
            result.to_csv(outdir / f"{contrast}_{label}_gsea_results.tsv", sep="\t", index=False)
            all_results.append(result)

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    combined = combined.sort_values(["gene_set", "contrast", "FDR q-val", "NOM p-val"], na_position="last")
    combined.to_csv(GSEA_DIR / "bulk_prerank_gsea_all_results.tsv", sep="\t", index=False)

    significant = combined[pd.to_numeric(combined.get("FDR q-val"), errors="coerce") < 0.25].copy()
    significant.to_csv(GSEA_DIR / "bulk_prerank_gsea_fdr_lt_0_25.tsv", sep="\t", index=False)

    top = (
        combined.assign(abs_NES=lambda x: pd.to_numeric(x["NES"], errors="coerce").abs())
        .sort_values(["gene_set", "contrast", "FDR q-val", "abs_NES"], ascending=[True, True, True, False])
        .groupby(["gene_set", "contrast"], observed=True)
        .head(25)
        .drop(columns=["abs_NES"])
    )
    top.to_csv(GSEA_DIR / "bulk_prerank_gsea_top25_per_contrast.tsv", sep="\t", index=False)
    write_gsea_heatmaps(combined)
    write_gsea_summary_plots(combined)
    return combined


def clean_term(term: str) -> str:
    term = str(term)
    if " (GO:" in term:
        term = term.split(" (GO:", 1)[0]
    return term


def wrap_label(label: str, width: int = 44) -> str:
    words = str(label).split()
    lines = []
    current = []
    for word in words:
        candidate = " ".join([*current, word])
        if current and len(candidate) > width:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines)


def write_gsea_heatmaps(combined: pd.DataFrame) -> None:
    if combined.empty or not {"gene_set", "contrast", "Term", "NES"}.issubset(combined.columns):
        return
    table = combined.copy()
    table["NES"] = pd.to_numeric(table["NES"], errors="coerce")
    table["FDR q-val"] = pd.to_numeric(table["FDR q-val"], errors="coerce")
    for gene_set, sub in table.groupby("gene_set", observed=True):
        selected_terms = (
            sub.assign(abs_NES=sub["NES"].abs())
            .sort_values(["FDR q-val", "abs_NES"], ascending=[True, False])
            .drop_duplicates("Term")
            .head(30)["Term"]
            .tolist()
        )
        if not selected_terms:
            continue
        heat = (
            sub[sub["Term"].isin(selected_terms)]
            .pivot_table(index="Term", columns="contrast", values="NES", aggfunc="first")
            .reindex(selected_terms)
        )
        height = max(4.8, 0.24 * len(heat) + 1.4)
        fig, ax = plt.subplots(figsize=(7.6, height))
        values = heat.to_numpy(dtype=float)
        vmax = np.nanmax(np.abs(values)) if np.isfinite(values).any() else 1
        im = ax.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(np.arange(heat.shape[1]))
        ax.set_xticklabels(heat.columns, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(np.arange(heat.shape[0]))
        ax.set_yticklabels(heat.index, fontsize=7)
        ax.set_title(f"Preranked GSEA NES: {gene_set}")
        cbar = fig.colorbar(im, ax=ax, shrink=0.75)
        cbar.set_label("NES")
        fig.tight_layout()
        label = safe_label(gene_set)
        fig.savefig(GSEA_DIR / f"bulk_prerank_gsea_{label}_nes_heatmap.png", dpi=300)
        fig.savefig(GSEA_DIR / f"bulk_prerank_gsea_{label}_nes_heatmap.pdf")
        plt.close(fig)


def write_gsea_summary_plots(combined: pd.DataFrame) -> None:
    if combined.empty or not {"gene_set", "contrast", "Term", "NES", "FDR q-val"}.issubset(combined.columns):
        return

    plot = combined.copy()
    plot["NES"] = pd.to_numeric(plot["NES"], errors="coerce")
    plot["FDR q-val"] = pd.to_numeric(plot["FDR q-val"], errors="coerce")
    plot["neg_log10_fdr"] = -np.log10(plot["FDR q-val"].clip(lower=1e-4))
    plot["contrast_label"] = plot["contrast"].map(CONTRAST_LABELS).fillna(plot["contrast"])
    plot["term_clean"] = plot["Term"].map(clean_term)
    plot["direction"] = np.where(plot["NES"] >= 0, "positive NES", "negative NES")

    write_gsea_significant_count_plot(plot)
    for gene_set, sub in plot.groupby("gene_set", observed=True):
        write_gsea_dotplot(sub, gene_set)
        write_gsea_top_barplots(sub, gene_set)


def write_gsea_significant_count_plot(plot: pd.DataFrame) -> None:
    sig = plot[plot["FDR q-val"] < 0.25].copy()
    if sig.empty:
        return
    counts = (
        sig.groupby(["gene_set", "contrast_label", "direction"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["positive NES", "negative NES"], fill_value=0)
    )
    labels = [f"{gene_set}\n{contrast}" for gene_set, contrast in counts.index]
    y = np.arange(len(counts))
    fig, ax = plt.subplots(figsize=(8.5, max(4.8, 0.45 * len(counts))))
    pos = counts["positive NES"].to_numpy()
    neg = -counts["negative NES"].to_numpy()
    ax.barh(y, pos, color="#b2182b", label="positive NES")
    ax.barh(y, neg, color="#2166ac", label="negative NES")
    ax.axvline(0, color="#444444", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Number of pathways with FDR q-val < 0.25")
    ax.set_title("Significant preranked GSEA pathways by contrast")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(GSEA_DIR / "bulk_prerank_gsea_significant_pathway_counts.png", dpi=300)
    fig.savefig(GSEA_DIR / "bulk_prerank_gsea_significant_pathway_counts.pdf")
    plt.close(fig)


def write_gsea_dotplot(sub: pd.DataFrame, gene_set: str) -> None:
    selected_terms = (
        sub.assign(abs_NES=sub["NES"].abs())
        .sort_values(["FDR q-val", "abs_NES"], ascending=[True, False])
        .groupby("contrast", observed=True)
        .head(8)["Term"]
        .drop_duplicates()
        .tolist()
    )
    if not selected_terms:
        return
    dot = sub[sub["Term"].isin(selected_terms)].copy()
    term_order = (
        dot.assign(abs_NES=dot["NES"].abs())
        .sort_values(["FDR q-val", "abs_NES"], ascending=[False, True])["Term"]
        .drop_duplicates()
        .tolist()
    )
    contrast_order = [label for key, label in CONTRAST_LABELS.items() if label in set(dot["contrast_label"])]
    y_lookup = {term: idx for idx, term in enumerate(term_order)}
    x_lookup = {contrast: idx for idx, contrast in enumerate(contrast_order)}
    dot["x"] = dot["contrast_label"].map(x_lookup)
    dot["y"] = dot["Term"].map(y_lookup)
    dot = dot.dropna(subset=["x", "y", "NES", "neg_log10_fdr"])
    if dot.empty:
        return

    sizes = 18 + 28 * dot["neg_log10_fdr"].clip(upper=4)
    vmax = max(1.0, float(np.nanmax(np.abs(dot["NES"]))))
    fig, ax = plt.subplots(figsize=(8.0, max(5.2, 0.32 * len(term_order) + 1.7)))
    scatter = ax.scatter(
        dot["x"],
        dot["y"],
        s=sizes,
        c=dot["NES"],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        edgecolor="#333333",
        linewidth=0.25,
    )
    ax.set_xticks(range(len(contrast_order)))
    ax.set_xticklabels(contrast_order, fontsize=9)
    ax.set_yticks(range(len(term_order)))
    ax.set_yticklabels([wrap_label(clean_term(term)) for term in term_order], fontsize=7)
    ax.set_title(f"Top preranked GSEA pathways: {gene_set}")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.6)
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.7)
    cbar.set_label("NES")
    for size_value in [1, 2, 3]:
        ax.scatter([], [], s=18 + 28 * size_value, color="#777777", alpha=0.7, label=f"-log10 FDR {size_value}")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    fig.tight_layout()
    label = safe_label(gene_set)
    fig.savefig(GSEA_DIR / f"bulk_prerank_gsea_{label}_top_pathway_dotplot.png", dpi=300)
    fig.savefig(GSEA_DIR / f"bulk_prerank_gsea_{label}_top_pathway_dotplot.pdf")
    plt.close(fig)


def write_gsea_top_barplots(sub: pd.DataFrame, gene_set: str) -> None:
    selected = (
        sub.assign(abs_NES=sub["NES"].abs())
        .sort_values(["contrast", "FDR q-val", "abs_NES"], ascending=[True, True, False])
        .groupby("contrast", observed=True)
        .head(10)
        .copy()
    )
    if selected.empty:
        return
    contrasts = [key for key in CONTRAST_LABELS if key in set(selected["contrast"])]
    fig, axes = plt.subplots(len(contrasts), 1, figsize=(8.0, max(6.0, 2.3 * len(contrasts))), sharex=False)
    if len(contrasts) == 1:
        axes = [axes]
    for ax, contrast in zip(axes, contrasts):
        current = selected[selected["contrast"] == contrast].sort_values("NES")
        colors = np.where(current["NES"] >= 0, "#b2182b", "#2166ac")
        ax.barh(np.arange(len(current)), current["NES"], color=colors)
        ax.axvline(0, color="#444444", linewidth=0.8)
        ax.set_yticks(np.arange(len(current)))
        ax.set_yticklabels([wrap_label(clean_term(term), width=50) for term in current["Term"]], fontsize=7)
        ax.set_title(CONTRAST_LABELS.get(contrast, contrast), fontsize=10)
        ax.set_xlabel("NES")
    fig.suptitle(f"Top GSEA NES per contrast: {gene_set}", y=0.995, fontsize=12)
    fig.tight_layout()
    label = safe_label(gene_set)
    fig.savefig(GSEA_DIR / f"bulk_prerank_gsea_{label}_top_nes_barplots.png", dpi=300)
    fig.savefig(GSEA_DIR / f"bulk_prerank_gsea_{label}_top_nes_barplots.pdf")
    plt.close(fig)


def write_readme() -> None:
    text = """# Bulk DE / PCA / GSEA Inputs

This folder contains bulk RNA-seq downstream outputs for d0 NR vs OVA
rechallenge comparisons.

Generated outputs:

- `bulk_pca_d0_nr_d2_d8_ova.png/pdf`: PCA on log2 median-ratio normalized counts.
- `bulk_pca_sample_coordinates.tsv`: PCA coordinates and sample metadata.
- `bulk_*_deseq_all_genes.tsv`: all-gene DESeq2 tables for each contrast.
- `bulk_*_volcano.png/pdf`: volcano plots using all genes.
- `bulk_*.rnk`: preranked GSEA input files.
- `bulk_deseq_d0_vs_rechallenge_all_contrasts.tsv`: combined all-gene DE table.
- `gsea_prerank/`: preranked GSEA result tables, top-pathway dot plots, bar plots,
  significant-count plots, and NES heatmaps.

Contrasts:

- DP d2 OVA vs d0 NR
- DP d8 OVA vs d0 NR
- DN d2 OVA vs d0 NR
- DN d8 OVA vs d0 NR

Ranked GSEA score:

```text
sign(log2FoldChange) * -log10(pvalue)
```

GSEA defaults to the GSEApy/Enrichr libraries:

- `MSigDB_Hallmark_2020`
- `GO_Biological_Process_2023`
- `KEGG_2019_Mouse`

Override with a comma-separated list of GSEApy library names or local GMT paths:

```bash
GSEA_GENE_SETS="/path/to/mouse_hallmark.gmt,/path/to/reactome.gmt" uv run python scripts/bulk_pca_de_gsea.py
```
"""
    (OUT_DIR / "README.md").write_text(text)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    GSEA_DIR.mkdir(parents=True, exist_ok=True)
    _, count_matrix, sample_meta = load_bulk_counts()
    write_pca(count_matrix, sample_meta)
    contrasts = write_deseq_outputs()
    run_preranked_gsea(contrasts)
    write_readme()


if __name__ == "__main__":
    main()

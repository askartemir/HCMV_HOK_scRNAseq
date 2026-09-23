#!/usr/bin/env python3
"""Plot SciViewer directional-analysis GSEA panels.

This script makes the manuscript-style normalized enrichment score (NES)
bar plot from a Table S5/Table S7-style GSEA workbook. The plotting step is
separate from the interactive SCIViewer selection step: it starts from
the saved GSEA result table generated from genes ranked by directional
correlation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_SHEET = "GSEA_MSigDB_Hallmark_2020"

DISPLAY_TERM_FIXES = {
    "Pperoxisome": "Peroxisome",
    "heme Metabolism": "Heme Metabolism",
    "E2F Targets": "E2F Target",
}


def clean_term(term: str) -> str:
    """Return a display label matching the manuscript figure style."""
    label = str(term).replace("_", " ").replace("-", " ").strip()
    return DISPLAY_TERM_FIXES.get(label, label)


def load_gsea_table(path: Path, sheet_name: str = DEFAULT_SHEET) -> pd.DataFrame:
    """Load and validate a GSEA result workbook."""
    df = pd.read_excel(path, sheet_name=sheet_name)
    required = {"Term", "NES", "FDR q-val"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    df = df.copy()
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df["FDR q-val"] = pd.to_numeric(df["FDR q-val"], errors="coerce")
    df = df.dropna(subset=["Term", "NES", "FDR q-val"])
    return df


def select_pathways(df: pd.DataFrame, fdr_cutoff: float, top_n: int | None) -> pd.DataFrame:
    """Filter pathways by FDR and optionally cap to strongest absolute NES values."""
    selected = df[df["FDR q-val"] <= fdr_cutoff].copy()
    selected = selected.sort_values("NES", ascending=True)
    if top_n is not None and len(selected) > top_n:
        selected = selected.reindex(selected["NES"].abs().sort_values(ascending=False).head(top_n).index)
        selected = selected.sort_values("NES", ascending=True)
    if selected.empty:
        raise ValueError(f"No pathways passed FDR <= {fdr_cutoff}")
    selected["display_term"] = selected["Term"].map(clean_term)
    return selected.reset_index(drop=True)


def plot_gsea_bar(
    selected: pd.DataFrame,
    output_path: Path,
    title: str | None = None,
    positive_color: str = "#e41a1c",
    negative_color: str = "#3b5aa9",
) -> None:
    """Save an NES horizontal bar plot."""
    colors = [positive_color if value > 0 else negative_color for value in selected["NES"]]
    fig_height = max(3.2, 0.32 * len(selected) + 1.2)
    fig, ax = plt.subplots(figsize=(6.2, fig_height))
    ax.barh(selected["display_term"], selected["NES"], color=colors, edgecolor="none")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Normalized Enrichment Score")
    ax.set_ylabel("")
    if title:
        ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=9)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.5)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_selected_pathways(selected: pd.DataFrame, output_path: Path) -> None:
    """Save the plotted pathway subset for traceability."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    keep = [col for col in ["Term", "NES", "FDR q-val", "NOM p-val", "FWER p-val", "Lead_genes"] if col in selected.columns]
    selected[keep].to_csv(output_path, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Table S5/S7-style GSEA workbook")
    parser.add_argument("--output", required=True, type=Path, help="Output figure path, usually .svg or .png")
    parser.add_argument("--sheet", default=DEFAULT_SHEET, help="Workbook sheet containing GSEA results")
    parser.add_argument("--fdr", default=0.1, type=float, help="FDR q-value cutoff")
    parser.add_argument("--top-n", default=None, type=int, help="Optional maximum number of pathways to plot")
    parser.add_argument("--title", default=None, help="Optional plot title")
    parser.add_argument("--selected-output", default=None, type=Path, help="Optional TSV of plotted pathways")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    table = load_gsea_table(args.input, sheet_name=args.sheet)
    selected = select_pathways(table, fdr_cutoff=args.fdr, top_n=args.top_n)
    plot_gsea_bar(selected, args.output, title=args.title)
    if args.selected_output is not None:
        write_selected_pathways(selected, args.selected_output)
    print(f"Wrote: {args.output}")
    if args.selected_output is not None:
        print(f"Wrote: {args.selected_output}")
    print(f"Plotted pathways: {len(selected)}")


if __name__ == "__main__":
    main()

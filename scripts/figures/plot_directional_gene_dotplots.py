#!/usr/bin/env python3
"""Plot manuscript-style directional-analysis gene dotplots.

The Figure 5 directional-analysis gene panels are dotplots, not heatmaps.
This script reproduces those panels from the all-cell AnnData object and a
Table S4/Table S6-style directional-correlation workbook.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


GROUP_ORDER = [
    "Mock",
    "Bystander",
    "Marginal Infection",
    "<10%",
    "10-20%",
    "20-30%",
    "30-40%",
    "40-50%",
    "50-60%",
    "60-70%",
    ">70%",
]

DA1_POSITIVE_FIGURE_GENES = [
    "CKB", "TCF15", "RHEBL1", "EIF1B", "UCHL1", "PLAAT3", "TUBB2B", "RNF157", "ZFR2",
    "HLA-B", "TIMP2", "RPP25", "MLEC", "EEF1A2", "UBB", "RBP7", "BEX2", "HAGH", "HES6",
    "ENO2", "PCDH1", "GABARAPL1", "SDF2L1", "SELENOK", "PRSS22", "DDIT3", "UBE2S",
    "ATP6V0B", "WFDC2", "CTSF", "UBE2T", "PGP", "PRADC1",
]

DA1_NEGATIVE_FIGURE_GENES = [
    "KLHL5", "SORL1", "RAC2", "FSTL1", "PYCARD", "TMEM256", "ACSS2", "PRSS12", "LGALS3",
    "SAT1", "EML2", "CCDC25", "HEBP2", "SNX6", "MRPL45", "ATP5ME", "OST4", "APLP2",
    "NDUFS4", "UQCRB", "PSME2", "SUCLG2", "DECR1", "FDPS", "ATOX1", "LGALS1", "PPFIBP1",
    "TAPBP", "IQGAP1", "CMTM6",
]

DA2_POSITIVE_FIGURE_GENES = [
    "SDF2L1", "HERC5", "EEF1A2", "MEIOC", "PLAAT3", "KIF1A", "RNF157", "TUBB2B",
    "TNFRSF18", "ZFR2", "HES6", "ENO2", "C1QL1", "NIBAN1", "HSPA13", "TCF15",
    "RHEBL1", "UCHL1", "EIF1B", "GCLM", "NRIP3", "CD55", "DDIT3", "ATF3", "TRIB3",
    "NMRK2", "C1QL2", "C1QL4", "PCK2", "STC2",
]

DA2_NEGATIVE_FIGURE_GENES = [
    "CDH1", "S100A9", "KRT6B", "IL1RN", "SERPINB7", "IL32", "C3", "CD24", "PYCARD",
    "TAPBP", "KRT8", "CD9", "LGALS3", "CD82", "MVP", "CD47", "TOLLIP", "S100A13",
    "IFI27L2", "CD81", "UQCRB", "VAMP8", "IRF2BP2", "CASP4", "PDCD4", "S100A8",
    "IRF6", "ISG20", "ISG15", "FOXQ1",
]

FIGURE_GENE_ORDERS = {
    ("directional_analysis_1", "positive"): DA1_POSITIVE_FIGURE_GENES,
    ("directional_analysis_1", "negative"): DA1_NEGATIVE_FIGURE_GENES,
    ("directional_analysis_2", "positive"): DA2_POSITIVE_FIGURE_GENES,
    ("directional_analysis_2", "negative"): DA2_NEGATIVE_FIGURE_GENES,
}

STATE_LABELS = {
    0: "Mock",
    1: "Bystander",
    2: "Marginal Infection",
}


@dataclass(frozen=True)
class DotplotResult:
    direction: str
    genes: pd.DataFrame
    dotplot_values: pd.DataFrame
    figure_path: Path
    table_path: Path


def load_correlation_table(path: Path, sheet_name: str | int = 0) -> pd.DataFrame:
    """Load directional-correlation results and standardize required columns."""
    df = pd.read_excel(path, sheet_name=sheet_name)
    required = {"Gene name", "R", "P"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    df = df[["Gene name", "R", "P"]].copy()
    df["Gene name"] = df["Gene name"].astype(str).str.strip()
    df["R"] = pd.to_numeric(df["R"], errors="coerce")
    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    return df.dropna(subset=["Gene name", "R", "P"])


def select_directional_genes(correlations: pd.DataFrame, var_names: pd.Index, direction: str, top_n: int) -> pd.DataFrame:
    """Select top correlated genes in the same order used by the figure."""
    if direction == "positive":
        selected = correlations[correlations["R"] > 0].sort_values("R", ascending=False)
    elif direction == "negative":
        selected = correlations[correlations["R"] < 0].sort_values("R", ascending=True)
    else:
        raise ValueError("direction must be 'positive' or 'negative'")
    selected = selected[selected["Gene name"].isin(var_names)].head(top_n).reset_index(drop=True)
    if selected.empty:
        raise ValueError(f"No {direction} correlated genes from the table were found in AnnData")
    return selected


def select_figure_ordered_genes(
    correlations: pd.DataFrame,
    var_names: pd.Index,
    direction: str,
    analysis: str,
) -> pd.DataFrame:
    """Select the final manuscript gene order visible in the figure PNG."""
    key = (analysis, direction)
    if key not in FIGURE_GENE_ORDERS:
        raise ValueError(f"No figure gene order is defined for {analysis!r} / {direction!r}")
    genes = FIGURE_GENE_ORDERS[key]
    missing_from_adata = [gene for gene in genes if gene not in var_names]
    if missing_from_adata:
        raise ValueError(f"Figure genes are missing from AnnData: {missing_from_adata}")
    ordered = pd.DataFrame({"Gene name": genes})
    merged = ordered.merge(correlations, on="Gene name", how="left")
    missing_from_table = merged[merged["R"].isna()]["Gene name"].tolist()
    if missing_from_table:
        raise ValueError(f"Figure genes are missing from the correlation table: {missing_from_table}")
    return merged


def make_dotplot_group(
    obs: pd.DataFrame,
    dpi: str = "3dpi",
    state_column: str = "Infection_state_bkgd",
    percent_column: str = "percent.cmv",
) -> pd.Series:
    """Create Mock/Bystander/Marginal plus high-infection viral-percent groups."""
    required = {state_column, percent_column, "dpi"}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(f"adata.obs is missing columns: {sorted(missing)}")

    out = pd.Series(pd.NA, index=obs.index, dtype="object")
    dpi_mask = obs["dpi"].astype(str).eq(dpi)
    states = pd.to_numeric(obs[state_column], errors="coerce")
    percent = pd.to_numeric(obs[percent_column], errors="coerce")

    for state, label in STATE_LABELS.items():
        out.loc[dpi_mask & states.eq(state)] = label

    high = dpi_mask & states.eq(3)
    out.loc[high & percent.lt(10)] = "<10%"
    out.loc[high & percent.ge(10) & percent.lt(20)] = "10-20%"
    out.loc[high & percent.ge(20) & percent.lt(30)] = "20-30%"
    out.loc[high & percent.ge(30) & percent.lt(40)] = "30-40%"
    out.loc[high & percent.ge(40) & percent.lt(50)] = "40-50%"
    out.loc[high & percent.ge(50) & percent.lt(60)] = "50-60%"
    out.loc[high & percent.ge(60) & percent.lt(70)] = "60-70%"
    out.loc[high & percent.ge(70)] = ">70%"
    return pd.Categorical(out, categories=GROUP_ORDER, ordered=True)


def expression_frame(adata: ad.AnnData, genes: list[str], use_raw: bool) -> pd.DataFrame:
    """Return cells x genes expression matrix."""
    source = adata.raw if use_raw else adata
    if source is None:
        raise ValueError("--use-raw was requested, but AnnData has no .raw matrix")

    missing = [gene for gene in genes if gene not in source.var_names]
    if missing:
        raise ValueError(f"Selected genes are missing from expression matrix: {missing[:10]}")

    matrix = source[:, genes].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(np.asarray(matrix), index=adata.obs_names, columns=genes)


def summarize_expression(sub: pd.DataFrame, average_method: str) -> pd.Series:
    """Summarize group expression with either direct or Seurat-style averaging."""
    if average_method == "direct":
        return sub.mean(axis=0)
    if average_method == "seurat":
        if (sub.to_numpy() < 0).any():
            raise ValueError("Seurat-style averaging requires nonnegative log-normalized values; try --use-raw")
        return pd.Series(np.log1p(np.expm1(sub).mean(axis=0)), index=sub.columns)
    raise ValueError("average_method must be 'direct' or 'seurat'")


def compute_dotplot_values(
    adata: ad.AnnData,
    genes: list[str],
    group_column: str,
    use_raw: bool,
    average_method: str,
) -> pd.DataFrame:
    """Compute average expression, z-scored average expression, and percent expressing."""
    expr = expression_frame(adata, genes, use_raw=use_raw)
    groups = adata.obs[group_column]

    records = []
    for group in GROUP_ORDER:
        mask = groups.astype(str).eq(group).to_numpy()
        n_cells = int(mask.sum())
        if n_cells == 0:
            mean_values = pd.Series(np.nan, index=genes)
            percent_values = pd.Series(0.0, index=genes)
        else:
            sub = expr.loc[mask, :]
            mean_values = summarize_expression(sub, average_method=average_method)
            percent_values = (sub > 0).sum(axis=0) / n_cells * 100
        for gene in genes:
            records.append(
                {
                    "group": group,
                    "gene": gene,
                    "n_cells": n_cells,
                    "average_expression": float(mean_values[gene]),
                    "percent_expressing": float(percent_values[gene]),
                }
            )

    values = pd.DataFrame.from_records(records)
    wide = values.pivot(index="gene", columns="group", values="average_expression").reindex(index=genes, columns=GROUP_ORDER)
    scale = wide.std(axis=1, ddof=0).replace(0, np.nan)
    scaled = wide.sub(wide.mean(axis=1), axis=0).div(scale, axis=0).fillna(0).clip(-2, 2)
    values = values.merge(
        scaled.stack().rename("average_normalized_expression").reset_index(),
        on=["gene", "group"],
        how="left",
    )
    return values


def plot_dotplot(
    values: pd.DataFrame,
    genes: list[str],
    output_path: Path,
    title: str,
    dot_scale: float = 90.0,
) -> None:
    """Save one manuscript-style dotplot."""
    x_lookup = {group: idx for idx, group in enumerate(GROUP_ORDER)}
    y_lookup = {gene: idx for idx, gene in enumerate(reversed(genes))}
    plot_df = values.copy()
    plot_df["x"] = plot_df["group"].map(x_lookup)
    plot_df["y"] = plot_df["gene"].map(y_lookup)
    plot_df["size"] = np.clip(plot_df["percent_expressing"], 0, 100) / 100 * dot_scale

    fig_height = max(5.0, 0.21 * len(genes) + 1.8)
    fig, ax = plt.subplots(figsize=(4.7, fig_height))
    scatter = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["size"],
        c=plot_df["average_normalized_expression"],
        cmap="Purples",
        vmin=-2,
        vmax=2,
        edgecolors="none",
    )

    ax.set_title(title, fontsize=10)
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(GROUP_ORDER, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(list(reversed(genes)), fontsize=8)
    ax.yaxis.tick_right()
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, len(GROUP_ORDER) - 0.5)
    ax.set_ylim(-0.5, len(genes) - 0.5)
    ax.grid(color="#e3e3e3", linewidth=0.6)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.32)
    cbar.set_label("Average\nNormalized\nExpression", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    handles = [
        ax.scatter([], [], s=size / 100 * dot_scale, color="black", edgecolors="none")
        for size in (10, 25, 50, 75, 100)
    ]
    ax.legend(
        handles,
        ["10", "25", "50", "75", "100"],
        title="Percent Expressing",
        frameon=False,
        fontsize=7,
        title_fontsize=8,
        loc="lower left",
        bbox_to_anchor=(1.55, 0.0),
        borderaxespad=0,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_dotplot_table(genes: pd.DataFrame, values: pd.DataFrame, output_path: Path) -> None:
    """Write selected genes and plotted dotplot values for traceability."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        genes.to_excel(writer, index=False, sheet_name="selected_genes")
        values.to_excel(writer, index=False, sheet_name="dotplot_values")
        values.pivot(index="gene", columns="group", values="average_expression").reindex(
            index=genes["Gene name"], columns=GROUP_ORDER
        ).to_excel(writer, sheet_name="average_expression")
        values.pivot(index="gene", columns="group", values="average_normalized_expression").reindex(
            index=genes["Gene name"], columns=GROUP_ORDER
        ).to_excel(writer, sheet_name="average_normalized_expression")
        values.pivot(index="gene", columns="group", values="percent_expressing").reindex(
            index=genes["Gene name"], columns=GROUP_ORDER
        ).to_excel(writer, sheet_name="percent_expressing")


def build_dotplot(
    adata: ad.AnnData,
    correlations: pd.DataFrame,
    direction: str,
    output_dir: Path,
    top_n: int,
    use_raw: bool,
    file_prefix: str,
    title: str,
    figure_order: bool,
    analysis: str = "directional_analysis_1",
    average_method: str = "seurat",
) -> DotplotResult:
    """Build and save one directional gene dotplot."""
    source = adata.raw if use_raw else adata
    if source is None:
        raise ValueError("--use-raw was requested, but AnnData has no .raw matrix")
    if figure_order:
        genes = select_figure_ordered_genes(
            correlations,
            pd.Index(source.var_names),
            direction=direction,
            analysis=analysis,
        )
    else:
        genes = select_directional_genes(correlations, pd.Index(source.var_names), direction=direction, top_n=top_n)
    values = compute_dotplot_values(
        adata,
        genes["Gene name"].tolist(),
        group_column="figure5_dotplot_group",
        use_raw=use_raw,
        average_method=average_method,
    )
    gene_count = len(genes)
    figure_path = output_dir / f"{file_prefix}_{direction}_top{gene_count}_genes_dotplot.png"
    table_path = output_dir / f"{file_prefix}_{direction}_top{gene_count}_genes_dotplot_values.xlsx"
    plot_dotplot(values, genes["Gene name"].tolist(), figure_path, title=title)
    write_dotplot_table(genes, values, table_path)
    return DotplotResult(direction, genes, values, figure_path, table_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True, type=Path, help="All-cell AnnData object with Infection_state_bkgd")
    parser.add_argument("--correlations", required=True, type=Path, help="Table S4/S6-style correlation workbook")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for dotplot outputs")
    parser.add_argument("--dpi", default="3dpi", help="DPI subset to plot")
    parser.add_argument("--positive-top-n", default=33, type=int, help="Number of positive genes to plot")
    parser.add_argument("--negative-top-n", default=30, type=int, help="Number of negative genes to plot")
    parser.add_argument(
        "--rank-order",
        action="store_true",
        help="Use strict correlation-rank order instead of the final Figure 5 display order",
    )
    parser.add_argument("--use-raw", action="store_true", help="Use adata.raw instead of adata.X")
    parser.add_argument(
        "--average-method",
        choices=["seurat", "direct"],
        default="seurat",
        help="Use Seurat DotPlot-style log-normalized averaging or direct averaging of stored values",
    )
    parser.add_argument("--file-prefix", default="directional_analysis_1", help="Output file prefix")
    parser.add_argument(
        "--analysis",
        choices=["directional_analysis_1", "directional_analysis_2"],
        default="directional_analysis_1",
        help="Manuscript figure gene-order set to use unless --rank-order is supplied",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.adata)
    adata.obs["figure5_dotplot_group"] = make_dotplot_group(adata.obs, dpi=args.dpi)
    keep = ~pd.isna(adata.obs["figure5_dotplot_group"])
    adata = adata[keep, :].copy()
    correlations = load_correlation_table(args.correlations)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    specs = [
        ("positive", args.positive_top_n, "sciViewer top positively correlated\ngene expression in dpi3 highly infected cells"),
        ("negative", args.negative_top_n, "sciViewer top negatively correlated\ngene expression in dpi3 highly infected cells"),
    ]
    for direction, top_n, title in specs:
        result = build_dotplot(
            adata,
            correlations,
            direction=direction,
            output_dir=args.output_dir,
            top_n=top_n,
            use_raw=args.use_raw,
            file_prefix=args.file_prefix,
            title=title,
            figure_order=not args.rank_order,
            analysis=args.analysis,
            average_method=args.average_method,
        )
        print(f"Wrote: {result.figure_path}")
        print(f"Wrote: {result.table_path}")


if __name__ == "__main__":
    main()

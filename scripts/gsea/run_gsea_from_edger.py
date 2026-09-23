#!/usr/bin/env python3
"""Reproduce Figure 4 GSEA panels from saved edgeR/GSEA outputs.

The manuscript Figure 4 volcano plots and upstream edgeR statistics are outside
the scope of this script. This script reproduces the GSEA NES panels from saved
edgeR/GSEA result workbooks.

Default use is `plot-saved`: read the saved GSEA sheets already present in the
edgeR workbooks and regenerate Figure 4E-H-style NES bar plots plus pathway TSVs.
The optional `run-prerank` command reruns GSEA from `Gene` and `logFC`.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DISPLAY_TERM_FIXES = {
    "Pperoxisome": "Peroxisome",
    "heme Metabolism": "Heme Metabolism",
    "E2F Targets": "E2F Target",
    "G2-M Checkpoint": "G2-M Checkpoint",
    "KRAS Signaling Dn": "KRAS Signaling Downregulation",
    "KRAS Signaling Up": "KRAS Signaling Upregulation",
    "UV Response Dn": "UV Response Downregulation",
    "UV Response Up": "UV Response Upregulation",
}


@dataclass(frozen=True)
class Figure4Panel:
    panel: str
    name: str
    title: str
    workbook_relative_path: Path
    gsea_sheet: str


@dataclass(frozen=True)
class GseaJob:
    name: str
    input_workbook: Path
    sheet_name: str = "Sheet1"
    gene_column: str = "Gene"
    score_column: str = "logFC"
    gene_sets: str = "MSigDB_Hallmark_2020"


FIGURE4_PANELS = [
    Figure4Panel(
        panel="E",
        name="mock_vs_highly_infected",
        title="Mock vs Highly infected cells",
        workbook_relative_path=Path("DE_mock_vs_high/6_Infected_edgeRDEresults.xlsx"),
        gsea_sheet="GSEA_mock_vs_high_cmv_transc_MSigDB_Hallmark_2024",
    ),
    Figure4Panel(
        panel="F",
        name="mock_vs_marginally_infected",
        title="Mock vs Marginally infected cells",
        workbook_relative_path=Path("DE_mock_vs_low/5_Bystander_bkgd_edgeRDEresults.xlsx"),
        gsea_sheet="GSEA_mock_vs_low_cmv_transc_cells_MSigDB_Hallmark_2024",
    ),
    Figure4Panel(
        panel="G",
        name="mock_vs_bystander",
        title="Mock vs Bystander cells",
        workbook_relative_path=Path("DE_mock_vs_bystander/4_Bystander_edgeRDEresults.xlsx"),
        gsea_sheet="GSEA_mock_vs_bystander_MSigDB_Hallmark_2024",
    ),
    Figure4Panel(
        panel="H",
        name="bystander_vs_highly_infected",
        title="Bystander vs Highly infected cells",
        workbook_relative_path=Path("7_Infected_edgeRDEresults.xlsx"),
        gsea_sheet="GSEA_highinf_vs_bystander_MSigDB_Hallmark_2024",
    ),
]


def clean_term(term: str) -> str:
    """Return a display label matching the manuscript figure style."""
    label = str(term).replace("_", " ").strip()
    return DISPLAY_TERM_FIXES.get(label, label)


def load_gsea_sheet(workbook: Path, sheet_name: str) -> pd.DataFrame:
    """Load a saved GSEA sheet and validate columns needed for plotting."""
    df = pd.read_excel(workbook, sheet_name=sheet_name)
    required = {"Term", "NES", "FDR q-val"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{workbook} / {sheet_name} is missing columns: {sorted(missing)}")
    df = df.copy()
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df["FDR q-val"] = pd.to_numeric(df["FDR q-val"], errors="coerce")
    df = df.dropna(subset=["Term", "NES", "FDR q-val"])
    return df


def select_pathways(df: pd.DataFrame, fdr_cutoff: float, top_n: int | None) -> pd.DataFrame:
    """Select FDR-significant pathways and keep strongest absolute NES values."""
    selected = df[df["FDR q-val"] <= fdr_cutoff].copy()
    if top_n is not None and len(selected) > top_n:
        keep = selected["NES"].abs().sort_values(ascending=False).head(top_n).index
        selected = selected.loc[keep].copy()
    selected = selected.sort_values("NES", ascending=True)
    if selected.empty:
        raise ValueError(f"No pathways passed FDR <= {fdr_cutoff}")
    selected["display_term"] = selected["Term"].map(clean_term)
    return selected.reset_index(drop=True)


def plot_nes_panel(
    selected: pd.DataFrame,
    output_path: Path,
    title: str,
    positive_color: str = "#e41a1c",
    negative_color: str = "#3b5aa9",
) -> None:
    """Save one Figure 4-style NES horizontal bar panel."""
    colors = [positive_color if value > 0 else negative_color for value in selected["NES"]]
    fig_height = max(3.2, 0.25 * len(selected) + 1.1)
    fig, ax = plt.subplots(figsize=(4.2, fig_height))
    ax.barh(selected["display_term"], selected["NES"], color=colors, edgecolor="none")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Normalized Enrichment Score")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=10)
    ax.tick_params(axis="both", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.5)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def write_selected_pathways(selected: pd.DataFrame, output_path: Path) -> None:
    """Save the plotted pathway subset for traceability."""
    keep = [col for col in ["Term", "NES", "FDR q-val", "NOM p-val", "FWER p-val", "Lead_genes"] if col in selected.columns]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    selected[keep].to_csv(output_path, sep="\t", index=False)


def write_panel_manifest(records: list[dict[str, object]], output_path: Path) -> None:
    """Write one TSV summarizing Figure 4 panel inputs/outputs."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame.from_records(records).to_csv(output_path, sep="\t", index=False)


def plot_saved_panels(input_root: Path, output_dir: Path, fdr_cutoff: float, top_n: int | None) -> None:
    """Plot Figure 4E-H from saved GSEA workbook sheets."""
    records: list[dict[str, object]] = []
    for panel in FIGURE4_PANELS:
        workbook = input_root / panel.workbook_relative_path
        if not workbook.exists():
            raise FileNotFoundError(f"Missing Figure 4{panel.panel} input workbook: {workbook}")
        gsea = load_gsea_sheet(workbook, panel.gsea_sheet)
        selected = select_pathways(gsea, fdr_cutoff=fdr_cutoff, top_n=top_n)

        prefix = f"figure4{panel.panel}_{panel.name}"
        svg_path = output_dir / f"{prefix}_gsea_fdr_{fdr_cutoff:g}.svg"
        png_path = output_dir / f"{prefix}_gsea_fdr_{fdr_cutoff:g}.png"
        tsv_path = output_dir / f"{prefix}_plotted_pathways.tsv"
        plot_title = f"Top Enriched Pathways (FDR <= {fdr_cutoff:g})\n{panel.title}"

        plot_nes_panel(selected, svg_path, title=plot_title)
        plot_nes_panel(selected, png_path, title=plot_title)
        write_selected_pathways(selected, tsv_path)

        records.append(
            {
                "panel": f"Figure 4{panel.panel}",
                "comparison": panel.title,
                "input_workbook": str(workbook),
                "gsea_sheet": panel.gsea_sheet,
                "fdr_cutoff": fdr_cutoff,
                "top_n": top_n if top_n is not None else "all",
                "plotted_pathways": len(selected),
                "svg_output": str(svg_path),
                "png_output": str(png_path),
                "pathway_tsv": str(tsv_path),
            }
        )
        print(f"Wrote: {svg_path}")
        print(f"Wrote: {png_path}")
        print(f"Wrote: {tsv_path}")
        print(f"Plotted pathways for Figure 4{panel.panel}: {len(selected)}")

    manifest_path = output_dir / f"figure4_gsea_panel_manifest_fdr_{fdr_cutoff:g}.tsv"
    write_panel_manifest(records, manifest_path)
    print(f"Wrote: {manifest_path}")


def load_ranked_genes(job: GseaJob) -> pd.DataFrame:
    """Load and sort preranked genes from an edgeR result workbook."""
    df = pd.read_excel(job.input_workbook, sheet_name=job.sheet_name)
    missing = {job.gene_column, job.score_column} - set(df.columns)
    if missing:
        raise ValueError(f"{job.input_workbook} is missing columns: {sorted(missing)}")

    ranked = df[[job.gene_column, job.score_column]].copy()
    ranked = ranked.dropna()
    ranked[job.score_column] = pd.to_numeric(ranked[job.score_column], errors="coerce")
    ranked = ranked.dropna()
    ranked = ranked.sort_values(job.score_column, ascending=False)
    ranked = ranked.drop_duplicates(subset=job.gene_column, keep="first")
    return ranked


def run_prerank(
    job: GseaJob,
    output_dir: Path,
    min_size: int = 5,
    max_size: int = 5000,
    permutations: int = 1000,
) -> pd.DataFrame:
    """Run gseapy prerank from Gene/logFC. Requires gseapy installed."""
    try:
        import gseapy as gp
    except ImportError as exc:
        raise ImportError("The run-prerank command requires gseapy. Use plot-saved for manuscript panel reproduction.") from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    ranked = load_ranked_genes(job)
    result = gp.prerank(
        rnk=ranked,
        gene_sets=job.gene_sets,
        outdir=str(output_dir / job.name),
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutations,
        seed=7,
        verbose=True,
    )
    result_table = result.res2d.copy()
    result_table.to_csv(output_dir / f"{job.name}_gsea_results.tsv", sep="\t", index=False)
    return result_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    plot_saved = subparsers.add_parser("plot-saved", help="Plot Figure 4E-H from saved workbook GSEA sheets")
    plot_saved.add_argument("--input-root", required=True, type=Path, help="Folder containing the edgeR workbooks")
    plot_saved.add_argument("--outdir", default=Path("outputs/figures/figure4_gsea"), type=Path)
    plot_saved.add_argument("--fdr", default=0.1, type=float)
    plot_saved.add_argument("--top-n", default=30, type=int)

    rerun = subparsers.add_parser("run-prerank", help="Rerun GSEA from one edgeR workbook's Gene/logFC columns")
    rerun.add_argument("--input", required=True, type=Path, help="edgeR result workbook (.xlsx)")
    rerun.add_argument("--name", required=True, help="Short output name for this comparison")
    rerun.add_argument("--outdir", default=Path("outputs/gsea"), type=Path)
    rerun.add_argument("--gene-sets", default="MSigDB_Hallmark_2020")
    rerun.add_argument("--min-size", default=5, type=int)
    rerun.add_argument("--max-size", default=5000, type=int)
    rerun.add_argument("--permutations", default=1000, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "plot-saved":
        plot_saved_panels(args.input_root, args.outdir, fdr_cutoff=args.fdr, top_n=args.top_n)
    elif args.command == "run-prerank":
        job = GseaJob(name=args.name, input_workbook=args.input, gene_sets=args.gene_sets)
        result = run_prerank(
            job,
            args.outdir,
            min_size=args.min_size,
            max_size=args.max_size,
            permutations=args.permutations,
        )
        print(f"Wrote: {args.outdir / f'{args.name}_gsea_results.tsv'}")
        print(f"GSEA result rows: {len(result)}")
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()

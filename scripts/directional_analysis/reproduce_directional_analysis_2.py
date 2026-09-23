#!/usr/bin/env python3
"""Reproduce Directional Analysis 2 supplemental outputs.

This script formats the saved SCIViewer Directional Analysis 2 export into the
submitted supplemental tables.

Main outputs:
- Table S6: raw Directional Analysis 2 projection-correlation results.
- Table S7: MSigDB Hallmark GSEA results for Directional Analysis 2.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Border, Side

TABLE_S6_SHEET = "Directional Analysis 2 Results"
TABLE_S7_SHEET = "GSEA_MSigDB_Hallmark_2020"
DA2_GSEA_SOURCE_SHEET = "GSEA_Unfiltered_Curated_MSigDB_Hallmark_2025"

TABLE_S6_LEGEND = [
    ("Column Name", "Description"),
    ("Gene name", "Official gene symbol standardized by HGNC (HUGO Gene Nomenclature Committee)"),
    ("R", "Pearson correlation coefficient: indicating relationship strength between gene expression among selected cells and selected direction"),
    ("P", "P-value: statistical significance of correlation (values < 0.05 typically considered statistically significant)"),
    (
        "Correlation Strength",
        "Correlation strength was defined based on the absolute Pearson correlation coefficient (|R|), independent of direction. Pearson correlation coefficients (R) were calculated between gene expression and projection coordinate derived from Directional Analysis #2. Direction (positive or negative) was assigned based on the sign of R.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Strength Categories (Magnitude). \n                                                                                                                                                                                                                                  Strong Correlation\n|R| ≥ 0.7\n\nModerate Correlation\n0.3 ≤ |R| < 0.7\n\nWeak Correlation\n|R| < 0.3\n\nDirectional Categories (Sign + Strength)\n\nStrong Positive Correlation\nR ≥ 0.7\n\nModerate Positive Correlation\n0.3 ≤ R < 0.7\n\nWeak Positive Correlation\n0 < R < 0.3\n\nStrong Negative Correlation\nR ≤ −0.7\n\nModerate Negative Correlation\n−0.7 < R ≤ −0.3\n\nWeak Negative Correlation\n−0.3 < R < 0",
    ),
]


def correlation_strength(r_value: float) -> str:
    """Return the qualitative correlation-strength label defined in the table legend."""
    if pd.isna(r_value):
        return "Other"
    if r_value >= 0.7:
        return "Strong Positive Correlation"
    if 0.3 <= r_value < 0.7:
        return "Moderate Positive Correlation"
    if r_value > 0:
        return "Weak Positive Correlation"
    if r_value <= -0.7:
        return "Strong Negative Correlation"
    if -0.7 < r_value <= -0.3:
        return "Moderate Negative Correlation"
    if r_value < 0:
        return "Weak Negative Correlation"
    return "Other"


def read_table(path: Path, sheet_name: str | int | None = None) -> pd.DataFrame:
    """Read a CSV/TSV/XLSX table based on file extension."""
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix == ".tsv":
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=0 if sheet_name is None else sheet_name)
    raise ValueError(f"Unsupported input table type: {path}")


def clean_projection_correlation(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize projection-correlation columns for Table S6."""
    df = df.copy()
    if "Gene name" not in df.columns:
        first = df.columns[0]
        if str(first).startswith("Unnamed") or first in {"index", "Gene", "gene"}:
            df = df.rename(columns={first: "Gene name"})
    required = {"Gene name", "R", "P"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Projection-correlation table is missing columns: {sorted(missing)}")

    out = df[["Gene name", "R", "P"]].copy()
    out["Gene name"] = out["Gene name"].astype(str).str.strip()
    out["R"] = pd.to_numeric(out["R"], errors="coerce")
    out["P"] = pd.to_numeric(out["P"], errors="coerce")
    out = out.dropna(subset=["Gene name", "R", "P"])
    out["Correlation Strength"] = out["R"].apply(correlation_strength)
    return out.reset_index(drop=True)


def add_table_s6_legend(table: pd.DataFrame) -> pd.DataFrame:
    """Add the side-by-side legend columns used in the submitted supplemental table."""
    out = table.copy()
    out["Unnamed: 4"] = np.nan
    out["Unnamed: 5"] = np.nan
    out["Table Legend"] = pd.Series([pd.NA] * len(out), dtype="object")
    out["Unnamed: 7"] = pd.Series([pd.NA] * len(out), dtype="object")
    for idx, (key, description) in enumerate(TABLE_S6_LEGEND):
        if idx < len(out):
            out.loc[idx, "Table Legend"] = key
            out.loc[idx, "Unnamed: 7"] = description
    return out


def format_table_s6(output_path: Path) -> None:
    """Apply small reviewer-facing formatting without changing table values."""
    workbook = load_workbook(output_path)
    worksheet = workbook[TABLE_S6_SHEET]
    for cell in ("E1", "F1", "H1"):
        worksheet[cell].value = None

    thin_border = Border(
        left=Side(style="thin"),
        right=Side(style="thin"),
        top=Side(style="thin"),
        bottom=Side(style="thin"),
    )
    for row in worksheet["G1:H6"]:
        for cell in row:
            cell.border = thin_border
    workbook.save(output_path)


def write_table_s6(input_path: Path, output_path: Path, sheet_name: str | None = None) -> None:
    raw = read_table(input_path, sheet_name=sheet_name)
    cleaned = clean_projection_correlation(raw)
    table = add_table_s6_legend(cleaned)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name=TABLE_S6_SHEET, index=False)
    format_table_s6(output_path)


def write_table_s7(gsea_workbook: Path, output_path: Path, source_sheet: str = DA2_GSEA_SOURCE_SHEET) -> None:
    df = pd.read_excel(gsea_workbook, sheet_name=source_sheet)
    drop_cols = [col for col in ["Name"] if col in df.columns]
    df = df.drop(columns=drop_cols)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name=TABLE_S7_SHEET, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    s6 = subparsers.add_parser("table-s6", help="Format projection-correlation export as Table S6")
    s6.add_argument("--input", required=True, type=Path, help="Projection-correlation CSV/XLSX")
    s6.add_argument("--output", required=True, type=Path, help="Output Table S6 XLSX")
    s6.add_argument("--sheet", default=None, help="Input sheet name for XLSX files")

    s7 = subparsers.add_parser("table-s7", help="Extract final DA2 GSEA sheet as Table S7")
    s7.add_argument("--input", required=True, type=Path, help="DA2 projection/GSEA workbook")
    s7.add_argument("--output", required=True, type=Path, help="Output Table S7 XLSX")
    s7.add_argument("--source-sheet", default=DA2_GSEA_SOURCE_SHEET)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "table-s6":
        write_table_s6(args.input, args.output, sheet_name=args.sheet)
    elif args.command == "table-s7":
        write_table_s7(args.input, args.output, source_sheet=args.source_sheet)
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()

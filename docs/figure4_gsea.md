# Figure 4 GSEA provenance

Figure 4A-D volcano plots and the upstream edgeR comparisons are outside the scope of this Python script. This repository reproduces the Figure 4E-H GSEA NES panels from saved edgeR/GSEA result workbooks. The required data workbooks are available through GEO accession `GSE348416`.

Cleaned script:

`../scripts/gsea/run_gsea_from_edger.py`

The manuscript-panel command is:

```bash
python scripts/gsea/run_gsea_from_edger.py plot-saved \
  --input-root "/path/to/h5ad and tsv files" \
  --outdir outputs/figures/figure4_gsea \
  --fdr 0.1 \
  --top-n 30
```

This reads saved GSEA sheets from the edgeR workbooks and regenerates Figure 4E-H-style NES bar panels plus TSVs containing the plotted pathways. It does not regenerate the upstream edgeR differential-expression results or volcano plots.

## Panel Mapping

| Panel | Comparison | Workbook | GSEA sheet |
|---|---|---|---|
| Figure 4E | Mock vs Highly infected cells | `DE_mock_vs_high/6_Infected_edgeRDEresults.xlsx` | `GSEA_mock_vs_high_cmv_transc_MSigDB_Hallmark_2024` |
| Figure 4F | Mock vs Marginally infected cells | `DE_mock_vs_low/5_Bystander_bkgd_edgeRDEresults.xlsx` | `GSEA_mock_vs_low_cmv_transc_cells_MSigDB_Hallmark_2024` |
| Figure 4G | Mock vs Bystander cells | `DE_mock_vs_bystander/4_Bystander_edgeRDEresults.xlsx` | `GSEA_mock_vs_bystander_MSigDB_Hallmark_2024` |
| Figure 4H | Bystander vs Highly infected cells | `7_Infected_edgeRDEresults.xlsx` | `GSEA_highinf_vs_bystander_MSigDB_Hallmark_2024` |

The script also has an optional `run-prerank` command that reruns GSEA from an edgeR workbook's `Gene` and `logFC` columns. The default manuscript-panel command uses saved GSEA sheets because that path is deterministic and reproduces the submitted analysis source.

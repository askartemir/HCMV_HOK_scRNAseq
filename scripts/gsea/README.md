# GSEA scripts

Cleaned script:

- `run_gsea_from_edger.py`: reproduces the Figure 4E-H GSEA NES panels from saved edgeR/GSEA result workbooks.

Default manuscript-panel command:

```bash
python scripts/gsea/run_gsea_from_edger.py plot-saved \
  --input-root "/path/to/h5ad and tsv files" \
  --outdir outputs/figures/figure4_gsea \
  --fdr 0.1 \
  --top-n 30
```

Panel mapping:

- Figure 4E: `DE_mock_vs_high/6_Infected_edgeRDEresults.xlsx`
- Figure 4F: `DE_mock_vs_low/5_Bystander_bkgd_edgeRDEresults.xlsx`
- Figure 4G: `DE_mock_vs_bystander/4_Bystander_edgeRDEresults.xlsx`
- Figure 4H: `7_Infected_edgeRDEresults.xlsx`

The `plot-saved` command is deterministic: it reads the GSEA sheets already saved in those workbooks and regenerates SVG/PNG panels plus pathway-audit TSVs.

Optional rerun from edgeR `Gene` and `logFC`:

```bash
python scripts/gsea/run_gsea_from_edger.py run-prerank \
  --input "/path/to/6_Infected_edgeRDEresults.xlsx" \
  --name mock_vs_highly_infected \
  --outdir outputs/gsea
```

# Figure scripts

Cleaned scripts:

- `plot_directional_gsea_panel.py`: plots the SciViewer directional-analysis GSEA NES bar panel from a Table S5/Table S7-style workbook.
- `plot_directional_gene_dotplots.py`: plots the top positively and negatively correlated gene dotplots from a Table S4/Table S6-style directional-correlation workbook and the all-cell AnnData object.

Example:

```bash
python scripts/figures/plot_directional_gsea_panel.py \
  --input outputs/tables/Table_S5_supplement_GSEA_directional_analysis_1.xlsx \
  --output outputs/figures/directional_analysis_1/directional_analysis_1_gsea_fdr_0.1.svg \
  --selected-output outputs/figures/directional_analysis_1/directional_analysis_1_gsea_fdr_0.1_plotted_pathways.tsv \
  --title "Directional Analysis 1 GSEA" \
  --fdr 0.1
```

```bash
python scripts/figures/plot_directional_gene_dotplots.py \
  --adata "/path/to/cmv.srt.soupx.filt.h5ad" \
  --correlations outputs/tables/Table_S4_supplement_directional_analysis_1_raw_results.xlsx \
  --output-dir outputs/figures/directional_analysis_1 \
  --file-prefix directional_analysis_1 \
  --analysis directional_analysis_1 \
  --use-raw \
  --average-method seurat
```

```bash
python scripts/figures/plot_directional_gene_dotplots.py \
  --adata "/path/to/cmv.srt.soupx.filt.h5ad" \
  --correlations outputs/tables/Table_S6_supplement_directional_analysis_2_raw_results.xlsx \
  --output-dir outputs/figures/directional_analysis_2 \
  --file-prefix directional_analysis_2 \
  --analysis directional_analysis_2 \
  --use-raw \
  --average-method seurat
```

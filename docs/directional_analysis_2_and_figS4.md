# Directional Analysis 2 and final Figure S4 provenance

The SCIViewer Directional Analysis 2 supplemental figure corresponds to final manuscript Figure S4. In the local manuscript figure folder, the associated image file was named:

`png_files/Fig_S3.png`

This note is included because earlier local figure filenames did not always match the final manuscript figure numbering.

Large input/intermediate files are available through GEO accession `GSE348416`.

## Exact supplemental-table sources

Directional Analysis 2 submitted supplemental tables were matched exactly to:

`results_12Jan25_dir1_3dpiinf_cells_hostonly.proj_correlation.xlsx`

Despite the source filename containing `dir1`, this workbook matches the submitted Directional Analysis 2 supplemental outputs:

- Table S6 matches the workbook's first sheet.
- Table S7 matches the workbook's `GSEA_Unfiltered_Curated_MSigDB_Hallmark_2025` sheet after dropping the source workbook's `Name` column.

Cleaned table script:

`../scripts/directional_analysis/reproduce_directional_analysis_2.py`

Example commands:

```bash
python scripts/directional_analysis/reproduce_directional_analysis_2.py table-s6 \
  --input /path/to/results_12Jan25_dir1_3dpiinf_cells_hostonly.proj_correlation.xlsx \
  --output outputs/tables/Table_S6_supplement_directional_analysis_2_raw_results.xlsx

python scripts/directional_analysis/reproduce_directional_analysis_2.py table-s7 \
  --input /path/to/results_12Jan25_dir1_3dpiinf_cells_hostonly.proj_correlation.xlsx \
  --output outputs/tables/Table_S7_supplement_GSEA_directional_analysis_2.xlsx
```

## Interactive SCIViewer notebook

Cleaned interactive notebook:

`../notebooks/reproduce_directional_analysis_2.ipynb`

An earlier version of the workflow loaded the all-cell object:

`cmv.srt.soupx.filt.updated.h5ad`

but the saved 12Jan25 selected-cell export and final Figure S4 are anchored to the 389-cell 3 dpi infected-cell subset. The saved selected-cell indices map to:

`infected_cells_3dpi_cmv_data.h5ad`

and also to:

`cmv.srt.soupx.filt.dpi3Infected.only_v3.h5ad`

Those selected-cell indices are specific to the infected-cell subset rather than the all-cell object. The notebook therefore opens `infected_cells_3dpi_cmv_data.h5ad` for the interactive SCIViewer step and uses:

```python
SCIViewer(adata, embedding_name="X_umap", use_raw=False)
```

The cleaned notebook preserves that interactive step and then lets the user choose whether to use the original saved export or a newly recreated export.

## Figure S4 C-D dotplots

The dotplots are generated with:

`../scripts/figures/plot_directional_gene_dotplots.py`

Example command:

```bash
python scripts/figures/plot_directional_gene_dotplots.py \
  --adata /path/to/cmv.srt.soupx.filt.h5ad \
  --correlations outputs/tables/Table_S6_supplement_directional_analysis_2_raw_results.xlsx \
  --output-dir outputs/figures/directional_analysis_2 \
  --file-prefix directional_analysis_2 \
  --analysis directional_analysis_2 \
  --use-raw \
  --average-method seurat
```

This uses the all-cell AnnData object because the manuscript dotplots are drawn across Mock, Bystander, Marginal Infection, and high-infection viral-transcript bins.

## Figure S4E GSEA panel

The Directional Analysis 2 GSEA NES panel is plotted with the same generic GSEA plotting script used for Directional Analysis 1. The final Figure S4 panel title says FDR <= 0.1, so use `--fdr 0.1` for the manuscript-matching version.

```bash
python scripts/figures/plot_directional_gsea_panel.py \
  --input outputs/tables/Table_S7_supplement_GSEA_directional_analysis_2.xlsx \
  --output outputs/figures/directional_analysis_2/directional_analysis_2_gsea_fdr_0.1.svg \
  --selected-output outputs/figures/directional_analysis_2/directional_analysis_2_gsea_fdr_0.1_plotted_pathways.tsv \
  --title "Directional Analysis 2 GSEA" \
  --fdr 0.1 \
  --top-n 30
```

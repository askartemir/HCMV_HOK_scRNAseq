# Directional Analysis 1 provenance

Directional Analysis 1 uses an interactive SCIViewer selection on the 3 dpi infected-cell AnnData object. The saved selected cells and projection-correlation output are treated as reproducibility inputs for regenerating the submitted supplemental tables and plotted panels.

Required input files:

- `infected_cells_3dpi_cmv_data.h5ad`: AnnData input; too large for GitHub, should be deposited in GEO/public data.
- `selected_cells_04Nov24_dir1_dpi3_infected.only.host.only.csv`: saved selected cells/projection for Directional Analysis 1.
- `results_04Nov24_dir1_dpi3_infected.only.host.only.proj_correlation.xlsx`: source projection-correlation export matching submitted Table S4 values.
- `04Nov24_dir1_dpi3_infected.only.host.only.filtgenes.fdr_10_to_genes_with_regulation_mapping.xlsx`: source workbook containing the GSEA sheet used for submitted Table S5.


Cleaned interactive notebook:

`../notebooks/reproduce_directional_analysis_1.ipynb`

Use this notebook to reopen SCIViewer, recreate the Directional Analysis 1 selection, export a new selected-cell/projection-correlation intermediate, and then run the deterministic downstream cells.

To run the interactive SCIViewer cell on this machine, select the Jupyter kernel named `Python (HCMV sciviewer)`, restart the notebook, and run from the top. The notebook sets `JAVA_HOME`, loads the `py5` extension, and enables the macOS GUI event loop before importing SCIViewer.

Cleaned script:

`../scripts/directional_analysis/reproduce_directional_analysis_1.py`

Example commands, assuming data files are placed locally outside Git-tracked paths:

```bash
python scripts/directional_analysis/reproduce_directional_analysis_1.py table-s4 \
  --input /path/to/results_04Nov24_dir1_dpi3_infected.only.host.only.proj_correlation.xlsx \
  --output outputs/tables/Table_S4_supplement_directional_analysis_1_raw_results.xlsx

python scripts/directional_analysis/reproduce_directional_analysis_1.py table-s5 \
  --input /path/to/04Nov24_dir1_dpi3_infected.only.host.only.filtgenes.fdr_10_to_genes_with_regulation_mapping.xlsx \
  --output outputs/tables/Table_S5_supplement_GSEA_directional_analysis_1.xlsx
```

The GSEA NES panel can then be plotted from Table S5:

```bash
python scripts/figures/plot_directional_gsea_panel.py \
  --input outputs/tables/Table_S5_supplement_GSEA_directional_analysis_1.xlsx \
  --output outputs/figures/directional_analysis_1/directional_analysis_1_gsea_fdr_0.1.svg \
  --selected-output outputs/figures/directional_analysis_1/directional_analysis_1_gsea_fdr_0.1_plotted_pathways.tsv \
  --title "Directional Analysis 1 GSEA" \
  --fdr 0.1
```

The top positively and negatively correlated gene dotplots for the Directional Analysis 1 manuscript panels can be plotted from the regenerated Table S4 workbook and the all-cell AnnData object:

```bash
python scripts/figures/plot_directional_gene_dotplots.py \
  --adata /path/to/cmv.srt.soupx.filt.h5ad \
  --correlations outputs/tables/Table_S4_supplement_directional_analysis_1_raw_results.xlsx \
  --output-dir outputs/figures/directional_analysis_1 \
  --file-prefix directional_analysis_1 \
  --use-raw \
  --average-method seurat
```

The dotplot script writes both PNG figures and audit workbooks containing the selected genes, group mean expression, z-scored average expression, and percent-expressing values used for plotting. The plotted groups are 3 dpi Mock, Bystander, Marginal Infection, and High Infection bins by percent HCMV transcripts. The `--use-raw --average-method seurat` settings mirror the Seurat/scCustomize dotplot convention used for the original R plotting code.

Large `.h5ad` and `.rds` files are expected to be supplied separately through GEO or an equivalent public data repository.

# HCMV HOK scRNA-seq reproducibility code

This repository contains Python code used to reproduce selected HCMV HOK scRNA-seq manuscript outputs: SCIViewer directional analyses, related supplemental tables, directional-analysis dotplots/GSEA panels, and Figure 4 GSEA panels derived from saved edgeR/GSEA workbooks.

## Scope

This repository is organized as a reproducibility package for manuscript-associated Python analyses. Large data objects are distributed through GEO accession `GSE348416`.

Reproduced outputs:

- Figure 4E-H GSEA panels
- Figure 5C-E directional-analysis panels
- Figure S4C-E directional-analysis panels
- Table S4
- Table S5
- Table S6
- Table S7

Upstream R analyses, edgeR outputs, and large R/Seurat objects are documented in `manifests/` and `docs/upstream_r_analysis.md`. The upstream R analysis code is maintained separately.

## Repository Layout

```text
notebooks/                      Clean interactive SCIViewer notebooks
scripts/directional_analysis/   Supplemental table scripts for Tables S4-S7
scripts/figures/                Directional-analysis dotplot and GSEA panel scripts
scripts/gsea/                   Figure 4 GSEA panel scripts
manifests/                      Data, code, and output provenance tables
docs/                           Method/provenance notes
data/                           No large data; see data/README.md
outputs/                        Generated outputs from local reruns
```

## Data Availability

Raw sequencing data, processed count matrices, Seurat objects, AnnData objects, and analysis intermediate files are available through NCBI GEO accession `GSE348416`. The GEO record remains private until publication; during peer review, access is provided through the private GEO reviewer token supplied in the journal submission system.

Large `.h5ad`, `.rds`, `.h5Seurat`, FASTQ, count matrix, and generated-output files are distributed through GEO or regenerated locally rather than stored in this repository. See:

- `manifests/h5ad_and_rds_manifest.tsv`
- `manifests/R_objects_manifest.tsv`
- `manifests/data_availability.tsv`

## Environment

Create the Python environment with either:

```bash
conda env create -f environment.yml
conda activate hcmv-hok-scrnaseq
```

or install from `requirements.txt` into a compatible Python environment.

The deterministic scripts above run without launching SCIViewer. They start from saved SCIViewer/edgeR outputs and regenerate the submitted tables and plotted panels.

SCIViewer is only required for manually recreating the interactive directional selections in `notebooks/reproduce_directional_analysis_1.ipynb` or `notebooks/reproduce_directional_analysis_2.ipynb`. For those notebooks, use a Jupyter kernel from this environment. SCIViewer uses `py5`, so Java 17 is required; the conda environment installs `openjdk=17`.

If registering a fresh kernel manually:

```bash
python -m ipykernel install --user --name hcmv-sciviewer --display-name "Python (HCMV sciviewer)"
```

SCIViewer/py5 platform notes:

- macOS: the notebooks enable the macOS GUI event loop before loading `py5`.
- Windows/Linux: the notebooks skip the macOS-only GUI command and load `py5` directly.
- If SCIViewer does not open, launch Jupyter from the activated conda environment and confirm that `JAVA_HOME` points to the Java installation in that environment.

## Reproduction Order

Run commands from the repository root. Replace `/path/to/...` with downloaded copies of the manuscript data files listed in `manifests/`.

1. Reproduce Directional Analysis 1 supplemental tables:

```bash
python scripts/directional_analysis/reproduce_directional_analysis_1.py table-s4 \
  --input /path/to/results_04Nov24_dir1_dpi3_infected.only.host.only.proj_correlation.xlsx \
  --output outputs/tables/Table_S4_supplement_directional_analysis_1_raw_results.xlsx

python scripts/directional_analysis/reproduce_directional_analysis_1.py table-s5 \
  --input /path/to/04Nov24_dir1_dpi3_infected.only.host.only.filtgenes.fdr_10_to_genes_with_regulation_mapping.xlsx \
  --output outputs/tables/Table_S5_supplement_GSEA_directional_analysis_1.xlsx
```

2. Reproduce Directional Analysis 2 supplemental tables:

```bash
python scripts/directional_analysis/reproduce_directional_analysis_2.py table-s6 \
  --input /path/to/results_12Jan25_dir1_3dpiinf_cells_hostonly.proj_correlation.xlsx \
  --output outputs/tables/Table_S6_supplement_directional_analysis_2_raw_results.xlsx

python scripts/directional_analysis/reproduce_directional_analysis_2.py table-s7 \
  --input /path/to/results_12Jan25_dir1_3dpiinf_cells_hostonly.proj_correlation.xlsx \
  --output outputs/tables/Table_S7_supplement_GSEA_directional_analysis_2.xlsx
```

3. Reproduce directional-analysis dotplots:

```bash
python scripts/figures/plot_directional_gene_dotplots.py \
  --adata /path/to/cmv.srt.soupx.filt.h5ad \
  --correlations outputs/tables/Table_S4_supplement_directional_analysis_1_raw_results.xlsx \
  --output-dir outputs/figures/directional_analysis_1 \
  --file-prefix directional_analysis_1 \
  --analysis directional_analysis_1 \
  --use-raw \
  --average-method seurat

python scripts/figures/plot_directional_gene_dotplots.py \
  --adata /path/to/cmv.srt.soupx.filt.h5ad \
  --correlations outputs/tables/Table_S6_supplement_directional_analysis_2_raw_results.xlsx \
  --output-dir outputs/figures/directional_analysis_2 \
  --file-prefix directional_analysis_2 \
  --analysis directional_analysis_2 \
  --use-raw \
  --average-method seurat
```

4. Reproduce directional-analysis GSEA panels:

```bash
python scripts/figures/plot_directional_gsea_panel.py \
  --input outputs/tables/Table_S5_supplement_GSEA_directional_analysis_1.xlsx \
  --output outputs/figures/directional_analysis_1/directional_analysis_1_gsea_fdr_0.1.svg \
  --selected-output outputs/figures/directional_analysis_1/directional_analysis_1_gsea_fdr_0.1_plotted_pathways.tsv \
  --title "Directional Analysis 1 GSEA" \
  --fdr 0.1

python scripts/figures/plot_directional_gsea_panel.py \
  --input outputs/tables/Table_S7_supplement_GSEA_directional_analysis_2.xlsx \
  --output outputs/figures/directional_analysis_2/directional_analysis_2_gsea_fdr_0.1.svg \
  --selected-output outputs/figures/directional_analysis_2/directional_analysis_2_gsea_fdr_0.1_plotted_pathways.tsv \
  --title "Directional Analysis 2 GSEA" \
  --fdr 0.1 \
  --top-n 30
```

5. Reproduce Figure 4 GSEA panels from saved edgeR/GSEA workbooks:

```bash
python scripts/gsea/run_gsea_from_edger.py plot-saved \
  --input-root "/path/to/h5ad and tsv files" \
  --outdir outputs/figures/figure4_gsea \
  --fdr 0.1 \
  --top-n 30
```

The interactive notebooks in `notebooks/` are included for the SCIViewer selection step. The deterministic scripts above start from saved SCIViewer/edgeR outputs and regenerate the manuscript tables and plotted panels.

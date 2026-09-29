# Data directory

This directory is a placeholder for optional local data copies.

Expected large inputs are available through GEO accession `GSE348416` and are listed in `../manifests/data_availability.tsv`, `../manifests/h5ad_and_rds_manifest.tsv`, and `../manifests/R_objects_manifest.tsv`. For reruns, place downloaded data files here or set `HCMV_DATA_DIR` to a directory containing the required files.

Files downloaded directly from GEO are prefixed with `GSE348416_`. The scripts and notebooks use the original project filenames without that prefix. After downloading GEO files into this directory, normalize the filenames from the repository root:

```bash
python scripts/prepare_geo_files.py data
```

Preview the renames first with:

```bash
python scripts/prepare_geo_files.py data --dry-run
```

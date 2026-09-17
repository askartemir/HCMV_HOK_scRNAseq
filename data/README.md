# Data directory

This directory is intentionally kept free of large data files.

Expected large inputs are listed in `../manifests/data_availability.tsv`, `../manifests/h5ad_and_rds_manifest.tsv`, and `../manifests/R_objects_manifest.tsv`. For local reruns, place downloaded data files here or set `HCMV_DATA_DIR` to a directory containing the required files. Do not commit large data files to GitHub.

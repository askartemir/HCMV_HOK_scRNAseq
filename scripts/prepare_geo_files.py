#!/usr/bin/env python3
"""Normalize GEO-downloaded filenames for local reproducibility runs.

NCBI GEO prepends the accession to supplementary file names, for example:

    GSE348416_cmv.srt.soupx.filt.h5ad

The analysis scripts and notebooks in this repository use the original project
filenames without that accession prefix. This helper removes the prefix in a
local data directory so the README commands work directly after download.
"""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_PREFIX = "GSE348416_"


def iter_targets(data_dir: Path, prefix: str, recursive: bool) -> list[Path]:
    pattern = f"{prefix}*"
    iterator = data_dir.rglob(pattern) if recursive else data_dir.glob(pattern)
    return sorted(path for path in iterator if path.name.startswith(prefix))


def strip_prefix(path: Path, prefix: str) -> Path:
    return path.with_name(path.name.removeprefix(prefix))


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Remove the GEO accession prefix from files downloaded for this "
            "repository's reproduction workflows."
        )
    )
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="data",
        type=Path,
        help="Directory containing GEO-downloaded files. Defaults to ./data.",
    )
    parser.add_argument(
        "--prefix",
        default=DEFAULT_PREFIX,
        help=f"Filename prefix to remove. Defaults to {DEFAULT_PREFIX!r}.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Also rename matching files in subdirectories.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned renames without changing files.",
    )
    args = parser.parse_args()

    data_dir = args.data_dir.expanduser().resolve()
    if not data_dir.exists():
        parser.error(f"Data directory does not exist: {data_dir}")
    if not data_dir.is_dir():
        parser.error(f"Data path is not a directory: {data_dir}")

    targets = iter_targets(data_dir, args.prefix, args.recursive)
    if not targets:
        print(f"No files found with prefix {args.prefix!r} in {data_dir}")
        return 0

    rename_pairs: list[tuple[Path, Path]] = []
    blocked: list[tuple[Path, Path]] = []
    for source in targets:
        destination = strip_prefix(source, args.prefix)
        if destination == source:
            continue
        if destination.exists():
            blocked.append((source, destination))
        else:
            rename_pairs.append((source, destination))

    if blocked:
        print("Refusing to overwrite existing files:")
        for source, destination in blocked:
            print(f"  {source} -> {destination}")
        print("Move or remove the destination files, then rerun this command.")
        return 1

    for source, destination in rename_pairs:
        print(f"{source} -> {destination}")
        if not args.dry_run:
            source.rename(destination)

    if args.dry_run:
        print(f"Dry run complete: {len(rename_pairs)} rename(s) would be applied.")
    else:
        print(f"Renamed {len(rename_pairs)} file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

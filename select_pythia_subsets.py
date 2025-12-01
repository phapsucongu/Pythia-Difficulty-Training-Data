import argparse
import os
from pathlib import Path
import shutil
import sys


def read_verbose_names(parquet_path: Path, limit: int | None = None):
    """Read the `verbose_name` column from a parquet file.

    Uses pandas if available, otherwise falls back to pyarrow.
    Returns a list of unique names as strings.
    """

    try:
        import pandas as pd  # type: ignore

        df = pd.read_parquet(parquet_path)
        if "verbose_name" not in df.columns:
            raise KeyError("Column 'verbose_name' not found in parquet file")
        series = df["verbose_name"].astype(str)
        if limit is not None:
            series = series.head(limit)
        return sorted(set(series.tolist()))
    except ImportError:
        try:
            import pyarrow.parquet as pq  # type: ignore

            table = pq.read_table(parquet_path, columns=["verbose_name"])
            col = table.column("verbose_name")
            values = [str(v) for v in col.to_pylist()]
            if limit is not None:
                values = values[:limit]
            return sorted(set(values))
        except ImportError:
            raise SystemExit(
                "Neither pandas nor pyarrow is installed. Please install one of them (e.g., 'pip install pandas' or 'pip install pyarrow')."
            )


def copy_subset_names(names, output_root: Path, pythia_root: Path, dry_run: bool = False):
    """Copy folders from output_root to pythia_root based on verbose_name list.

    For each name N in names, copies directory `output_root / N` to `pythia_root / N`.
    Returns statistics dict.
    """

    pythia_root.mkdir(parents=True, exist_ok=True)

    stats = {"total": 0, "copied": 0, "missing": 0, "skipped": 0}
    for name in names:
        stats["total"] += 1
        src = output_root / name
        dst = pythia_root / name

        if not src.exists() or not src.is_dir():
            print(f"MISSING: {src}")
            stats["missing"] += 1
            continue

        if dst.exists():
            print(f"SKIP   (already exists): {dst}")
            stats["skipped"] += 1
            continue

        if dry_run:
            print(f"DRY    would copy {src} -> {dst}")
            stats["copied"] += 1
        else:
            print(f"COPY   {src} -> {dst}")
            shutil.copytree(src, dst)
            stats["copied"] += 1

    return stats


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Read 'training_data.parquet', take the 'verbose_name' column, and "
            "copy matching subfolders from an 'output' directory into a new 'pythia' directory."
        )
    )
    parser.add_argument(
        "--parquet",
        default="training_data.parquet",
        help="Path to training_data.parquet (default: training_data.parquet in current directory)",
    )
    parser.add_argument(
        "--output-root",
        default="output",
        help="Root directory containing extracted subfolders (default: output)",
    )
    parser.add_argument(
        "--pythia-root",
        default="pythia",
        help="Destination root directory for selected subfolders (default: pythia)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optionally limit to the first N verbose_name entries (useful for testing)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be copied without actually copying",
    )

    args = parser.parse_args(argv)

    parquet_path = Path(args.parquet).resolve()
    if not parquet_path.exists():
        print(f"Parquet file not found: {parquet_path}")
        return 2

    output_root = Path(args.output_root).resolve()
    if not output_root.exists() or not output_root.is_dir():
        print(f"Output root directory not found or not a directory: {output_root}")
        return 2

    pythia_root = Path(args.pythia_root).resolve()

    print(f"Reading verbose_name from: {parquet_path}")
    names = read_verbose_names(parquet_path, limit=args.limit)
    print(f"Found {len(names)} unique verbose_name values")

    stats = copy_subset_names(names, output_root, pythia_root, dry_run=args.dry_run)
    print(
        "Summary: total_names={total} copied_or_planned={copied} "
        "missing={missing} skipped_existing={skipped}".format(**stats)
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

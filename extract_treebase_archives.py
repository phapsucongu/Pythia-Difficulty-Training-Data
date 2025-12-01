import argparse
import os
from pathlib import Path
import tarfile
import sys


def find_archives(input_dir: Path):
    """
    Yield all .tar.gz files under the input directory recursively.
    """
    yield from input_dir.rglob("*.tar.gz")


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def is_within_directory(directory: Path, target: Path) -> bool:
    try:
        directory_resolved = directory.resolve(strict=False)
        target_resolved = target.resolve(strict=False)
    except Exception:
        # Fallback to absolute path compare if resolution fails
        directory_resolved = Path(os.path.abspath(str(directory)))
        target_resolved = Path(os.path.abspath(str(target)))
    dir_str = str(directory_resolved)
    tgt_str = str(target_resolved)
    if not dir_str.endswith(os.sep):
        dir_str = dir_str + os.sep
    return tgt_str.startswith(dir_str)


def safe_extract(tar: tarfile.TarFile, path: Path):
    """
    Safely extract tar members, preventing path traversal outside of `path`.
    """
    for member in tar.getmembers():
        member_path = path / member.name
        if not is_within_directory(path, member_path):
            raise RuntimeError(
                f"Unsafe path detected in archive: {member.name} -> {member_path}"
            )
    tar.extractall(str(path))


def extract_archive(archive: Path, dest_root: Path, flatten: bool, force: bool, dry_run: bool) -> str:
    """
    Extract a single .tar.gz archive.
    - When not flattening, extracts into dest_root/<archive.parent.name>
    - When flattening, extracts directly into dest_root
    Returns a status string.
    """
    if flatten:
        dest = dest_root
    else:
        # Put each archive into its own folder, named by its containing folder (e.g., 115_0.phy)
        dest = dest_root / archive.parent.name

    ensure_dir(dest)

    # If not forcing, and destination already contains files, skip to avoid duplicates
    if not force and any(dest.iterdir()):
        return f"SKIP  already extracted? {archive} -> {dest}"

    if dry_run:
        return f"DRY   would extract {archive} -> {dest}"

    try:
        with tarfile.open(archive, mode="r:gz") as tf:
            safe_extract(tf, dest)
        return f"DONE  extracted {archive} -> {dest}"
    except tarfile.ReadError as e:
        return f"FAIL  not a valid tar.gz: {archive} ({e})"
    except Exception as e:
        return f"FAIL  error extracting {archive}: {e}"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract all .tar.gz archives under a TreeBASE 'trees' folder into an 'output' folder."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="trees",
        help="Input directory to scan for .tar.gz (default: trees)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="Destination output directory (default: output)",
    )
    parser.add_argument(
        "--flatten",
        action="store_true",
        help="Extract all files directly into output (no per-archive subfolders).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-extraction even if destination contains files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print planned actions, do not extract.",
    )

    args = parser.parse_args(argv)

    input_dir = Path(args.input).resolve()
    output_dir = Path(args.output).resolve()
    ensure_dir(output_dir)

    if not input_dir.exists() or not input_dir.is_dir():
        print(f"Input directory not found or not a directory: {input_dir}")
        return 2

    archives = sorted(find_archives(input_dir))
    if not archives:
        print(f"No .tar.gz archives found under: {input_dir}")
        return 0

    print(f"Found {len(archives)} archives under {input_dir}")
    print(f"Output directory: {output_dir}")
    if args.flatten:
        print("Mode: flatten (extract all files directly into output)")
    else:
        print("Mode: per-archive subfolders (output/<parent.name>)")
    if args.force:
        print("Option: force re-extraction enabled")
    if args.dry_run:
        print("Option: dry-run (no changes)")

    count_done = count_skip = count_fail = 0
    for arc in archives:
        status = extract_archive(arc, output_dir, args.flatten, args.force, args.dry_run)
        print(status)
        if status.startswith("DONE") or status.startswith("DRY"):
            count_done += 1
        elif status.startswith("SKIP"):
            count_skip += 1
        else:
            count_fail += 1

    print(
        f"Summary: processed={len(archives)} ok_or_planned={count_done} skipped={count_skip} failed={count_fail}"
    )
    return 0 if count_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())

import argparse
from pathlib import Path
import shutil
import sys


def folder_size(path: Path) -> int:
    """Compute total size (in bytes) of all files under path."""

    total = 0
    for p in path.rglob("*"):
        if p.is_file():
            try:
                total += p.stat().st_size
            except OSError:
                continue
    return total


def list_pythia_folders(pythia_root: Path):
    """Return list of immediate subdirectories under pythia_root."""

    return [p for p in pythia_root.iterdir() if p.is_dir()]


def find_one_msa_file(folder: Path):
    """Find a single MSA file in folder.

    Heuristic: prefer files with extensions commonly used for MSAs.
    Falls back to the first regular file if none match the preferred list.
    """

    preferred_exts = {".phy", ".phylip", ".fasta", ".fa", ".aln"}
    candidates_pref = []
    candidates_any = []

    for p in folder.rglob("*"):
        if not p.is_file():
            continue
        candidates_any.append(p)
        if p.suffix.lower() in preferred_exts:
            candidates_pref.append(p)

    if candidates_pref:
        return sorted(candidates_pref)[0]
    if candidates_any:
        return sorted(candidates_any)[0]
    return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Sort folders in 'pythia' by size, then with a given step "
            "(e.g., 100) select folders and copy one MSA file from each "
            "into 'test_input'."
        )
    )
    parser.add_argument(
        "--pythia-root",
        default="pythia",
        help="Root directory containing pythia subfolders (default: pythia)",
    )
    parser.add_argument(
        "--test-root",
        default="test_input",
        help="Directory where selected MSA files will be copied (default: test_input)",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=100,
        help="Step size when picking folders from sorted list (default: 100)",
    )
    parser.add_argument(
        "--max-folders",
        type=int,
        default=None,
        help="Optional cap on number of folders to process (after stepping)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be copied without actually copying",
    )

    args = parser.parse_args(argv)

    pythia_root = Path(args.pythia_root).resolve()
    if not pythia_root.exists() or not pythia_root.is_dir():
        print(f"Pythia root directory not found or not a directory: {pythia_root}")
        return 2

    test_root = Path(args.test_root).resolve()
    test_root.mkdir(parents=True, exist_ok=True)

    folders = list_pythia_folders(pythia_root)
    if not folders:
        print(f"No subfolders found under: {pythia_root}")
        return 0

    # Compute sizes
    print(f"Computing sizes for {len(folders)} folders under {pythia_root} ...")
    sized = []
    for f in folders:
        sz = folder_size(f)
        sized.append((sz, f))
    sized.sort(key=lambda x: x[0])  # ascending by size

    print("Example sizes (smallest 5):")
    for sz, f in sized[:5]:
        print(f"  {f.name}: {sz} bytes")

    # Step selection
    step = max(1, args.step)
    selected = sized[::step]
    if args.max_folders is not None:
        selected = selected[: args.max_folders]

    print(
        f"Selected {len(selected)} folders with step={step} "
        f"(from total {len(sized)})"
    )

    copied = 0
    missing_msa = 0
    for sz, folder in selected:
        msa = find_one_msa_file(folder)
        if msa is None:
            print(f"MISSING_MSA in {folder}")
            missing_msa += 1
            continue

        dst = test_root / msa.name
        if dst.exists():
            # Avoid overwriting: prepend folder name
            dst = test_root / f"{folder.name}__{msa.name}"

        if args.dry_run:
            print(f"DRY  would copy {msa} -> {dst}")
        else:
            print(f"COPY {msa} -> {dst}")
            shutil.copy2(msa, dst)
        copied += 1

    print(
        f"Summary: folders_selected={len(selected)} copied_files={copied} "
        f"folders_without_msa={missing_msa}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

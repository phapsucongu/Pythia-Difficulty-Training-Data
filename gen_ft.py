import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from collections import Counter

# Load MSA file (FASTA format)
def load_msa(filepath, format=None):
    if format is None:
        ext = filepath.lower().rsplit('.', 1)[-1]
        if ext in {"fasta", "fa"}:
            format = "fasta"
        elif ext in {"phy", "phylip"}:
            format = "phylip-relaxed"
        elif ext == "aln":
            format = "clustal"
        else:
            format = "fasta"
    try:
        alignment = AlignIO.read(filepath, format)
        return alignment
    except Exception as e:
        print(f"Failed to load {filepath} (format={format}): {e}")
        return None

# % parsimony-informative and variable sites
def compute_variable_and_parsimony_sites(alignment):
    n_sites = alignment.get_alignment_length()
    msa_array = np.array([list(rec.seq) for rec in alignment], dtype=str)
    pis = 0
    variable = 0
    for col in msa_array.T:
        states, counts = np.unique(col, return_counts=True)
        if len(states) > 1:
            variable += 1
            if np.sum(counts >= 2) >= 2:
                pis += 1
    return 100 * pis / n_sites, 100 * variable / n_sites

# RCV
def compute_rcv(msa_array, alphabet):
    total_taxa, total_sites = msa_array.shape
    counts = {base: np.sum(msa_array == base) for base in alphabet}
    total_counts = sum(counts.values())
    freqs = {base: counts[base] / total_counts for base in alphabet}
    rcv = 0
    for row in msa_array:
        counts_t = {base: np.sum(row == base) for base in alphabet}
        freqs_t = {base: counts_t[base] / total_sites for base in alphabet}
        rcv += sum(abs(freqs_t[b] - freqs[b]) for b in alphabet)
    return rcv / total_taxa * 100

# RCVT
def compute_rcvt(msa_array, alphabet):
    total_sites = msa_array.shape[1]
    freqs_all = {b: np.sum(msa_array == b) / msa_array.size for b in alphabet}
    rcvt_list = []
    for row in msa_array:
        freqs_t = {b: np.sum(row == b) / total_sites for b in alphabet}
        rcvt = sum(abs(freqs_t[b] - freqs_all[b]) for b in alphabet)
        rcvt_list.append(rcvt)
    return np.mean(rcvt_list) * 100

# Treeness from UPGMA tree
def compute_treeness(alignment):
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    total_branch_length = tree.total_branch_length()
    internal_branch_length = sum(
        clade.branch_length for clade in tree.get_nonterminals() if clade.branch_length
    )
    if total_branch_length == 0:
        return 0.0
    return internal_branch_length / total_branch_length * 100

# Long Branch Score: std of terminal branch lengths
def compute_lb_score(tree):
    tips = tree.get_terminals()
    lengths = [t.branch_length for t in tips if t.branch_length is not None]
    return np.std(lengths) if lengths else np.nan

# Column entropy standard deviation
def compute_column_entropy_std(msa_array):
    from scipy.stats import entropy
    ce_list = []
    for col in msa_array.T:
        vals, counts = np.unique(col, return_counts=True)
        probs = counts / counts.sum()
        ce_list.append(entropy(probs, base=2))
    return float(np.std(ce_list)) if ce_list else np.nan

# Main function to process a list of MSA files
def analyze_msa_list(msa_files, output_csv="msa_features_output.csv", alphabet="ACGT"):
    import csv

    fieldnames = [
        "filename",
        "% parsimony-informative",
        "% variable",
        "RCV",
        "RCVT",
        "Treeness",
        "LB score (std tip lengths)",
        "degree_dispersion_missing",
        "average_amount_actual",
        "freq_A",
        "freq_C",
        "freq_G",
        "freq_T",
        "ce_std",
    ]

    write_header = True
    try:
        if os.path.exists(output_csv) and os.path.getsize(output_csv) > 0:
            write_header = False
    except Exception:
        write_header = True

    with open(output_csv, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()

        rows = []
        for filepath in msa_files:
            alignment = load_msa(filepath)
            if alignment is None:
                continue
            msa_array = np.array([list(rec.seq.upper()) for rec in alignment], dtype=str)
            pis_pct, var_pct = compute_variable_and_parsimony_sites(alignment)
            rcv = compute_rcv(msa_array, alphabet)
            rcvt = compute_rcvt(msa_array, alphabet)
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            treeness = compute_treeness(alignment)
            lb_score = compute_lb_score(tree)

            num_sites = msa_array.shape[1]
            num_taxa = msa_array.shape[0]
            prop_gaps = np.sum(msa_array == '-') / msa_array.size
            prop_inv = 1 - (var_pct / 100)
            degree_dispersion_missing = prop_gaps * num_taxa
            average_amount_actual = (1 - prop_gaps) * (1 - prop_inv) * num_sites

            flat = msa_array.flatten()
            freq_a = np.sum(flat == 'A') / flat.size * 100
            freq_c = np.sum(flat == 'C') / flat.size * 100
            freq_g = np.sum(flat == 'G') / flat.size * 100
            freq_t = np.sum(flat == 'T') / flat.size * 100

            ce_std = compute_column_entropy_std(msa_array)

            print(
                f"Processed {filepath}: PI%={pis_pct:.2f}, Var%={var_pct:.2f}, RCV={rcv:.2f}, RCVT={rcvt:.2f}, Treeness={treeness:.2f}, LB Score={lb_score:.2f}, ce_std={ce_std:.4f}"
            )
            row = {
                "filename": os.path.basename(filepath),
                "% parsimony-informative": pis_pct,
                "% variable": var_pct,
                "RCV": rcv,
                "RCVT": rcvt,
                "Treeness": treeness,
                "LB score (std tip lengths)": lb_score,
                "degree_dispersion_missing": degree_dispersion_missing,
                "average_amount_actual": average_amount_actual,
                "freq_A": freq_a,
                "freq_C": freq_c,
                "freq_G": freq_g,
                "freq_T": freq_t,
                "ce_std": ce_std,
            }
            writer.writerow(row)
            f.flush()
            rows.append(row)
    df = pd.DataFrame(rows)
    return df



def find_msa_files(root_dir: Path, extensions):
    """Recursively find MSA files under root_dir matching given extensions set."""
    for p in root_dir.rglob("*"):
        if p.is_file() and p.suffix.lower().lstrip('.') in extensions:
            yield str(p)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Batch compute MSA-derived features for all alignment files under a directory."
    )
    # Re-added missing arguments
    parser.add_argument(
        "--input-dir",
        default="output",
        help="Root directory containing MSA files (default: output)",
    )
    parser.add_argument(
        "--extensions",
        default="fasta,fa,phy,phylip,aln",
        help="Comma-separated list of allowed file extensions (default: fasta,fa,phy,phylip,aln)",
    )
    parser.add_argument(
        "--alphabet",
        default="ACGT",
        help="Alphabet characters relevant for RCV/RCVT (default: ACGT)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of MSA files processed (for testing)",
    )
    parser.add_argument(
        "--output-csv",
        default="msa_features_output.csv",
        help="Output CSV file path (default: msa_features_output.csv)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List files that would be processed without computing features",
    )
    args = parser.parse_args(argv)

    root = Path(args.input_dir).resolve()
    if not root.exists() or not root.is_dir():
        print(f"Input directory not found or not a directory: {root}")
        return 2

    exts = {e.strip().lower() for e in args.extensions.split(',') if e.strip()}
    files = list(find_msa_files(root, exts))
    if args.limit is not None:
        files = files[: args.limit]

    if not files:
        print(f"No MSA files found under {root} with extensions: {', '.join(sorted(exts))}")
        return 0

    print(f"Found {len(files)} MSA files under {root}")
    if args.dry_run:
        for f in files:
            print(f"DRY {f}")
        print("Dry-run complete. No computations performed.")
        return 0

    df = analyze_msa_list(files, output_csv=args.output_csv, alphabet=args.alphabet)
    print(f"Wrote features for {len(df)} alignments to {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

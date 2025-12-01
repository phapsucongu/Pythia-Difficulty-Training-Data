import os
import argparse
from pathlib import Path
from collections import Counter, defaultdict
from typing import List

import numpy as np
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

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

def generate_parsimony_trees(alignment, n_trees=24):
    calculator = DistanceCalculator("identity")
    trees = []
    for _ in range(n_trees):
        seqs = list(alignment)
        import random
        random.shuffle(seqs)
        sub_align = AlignIO.MultipleSeqAlignment(seqs)
        dm = calculator.get_distance(sub_align)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        trees.append(tree)
    return trees

def compute_leaf_instability(trees, k=3):
    taxa = [tip.name for tip in trees[0].get_terminals()]
    neighborhoods = defaultdict(list)
    not_found = set()
    for tree in trees:
        for taxon in taxa:
            node = tree.find_any(name=taxon)
            if node is None:
                not_found.add(taxon)
                neighborhoods[taxon].append(())
                continue
            nearby = node.get_terminals()
            close = sorted([t.name for t in nearby if t.name != taxon])[:k]
            neighborhoods[taxon].append(tuple(close))
    if not_found:
        print(f"Warning: could not find taxon(s) in some trees: {sorted(not_found)}")
    instability_scores = {}
    for taxon, neigh_list in neighborhoods.items():
        counts = Counter(neigh_list)
        most_common = counts.most_common(1)[0][1] if counts else 0
        instability = 1 - (most_common / len(neigh_list))
        instability_scores[taxon] = instability
    all_scores = list(instability_scores.values())
    return {
        "rogue_mean": float(np.mean(all_scores)) if all_scores else np.nan,
        "rogue_max": float(np.max(all_scores)) if all_scores else np.nan,
        "rogue_taxa_>0.5": int(sum(s > 0.5 for s in all_scores))
    }

def compute_consistency_index(alignment, tree):
    n_sites = alignment.get_alignment_length()
    if n_sites == 0:
        return np.nan
    seq_dict = {rec.id: str(rec.seq) for rec in alignment}
    import io
    s = tree.format("newick")
    root_tree = Phylo.read(io.StringIO(s), "newick")

    def fitch_pass(node, site):
        if node.is_terminal():
            name = node.name
            char = seq_dict.get(name)
            if char is None:
                for key in seq_dict:
                    if key.startswith(name):
                        char = seq_dict[key]
                        break
                else:
                    return {"?"}, 0
            return {char[site]}, 0
        else:
            sets = []
            changes = 0
            for child in node.clades:
                s_child, c_child = fitch_pass(child, site)
                sets.append(s_child)
                changes += c_child
            inter = set.intersection(*sets) if sets else set()
            if inter:
                return inter, changes
            else:
                union = set.union(*sets)
                return union, changes + 1

    total_min = 0
    total_actual = 0
    for site in range(n_sites):
        site_chars = [seq_dict[t][site] for t in seq_dict]
        unique_states = set(site_chars)
        total_min += max(len(unique_states) - 1, 0)
        _, changes = fitch_pass(root_tree.clade, site)
        total_actual += changes
    if total_actual == 0:
        return np.nan
    return float(total_min) / float(total_actual)

def write_trees_to_file(trees: List[Phylo.BaseTree.Tree], out_path: str):
    # Mỗi cây trên 1 dòng, dòng tiếp theo là metadata
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        for i, tree in enumerate(trees):
            newick_str = tree.format("newick").replace("\n", "")
            f.write(newick_str + "\n")
            f.write(f"[&R] TREE_{i}\n")

def find_msa_files(root_dir: Path, extensions: set) -> List[str]:
    files = []
    for p in root_dir.rglob("*"):
        if p.is_file() and p.suffix.lower().lstrip('.') in extensions:
            files.append(str(p))
    return files

def process_msa_for_rogue_ci(msa_files: List[str], input_dir: Path, tree_dir: str, csv_out: str):
    import csv
    os.makedirs(tree_dir, exist_ok=True)
    fieldnames = ["filename", "rogue_mean", "rogue_max", "rogue_taxa_>0.5", "CI"]
    write_header = not os.path.exists(csv_out) or os.path.getsize(csv_out) == 0
    with open(csv_out, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        for filepath in msa_files:
            alignment = load_msa(filepath)
            if alignment is None:
                continue
            trees = generate_parsimony_trees(alignment, n_trees=24)
            rogue_stats = compute_leaf_instability(trees)
            ci_value = compute_consistency_index(alignment, trees[0])
            relpath = Path(filepath).relative_to(input_dir)
            base_name = relpath.as_posix().replace("/", "_").replace(".", "_")
            tree_file_name = f"{base_name}_trees.nwk"
            tree_file = os.path.join(tree_dir, tree_file_name)
            write_trees_to_file(trees, tree_file)
            writer.writerow({
                "filename": str(relpath),
                "rogue_mean": rogue_stats["rogue_mean"],
                "rogue_max": rogue_stats["rogue_max"],
                "rogue_taxa_>0.5": rogue_stats["rogue_taxa_>0.5"],
                "CI": ci_value
            })
            print(f"Processed {filepath} → CI={ci_value:.3f}, rogue_mean={rogue_stats['rogue_mean']:.3f}")

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Compute rogue scores and CI from alignments.")
    parser.add_argument("--input-dir", default="output", help="Directory containing MSA files")
    parser.add_argument("--extensions", default="fasta,fa,phy,phylip,aln", help="Allowed extensions")
    parser.add_argument("--output-csv", default="rogue_ci_summary.csv", help="Output CSV file")
    parser.add_argument("--tree-dir", default="output_pars_tree", help="Directory to save tree files")
    parser.add_argument("--limit", type=int, default=None, help="Optional limit on number of files")
    parser.add_argument("--dry-run", action="store_true", help="List files without running")
    args = parser.parse_args(argv)
    root = Path(args.input_dir).resolve()
    if not root.exists() or not root.is_dir():
        print(f"Input directory not found: {root}")
        return 2
    exts = {e.strip().lower() for e in args.extensions.split(',') if e.strip()}
    files = find_msa_files(root, exts)
    if args.limit is not None:
        files = files[: args.limit]
    if not files:
        print(f"No MSA files found under {root}")
        return 0
    print(f"Found {len(files)} MSA files under {root}")
    if args.dry_run:
        for f in files:
            print(f"DRY {f}")
        print("Dry-run complete.")
        return 0
    process_msa_for_rogue_ci(files, input_dir=root, tree_dir=args.tree_dir, csv_out=args.output_csv)
    print(f"Wrote output to {args.output_csv}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

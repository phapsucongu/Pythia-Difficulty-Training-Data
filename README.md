# TreeBASEMirror

This repo contains MSAs from [TreeBASE](https://www.treebase.org/treebase-web/home.html) which we downloaded for our experiments, 
and - for most MSAs - the best trees we found after 100 [RAxML-NG](https://github.com/amkozlov/raxml-ng) tree searches under the "GTR+G" model for DNA and "LG" model for AA MSAs (with their respective inferred 
model parameters). We substituted the original taxon names since some of them were quite creative and made the parsers in the software we had to use
for our experiments cry. 

There is also a [SQLite](https://www.sqlite.org/index.html) database with our example scripts [here](https://github.com/angtft/RAxMLGroveScripts), that can be used for easier lookup of MSAs.


## Quick Start

- Install Python 3.10+ and dependencies:

```cmd
cd /d d:\test\TreeBASEMirror-main
python -m pip install -r requirements.txt
```

- Extract any `.tar.gz` archives from `trees/` into `output/` (optional):

```cmd
python extract_treebase_archives.py --input trees --output output
```

## Compute MSA Features (gen_ft.py)

`gen_ft.py` scans a directory for alignments and computes per-alignment features, appending results to a CSV as each file is processed.

- Supported formats: `fasta, fa, phy, phylip, aln`
- Output CSV default: `msa_features_output.csv`
- Writes each row immediately (safe on long runs)

Run examples:

```cmd
python gen_ft.py --input-dir output
```

Specify options:

```cmd
python gen_ft.py --input-dir output ^
	--extensions fasta,fa,phy,phylip,aln ^
	--alphabet ACGT ^
	--limit 100 ^
	--output-csv msa_features_output.csv
```


## Rogue Scores and CI (gen_tree.py)

`gen_tree.py` generates a set of proxy parsimony trees per alignment, computes simple rogue/instability statistics and a Fitch-consistency index (CI), and saves the trees.

- Input directory default: `output`
- Output CSV default: `rogue_ci_summary.csv`
- Trees directory default: `output_pars_tree`

Run examples:

```cmd
python gen_tree.py --input-dir output
```

With options:

```cmd
python gen_tree.py --input-dir output ^
	--extensions fasta,fa,phy,phylip,aln ^
	--output-csv rogue_ci_summary.csv ^
	--tree-dir output_pars_tree ^
	--limit 50


## Notes

- Outputs and large data files are ignored via `.gitignore` (e.g., `output/`, `trees/`, `*.csv`, alignment formats). Commit only source code.
- If you run into missing packages, ensure `requirements.txt` is installed.

## References
* Piel, W. H., Chan, L., Dominus, M. J., Ruan, J., Vos, R. A., and V. Tannen (2009)
**TreeBASE v. 2: A Database of Phylogenetic Knowledge.**
*e-BioSphere 2009*.

* Vos, R. A., Balhoff, J. P., Caravas, J. A., Holder, M. T., Lapp, H. Maddison, W. P., Midford, P. E., Priyam, A., Sukumaran, J., Xia, X. and A. Stoltzfus (2012)
**NeXML: rich, extensible, and verifiable representation of comparative data and metadata.**
*Systematic Biology* 61(4): 675-689. 

* Alexey M. Kozlov, Diego Darriba, Tom&aacute;&scaron; Flouri, Benoit Morel, and Alexandros Stamatakis (2019)
**RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference.** 
*Bioinformatics*, 35 (21), 4453-4455 
doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)

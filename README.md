# IS–determinant proximity analysis

A reproducible command-line workflow for describing the genomic proximity of insertion sequences (ISs) to antimicrobial resistance (AMR), virulence-factor (VF), and metal-resistance (METRES) determinants on plasmids.

> **Research software / pre-release.** This repository is a publication-preparation package.

## Features

- Nearest annotated IS for each determinant; distances of genomic intervals
- Fractions within user-defined proximity thresholds, with both all-determinant and IS-positive-plasmid denominators
- All determinant–IS pairs within a specified window, recurrent pairs and exact distances
- Candidate IS–determinant–IS structures, including orientation-independent gene-order labels (not strand-aware)
- Optional circular-origin-aware **proximity** and within-plasmid positional permutations
- CSV outputs and JSON run metadata

## Requirements

Python 3.9+; install with:

```bash
python -m pip install -r requirements.txt
```

## Quick start

Run the small **synthetic** example bundled with this repository:

```bash
python IS_proximity_publication.py \
  --analysis AMR \
  --target examples/amr.xlsx \
  --is-file examples/is.csv \
  --lengths examples/lengths.csv \
  --circular --permutations 1000 --seed 42 \
  --replicon Demo --output results/demo
```

To analyze real datasets, replace the example files with your own. Use `python IS_proximity_publication.py --help` for all parameters.

## Input formats

**All coordinates are 1-based, inclusive, and non-wrapping.** Sequence IDs must match exactly across files.

| Analysis | Input | Expected columns |
|---|---|---|
| AMR | Excel/CSV with header | `Contig id`, `Gene symbol`, `Start`, `Stop` |
| VF | Excel/CSV with header | `ID`, `Unnamed: 7`, `Unnamed: 8`, `Unnamed: 14` (description includes a gene name in parentheses) |
| METRES | Headerless CSV/Excel | Positional columns 0, 1, 7, 8, 14; gene label derived from hit name before first underscore |
| IS | Headerless CSV/Excel | Positional columns 0, 1, 7, 8; IS label is before first `|`; names filtered with `^IS` by default |
| Plasmid lengths | CSV with header | `Contig id`, `Length_bp` |

`--lengths` is required for `--circular` or `--permutations`. These legacy input schemas reflect the source annotation exports; confirm that your files conform before analysis. See [data dictionary](docs/DATA_DICTIONARY.md).

## Outputs

Nine CSV tables (`global_summary`, `nearest_IS`, `IS_within_window`, `IS_pairs`, `IS_recurrent_distances`, `candidate_modules_all`, `candidate_module_families`, `candidate_same_IS_families`, `permutation_null`) plus a `run_metadata.json` file. Empty tables are still written. Output filenames are prefixed by the replicon label and analysis type.

## Interpretation and limitations

- **Proximity is not proof of IS-mediated mobilization.** An IS-flanked arrangement is only a *candidate* module, not an experimentally demonstrated transposon.
- Candidate modules must be bounded by consecutive, nonoverlapping IS annotations, with fully enclosed target determinants. Other unannotated genes may occur between them.
- Circular correction applies to **distance**, not module calling; modules spanning the origin will be missed.
- The permutation test randomly relocates each determinant within its original plasmid, preserving gene length and IS positions. It does **not** control for gene density, shared ancestry, functional genomic compartments, or host/Inc composition. P-values are exploratory.
- Repeated or closely related plasmids may inflate recurrent patterns. Consider dereplication, stratification and sensitivity analyses.
- Report annotation database/software versions, detection cutoffs, source sampling, exclusion criteria, denominators, random seed and permutations.

# IS–determinant proximity pipeline

This is a revised, command-line version of `IS_proximity.py`, intended to support reproducible comparative analysis of AMR, virulence-factor (VF), and metal-resistance (METRES) annotations.

## Install

```bash
python -m pip install pandas numpy openpyxl
```

## Example: AMR

```bash
python IS_proximity_publication.py \
  --analysis AMR \
  --target amr_annotations.xlsx \
  --is-file is_annotations.csv \
  --replicon IncFIB \
  --output results/IncFIB_AMR
```

## Example: circular plasmids and positional permutation test

```bash
python IS_proximity_publication.py \
  --analysis VF \
  --target vf_annotations.xlsx \
  --is-file is_annotations.csv \
  --lengths plasmid_lengths.csv \
  --circular \
  --permutations 1000 \
  --seed 42 \
  --replicon IncFIB \
  --output results/IncFIB_VF
```

`plasmid_lengths.csv` must contain `Contig id,Length_bp` and a row for every contig in the analysis. Sequence identifiers must match exactly across files.

For metal resistance, use `--analysis METRES` and the original headerless BLAST-like annotation format.

## Output

The output directory contains ten files: summary, nearest IS, IS within window, IS pairs, exact recurrent distances, candidate modules, candidate module families, same-IS-label families, permutation null results, and run metadata JSON.

## Important interpretation

- Input coordinates are **1-based, inclusive**, and intervals are non-wrapping.
- Proximity percentages are reported with two denominators: all target determinants and determinants on IS-positive plasmids.
- Circular-origin correction is available for proximity **only**. Candidate IS–gene–IS modules crossing the origin are not identified.
- A module is called only when target genes lie fully between two consecutive non-overlapping annotated ISs. This is a *candidate structure*, not a validated transposon. Strand orientation is not considered.
- Exact recurrence can be inflated by clonally related plasmids. Dereplicate or stratify plasmids before interpreting recurrent structures.
- The positional permutation null relocates target intervals uniformly on the same plasmid, preserving their lengths and observed IS annotations. It does not control for local gene density, phylogeny, or the number of related plasmids.
- Inspect excluded rows, annotation naming, coordinates, and matching contig identifiers before drawing conclusions.
- For publication, report annotation software/database versions, thresholds, the number of permutations, seed, whether circular correction was used, the proximity denominator, and limitations of candidate module calling.

## Validation

The revised pipeline was smoke-tested using a small synthetic AMR dataset with an IS-free plasmid, an IS excluded by name, and a gene near a circular origin. This does **not** validate biological performance on your real datasets.

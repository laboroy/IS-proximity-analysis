# Data dictionary

## Source input

- `Contig id`: exact identifier of the plasmid sequence (must agree across files)
- `Gene symbol` / parsed `Gene`: annotation label of the determinant
- `Start`, `Stop`: 1-based inclusive genomic coordinates, regardless of strand
- `Length_bp`: positive integer sequence length
- IS annotations: original name before `|`, filtered by the configured regex

## Main outputs

- `Target_ID`, `IS_ID`: run-specific feature identifiers, not stable external accessions
- `Gene_start`, `Gene_end`, `IS_start`, `IS_end`: normalized ascending coordinates
- `Distance_bp`: minimum gap between intervals; zero if overlapping; adjacent nonoverlapping intervals have distance 1 bp
- `Within_<N>bp`: nearest annotated IS is within N bp
- `occurrences`: number of annotated determinant occurrences / structural modules
- `plasmids` / `n_plasmids`: distinct sequence IDs containing the association
- `Canonical_pattern`: lexicographically normalized label and target-gene order; **not** nucleotide sequence identity or strand-aware orientation
- `Same_IS_label`: same IS label at both boundaries; does not prove transposition
- `Observed_fraction_all`: proportion among **all** target occurrences, including IS-negative plasmids
- `Null_mean`, `Null_sd`: permutation-based positional-null summary
- `Empirical_p_enrichment`: one-sided `(1 + number of null fractions >= observed)/(permutations + 1)`; exploratory, unadjusted

## Missing data and filtering

Malformed coordinates, missing labels and exact duplicate coordinate/label combinations are excluded. The pipeline logs total excluded/duplicate records together; it does not currently export a row-level rejection audit. Verify data quality upstream. Non-IS names are excluded by the IS regex. The analysis counts determinant *occurrences*, not necessarily unique genes or independent plasmid lineages.

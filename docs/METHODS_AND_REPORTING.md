# Methods and reporting guidance

## Adaptable Methods text (verify all settings before using)

> Insertion sequence (IS) proximity to antimicrobial resistance, virulence, and metal-resistance determinants was assessed using a custom Python workflow (repository URL, release tag, and access date to be inserted). Annotation coordinates were interpreted as 1-based inclusive intervals. For each determinant, the minimum interval distance to an annotated IS on the same plasmid was calculated, with overlapping intervals assigned a distance of zero. Where complete plasmid lengths were available, distances were corrected for circular sequence boundaries. The proportions of determinants within [thresholds] of an IS were calculated both across all annotated determinants and among determinants on plasmids containing at least one annotated IS. Candidate IS-flanked structures were identified as determinants fully enclosed by two consecutive nonoverlapping IS annotations within [maximum span] bp. These structures were treated as putative genomic arrangements rather than experimentally confirmed transposable elements. [If used:] A within-plasmid positional permutation test ([number] permutations, seed [seed]) randomly reassigned determinant positions while preserving their lengths and observed IS positions to assess proximity enrichment under a simple positional null model.

## Before publication, document

1. Source assemblies, sequence completeness/circularity and sample selection criteria.
2. Plasmid typing approach and handling of multireplicon plasmids.
3. Annotation tool versions, databases and similarity thresholds for IS, ARG, VF and metal resistance.
4. Coordinate conversion and ID matching across files; excluded records.
5. Exact denominator for every reported percentage; numbers of plasmids and gene occurrences.
6. IS name filtering, proximity thresholds, maximum module span and circular-distance setting.
7. Permutation count, seed, null-model assumptions and whether multiple comparisons were adjusted.
8. Effects of plasmid length, IS density, host genus, Inc group and shared ancestry.
9. Repository release/commit, software environment, data accession numbers and any access restrictions.

## Scientific caution

IS proximity alone does not establish IS-mediated acquisition. IS-flanked regions require structural and mechanistic validation. Host-associated virulence cargo is consistent with coadaptation but does not demonstrate reciprocal coevolution. A positional null does not establish evolutionary causation.

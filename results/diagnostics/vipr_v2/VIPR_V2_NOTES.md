# VIPR copy-dependence v2 — scope notes

## Why the 94-nt vrRNA test is NOT run

The instruction permits testing whole 94-nt vrRNA units only if additional
homologous vrRNAs are confidently located. They are not.

A full-length scan of NC_028808.2 (both strands, all offsets) against the labelled
94-nt vrRNA gives exactly one perfect match (the label itself, 29135). The best
non-self match is 35/94 mismatches = 63% identity at 27855. Random DNA gives ~25%
identity, so 63% is elevated but far short of a confident homolog call, and that
position lies inside the off-label candidate NC_028808.2:[27682,27906) which
remains explicitly unconfirmed. Using it as a "homologous vrRNA" would promote an
unconfirmed candidate to ground truth.

Recorded as `missing_reference` for the 94-nt arm. The off-label candidate's status
is unchanged: sequence-supported repeat candidate, vrRNA identity unconfirmed.

## v1 defect this replaces

v1 used 5-nt units with a 10-nt flank, so the nearest perturbed unit sat 10 nt
(2 periods) from the target and the copy arm necessarily rewrote the target's
immediate context. Its three targets disagreed in sign (+0.189, -0.021, -0.413).

## v2 design

Primary region: the mapped 60-nt YNNGG segment, window w0003 [1487,1547)
(genomic NC_028808.2:[29169,29229), minus strand).

- flank 15 nt, so only units >= 3 periods from the target are perturbed
- two unit definitions: 12 single 5-nt units, and 6 paired 10-nt blocks
- targets at both ends and the middle of the region
- measured separation of the nearest perturbed unit: 15 nt (single), 20 nt (pair)
- mutation, budget-matched and distance-matched non-copy arm, 5 fixed seeds and
  target-only re-scoring are reused unchanged from copy_dependence.py

Known weakness: the middle `pair` target leaves only one distant block (10 nt) to
perturb, so its copy arm has the smallest budget in the set. Read that row with
that in mind.

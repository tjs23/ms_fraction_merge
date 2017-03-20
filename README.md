
** Merging MS Fraction Abundances **

This software merges fractionated mass specrometry abundance profiles from
replicate experiments where fraction data is not immediately comparable between
replicates and must first be aligned to combine equivalent fractions.

Merging is achieved by progressive, pairwise aggregation of abundance profiles,
using a scheme that maximises the sum of spectral count Pearson correlation
coefficients of components (typically proteins) shared between aligned
fractions. Fraction alignment involves an exhaustive search of relative offset
and linear width scaling.

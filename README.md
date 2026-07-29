# Hypercube probabilistic estimates

This repository records small, proved refinements of probabilistic estimates on
the Boolean hypercube.  It is separate from the
[A183068 supercongruence repository](https://github.com/rbajaj5/a183068-supercongruence):
the results here may help analyze spectra arising there, but they are not
supercongruences.

## Current result

| Result | Mathematical status | What improves | Novelty status |
|---|---|---|---|
| [Exact hashing of affine Fourier spectra](notes/AffineSpectrumHashing.md) | Complete elementary proof; exhaustive exact checks in the stated finite ranges | Replaces a generic `2 log₂(s)` collision budget by an exact formula and a `log₂(s)`-scale sufficient budget when the support is an affine subspace | The random-matrix rank formula is classical.  The contribution here is its explicit specialization to matching-phase/dyadic-defect spectra; no literature-priority claim is made. |

## Economist-style assessment

Scores are out of 10 and measure the result as it stands, not a hoped-for
future theorem.

| Criterion | Score | Reason |
|---|---:|---|
| Mathematical impact | 3 | A useful parameter sharpening, not a deep new theorem |
| Reusability | 8 | Applies to any affine Fourier support over `F₂` |
| Deployment value | 6 | Can reduce the number of hash buckets and measurements in structured sparse-spectrum tests |
| Proof confidence | 9 | Short proof plus exact exhaustive checks |
| Literature novelty | 2 | The rank formula is standard; the application is the new organizational point |

## Verification

Run:

```text
python verification/verify_affine_spectrum_hashing.py
```

The checker uses exact integer and rational arithmetic.  It has no third-party
dependencies.

## Source boundary

The motivating reading report,
*Testing of Hypercube Functions: Booleanity and Sparsity* (Jiang Yu, 2025),
summarizes earlier work and does not claim the refinement proved here.  The
generic hashing estimate is Proposition 3 of Gopalan--O'Donnell--Servedio--
Shpilka--Wimmer.  This repository does not redistribute the report.

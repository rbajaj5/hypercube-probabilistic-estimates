# Hypercube probabilistic estimates

This repository records small, proved refinements of probabilistic estimates.
It began with estimates on the Boolean hypercube and now also includes a
separate module on small deviations of sums of independent nonnegative random
variables. It is separate from the
[A183068 supercongruence repository](https://github.com/rbajaj5/a183068-supercongruence):
the hypercube results here may help analyze spectra arising there, but none of
the results in this repository are supercongruences. The hypercube and
small-deviation modules are mathematically independent.

## Current results

| Result | Mathematical status | What improves | Novelty status |
|---|---|---|---|
| [Exact hashing of affine Fourier spectra](notes/AffineSpectrumHashing.md) | Complete elementary proof; exhaustive exact checks in the stated finite ranges | Replaces a generic `2 log₂(s)` collision budget by an exact formula and a `log₂(s)`-scale sufficient budget when the support is an affine subspace | The random-matrix rank formula is classical.  The contribution here is its explicit specialization to matching-phase/dyadic-defect spectra; no literature-priority claim is made. |
| [Matroid law for exact linear hashing](notes/MatroidHashingLaw.md) | Complete proof from inclusion--exclusion and the Crapo--Rota Critical Theorem; exact finite checks | Gives the exact probability for every finite support through the characteristic polynomial of its difference matroid | Classical theorem in a new hashing/Fourier presentation; no priority claim for the probability law |
| [Exact projector matrix-hypercube law](notes/ProjectorMatrixHypercubeLaw.md) | Complete elementary enumeration for two real orthonormal projector bases | Replaces a variance-only matrix-Rademacher estimate by the complete angle-sensitive operator-norm distribution; the Zhang projector frame is the extremal mutually unbiased case | Elementary structured calculation; no standalone priority claim |
| [Sharp small-deviation bounds in dimensions at most four](notes/LowDimensionalSmallDeviations.md) | Complete corollary of Samuels' theorem and the finite-dimensional inequality in Paulin's proof; explicit extremizers | For \(2\leq n\leq4\) and \(0<\delta<1\), strictly improves the \(\delta(n/(n+\delta))^n\) branch of [arXiv:2607.23980](https://arxiv.org/abs/2607.23980) | The sharp low-dimensional theorem is classical; no priority claim |
| [Concurrent small-deviation work](CONCURRENT_WORK.md) | Formula, chronology, and dependency audit of three July 2026 manuscripts | Separates their common Vlassis--Thomas proof input from the independent Samuels--Paulin route for \(n\leq4\) | Chronology only; no priority inference |

## Assessment of the hypercube hashing module

Scores are out of 10 and measure the result as it stands, not a hoped-for
future theorem.

| Criterion | Score | Reason |
|---|---:|---|
| Mathematical impact | 4 | A useful exact synthesis and parameter sharpening, not a new foundational theorem |
| Reusability | 9 | The matroid law applies to every finite support over `F₂` |
| Deployment value | 6 | Gives exact measurement budgets for structured sparse spectra |
| Proof confidence | 9 | Short proof plus exact exhaustive checks |
| Literature novelty | 2 | The rank formula is standard; the application is the new organizational point |

## Verification

Run:

```text
python verification/verify_affine_spectrum_hashing.py
python verification/verify_projector_matrix_hypercube.py
python verification/verify_low_dimensional_bounds.py
```

All three checkers use exact integer and rational arithmetic and have no
third-party dependencies. The small-deviation checker audits the displayed
formulas, the two sharpness families, the concurrent-work comparisons, and
the reduction of the Samuels candidates on a large exact grid. It does not
replace the cited extremal theorems.

See [LITERATURE_PRIORITY.md](LITERATURE_PRIORITY.md) for the documented
priority boundaries and the hypercube module's supercongruence-defect
application.

## Hypercube source boundary

The motivating reading report,
*Testing of Hypercube Functions: Booleanity and Sparsity* (Jiang Yu, 2025),
summarizes earlier work and does not claim the refinement proved here.  The
generic hashing estimate is Proposition 3 of Gopalan--O'Donnell--Servedio--
Shpilka--Wimmer.  This repository does not redistribute the report.

## Small-deviation source boundary

The \(n\leq4\) theorem follows from classical results of Samuels and the
finite-dimensional inequality in Paulin's proof. Its proof is independent of
the new Vlassis--Thomas calibration theorem used by the three concurrent July
2026 manuscripts. See [CONCURRENT_WORK.md](CONCURRENT_WORK.md) for the exact
comparison and [LITERATURE_PRIORITY.md](LITERATURE_PRIORITY.md) for the
documented priority boundary.

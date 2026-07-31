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
| [Heavy-tail projector mixtures](notes/ProjectorHeavyTailMixtures.md) | Exact Pareto scale-mixture and quartic-shock calculations | Gives the exact angular tail coefficient, its phase transition at exponent \(2\), common-shock exponent quartering, and the independent-shock logarithmic correction | Elementary Pareto/Erlang consequences of the projector law; no general regular-variation or systemic-risk priority claim |
| [Knot-volume hypercube mixtures](notes/KnotVolumeHypercubeMixtures.md) | Reproducible numerical experiment plus exact finite pushforward law | Turns fractional hyperbolic volumes of 16 seeded census knots into angle-labeled projector cubes, a normalized 256-configuration law, and exact Pareto tail coefficients | Experimental application of the projector law; no knot-volume equidistribution or random-knot claim |
| [Collatz parity cubes](notes/CollatzParityHypercube.md) | Exact finite-horizon parity cube plus classical Kolmogorov and Berry--Esseen consequences; exact affine-residue transfer law | Adds simultaneous prefix control, an \(O(m^{-1/2})\) endpoint normal approximation, fixed-mode circle decay, and a complete exceptional-prime-\(3\) residue law | Classical probability applied to the Bernstein--Lagarias parity coding, plus an elementary finite-field transfer calculation; no Collatz-conjecture claim |
| [Hypercube walks and matrix mixing](notes/HypercubeWalkMatrixMixing.md) | Exact Walsh spectrum for the coordinate-refresh walk, observable-specific relaxation, and beta-overlap matrix estimates | Separates full-cube \(d\log d\) mixing from degree-two matrix-energy \(d\)-scale relaxation; quantifies almost-orthogonal random projector frames using the GOE angular law | Classical hypercube mixing plus an exact matrix-observable specialization; no GOE spectral or Collatz-orbit claim |
| [Hex winner noise and mixing](notes/HexWinnerNoiseMixing.md) | Exact topology, spectra through \(3\times3\), influence through \(4\times4\), and a seeded batched noise/pivotality simulation through \(31\times31\) | Separates the raw tie-free winner from the complement-odd crossing contrast and supplies finite-size winner-noise and influence diagnostics | Classical Hex topology and Boolean noise analysis plus a reproducible experiment; no efficient winning-strategy or scaling-limit claim |
| [Reflection-positive hypercube pictures](notes/ReflectionPositiveHypercubePictures.md) | Realizes mirror gluing as a positive autocorrelation Hankel kernel; proves complete monotonicity, log-convexity, moment reconstruction, and influence-only bounds | Turns a few refresh/noise lags into spectral consistency constraints and adds a certified \(4\times4\) Hex mixing interval plus an exact two-step anchor | Classical finite moment and Walsh identities inspired by Jaffe--Liu's picture-language program; no Ising-duality, scaling-limit, or priority claim |
| [Noisy reciprocity on a two-bit cube](notes/NoisyReciprocityHypercube.md) | Exact four-state Markov law for memory-one Prisoner's Dilemma strategies; stationary payoffs and the signed spectrum of noisy generous Tit-for-Tat | Identifies the negative retaliation mode, quantifies how generosity shortens error cascades, and recovers reflection positivity after two-step subsampling | Classical Markov/game calculation motivated by Veritasium and Wu--Axelrod; no universal optimal-strategy or evolutionary-scaling claim |
| [Novelty perturbations of an optimized hypercube](notes/NoveltyRefreshHypercube.md) | Exact nonuniform Bernoulli heat-bath law, centered subset spectrum, escape time, entropy, Poisson-binomial prevalence, and chi-squared mixing formula | Separates rare escape, stationary diversity, old-objective cost, and relaxation; novelty changes the first three while leaving the refresh spectrum fixed | Classical product-measure heat-bath calculation motivated by Tagore's *The Fugitive III*, poem 26; no claim that entropy models aesthetic value or that the formulas are new |
| [Latent environments on a hypercube](notes/LatentEnvironmentHypercube.md) | Exact two-environment Bernoulli mixture, count law, overdispersion, centered Fourier hierarchy, and posterior odds | Shows how a shared hidden environment creates order-\(d^2\) excess count variance and interactions invisible to coordinate means | Elementary finite-mixture calculation motivated by the supplied *Guns, Germs, and Steel* reading; no historical-causation or priority claim |
| [Fixed-drawing plane saturation as a conflict hypercube](notes/PlaneSaturationConflictHypercube.md) | Exact reduction to independent/maximal independent sets, Bernoulli and hard-core laws, Boolean indicators, and random-greedy bounds | Converts prescribed-edge crossings into a conflict graph and gives exact planarity/saturation probabilities plus a certified saturated sampler | Classical independent-set probability inspired by Clifton--Salia; strictly narrower than their unlabeled plane-saturation model and no improvement of their \(1/16\) theorem |
| [Hard-core heat-bath mixing on conflict hypercubes](notes/HardCoreConflictWalkMixing.md) | Exact reversible hard-core chain, occupancy drift, first-edge escape, adjacent-state contraction, and a finite-time total-variation bound | Adds a tractable sampler for Bernoulli edge sets conditioned on fixed-drawing planarity and quantifies low-activity mixing through crossing degree | Classical Glauber/path-coupling specialization; sufficient rather than optimal threshold, with no unlabeled plane-saturation or new general mixing claim |
| [Exact discrete loop equations for projector hypercubes](notes/DiscreteLoopEquationsProjectorCubes.md) | Exact Rademacher integration by parts plus full finite-rank resolvent difference; matrix, trace, rank-one cavity, and multi-resolvent cavity identities | Gives a zero-remainder loop diagnostic, retains higher discrete product terms, exposes perturbations beyond the second moment, and adds a finite-sample residual certificate | Elementary finite-cube specialization motivated by Bourgade--Huang's cumulant and switching calculus; not their closed hierarchy, and no local-law or random-matrix-universality claim |
| [Hopf-fiber projective kernel](notes/HopfFiberHypercubeKernel.md) | Embeds Boolean vertices in \(\mathbb{CP}^{d-1}\), derives the exact Fubini--Study/Hamming kernel, and transfers it to product noise, refresh walks, and Hex Walsh masses | Adds an exact projectively weighted Hex correlation computable from existing spectra, plus a finite Gaussian determinant for the cube Laplacian | Elementary consequences of standard Hopf-fiber Fourier and Fubini--Study identities; no endorsement of the source preprint's physical claims or priority claim |
| [KAN baselines on the Boolean cube](notes/KANHypercubeBaselines.md) | Exact collapse of every scalar one-layer KAN on binary inputs to Walsh degree at most one; optimal \(L^2\) error and noise/mixing remainder certificates | Replaces unidentifiable spline resolution by controlled Walsh interactions and quantifies the gain exactly for Hex through \(3\times3\) | Elementary finite-cube projection identities motivated by arXiv:2407.11075's evaluation principles; no KAN-performance or priority claim |
| [State coverage, supervision, and coordination on a hypercube](notes/StateCoverageCoordinationHypercube.md) | Sharp total-variation certificate, exact finite-buffer and sparse-signal laws, Walsh noise smoothing, oracle-loss guarantees, and a hidden-protocol lower ceiling | Separates distribution-shift gains, move-level oracle consistency, trace perfection, and protocol information | Classical finite-probability and Fourier identities motivated by OvercookedV2 and arXiv:2607.08984; no reinforcement-learning or priority claim |
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
python verification/verify_projector_heavy_tails.py
python verification/verify_knot_volume_hypercubes.py
python verification/verify_collatz_parity_estimates.py
python verification/verify_hypercube_walk_matrix_mixing.py
python verification/verify_hex_winner_noise.py
python verification/verify_hex_noise_simulation.py
python verification/verify_reflection_positive_hypercube.py
python verification/verify_noisy_reciprocity_hypercube.py
python verification/verify_novelty_refresh_hypercube.py
python verification/verify_latent_environment_hypercube.py
python verification/verify_plane_saturation_conflict_hypercube.py
python verification/verify_hard_core_conflict_walk.py
python verification/verify_discrete_loop_projector_cubes.py
python verification/verify_hopf_fiber_hypercube.py
python verification/verify_kan_hypercube_baselines.py
python verification/verify_state_coverage_coordination.py
python verification/verify_low_dimensional_bounds.py
```

All nineteen checkers have no third-party dependencies. The knot-volume checker
uses floating-point arithmetic to audit the fixed numerical dataset and exact
finite-law formulas; the other checkers use exact integer and rational
arithmetic except for the Hex simulation-data audit, which recomputes the
reported floating-point estimates from fixed integer counts. The optional
Hex data generator uses PyTorch and selects CUDA automatically when available.
The small-deviation checker audits the displayed formulas, the
two sharpness families, the concurrent-work comparisons, and the reduction
of the Samuels candidates on a large exact grid. It does not replace the
cited extremal theorems.

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

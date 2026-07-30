# Literature-priority audit

## Verdict

The exact injectivity probability for a random linear hash is not a new
standalone theorem.  In full generality it is the binary hashing formulation
of the Crapo--Rota Critical Theorem: the answer is an evaluation of the
characteristic polynomial of the matroid represented by the distinct nonzero
differences of the hashed set.

The affine-subspace product formula is the projective-geometry specialization
of that theorem, equivalently the standard full-column-rank probability for a
random binary matrix.

The repository-specific contribution is therefore limited but useful:

1. identify the difference matroid as the exact invariant for linear hashing;
2. translate the general Critical Theorem into Fourier-bucket language;
3. apply it to the affine Walsh supports supplied by the dyadic
   supercongruence-defect theorem; and
4. state both internal-spectrum and cross-spectrum measurement laws.

No literature-priority claim is made for the underlying probability formula.

## Source matrix

| Source | Relevant content | Relationship to this repository |
|---|---|---|
| Crapo--Rota (1970), *Combinatorial Geometries*, as summarized in [Alfarano--Byrne (2023)](https://arxiv.org/abs/2305.07567) | Origin of the Critical Theorem and its characteristic-polynomial count for representable matroids | Contains the general counting mechanism |
| [Handbook of the Tutte Polynomial, Theorem 16.12](https://www.routledge.com/Handbook-of-the-Tutte-Polynomial-and-Related-Topics/Ellis-Monaghan-Moffatt/p/book/9780367659682) | Counts ordered tuples of codewords whose union has prescribed support by a characteristic-polynomial evaluation | Gives the exact theorem used in the proof |
| [Alon--Dietzfelbinger--Miltersen--Petrank--Tardos (1997), *Linear Hashing*](https://www.brics.dk/RS/97/16/) | Studies random linear/affine transformations as universal hash families | Establishes that structured linear hashing is an old subject |
| [Gopalan--O'Donnell--Servedio--Shpilka--Wimmer (2011)](https://doi.org/10.1137/100785429) | Random coset hashing of Fourier characters; pair-collision and `2 log(s)` union bound | The immediate Fourier-testing baseline |
| [Haviv--Regev (2015)](https://arxiv.org/abs/1504.01649) | Improves Booleanity testing and learning bounds for Fourier-sparse Boolean functions | Prevents an overclaim that the present specialization improves the best general algorithms |
| [Chakraborty--Datta--Dutta--Ghosh--Sanyal (2024)](https://doi.org/10.4230/LIPIcs.MFCS.2024.40) | Extends Fourier-sparsity testing and coset hashing to other finite abelian groups | Shows the broader hashing program is active |
| [*Price of Parsimony* (NeurIPS 2025)](https://papers.nips.cc/paper_files/paper/2025/hash/f17376c941d5882050e2e366bb74dffa-Abstract-Conference.html) | Gives a tolerant Fourier-sparsity estimator with near-linear dependence on sparsity | Confirms that our result is a structured certificate, not a new general query-complexity record |
| [Ferrere et al. (2026)](https://arxiv.org/abs/2510.07088), *Fourier Analysis on the Boolean Hypercube via Hoeffding Functional Decomposition* | Builds a generalized Fourier expansion for arbitrary full-support input measures | A complementary nonuniform-measure direction; it does not supply or supersede the difference-matroid hashing law |

## Searches performed

The audit used combinations of:

- `"linear hashing" finite field`;
- `"random linear map" injective affine subspace`;
- `"Crapo Rota critical theorem" codeword support`;
- `"characteristic polynomial" matroid random linear map kernel`;
- `"Fourier sparsity" affine support coset hashing`;
- `"binary matroid" critical exponent`; and
- recent Fourier-sparsity testing literature through 2025.

The absence of a search hit is not evidence of priority.  The positive
Crapo--Rota identification is decisive: the exact general law is classical.

## Consequence for the supercongruence program

The audit produces a genuine but auxiliary improvement.  The dyadic
joint-spectrum theorem already proves that every scalar Walsh support is an
affine space with dimension determined by a rank $2s_\lambda$.  The
Critical-Theorem translation turns that rank data into:

- exact success probabilities for random linear measurements of each scalar
  spectrum;
- a simultaneous internal-separation bound over many scalar tests; and
- an exact characteristic-polynomial law when all cross-spectrum collisions
  must also be excluded.

This improves experimental design and finite certification.  It does not
increase a $p$-adic exponent or prove an additional supercongruence.

## Small-deviation module

### Verdict

The sharp bound recorded in
[notes/LowDimensionalSmallDeviations.md](notes/LowDimensionalSmallDeviations.md)
for \(n\leq4\) is not a new standalone theorem. It follows from:

1. Samuels' extremal theorem for sums of at most four independent
   nonnegative random variables;
2. the stronger finite-dimensional inequality inside Paulin's proof; and
3. the elementary mean-completion transformation
   \(Y_i=X_i+1-\mathbb EX_i\).

The useful contribution of this module is organizational: it states the exact
low-dimensional consequence in the notation of the July 2026 work, proves
the strict comparison with their shared \(0<\delta<1\) bound, records both
sharpness families, and separates their common new-preprint dependency from
the independent classical proof path.

The exact \(n=2\) admissible-merger proof is explicitly recorded by
Nie--Wei. No priority claim is made for either that proof route or the
low-dimensional constant.

### Source matrix

| Source | Relevant content | Relationship to this module |
|---|---|---|
| [Samuels (1966)](https://doi.org/10.1214/aoms/1177699614) and Samuels (1968) | Extremal lower-tail theorem proved through four summands | Supplies the classical low-dimensional result |
| [Feige (2004/2005 manuscript)](https://www.wisdom.weizmann.ac.il/~feige/Others/newmarkov.pdf) | States the endpoint extremizers and reports Samuels' theorem for \(n\leq4\) | Supplies the centered-small-deviation formulation and historical boundary |
| [Paulin (2017)](https://arxiv.org/abs/1703.05152) | Gives a stronger finite-dimensional two-endpoint inequality inside the proof that Samuels' conjecture implies Feige's conjecture | Reduces the finite Samuels family to the two endpoint candidates |
| [Vlassis--Thomas (2026)](https://arxiv.org/abs/2607.08415) | Proves validity of Gaffke's Dirichlet merger | Shared input to the concurrent geometric proofs, but not to the \(n\leq4\) Samuels--Paulin proof |
| [Ming--Ramdas--Shen--Wang--Waudby-Smith (2026)](https://arxiv.org/abs/2607.18661) | Constructs the exact admissible bivariate dominator \(K_2^{\mathrm{ad}}\) | Supplies the merger used in the short \(n=2\) proof |
| [Fu--Han--Wang--Yan--Zhang--Zhou (2026)](https://arxiv.org/abs/2607.23980) | Gives the dimension-dependent bound for every \(n,\delta\), sharp for \(\delta\geq1\) | Immediate baseline sharpened here in dimensions \(2,3,4\) |
| [Stander (2026)](https://doi.org/10.5281/zenodo.21622951) | Conditionally derives the shared \(0<\delta\leq1\) formula from Vlassis--Thomas | Makes the shared recent-preprint dependency explicit |
| [Nie--Wei (2026)](https://arxiv.org/abs/2607.24528) | Gives the shared \(0<\delta\leq1\) bound and the exact all-slack result for \(n=2\) | Supplies prior art for the modern bivariate proof and confirms the general scope boundary |

The three July manuscripts and their exact overlap are compared in
[CONCURRENT_WORK.md](CONCURRENT_WORK.md). Their public versions appeared
within hours of one another, so this repository records chronology without
drawing a priority conclusion.

A separate [2025 preprint](https://arxiv.org/abs/2508.07316) claims the full
arbitrary-\(\delta\) conjecture, while the later July 2026 papers describe the
general small-slack regime as open. This repository does not adjudicate or
rely on that claim.

## Projector matrix-hypercube module

### Verdict

The
[exact two-basis operator-norm law](notes/ProjectorMatrixHypercubeLaw.md)
is an elementary sixteen-vertex calculation, not a claimed new general
matrix-concentration theorem. Its value is that it extracts a valid
probabilistic consequence from the rank-one projector geometry behind the
counterexample to Zhang's Conjecture 4.1:

1. it gives the complete angle-sensitive norm distribution and every
   positive norm moment;
2. it identifies the mutually unbiased Zhang frame as the maximizer of the
   expected norm within the two-basis family; and
3. it shows exactly what information is lost by the common variance proxy
   \(\|\sum_iP_i^2\|=2\), which is angle-independent.

The module neither uses the false conjectured inequality nor changes the
scalar small-deviation and linear-hashing results.

### Source matrix

| Source | Relevant content | Relationship to this module |
|---|---|---|
| [Zhang (2014/2018)](https://arxiv.org/abs/1411.5058) | States the four-matrix auxiliary conjecture in the noncommutative AM--GM program | Supplies the conjectural setting; the exact projector counterexample supplies this module's distinguished frame |
| [Tropp (2010/2012)](https://arxiv.org/abs/1004.4389) | Gives dimension-dependent subgaussian bounds for matrix Gaussian and Rademacher series in terms of a matrix variance proxy | Supplies the general concentration baseline that the exact structured law calibrates |

### Heavy-tail companion boundary

The [Pareto-mixture companion](notes/ProjectorHeavyTailMixtures.md) contains
exact consequences of the finite projector law:

- conditioning a bounded angular variable against a standard Pareto radial
  scale;
- the elementary concavity/convexity transition of its \(\alpha\)-moment at
  \(\alpha=2\);
- degree-four homogeneity under a shared scale; and
- the four-stage Erlang survival formula for a product of independent Pareto
  variables.

These are self-contained calculations. The note does not claim a new general
result in regular variation, multivariate extremes, or systemic-risk theory.

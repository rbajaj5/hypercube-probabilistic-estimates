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

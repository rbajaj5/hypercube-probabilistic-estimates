# Exact hashing of affine Fourier spectra

## 1. Motivation

Random coset hashing is a standard way to separate a small set of Fourier
frequencies.  Let

```math
L:\mathbb F_2^n\longrightarrow\mathbb F_2^t
```

be a uniformly random linear map.  Frequencies in the same fiber of $L$ land
in the same bucket.  For two distinct frequencies $a,b$,

```math
\Pr[L(a)=L(b)]=2^{-t}.
```

For an arbitrary set of at most $s+1$ frequencies, a union bound over pairs
therefore gives the familiar sufficient scale

```math
t\ \geq\ 2\log_2 s+\log_2(1/\delta)
```

for collision probability at most $\delta$.  This is the estimate used in
Proposition 3 of Gopalan--O'Donnell--Servedio--Shpilka--Wimmer.

The pair count is unnecessarily large when the support has additive
structure.  This note gives the difference-set refinement and the exact answer
for an affine subspace.

## 2. A difference-set bound

For a finite nonempty set $A\subseteq\mathbb F_2^n$, define

```math
D(A)=A+A=\{a+b:a,b\in A\}.
```

### Theorem 1

For a uniformly random linear map
$L:\mathbb F_2^n\to\mathbb F_2^t$,

```math
\Pr[L|_A\text{ is not injective}]
\ \leq\
\frac{|D(A)|-1}{2^t}.
```

### Proof

The restriction $L|_A$ is not injective exactly when there are distinct
$a,b\in A$ with $L(a+b)=0$.  Thus some nonzero member of $D(A)$ lies in
$\ker L$.  For each fixed nonzero $d$, the vector $L(d)$ is uniform in
$\mathbb F_2^t$, so

```math
\Pr[L(d)=0]=2^{-t}.
```

A union bound over the distinct nonzero elements of $D(A)$ proves the
claim. $\square$

This can be much smaller than the pair bound: it counts distinct differences,
not pairs that produce them.

## 3. Exact affine-subspace formula

### Theorem 2

Let $A=a+V\subseteq\mathbb F_2^n$, where $V$ is a linear subspace of
dimension $r$.  Then

```math
\Pr[L|_A\text{ is injective}]
=
\prod_{j=0}^{r-1}\left(1-2^{j-t}\right)
\quad\text{if }t\geq r,
```

whereas the probability is zero if $t<r$. Consequently, for $t\geq r$,

```math
\Pr[L|_A\text{ is not injective}]
\leq
\sum_{j=0}^{r-1}2^{j-t}
=
\frac{2^r-1}{2^t}
<
2^{r-t}.
```

In particular, for $0<\delta<1$, it is sufficient to take

```math
t\geq r+\left\lceil\log_2(1/\delta)\right\rceil.
```

### Proof

Translation by $a$ does not affect collisions.  Hence $L|_A$ is injective
if and only if $L|_V$ has rank $r$.  After choosing a basis of $V$, this
restriction is a uniformly random $t\times r$ binary matrix.

If $t<r$, full column rank is impossible.  If $t\geq r$, expose the
columns in order.  The first column has $2^t-1$ choices outside the zero
span, the second has $2^t-2$ choices outside the first column's span, and
the $j$-th zero-indexed column has $2^t-2^j$ choices outside the span of
the preceding columns.  Dividing by the $2^{tr}$ possible matrices gives

```math
\prod_{j=0}^{r-1}\frac{2^t-2^j}{2^t}
=
\prod_{j=0}^{r-1}(1-2^{j-t}).
```

Finally,

```math
1-\prod_j(1-x_j)\leq\sum_jx_j
\qquad(0\leq x_j\leq1)
```

gives the displayed failure bound. $\square$

If $s=|A|=2^r$, the sufficient number of hash bits is therefore

```math
t\geq \log_2 s+\left\lceil\log_2(1/\delta)\right\rceil,
```

rather than the generic pair-union scale
$2\log_2s+\log_2(1/\delta)$.

## 4. Matching phases and dyadic defect spectra

Consider the Boolean phase

```math
g_{\delta,k}(x_1,y_1,\ldots,x_k,y_k,z)
=
(-1)^{\delta z+\sum_{j=1}^k x_jy_j},
\qquad \delta\in\{0,1\}.
```

For one pair,

```math
(-1)^{xy}
```

has four nonzero Walsh coefficients, with values
$\tfrac12,\tfrac12,\tfrac12,-\tfrac12$.  Products on disjoint pairs tensor
their Fourier supports, while multiplication by the character
$(-1)^{\delta z}$ translates the support.  Therefore

```math
\mathrm{supp}\,\widehat g_{\delta,k}
=
\delta e_z+
\mathrm{span}\{e_{x_1},e_{y_1},\ldots,e_{x_k},e_{y_k}\}.
```

This is an affine subspace of dimension $2k$ and cardinality

```math
s=2^{2k}=4^k.
```

Theorem 2 now gives the exact separation probability

```math
\Pr[L\text{ separates }\mathrm{supp}\,\widehat g_{\delta,k}]
=
\prod_{j=0}^{2k-1}(1-2^{j-t})
\quad\text{if }t\geq2k,
```

and the probability is zero if $t<2k$. Thus the hashing budget for this
spectrum is

```math
t\geq 2k+\left\lceil\log_2(1/\delta)\right\rceil
=
\log_2(4^k)+\left\lceil\log_2(1/\delta)\right\rceil.
```

The same matching phases occur in the exact Walsh spectra of dyadic
supercongruence defect coordinates.  This observation converts a structural
description of those spectra into an exact measurement/hash-success law.

## 5. What is and is not claimed

The full-rank probability for a random binary matrix is classical.  The
difference-set union bound is elementary.  The useful point here is their
explicit application to affine Walsh supports, including the matching-phase
spectra above.

No claim is made that the general sparse-Fourier testing algorithm can simply
replace its $O(s^2)$ buckets by $O(s)$ buckets for arbitrary inputs.  The
improvement applies when one already knows, or separately verifies, that the
relevant support is an affine subspace of dimension $\log_2s$.

## References

1. Parikshit Gopalan, Ryan O'Donnell, Rocco A. Servedio, Amir Shpilka, and
   Karl Wimmer, “Testing Fourier Dimensionality and Sparsity,” *SIAM Journal
   on Computing* 40(4), 1075–1100 (2011),
   [doi:10.1137/100785429](https://doi.org/10.1137/100785429).
2. Tom Gur and Omer Tamuz, “Testing Booleanity and the Uncertainty
   Principle,” [arXiv:1204.0944](https://arxiv.org/abs/1204.0944).
3. Jiang Yu, “Testing of Hypercube Functions: Booleanity and Sparsity,”
   harmonic-analysis reading report (2025),
   [source PDF](https://acm.sjtu.edu.cn/~jy_15924374500/HA_final_Project.pdf).
4. “Joint defect spectrum: exact Walsh support and effective dimension,”
   [related dyadic-supercongruence note](https://github.com/rbajaj5/a183068-supercongruence/blob/main/related-results/DyadicHypercubeJointSpectrum.md).

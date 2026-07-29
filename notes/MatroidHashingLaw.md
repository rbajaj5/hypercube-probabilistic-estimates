# The matroid law for exact linear hashing

## 1. Statement

Let $A$ be a finite nonempty subset of $\mathbb F_2^n$, and let

```math
E_A=(A+A)\setminus\{0\}.
```

Regard the distinct vectors in $E_A$ as the columns of a binary matrix.
They represent a simple binary matroid $M_A$. Write

```math
\rho=\mathrm{rank}(M_A)
=
\dim\mathrm{span}(E_A),
```

and let its characteristic polynomial be

```math
\chi_{M_A}(q)
=
\sum_{X\subseteq E_A}
(-1)^{|X|}q^{\rho-\mathrm{rank}(X)}.
```

Choose a linear map

```math
L:\mathbb F_2^n\longrightarrow\mathbb F_2^t
```

uniformly at random.

### Theorem 1 (exact hashing law)

The probability that $L$ is injective on $A$ is

```math
\Pr[L|_A\text{ is injective}]
=
\frac{\chi_{M_A}(2^t)}{2^{t\rho}}.
\qquad\text{(1)}
```

Thus the injectivity probability profile, as $t$ varies, is determined by the
binary matroid of the distinct nonzero differences of $A$, not merely by
$|A|$.

## 2. Direct proof

For $d\in E_A$, let $B_d$ be the bad event $L(d)=0$. The restriction
$L|_A$ is injective exactly when none of the events $B_d$ occurs.

For $X\subseteq E_A$, all events $B_d$, $d\in X$, occur exactly when
every row of $L$ annihilates $\mathrm{span}(X)$. A single uniformly random
row does this with probability $2^{-\mathrm{rank}(X)}$. Since the $t$
rows are independent,

```math
\Pr\left[\bigcap_{d\in X}B_d\right]
=
2^{-t\,\mathrm{rank}(X)}.
```

Inclusion--exclusion gives

```math
\begin{aligned}
\Pr[L|_A\text{ is injective}]
&=
\sum_{X\subseteq E_A}
(-1)^{|X|}2^{-t\,\mathrm{rank}(X)}\\
&=
2^{-t\rho}
\sum_{X\subseteq E_A}
(-1)^{|X|}(2^t)^{\rho-\mathrm{rank}(X)}\\
&=
\frac{\chi_{M_A}(2^t)}{2^{t\rho}}.
\end{aligned}
```

This proves (1). $\square$

## 3. Critical-Theorem interpretation

The same calculation is the binary case of the Crapo--Rota Critical
Theorem.

Let $D$ be the matrix whose columns are the elements of $E_A$, and let
$C$ be its row-space code. Each row $\ell$ of $L$ gives a codeword

```math
(\ell(d))_{d\in E_A}\in C.
```

The restriction $L|_A$ is injective precisely when the union of the supports
of the resulting $t$ codewords is all of $E_A$. The Critical Theorem says
that the number of ordered $t$-tuples of codewords with full union support
is $\chi_{M_A}(2^t)$. The ambient rows induce independent uniform codewords
in $C$, so division by $|C|^t=2^{t\rho}$ recovers (1).

This establishes the correct priority boundary: equation (1) is a
specialization of a classical theorem, while its use as an exact certificate
for the structured Fourier spectra below is an application.

## 4. Recovery of the affine formula

Suppose $A=a+V$, with $\dim V=r$. Then

```math
E_A=V\setminus\{0\},
```

whose binary matroid is the projective geometry
$\mathrm{PG}(r-1,2)$. Its characteristic polynomial is

```math
\chi_{M_A}(q)
=
\prod_{j=0}^{r-1}(q-2^j).
```

Substituting $q=2^t$ in (1) gives

```math
\Pr[L|_A\text{ is injective}]
=
\prod_{j=0}^{r-1}(1-2^{j-t})
```

when $t\geq r$, and zero otherwise. This recovers the formula in
[AffineSpectrumHashing.md](AffineSpectrumHashing.md).

## 5. Dyadic defect spectra

For the finite dyadic defect map in the
[joint-spectrum theorem](https://github.com/rbajaj5/a183068-supercongruence/blob/main/related-results/DyadicHypercubeJointSpectrum.md),
every scalar output test $\lambda$ has Walsh support

```math
\mathcal A_\lambda
=
\{\xi:\xi|_{R_\lambda}=r_\lambda\}.
```

It is an affine space of dimension $2s_\lambda$. Therefore a random
$t$-bit linear hash separates every frequency inside this scalar spectrum
with exact probability

```math
p_\lambda(t)
=
\prod_{j=0}^{2s_\lambda-1}(1-2^{j-t}),
\qquad t\geq2s_\lambda.
\qquad\text{(2)}
```

For a finite collection $\Lambda$ of scalar tests, a single random map
separates each spectrum internally with probability at least

```math
1-\sum_{\lambda\in\Lambda}\left(1-p_\lambda(t)\right).
\qquad\text{(3)}
```

This union bound permits collisions between different scalar spectra. If
cross-spectrum separation is also required, set

```math
\mathcal A=\bigcup_{\lambda\in\Lambda}\mathcal A_\lambda
```

and apply Theorem 1 to $\mathcal A$. The resulting exact probability is

```math
\frac{\chi_{M_{\mathcal A}}(2^t)}
{2^{t\,\mathrm{rank}(M_{\mathcal A})}}.
\qquad\text{(4)}
```

Equations (2)--(4) turn the dyadic rank/radical data into an exact
measurement-design certificate. They do not strengthen the underlying
supercongruence modulus and do not assert pseudorandomness of arithmetic
coefficient vectors.

## 6. References

1. Henry H. Crapo and Gian-Carlo Rota, *On the Foundations of
   Combinatorial Theory: Combinatorial Geometries*, MIT Press, 1970.
2. Joanna A. Ellis-Monaghan and Iain Moffatt (eds.), *Handbook of the Tutte
   Polynomial and Related Topics*, Section 16.2.3, CRC Press, 2022.
3. Gianira N. Alfarano and Eimear Byrne, “The Critical Theorem for
   q-Polymatroids,”
   [arXiv:2305.07567](https://arxiv.org/abs/2305.07567).
4. Noga Alon, Martin Dietzfelbinger, Peter Bro Miltersen, Erez Petrank, and
   Gábor Tardos,
   [*Linear Hashing*](https://www.brics.dk/RS/97/16/), BRICS Report
   RS-97-16, 1997.
5. Parikshit Gopalan, Ryan O'Donnell, Rocco A. Servedio, Amir Shpilka, and
   Karl Wimmer, “Testing Fourier Dimensionality and Sparsity,”
   [doi:10.1137/100785429](https://doi.org/10.1137/100785429).

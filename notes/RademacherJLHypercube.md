# Exact Rademacher Johnson--Lindenstrauss laws for cube pairs

## Scope and source boundary

[Li's unified Johnson--Lindenstrauss analysis](https://arxiv.org/abs/2402.10232)
covers spherical, binary-coin, sparse, Gaussian, and sub-Gaussian projection
models.  This note specializes the binary-coin construction to differences of
Boolean-cube vertices.  Because such a difference has exactly \(h\) nonzero
unit coordinates, its finite-dimensional distortion law is elementary and
exact.

The formulas below are classical Rademacher/binomial calculations.  They do
not improve the general Johnson--Lindenstrauss dimension bound or the paper's
Hanson--Wright analysis.  They provide distance-sensitive finite calibration
that a distribution-free tail estimate does not display.

## Projection of a Hamming-\(h\) pair

Let \(x,y\in\{0,1\}^d\) have Hamming distance \(h\geq1\).  Let
\(A\in\mathbb R^{m\times d}\) have independent entries

\[
A_{ri}=\frac{\varepsilon_{ri}}{\sqrt m},
\qquad
\Pr(\varepsilon_{ri}=1)=\Pr(\varepsilon_{ri}=-1)=\frac12.
\]

Signs in \(x-y\) can be absorbed into the Rademacher entries.  Therefore

\[
\frac{\|A(x-y)\|_2^2}{\|x-y\|_2^2}
=R_{m,h}
:=\frac1{mh}\sum_{r=1}^m S_{h,r}^2,
\tag{1}
\]

where the \(S_{h,r}\) are independent copies of

\[
S_h=\varepsilon_1+\cdots+\varepsilon_h.
\]

In particular, the law depends on the pair only through \(h\), not through
the ambient dimension or the locations of the differing coordinates.

## Complete finite distribution

Define the probability generating polynomial

\[
P_h(z)
=2^{-h}\sum_{k=0}^h\binom hk z^{(h-2k)^2}.
\tag{2}
\]

Then

\[
\Pr\!\left(R_{m,h}=\frac q{mh}\right)
=[z^q]P_h(z)^m.
\tag{3}
\]

Consequently, for every \(\delta>0\), the exact distortion probability is

\[
p_{m,h}(\delta)
=
\sum_{\substack{q\geq0\\|q/(mh)-1|\geq\delta}}
[z^q]P_h(z)^m.
\tag{4}
\]

This finite convolution is inexpensive when \(m\) and \(h\) are moderate and
has no asymptotic or unspecified constant.

Two low-distance cases reduce to one binomial variable:

\[
R_{m,1}=1
\quad\text{almost surely},
\tag{5}
\]

\[
R_{m,2}=\frac{2B}{m},
\qquad B\sim\operatorname{Bin}(m,\tfrac12),
\tag{6}
\]

and

\[
R_{m,3}=\frac{m+8B}{3m},
\qquad B\sim\operatorname{Bin}(m,\tfrac14).
\tag{7}
\]

Thus Rademacher projection preserves every cube edge exactly, a property
hidden by a generic JL bound.

## Exact first two moments

The Rademacher identities

\[
\mathbb E S_h^2=h,\qquad
\mathbb E S_h^4=3h^2-2h
\]

give

\[
\mathbb E R_{m,h}=1,
\qquad
\operatorname{Var}(R_{m,h})
=\frac{2(h-1)}{mh}
=\frac2m\left(1-\frac1h\right).
\tag{8}
\]

The usual \(2/m\) variance scale is approached as \(h\to\infty\), but it is
strictly loose at every finite Hamming distance and vanishes at \(h=1\).
Chebyshev immediately yields the distance-sensitive certificate

\[
\Pr(|R_{m,h}-1|\geq\delta)
\leq
\frac{2(h-1)}{mh\delta^2}.
\tag{9}
\]

Equation (4) is normally sharper; (9) is useful as a closed-form check.

## Exact collapse probability

A projected pair collides precisely when every row sum is zero.  Hence

\[
\Pr(Ax=Ay)
=
\begin{cases}
0,&h\text{ odd},\\[3pt]
\left(2^{-h}\binom h{h/2}\right)^m,&h\text{ even}.
\end{cases}
\tag{10}
\]

Parity therefore rules out collisions at all odd distances.  This statement
is special to the binary-coin matrix; it is not true for an arbitrary JL
distribution.

## The complete cube

The \(d\)-cube has

\[
N_{d,h}=2^{d-1}\binom dh
\tag{11}
\]

unordered vertex pairs at distance \(h\).  If \(D_\delta\) denotes the number
of pairs whose distortion is at least \(\delta\), linearity of expectation
and (4) give

\[
\mathbb E D_\delta
=2^{d-1}\sum_{h=1}^d\binom dh\,p_{m,h}(\delta).
\tag{12}
\]

For a uniformly selected distinct unordered pair,

\[
\Pr(\text{distorted})
=\frac1{2^d-1}\sum_{h=1}^d\binom dh\,p_{m,h}(\delta).
\tag{13}
\]

Likewise, the expected number of collapsed pairs is

\[
\mathbb E C
=2^{d-1}
\sum_{\substack{2\leq h\leq d\\h\ {\rm even}}}
\binom dh
\left(2^{-h}\binom h{h/2}\right)^m.
\tag{14}
\]

Markov's inequality gives
\(\Pr(A\text{ is noninjective on the cube})\leq\mathbb EC\).
This union bound can exceed one and is not asserted to be the exact
noninjectivity probability, because different pair-collision events are
dependent.

## What is verified

`verification/verify_rademacher_jl_hypercube.py` uses exact rational
arithmetic to verify:

1. the complete convolution law (2)--(4);
2. the exact mean and distance-sensitive variance (8);
3. the \(h=1,2,3\) special laws;
4. the collapse probability and Chebyshev certificate;
5. the Hamming-shell pair count (11); and
6. the expected full-cube collision count, independently audited by
   exhaustive enumeration of small Rademacher matrices.

## Relation to the other modules

- [Exact linear hashing](AffineSpectrumHashing.md) concerns maps over
  \(\mathbb F_2\).  The real-valued Rademacher projection here has a different
  collision law.
- [Hypercube walks and matrix mixing](HypercubeWalkMatrixMixing.md) also
  separates behavior by Walsh/Hamming scale.
- [Rare-event mass separation](RareEventMassSeparationHypercube.md) warns
  about exponentially rare sets.  Equations (12)--(14) aggregate similarly
  rare failures over exponentially many cube pairs.

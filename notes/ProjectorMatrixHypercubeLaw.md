# An exact matrix-hypercube law from two projector bases

## 1. Motivation

Teng Zhang's Conjecture 4.1 proposed a quartic operator-norm inequality for
four positive semidefinite matrices as an auxiliary route toward a
noncommutative arithmetic--geometric mean problem
([arXiv:1411.5058](https://arxiv.org/abs/1411.5058)). An exact counterexample
uses the four rank-one projectors from two mutually unbiased bases of
\(\mathbb R^2\).

The false inequality cannot be used as an estimate. Its projector geometry
does, however, give a useful matrix-valued Boolean-hypercube model. This note
computes the complete operator-norm law for a Rademacher sum of any two real
orthonormal projector bases in dimension two. The Zhang configuration is the
mutually unbiased member of the family and maximizes the expected norm.

This is a structured exact calculation, not a claimed improvement to general
matrix concentration or to the scalar small-deviation theorem elsewhere in
this repository.

## 2. The two-basis family

Let \(A,B,C,D\) be real rank-one orthogonal projectors on \(\mathbb R^2\)
such that

```math
A+B=C+D=I_2.
```

Set

```math
Z=A-B,\qquad Q=D-C.
```

Then \(Z\) and \(Q\) are real symmetric traceless involutions:

```math
Z^2=Q^2=I_2.
```

There is a parameter \(q\in[-1,1]\) such that

```math
ZQ+QZ=2qI_2,
\qquad
q=\frac12\operatorname{tr}(ZQ).
\tag{1}
```

After choosing the basis of \(A,B\), one may write

```math
Z=
\begin{pmatrix}1&0\\0&-1\end{pmatrix},
\qquad
Q=
\begin{pmatrix}q&r\\r&-q\end{pmatrix},
\qquad
q^2+r^2=1.
```

If the second basis is rotated by an angle \(\theta\), then
\(q=\cos(2\theta)\). The bases are mutually unbiased exactly when \(q=0\).

Let \(\varepsilon_A,\varepsilon_B,\varepsilon_C,\varepsilon_D\) be
independent uniform signs and define the matrix-valued hypercube observable

```math
S
=
\varepsilon_AA+\varepsilon_BB+\varepsilon_CC+\varepsilon_DD.
\tag{2}
```

## 3. Exact operator-norm distribution

### Theorem

For the uniform measure on \(\{-1,1\}^4\), the operator norm of \(S\) has
the atomic law

```math
\boxed{
\mathcal L(\|S\|)
=
\frac18\delta_0
+\frac58\delta_2
+\frac18\delta_{\sqrt{2+2q}}
+\frac18\delta_{\sqrt{2-2q}}.}
\tag{3}
```

Atoms at the same point are combined. Consequently, for every \(p>0\),

```math
\mathbb E\|S\|^p
=
\frac58\,2^p
+\frac18(2+2q)^{p/2}
+\frac18(2-2q)^{p/2}.
\tag{4}
```

In particular,

```math
\mathbb E\|S\|^2=3
\tag{5}
```

for every relative angle, while

```math
\frac32
\leq
\mathbb E\|S\|
=
\frac54+
\frac{\sqrt{2+2q}+\sqrt{2-2q}}8
\leq
\frac{5+\sqrt2}{4}.
\tag{6}
```

The upper equality in (6) holds exactly for mutually unbiased bases
\((q=0)\); the lower equality holds when the two bases coincide, up to
permutation \((|q|=1)\).

### Proof

Put

```math
m=
\frac{\varepsilon_A+\varepsilon_B
      +\varepsilon_C+\varepsilon_D}{2},
\qquad
x=\frac{\varepsilon_A-\varepsilon_B}{2},
\qquad
y=\frac{\varepsilon_D-\varepsilon_C}{2}.
```

Using \(A=(I_2+Z)/2\), \(B=(I_2-Z)/2\),
\(C=(I_2-Q)/2\), and \(D=(I_2+Q)/2\), equation (2) becomes

```math
S=mI_2+xZ+yQ.
\tag{7}
```

Classify the sixteen sign vectors by whether the \(A,B\) pair and the
\(C,D\) pair are aligned or opposite.

1. If both pairs are aligned, then \(x=y=0\). Two sign vectors give
   \(S=0\), and two give \(\|S\|=2\).
2. If exactly one pair is aligned, then \(S\) is a signed identity plus a
   signed involution. All eight sign vectors give eigenvalues \(0,\pm2\),
   so \(\|S\|=2\).
3. If both pairs are opposite, then \(m=0\) and \(x,y\in\{-1,1\}\).
   Equations (1) and (7) give

   ```math
   S^2
   =
   (xZ+yQ)^2
   =
   (2+2xyq)I_2.
   ```

   Two sign vectors have \(xy=1\) and two have \(xy=-1\), producing the
   last two atoms in (3).

This proves (3), and (4)--(5) follow directly. For (6), observe that

```math
\left(\sqrt{2+2q}+\sqrt{2-2q}\right)^2
=
4+4\sqrt{1-q^2}
\leq8.
```

The maximum occurs exactly at \(q=0\), and the minimum at \(|q|=1\).
\(\square\)

## 4. The Zhang projector hypercube

The counterexample projectors are

```math
A=
\begin{pmatrix}1&0\\0&0\end{pmatrix},
\quad
B=
\begin{pmatrix}0&0\\0&1\end{pmatrix},
\quad
C=\frac12
\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
\quad
D=\frac12
\begin{pmatrix}1&1\\1&1\end{pmatrix}.
\tag{8}
```

Here \(q=0\), so (3) specializes to

```math
\boxed{
\mathcal L(\|S\|)
=
\frac18\delta_0
+\frac14\delta_{\sqrt2}
+\frac58\delta_2.}
\tag{9}
```

Thus the exact upper-tail function is

```math
\mathbb P(\|S\|\geq t)
=
\begin{cases}
1,&t\leq0,\\
\frac78,&0<t\leq\sqrt2,\\
\frac58,&\sqrt2<t\leq2,\\
0,&t>2.
\end{cases}
\tag{10}
```

The moments and norm Laplace transform are

```math
\mathbb E\|S\|^p
=
\frac58\,2^p+\frac14\,2^{p/2},
\qquad p>0,
\tag{11}
```

and

```math
\mathbb E e^{\lambda\|S\|}
=
\frac18+\frac14e^{\sqrt2\lambda}+\frac58e^{2\lambda}.
\tag{12}
```

This is a strict calibration improvement over applying a generic
matrix-Rademacher tail inequality to this particular four-point frame.
Indeed,

```math
\sigma^2
:=
\left\|\sum_{P\in\{A,B,C,D\}}P^2\right\|
=
\|2I_2\|
=2
\tag{13}
```

for every \(q\). The variance proxy therefore cannot distinguish the
relative angle. The standard dimension-dependent subgaussian estimate sees
only \(\sigma^2=2\), whereas (3) gives the exact angle-sensitive law and
vanishes identically beyond the true endpoint \(2\). The general
matrix-Rademacher inequality used for this comparison is Theorem 4.1 of
[Tropp](https://arxiv.org/abs/1004.4389).

## 5. What this does and does not improve

- It adds an exact matrix-valued estimate on the Boolean hypercube, including
  all norm moments and tails, for the projector frame extracted from the
  Zhang counterexample.
- It strengthens a variance-only analysis of this structured frame because
  the exact law detects the relative basis angle while the variance proxy
  does not.
- It does not improve the repository's exact affine/matroid hashing laws;
  those concern collisions of Fourier frequencies under random linear maps.
- It does not improve the scalar independent-sum bound from
  arXiv:2607.23980; that theorem has different hypotheses and no
  operator-valued four-factor step.
- It does not use Zhang's false conjectured inequality. It uses only the
  valid rank-one projector configuration that exposed the failure.

## Verification

Run

```text
python verification/verify_projector_matrix_hypercube.py
```

The checker enumerates all sixteen sign vectors for several exact rational
two-basis frames, verifies the general atom formula and angle-independent
second moment, and separately recovers the Zhang matrices and their
\(0,\sqrt2,2\) norm distribution.

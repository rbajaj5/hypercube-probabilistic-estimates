# A physics-informed kernel filter on the Boolean cube

## Scope and source boundary

[PIKS](https://arxiv.org/abs/2607.27062) studies kernel regression augmented
by observations of a linear differential operator and proves universal
consistency in continuous, possibly misspecified settings.  This note records
a finite Boolean-cube analogue in which the operator is the graph Laplacian.
The simultaneous Walsh diagonalization makes the population estimator and its
risk completely explicit.

This is not a discretization of the paper's PDE experiments and does not
extend its consistency theorem.  On a finite cube, every strictly
positive-definite kernel already spans every function, so the paper's central
misspecification issue disappears.  The contribution here is an exact
spectral calibration of structural supervision.

## Kernel and operator

Let \(\Omega=\{0,1\}^d\), with uniform probability, and write

\[
\chi_S(x)=(-1)^{\sum_{i\in S}x_i},
\qquad
\widehat f(S)=2^{-d}\sum_x f(x)\chi_S(x).
\]

For \(0<\rho\leq1\), consider the Hamming product kernel

\[
K_\rho(x,y)
=\prod_{i=1}^d\left(1+\rho(-1)^{x_i+y_i}\right)
=\sum_{S\subseteq[d]}\rho^{|S|}\chi_S(x)\chi_S(y).
\tag{1}
\]

Its uniform integral operator has eigenvalue
\(\kappa_S=\rho^{|S|}>0\) on \(\chi_S\).  It is therefore strictly positive
definite and its RKHS is the entire \(2^d\)-dimensional function space, with

\[
\|f\|_{\mathcal H_\rho}^2
=\sum_S\frac{\widehat f(S)^2}{\kappa_S}.
\tag{2}
\]

Use the unnormalized cube Laplacian

\[
(Lf)(x)=\sum_{i=1}^d\bigl(f(x)-f(x\oplus e_i)\bigr).
\tag{3}
\]

Then

\[
L\chi_S=\mu_S\chi_S,\qquad \mu_S=2|S|.
\tag{4}
\]

Thus kernel smoothness and the structural operator are aligned by Walsh
degree.

## Exact population estimator

Let \(y\) be a value target and \(w\) an operator target.  For
\(\gamma\geq0\) and \(\lambda>0\), minimize

\[
\mathcal J(u)
=\|u-y\|_2^2
+\gamma\|Lu-w\|_2^2
+\lambda\|u\|_{\mathcal H_\rho}^2.
\tag{5}
\]

This is the full-population finite-cube analogue of a regularized
physics-informed objective.  Parseval separates (5) into scalar quadratics,
giving the unique solution

\[
\widehat u_{\lambda,\gamma}(S)
=
\frac{\widehat y(S)+\gamma\mu_S\widehat w(S)}
{1+\gamma\mu_S^2+\lambda/\kappa_S}.
\tag{6}
\]

No matrix inversion is required once the Walsh transform is available.

If \(y=f\) and \(w=Lf\) are noiseless, put

\[
a_S(\gamma)
=
\frac{1+\gamma\mu_S^2}
{1+\gamma\mu_S^2+\lambda/\kappa_S}.
\tag{7}
\]

Then \(\widehat u(S)=a_S(\gamma)\widehat f(S)\), and

\[
\|u-f\|_2^2
=
\sum_S
\left(
\frac{\lambda/\kappa_S}
{1+\gamma\mu_S^2+\lambda/\kappa_S}
\right)^2\widehat f(S)^2,
\tag{8}
\]

\[
\|L(u-f)\|_2^2
=
\sum_S\mu_S^2
\left(
\frac{\lambda/\kappa_S}
{1+\gamma\mu_S^2+\lambda/\kappa_S}
\right)^2\widehat f(S)^2.
\tag{9}
\]

Increasing \(\gamma\) strictly reduces the bias of every nonconstant active
mode.  It cannot help the constant mode, which lies in the nullspace of \(L\).
At fixed degree \(k>0\), the physics term changes the effective ridge penalty
from \(\lambda/\rho^k\) to

\[
\frac{\lambda/\rho^k}{1+4\gamma k^2}.
\tag{10}
\]

## Exact noisy-supervision tradeoff

Suppose, in one Walsh mode,

\[
\widehat y=\widehat f+\xi,\qquad
\widehat w=\mu\widehat f+\zeta,
\]

where the centered noises are independent with variances
\(\sigma_y^2,\sigma_w^2\).  Put \(r=\lambda/\kappa\).  Equation (6) gives

\[
\mathbb E(\widehat u-\widehat f)^2
=
\frac{r^2\widehat f^2+\sigma_y^2
+\gamma^2\mu^2\sigma_w^2}
{(1+\gamma\mu^2+r)^2}.
\tag{11}
\]

For a nonconstant mode and \(\sigma_w^2>0\), the exact risk-minimizing weight
is

\[
\gamma_*
=
\frac{r^2\widehat f^2+\sigma_y^2}
{\sigma_w^2(1+r)}.
\tag{12}
\]

The eigenvalue \(\mu\) controls how quickly the gain is realized but cancels
from the interior optimum.  With noiseless operator observations, risk
decreases toward zero as \(\gamma\to\infty\) for every nonconstant mode.  With
noisy observations, taking \(\gamma\) arbitrarily large leaves the noise floor
\(\sigma_w^2/\mu^2\); structural information should not be weighted as if it
were exact.

## What this improves

The general PIKS representer theorem produces a block linear system.  In this
specific full-cube, translation-invariant setting, (6) replaces that solve by
one scalar filter per Walsh set.  Equations (8)--(12) then distinguish:

- regularization bias from insufficient structural weight;
- the Laplacian's unavoidable constant-mode blind spot;
- exponential kernel suppression \(\rho^{|S|}\) at high degree; and
- variance injected by an imperfect structural oracle.

These are exact finite identities, not new general RKHS or PDE theorems.

## What is verified

`verification/verify_physics_informed_hypercube_kernel.py` uses exact rational
arithmetic through dimension five to verify:

1. the kernel expansion and every Walsh kernel eigenvalue;
2. every cube-Laplacian eigenvalue;
3. the coefficientwise normal equations for (6);
4. the noiseless bias and operator-residual formulas;
5. the noisy risk formula by exhaustive symmetric-noise enumeration; and
6. the stationary point (12) and its global minimality.

## Relation to the other modules

- [Hypercube walks and matrix mixing](HypercubeWalkMatrixMixing.md) use the
  same degree spectrum for a Markov generator.
- [Hopf-fiber projective kernels](HopfFiberHypercubeKernel.md) weight Hamming
  distances geometrically; the kernel here is chosen for strictly positive
  Walsh eigenvalues and regression.
- [KAN baselines](KANHypercubeBaselines.md) quantify approximation by Walsh
  degree.  Here all degrees are representable, but the kernel and operator
  impose different degree-dependent regularization.

# Reflection-positive pictures for hypercube correlations

## 1. Source and scope

Arthur Jaffe and Zhengwei Liu's
[*A Mathematical Picture Language
Program*](https://arthurjaffe.com/Assets/pdf/PictureLanguage.pdf)
advocates two ideas that fit the present repository:

1. Fourier transform turns one pictorial gluing operation into another; and
2. gluing a picture to its reflected copy should produce a positive
   quantity.

The article is a research program, not a source of a new hypercube or Hex
estimate.  It also presents lattice duality, the Ising model, and
Kramers--Wannier duality as directions and questions, not as new estimates
proved there.

This note gives a finite, exact simulation of the gluing principle.  The
"pictures" are Boolean observables, a seam is a coordinate-refresh
transition, and reflection positivity becomes positivity of an
autocorrelation Hankel matrix.  All results below are elementary consequences
of the Walsh spectrum; no literature-priority claim is made.

## 2. A doubled Boolean picture

Let \(P\) be the coordinate-refresh operator on
\(\Omega_d=\{-1,1\}^d\): choose a coordinate uniformly and replace it by a
fresh fair sign.  Under the uniform measure, \(P\) is self-adjoint and

\[
 P\chi_S=\lambda_{|S|}\chi_S,
 \qquad
 \lambda_k=1-\frac{k}{d}\in[0,1].
 \tag{1}
\]

For a real observable \(f\), put

\[
 C_t=\langle f,P^t f\rangle
 =\mathbb E[f(X_0)f(X_t)].
 \tag{2}
\]

Interpret \(P^i f\) and \(P^j f\) as two halves of a picture, with one half
reflected before gluing.  Their glued value is

\[
 C_{i+j}
 =
 \langle P^i f,P^j f\rangle.
 \tag{3}
\]

Therefore every finite Hankel matrix

\[
 H^{(m)}_a
 =
 \bigl[C_{a+i+j}\bigr]_{0\leq i,j<m}
 \tag{4}
\]

is positive semidefinite.  Indeed, when \(a=2b\),

\[
 \sum_{i,j}u_i u_jC_{2b+i+j}
 =
 \left\|\sum_i u_iP^{b+i}f\right\|_2^2\geq0.
\]

The Walsh formula below proves (4) for every integer \(a\geq0\), including
odd \(a\).

## 3. The exact Hausdorff-moment law

Use normalized Walsh masses

\[
 W_k(f)=\sum_{|S|=k}\widehat f(S)^2.
\]

Then

\[
 \boxed{
 C_t=\sum_{k=0}^d W_k(f)\lambda_k^t.}
 \tag{5}
\]

Thus \((C_t)\) is the moment sequence of the positive measure

\[
 \nu_f=\sum_{k=0}^d W_k(f)\,\delta_{\lambda_k}
 \quad\text{on }[0,1].
 \tag{6}
\]

This immediately yields the following exact constraints.

### Complete monotonicity

With \(\Delta C_t=C_{t+1}-C_t\),

\[
 \boxed{
 (-1)^r\Delta^rC_t
 =
 \sum_{k=0}^d
 W_k(f)\lambda_k^t(1-\lambda_k)^r
 \geq0.}
 \tag{7}
\]

In particular, \(C_t\geq C_{t+1}\geq0\), and every shifted difference
Hankel matrix \(H_a^{(m)}-H_{a+1}^{(m)}\) is also positive semidefinite.

### Log-convexity

Every \(2\times2\) Hankel minor is nonnegative:

\[
 \boxed{
 C_tC_{t+2}\geq C_{t+1}^2.}
 \tag{8}
\]

Whenever \(C_t>0\), the ratios \(C_{t+1}/C_t\) are nondecreasing.

### Spectral-complexity certificate

The rank of the infinite Hankel matrix is the number of Walsh levels with
nonzero mass.  More precisely,

\[
 \operatorname{rank}H_0^{(m)}
 =
 \min\!\left(
 m,\#\{k:W_k(f)>0\}
 \right).
 \tag{9}
\]

This follows from the Vandermonde factorization

\[
 H_0^{(m)}
 =
 \sum_{k:W_k>0}
 W_k
 (1,\lambda_k,\ldots,\lambda_k^{m-1})^{\mathsf T}
 (1,\lambda_k,\ldots,\lambda_k^{m-1}).
 \tag{10}
\]

Since the possible nodes \(\lambda_k=1-k/d\) are known, the first \(d+1\)
exact correlations \(C_0,\ldots,C_d\) determine every \(W_k\) by a
Vandermonde solve.  Sparse spectra need fewer moments through the usual
finite-moment/Prony reduction.  For noisy simulations, fitting nonnegative
\(W_k\)'s to (5) enforces all constraints (7)--(8) automatically.

## 4. Influence-only correlation bounds

Assume \(f\in\{-1,1\}\), let

\[
 \mu=\mathbb E f,\qquad
 V=1-\mu^2,\qquad
 I=\operatorname{Inf}(f)=\sum_k kW_k.
\]

The centered autocovariance

\[
 D_t=C_t-\mu^2=\sum_{k\geq1}W_k\lambda_k^t
\tag{11}
\]

satisfies

\[
 D_0=V,\qquad D_1=V-\frac{I}{d}.
 \tag{12}
\]

If \(V>0\), set \(a=D_1/V\).  If the first active Walsh level is \(q\),
set \(b=1-q/d\).  The moment measure for \(D_t/V\) is supported on
\([0,b]\) and has mean \(a\).  Jensen's inequality and
\(\lambda^t\leq b^{t-1}\lambda\) give, for every integer \(t\geq1\),

\[
 \boxed{
 V a^t
 \leq D_t
 \leq D_1 b^{t-1}.}
 \tag{13}
\]

For a balanced Boolean observable with \(q=1\), this becomes

\[
 \boxed{
 \left(1-\frac{I}{d}\right)^t
 \leq C_t
 \leq
 \left(1-\frac{I}{d}\right)
 \left(1-\frac1d\right)^{t-1}.}
 \tag{14}
\]

Consequently,

\[
 \boxed{
 \frac{1-
 (1-I/d)(1-1/d)^{t-1}}2
 \leq
 \Pr(f(X_t)\ne f(X_0))
 \leq
 \frac{1-(1-I/d)^t}{2}.}
 \tag{15}
\]

This is useful when total influence is known or estimable but the full
Walsh spectrum is not.

Two consecutive centered autocovariances sharpen the lower bound.  Put
\(r_t=D_t/D_{t-1}\) when \(D_{t-1}>0\).  Log-convexity and the spectral
support imply, for \(s\geq0\),

\[
\boxed{
 D_t r_t^s
 \leq D_{t+s}
 \leq D_t b^s.}
\tag{16}
\]

For a balanced observable, \(D_t=C_t\).

## 5. Exact Hex consequences

The Hex winner \(G_n\) from
[the Hex note](HexWinnerNoiseMixing.md) is balanced, and its first active
level is \(1\).  The exact influences already computed there therefore give
rigorous correlation intervals without needing the full spectrum.

| board | \(d\) | \(I\) | \(a=1-I/d\) |
|---|---:|---:|---:|
| \(2\times2\) | 4 | \(3/2\) | \(5/8\) |
| \(3\times3\) | 9 | \(249/128\) | \(301/384\) |
| \(4\times4\) | 16 | \(2405/1024\) | \(13979/16384\) |

At two refresh steps, (14) and exact enumeration give:

| board | lower bound \(a^2\) | exact \(C_2\) | upper bound \(a(1-1/d)\) |
|---|---:|---:|---:|
| \(2\times2\) | \(25/64\) | \(27/64\) | \(15/32\) |
| \(3\times3\) | \(90601/147456\) | \(2177/3456\) | \(301/432\) |
| \(4\times4\) | \(195412441/268435456\) | \(1544691/2097152\) | \(209685/262144\) |

For \(4\times4\), this adds a certified mixing interval based only on the
previously recorded influence, while the exact \(C_2\) supplies a new anchor
for simulations.  For the larger seeded boards, an influence confidence
interval can be propagated through (14)--(15); a point estimate alone should
not be presented as a rigorous bound.

## 6. Continuous-time/product-noise picture

For product noise, write \(\rho=e^{-\tau}\) and

\[
 R_f(\tau)
 =
 \operatorname{Stab}_{e^{-\tau}}(f)
 =
 \sum_kW_k e^{-k\tau}.
 \tag{17}
\]

The semigroup identity \(T_\rho T_\sigma=T_{\rho\sigma}\) is the Boolean
version of gluing two seams.  Walsh transform changes this convolution into
pointwise multiplication.  Reflection positivity gives, for arbitrary
\(\tau_1,\ldots,\tau_m\geq0\),

\[
 \bigl[R_f(\tau_i+\tau_j)\bigr]_{i,j=1}^m\succeq0.
 \tag{18}
\]

It also gives

\[
 (-1)^rR_f^{(r)}(\tau)\geq0,
 \qquad
 R_f\!\left(\frac{s+t}{2}\right)^2
 \leq R_f(s)R_f(t).
 \tag{19}
\]

These are exact consistency tests for every product-noise curve in the
repository.

## 7. What this does not prove about Ising

The picture-language article points toward lattice symmetry,
Kramers--Wannier duality, and reflection-positive continuum limits.  Those
remarks do not by themselves imply a new critical temperature, an Ising
scaling limit for Hex, or a crossing-probability exponent.  Establishing any
of those would require a specified Gibbs measure, boundary conditions, and
an actual duality or comparison theorem.

The finite moment law above is narrower but complete: it gives exact
positivity and interpolation constraints for the refresh and product-noise
models already used here.

## 8. Verification

Run

```text
python verification/verify_reflection_positive_hypercube.py
```

The checker uses only exact integer and rational arithmetic.  It exhausts
every Boolean observable through dimension \(3\), verifies finite
differences, Hankel and shifted-Hankel positivity, log-convex and
influence-only bounds, reconstructs the small-Hex Walsh masses from their
moments, and checks the displayed Hex values through \(4\times4\).

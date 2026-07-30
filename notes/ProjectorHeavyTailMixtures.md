# Heavy-tail scale mixtures of the projector hypercube

## 1. Why the finite cube is not itself heavy-tailed

The
[two-basis projector law](ProjectorMatrixHypercubeLaw.md)
computes the operator norm of

```math
S_q
=
\varepsilon_AA+\varepsilon_BB+\varepsilon_CC+\varepsilon_DD,
```

where the four signs are independent and \(q\in[-1,1]\) records the
relative angle between two orthonormal projector bases. Its norm

```math
W_q:=\|S_q\|
```

is supported in \([0,2]\), so the finite hypercube cannot create an
asymptotic heavy tail by itself.

Heavy tails appear after a radial scale mixture. The exact angular law then
determines the power-law constant rather than merely supplying a generic
moment bound. This note gives two such calculations:

1. a Pareto mixture of the matrix-hypercube observable; and
2. common versus independent Pareto shocks in Zhang's quartic matrix
   functional.

Every formula below is exact for a standard Pareto scale. No asymptotic
regular-variation theorem is required.

## 2. Exact Pareto mixture of the hypercube law

Let \(R\) be independent of the signs and have the standard Pareto tail

```math
\mathbb P(R>r)=r^{-\alpha},
\qquad r\geq1,
\qquad \alpha>0.
\tag{1}
```

The exact projector-hypercube law is

```math
\mathcal L(W_q)
=
\frac18\delta_0
+\frac58\delta_2
+\frac18\delta_{\sqrt{2+2q}}
+\frac18\delta_{\sqrt{2-2q}}.
\tag{2}
```

### Theorem 1 (exact radial black-swan coefficient)

For every \(t\geq2\),

```math
\boxed{
\mathbb P(RW_q>t)
=
C_\alpha(q)t^{-\alpha},}
\tag{3}
```

where

```math
C_\alpha(q)
=
\frac58\,2^\alpha
+\frac18(2+2q)^{\alpha/2}
+\frac18(2-2q)^{\alpha/2}.
\tag{4}
```

Thus the projector geometry does not change the Pareto exponent
\(\alpha\), but it determines the extreme-event coefficient exactly.

### Proof

Condition on \(W_q=w\). Because \(0\leq w\leq2\) and \(t\geq2\), equation
(1) applies at \(t/w\) whenever \(w>0\). Hence

```math
\mathbb P(RW_q>t\mid W_q=w)
=
\left(\frac wt\right)^\alpha.
```

The conditional probability is zero when \(w=0\). Averaging and using (2)
gives

```math
\mathbb P(RW_q>t)
=
t^{-\alpha}\mathbb E W_q^\alpha
=
C_\alpha(q)t^{-\alpha}.
```

\(\square\)

## 3. A phase transition at the second-moment boundary

The angle-dependent part of (4) is a positive multiple of

```math
g_\beta(q)
=(1+q)^\beta+(1-q)^\beta,
\qquad
\beta=\frac{\alpha}{2}.
```

On \((-1,1)\),

```math
g_\beta''(q)
=
\beta(\beta-1)
\left((1+q)^{\beta-2}+(1-q)^{\beta-2}\right).
```

Therefore:

- if \(0<\alpha<2\), then \(g_\beta\) is strictly concave, and
  \(C_\alpha(q)\) is maximized at \(q=0\);
- if \(\alpha=2\), then \(C_2(q)=3\) for every \(q\); and
- if \(\alpha>2\), then \(g_\beta\) is strictly convex, and
  \(C_\alpha(q)\) is minimized at \(q=0\).

The mutually unbiased Zhang frame is exactly \(q=0\). Coincident bases,
up to permutation, have \(|q|=1\). Their constants are

```math
C_\alpha(0)
=
\frac58\,2^\alpha+\frac14\,2^{\alpha/2},
\tag{5}
```

and

```math
C_\alpha(\pm1)
=
\frac34\,2^\alpha.
\tag{6}
```

In the strictly sub-quadratic regime \(0<\alpha<2\), where the Pareto
variance is infinite, the exact angular
amplification of the Zhang frame relative to coincident bases is

```math
\frac{C_\alpha(0)}{C_\alpha(\pm1)}
=
\frac56+\frac13\,2^{-\alpha/2},
\tag{7}
```

which has limit \(7/6\) as \(\alpha\downarrow0\) and decreases to \(1\) as
\(\alpha\uparrow2\). This is the precise black-swan contribution of the
mutually unbiased geometry to the radial mixture: the increase in the tail
constant approaches \(16\frac23\%\), but the geometry does not manufacture
the power law.

The reversal for \(\alpha>2\) is also informative. High positive moments
favor the coincident frame because one of its two angle-sensitive atoms
merges with the endpoint \(2\), while the other collapses to zero.

## 4. Zhang's quartic functional under a common shock

For positive semidefinite matrices define

```math
\Phi(A,B,C,D)
=
A(BC+CB)D+D(BC+CB)A.
\tag{8}
```

For the Zhang projectors \(P_A,P_B,P_C,P_D\),

```math
\|\Phi(P_A,P_B,P_C,P_D)\|
=
\kappa,
\qquad
\kappa:=\frac{1+\sqrt2}{4}.
\tag{9}
```

Apply one common Pareto shock to all four matrices:

```math
A_R=RP_A,\quad B_R=RP_B,\quad
C_R=RP_C,\quad D_R=RP_D.
```

The functional is homogeneous of degree four, so

```math
\|\Phi(A_R,B_R,C_R,D_R)\|
=
\kappa R^4.
\tag{10}
```

### Theorem 2 (common-shock exponent quartering)

For \(t\geq\kappa\),

```math
\boxed{
\mathbb P\!\left(
\|\Phi(A_R,B_R,C_R,D_R)\|>t
\right)
=
\kappa^{\alpha/4}t^{-\alpha/4}.}
\tag{11}
```

Thus a systemic radial shock with exponent \(\alpha\) produces a quartic
observable with exponent \(\alpha/4\).

The right-hand side proposed in Zhang's conjecture becomes

```math
G_R
:=
\frac1{64}\|A_R+B_R+C_R+D_R\|^4
=
\frac14R^4.
\tag{12}
```

For \(t\geq\kappa\), equations (11)--(12) give the exact exceedance ratio

```math
\boxed{
\frac{
\mathbb P(\|\Phi(A_R,B_R,C_R,D_R)\|>t)
}{
\mathbb P(G_R>t)
}
=(1+\sqrt2)^{\alpha/4}.}
\tag{13}
```

The tail-index change in (11) comes from quartic homogeneity and the shared
shock. The noncommutative counterexample contributes the exact multiplicative
amplification in (13).

## 5. Four independent shocks behave differently

Now let \(R_A,R_B,R_C,R_D\) be independent standard Pareto variables with
the same exponent \(\alpha\), and scale the four Zhang projectors
separately. Multilinearity gives

```math
\left\|
\Phi(R_AP_A,R_BP_B,R_CP_C,R_DP_D)
\right\|
=
\kappa R_AR_BR_CR_D.
\tag{14}
```

For \(x\geq1\), the logarithms \(\log R_i\) are independent exponential
variables with rate \(\alpha\). Their sum has the four-stage Erlang
survival function

```math
\mathbb P(R_AR_BR_CR_D>x)
=
x^{-\alpha}
\sum_{j=0}^{3}
\frac{(\alpha\log x)^j}{j!}.
\tag{15}
```

Consequently, for \(t\geq\kappa\),

```math
\boxed{
\mathbb P\!\left(
\left\|\Phi(R_AP_A,R_BP_B,R_CP_C,R_DP_D)\right\|>t
\right)
=
\left(\frac{\kappa}{t}\right)^\alpha
\sum_{j=0}^{3}
\frac{\left(\alpha\log(t/\kappa)\right)^j}{j!}.}
\tag{16}
```

The independent-shock tail retains exponent \(\alpha\), with a cubic
logarithmic enhancement. This contrasts with the exponent \(\alpha/4\) in
the common-shock model. In this example, the strongest black-swan
amplification is therefore a dependence phenomenon: four simultaneous
exposures to one systemic shock are much heavier-tailed than four
independent shocks.

## 6. Scope

- The bounded sign cube supplies an exact angular multiplier, not a heavy
  tail by itself.
- The Pareto mixture turns the previously computed norm moments into exact
  tail coefficients.
- The threshold \(\alpha=2\) is a genuine geometric phase transition:
  the Zhang frame maximizes the coefficient below it, is angle-neutral at
  it, and minimizes the coefficient above it.
- The \(\alpha/4\) common-shock exponent is caused by degree-four
  homogeneity and perfect radial dependence. It should not be attributed to
  noncommutativity alone.
- No claim of a new general theorem in regular variation, matrix
  concentration, or systemic-risk theory is made.

## Verification

Run

```text
python verification/verify_projector_heavy_tails.py
```

The checker verifies the exact hypercube atom counts, the angle-independent
second moment, the even-moment phase transition on rational projector
frames, common- and independent-scale homogeneity of the quartic functional,
and the polynomial identity behind the four-stage Erlang survival law.

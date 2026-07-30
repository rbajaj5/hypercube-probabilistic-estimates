# Novelty perturbations of an optimized hypercube

## 1. Motivation and boundary

Rabindranath Tagore's *The Fugitive III*, poem 26, is often presented as
[*A Wrong Man in Workers'
Paradise*](https://www.gutenberg.org/ebooks/7971).  An artist introduces
objects with no instrumental purpose into a society organized entirely
around useful work.  The objects alter attention and preference even though
the old objective function assigns them no value.

The poem is motivation, not mathematical evidence.  This note isolates one
testable abstraction: start with an absorbing all-work state, inject rare
novelty through coordinate refreshes, and quantify the resulting escape,
entropy, cost, and mixing.

## 2. The novelty-refresh chain

Let

```math
\Omega_d=\{0,1\}^d,
```

where \(x_i=0\) means that coordinate \(i\) is in its productive baseline
and \(x_i=1\) marks a novelty.  Fix selection probabilities

```math
q_i>0,\qquad \sum_{i=1}^d q_i=1,
```

and novelty probabilities \(0<\nu_i<1\).  One step of the chain is:

1. choose coordinate \(I=i\) with probability \(q_i\);
2. replace \(X_i\) by an independent
   \(\operatorname{Bernoulli}(\nu_i)\); and
3. leave all other coordinates unchanged.

At \(\nu_i=0\) for every \(i\), the all-zero state is absorbing.  Positive
\(\nu_i\)'s perturb this optimized but degenerate system into an irreducible
heat-bath chain.

## 3. Stationary law and exact spectrum

The stationary law is the product measure

```math
\boxed{
\pi_\nu(x)
=
\prod_{i=1}^d
\nu_i^{x_i}(1-\nu_i)^{1-x_i}.}
\tag{1}
```

Each single-coordinate refresh is reversible with respect to its Bernoulli
factor, so their \(q_i\)-weighted mixture is reversible with respect to
\(\pi_\nu\).

Define centered coordinate functions

```math
\phi_i(x)=
\frac{x_i-\nu_i}{\sqrt{\nu_i(1-\nu_i)}},
\qquad
\phi_S(x)=\prod_{i\in S}\phi_i(x).
\tag{2}
```

If the refreshed coordinate lies in \(S\), its conditional centered mean is
zero; otherwise \(\phi_S\) is unchanged.  Hence

```math
\boxed{
P\phi_S=
\left(1-\sum_{i\in S}q_i\right)\phi_S.}
\tag{3}
```

The \(2^d\) functions \(\phi_S\) form an orthonormal basis of
\(L^2(\pi_\nu)\).  Therefore the complete spectrum is

```math
\left\{
1-\sum_{i\in S}q_i:S\subseteq[d]
\right\},
\tag{4}
```

and the spectral gap is

```math
\boxed{\operatorname{gap}(P)=\min_iq_i.}
\tag{5}
```

The novelty probabilities change the stationary law and its entropy, but
not these eigenvalues.  In the uniform case \(q_i=1/d\), level \(k\) still
has eigenvalue \(1-k/d\), exactly as for the fair coordinate-refresh walk.

## 4. Escape and finite-time novelty

Let \(0^d\) be the all-work state and

```math
\tau=\min\{t\geq1:X_t\ne0^d\}.
```

While the chain remains at \(0^d\), one step creates a novelty with
probability

```math
r=\sum_{i=1}^dq_i\nu_i.
```

Thus

```math
\boxed{
\Pr_{0^d}(\tau>t)=(1-r)^t,
\qquad
\mathbb E_{0^d}\tau=\frac1r.}
\tag{6}
```

For coordinate \(i\), the probability that it has not yet been refreshed
after \(t\) steps is \((1-q_i)^t\).  Starting from \(0^d\),

```math
\boxed{
\mathbb E_{0^d}|X_t|
=
\sum_{i=1}^d
\nu_i\left[1-(1-q_i)^t\right].}
\tag{7}
```

At stationarity, \(K=|X|\) has the Poisson-binomial generating function

```math
\mathbb E_{\pi_\nu}z^K
=
\prod_{i=1}^d(1-\nu_i+\nu_i z),
\tag{8}
```

so

```math
\mathbb E_{\pi_\nu}K=\sum_i\nu_i,\qquad
\Pr_{\pi_\nu}(K>0)=1-\prod_i(1-\nu_i).
\tag{9}
```

Equation (6) measures the first disruptive event; (9) measures the
long-run prevalence of any novelty.  They are different quantities.

## 5. Entropy versus productive cost

The Shannon entropy of the stationary law is

```math
\boxed{
H(\pi_\nu)
=
\sum_{i=1}^d h(\nu_i),\qquad
h(u)=-u\log u-(1-u)\log(1-u).}
\tag{10}
```

If the old objective charges one unit for every novelty coordinate, its
stationary expected cost is

```math
C(\pi_\nu)=\sum_i\nu_i.
\tag{11}
```

For the uniform choice \(\nu_i=\nu\),

```math
H(\pi_\nu)=d\,h(\nu),\qquad C(\pi_\nu)=d\nu,
\tag{12}
```

and

```math
\frac{H(\pi_\nu)}{C(\pi_\nu)}
=
\frac{h(\nu)}{\nu}
=
\log\frac1\nu+1+O(\nu)
\quad(\nu\downarrow0).
\tag{13}
```

Rare novelty therefore produces a large entropy gain per unit of expected
old-objective cost.  This is an entropy statement, not a claim that entropy
is identical to aesthetic or human value.

## 6. Exact mixing from the all-work state

By the reversible spectral expansion and (3),

```math
\boxed{
\chi^2\!\left(P^t(0^d,\cdot)\Vert\pi_\nu\right)
=
\sum_{\varnothing\ne S\subseteq[d]}
\left[
\prod_{i\in S}\frac{\nu_i}{1-\nu_i}
\right]
\left(1-\sum_{i\in S}q_i\right)^{2t}.}
\tag{14}
```

Consequently,

```math
\left\|P^t(0^d,\cdot)-\pi_\nu\right\|_{\rm TV}
\leq
\frac12\sqrt{\chi^2(P^t(0^d,\cdot)\Vert\pi_\nu)}.
\tag{15}
```

There is also a coupling bound.  Once every coordinate has been refreshed,
the state is exactly distributed as \(\pi_\nu\).  A union bound gives

```math
\boxed{
\left\|P^t(x,\cdot)-\pi_\nu\right\|_{\rm TV}
\leq
\sum_{i=1}^d(1-q_i)^t.}
\tag{16}
```

For uniform selection this becomes \(d(1-1/d)^t\).  The coupon-collector
bound ignores how rare novelty is, while (14) retains the \(\nu_i\)
dependence and can be much sharper near the all-work state.

## 7. Interpretation

This model separates four effects that can otherwise be conflated:

1. **escape:** rare novelty breaks the absorbing baseline on the
   \(1/r\) time scale;
2. **prevalence:** the invariant number of novelties is Poisson-binomial;
3. **diversity:** stationary entropy increases by \(\sum_i h(\nu_i)\);
4. **relaxation:** the spectral gap remains \(\min_iq_i\).

The last point is the useful surprise.  For heat-bath refreshes, changing
what a coordinate prefers does not change how fast coordinate information
is forgotten.  Different update mechanisms—copying, strategic response,
or reinforcement—need not share this invariance.

## Verification

Run

```text
python verification/verify_novelty_refresh_hypercube.py
```

The checker uses exact rational arithmetic.  It exhausts five nonuniform
chains through dimension five and verifies stochasticity, detailed
balance, stationarity, every eigenfunction in (3), escape and finite-time
means, the Poisson-binomial law, the chi-squared identity, total-variation
comparison, and the coupon-collector bound.

If the novelty rates share a hidden common cause rather than remaining
independent, see [Latent environments on a hypercube](LatentEnvironmentHypercube.md)
for the resulting overdispersion and higher-order centered interactions.

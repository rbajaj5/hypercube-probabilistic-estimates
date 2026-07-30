# Hard-core heat-bath mixing on conflict hypercubes

## Scope

The [fixed-drawing plane-saturation module](PlaneSaturationConflictHypercube.md)
turns prescribed geometric edges into vertices of a crossing-conflict graph
\(C\).  Plane drawings are independent sets of \(C\), and fixed-drawing
saturated states are maximal independent sets.

This note adds dynamics.  It studies the standard single-site heat-bath walk
on the independent sets of an arbitrary finite conflict graph.  The result is
a rigorous sampler for the Bernoulli edge cube conditioned on planarity, with
an explicit path-coupling mixing estimate in the low-activity regime.

The walk samples all plane states, not only saturated ones.  The probability
that a stationary sample is saturated remains the exact ratio
\(M_C(\lambda)/Z_C(\lambda)\) from the preceding module.

## 1. The heat-bath chain

Let \(C\) have \(m\geq1\) vertices and maximum degree \(\Delta\).  Write
\(\mathcal I(C)\) for its independent sets.  Fix an activity \(\lambda>0\)
and put

\[
\rho=\frac{\lambda}{1+\lambda}.
\tag{1}
\]

From \(I\in\mathcal I(C)\):

1. choose \(v\in V(C)\) uniformly;
2. if \(N_C(v)\cap I\ne\varnothing\), set \(v\) absent;
3. otherwise set \(v\) present with probability \(\rho\) and absent with
   probability \(1-\rho\).

Only one Boolean coordinate is refreshed, and every intermediate state
remains independent.

## 2. Exact stationary law

Let

\[
Z_C(\lambda)=\sum_{I\in\mathcal I(C)}\lambda^{|I|}.
\]

### Proposition 1

The chain is irreducible, aperiodic, reversible, and has stationary law

\[
\boxed{\mu_\lambda(I)=\frac{\lambda^{|I|}}{Z_C(\lambda)}.}
\tag{2}
\]

Irreducibility follows by deleting the current set and then adding the target
set.  The chain has a self-loop at every state.  For adjacent independent
sets \(I\) and \(I\cup\{v\}\),

\[
\mu_\lambda(I)\frac{1}{m}\rho
=
\mu_\lambda(I\cup\{v\})\frac{1}{m}(1-\rho),
\tag{3}
\]

because \(\rho/(1-\rho)=\lambda\).  All other detailed-balance identities
are trivial.

Equation (2) is precisely the product Bernoulli-\(\rho\) edge law conditioned
on having no conflicts.

## 3. A local drift identity

For an independent set \(I\), define the number of eligible coordinates

\[
a(I)=
\left|
\{v\in V(C):N_C(v)\cap I=\varnothing\}
\right|.
\tag{4}
\]

This includes every currently occupied vertex, since \(I\) is independent.
The one-step occupancy drift is

\[
\boxed{
\mathbb E\!\left[|X_{t+1}|-|X_t|\mid X_t=I\right]
=\frac{\rho\,a(I)-|I|}{m}.}
\tag{5}
\]

At stationarity the mean drift vanishes, giving

\[
\boxed{\mathbb E_{\mu_\lambda}|I|
=\rho\,\mathbb E_{\mu_\lambda}a(I).}
\tag{6}
\]

This balances deletion pressure against the number of currently unblocked
edges.

Starting from the empty drawing, every coordinate is eligible.  Therefore

\[
\Pr(\varnothing\to\varnothing)=1-\rho,
\]

and the first time at which any edge appears is geometric with success
probability \(\rho\):

\[
\boxed{\Pr(T_{\rm first}>t)=(1-\rho)^t,\qquad
\mathbb ET_{\rm first}=\frac1\rho.}
\tag{7}
\]

## 4. Adjacent-state coupling

Couple two copies by choosing the same coordinate and using the same uniform
random number for their heat-bath decisions.  Suppose two independent sets
differ at exactly one coordinate \(v\).

- If \(v\) is selected, the two copies coalesce.
- If a neighbor of \(v\) is selected, at most one new disagreement is
  created, with probability at most \(\rho\).
- Any other update preserves the single disagreement.

Hence the expected Hamming distance after one step is at most

\[
\left(
1-\frac{1-\Delta\rho}{m}
\right)d_{\rm H}(X_0,Y_0).
\tag{8}
\]

The independent-set state graph is connected by single-coordinate moves,
and its graph distance is Hamming distance: remove \(I\setminus J\), then add
\(J\setminus I\).  Path coupling therefore extends (8) to arbitrary initial
states.

## 5. Mixing estimate

Assume the Dobrushin condition

\[
\boxed{\Delta\rho
=\Delta\frac{\lambda}{1+\lambda}<1.}
\tag{9}
\]

Put

\[
\kappa=1-\frac{1-\Delta\rho}{m}<1.
\]

Coupling a chain started at \(x\) with a stationary chain and using the
diameter bound \(d_{\rm H}\leq m\) gives

\[
\boxed{
\max_x\|P^t(x,\cdot)-\mu_\lambda\|_{\rm TV}
\leq m\kappa^t
\leq
m\exp\!\left(
-\frac{1-\Delta\rho}{m}t
\right).}
\tag{10}
\]

In particular,

\[
\boxed{
t_{\rm mix}(\varepsilon)
\leq
\left\lceil
\frac{m}{1-\Delta\rho}
\left(\log m+\log\frac1\varepsilon\right)
\right\rceil.}
\tag{11}
\]

For \(\Delta>1\), condition (9) is equivalent to
\(\lambda<1/(\Delta-1)\).  It is sufficient, not necessary: many conflict
graphs mix rapidly at larger activity, but (10) makes no such claim.

## 6. Geometric interpretation

For a fixed prescribed drawing:

- a heat-bath update proposes adding or deleting one geometric edge;
- a proposed addition is suppressed exactly when it would cross a selected
  prescribed edge;
- the invariant distribution is Bernoulli edge sampling conditioned on
  planarity;
- increasing \(\lambda\) favors denser plane drawings but eventually leaves
  the elementary contraction regime (9).

This separates two sampling objectives:

1. the random-order greedy process from the preceding note always returns a
   saturated drawing;
2. the heat-bath walk has a tractable reversible law and returns a saturated
   drawing with stationary probability \(M_C(\lambda)/Z_C(\lambda)\).

## 7. Verification and novelty boundary

Run

```text
python verification/verify_hard_core_conflict_walk.py
```

The checker uses exact rational arithmetic.  On abstract conflict graphs and
the crossing graphs of convex \(K_4,K_5,K_6\), it verifies:

1. stochasticity, irreducibility, stationarity, and every detailed-balance
   identity;
2. the drift and first-escape formulas (5)--(7);
3. every adjacent-state common-uniform coupling and bound (8);
4. exact finite-time total variation against (10);
5. the stationary hard-core law and saturation probability inherited from
   the static module.

These are standard finite hard-core Glauber/path-coupling facts specialized
to crossing-conflict hypercubes.  No new general mixing theorem, optimal
threshold, or result about the unlabeled plane-saturation ratio is claimed.

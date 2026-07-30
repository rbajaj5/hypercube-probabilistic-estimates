# Fixed-drawing plane saturation as a conflict hypercube

## Scope and source boundary

Clifton and Salia introduced the plane-saturation ratio for a plane subgraph
of a planar host graph in
[arXiv:2403.02458](https://arxiv.org/abs/2403.02458).  Their definition permits
the unlabeled subgraph to be embedded into the host in different ways.  That
freedom is essential to their sharp \(1/16\) theorem.

This note studies a different, deliberately simpler model:

- the vertices are labeled and placed in advance;
- every allowed edge has a prescribed arc (straight chords are the main
  example);
- an edge may be accepted only with that prescribed arc.

Assume the prescribed arcs are simple, do not overlap, and meet in their
interiors only through crossings, so that planarity is determined pairwise.

In this fixed-drawing model, crossing is a pairwise conflict and the state
space is an ordinary Boolean edge cube.  The resulting identities do not
prove, reprove, or improve the \(1/16\) theorem.

## 1. The crossing-conflict graph

Let \(\mathcal E=\{e_1,\ldots,e_m\}\) be the allowed prescribed edges.  Define
the crossing-conflict graph \(C\) by

\[
V(C)=\mathcal E,\qquad
e_i e_j\in E(C)
\quad\Longleftrightarrow\quad
e_i\text{ and }e_j\text{ cross in their interiors}.
\tag{1}
\]

A drawing is encoded by \(x\in\{0,1\}^{m}\), or equivalently by
\(A=\{e_i:x_i=1\}\).

### Proposition 1

The selected drawing is plane if and only if \(A\) is an independent set of
\(C\).  It is plane-saturated with respect to the fixed drawing if and only if
\(A\) is a maximal independent set of \(C\).

The proof is immediate: a crossing is exactly a conflict edge of \(C\), and a
missing prescribed edge can be added exactly when it has no selected
neighbor in \(C\).

This reduction is useful because it separates two issues:

1. geometry constructs the conflict graph \(C\);
2. all subsequent probability calculations take place on
   \(\{0,1\}^{V(C)}\).

## 2. Bernoulli edge sampling

Select every prescribed edge independently with probability \(p\), where
\(0<p<1\), and put

\[
\lambda=\frac{p}{1-p}.
\]

Let

\[
Z_C(\lambda)
=\sum_{\substack{I\subseteq V(C)\\I\ {\rm independent}}}
\lambda^{|I|}
\tag{2}
\]

be the independence polynomial, and let

\[
M_C(\lambda)
=\sum_{\substack{I\subseteq V(C)\\I\ {\rm maximal\ independent}}}
\lambda^{|I|}
\tag{3}
\]

be the maximal-independence polynomial.

### Theorem 2

For the Bernoulli-\(p\) edge cube,

\[
\boxed{\Pr(\text{plane})=(1-p)^m Z_C(\lambda),}
\tag{4}
\]

\[
\boxed{\Pr(\text{plane and saturated})=(1-p)^m M_C(\lambda),}
\tag{5}
\]

and therefore

\[
\boxed{\Pr(\text{saturated}\mid\text{plane})
=\frac{M_C(\lambda)}{Z_C(\lambda)}.}
\tag{6}
\]

Conditioned on planarity, the selected set has the hard-core law

\[
\Pr(A=I\mid A\text{ plane})
=\frac{\lambda^{|I|}}{Z_C(\lambda)}
\qquad(I\text{ independent}).
\tag{7}
\]

Consequently,

\[
\boxed{
\mathbb E(|A|\mid A\text{ plane})
=\lambda\frac{Z_C'(\lambda)}{Z_C(\lambda)}.}
\tag{8}
\]

These follow by factoring
\(p^{|I|}(1-p)^{m-|I|}=(1-p)^m\lambda^{|I|}\).

At \(p=1/2\), (4)--(6) simply count independent and maximal independent
sets and divide by \(2^m\).

## 3. Boolean indicator formulas

The planarity indicator on the full edge cube is

\[
\mathbf 1_{\rm plane}(x)
=
\prod_{\{i,j\}\in E(C)}(1-x_ix_j).
\tag{9}
\]

For a plane state, a missing edge \(i\) is blocked exactly when at least one
of its conflict neighbors is present.  Hence the saturation indicator is

\[
\mathbf 1_{\rm sat}(x)
=
\mathbf 1_{\rm plane}(x)
\prod_{i=1}^{m}
\left[
x_i+(1-x_i)
\left(
1-\prod_{j\in N_C(i)}(1-x_j)
\right)
\right].
\tag{10}
\]

For an isolated vertex of \(C\), the corresponding factor in (10) is
\(x_i\): a prescribed edge that crosses nothing must be present in every
saturated drawing.

Equations (9)--(10) are exact Boolean polynomials.  They make the model
available to the Walsh/noise machinery used elsewhere in this repository,
although no low-degree concentration is asserted here.

## 4. Random-order greedy saturation

Give the \(m\) prescribed edges a uniformly random order.  Visit them in
that order and accept an edge if it conflicts with none of the previously
accepted edges.  The output \(G_C\) is always a maximal independent set, so
it is always a fixed-drawing saturated state.

For \(i\in V(C)\), if \(i\) is earliest in its closed neighborhood
\(N_C[i]\), it is certainly accepted.  Thus

\[
\Pr(i\in G_C)\geq\frac{1}{d_C(i)+1}
\]

and

\[
\boxed{
\mathbb E|G_C|
\geq
\sum_{i\in V(C)}\frac{1}{d_C(i)+1}
\geq
\frac{m^2}{m+2|E(C)|}
\geq
\frac{m}{\Delta(C)+1}.}
\tag{11}
\]

The middle inequality is Cauchy--Schwarz.  Independently of the random
order, every maximal independent set \(I\) satisfies

\[
\boxed{|I|\geq\frac{m}{\Delta(C)+1},}
\tag{12}
\]

because the closed neighborhoods of the vertices of \(I\) cover \(V(C)\).

There is also an exact finite recursion for the greedy output distribution.
If \(R_S(y)=\mathbb E y^{|G_{C[S]}|}\), then

\[
R_\varnothing(y)=1,\qquad
R_S(y)
=\frac{y}{|S|}
\sum_{i\in S}R_{S\setminus N_C[i]}(y).
\tag{13}
\]

Indeed, the first edge in the random order is uniform in \(S\), is accepted,
and deletes its closed conflict neighborhood.  Equation (13) permits exact
dynamic programming without enumerating all \(m!\) orders.

## 5. Component factorization

If \(C=C_1\sqcup\cdots\sqcup C_r\), then

\[
Z_C(\lambda)=\prod_{a=1}^{r}Z_{C_a}(\lambda),
\qquad
M_C(\lambda)=\prod_{a=1}^{r}M_{C_a}(\lambda).
\tag{14}
\]

The greedy generating polynomial also factors:

\[
R_C(y)=\prod_{a=1}^{r}R_{C_a}(y).
\tag{15}
\]

For example, if \(C\) consists of \(r\) disjoint conflict pairs and \(s\)
isolated vertices, then

\[
Z_C(\lambda)=(1+2\lambda)^r(1+\lambda)^s,
\qquad
M_C(\lambda)=(2\lambda)^r\lambda^s.
\tag{16}
\]

Every greedy drawing then contains exactly \(r+s\) edges.

## 6. Convex chord examples

Place \(n\) labeled points on a circle and prescribe straight chords.  Two
chords with four distinct endpoints conflict exactly when their endpoints
alternate cyclically.

- For \(K_4\) drawn on a convex quadrilateral, the two diagonals form one
  conflict pair and the four boundary edges are isolated in \(C\).
- For \(K_5\) drawn on a convex pentagon, the five diagonals induce a
  \(5\)-cycle in \(C\), while the five boundary edges are isolated.

Thus the geometry can be compressed into a sparse conflict graph before any
hypercube enumeration is performed.

## 7. Verification and limitations

Run

```text
python verification/verify_plane_saturation_conflict_hypercube.py
```

The checker uses exact integer and rational arithmetic.  It constructs
conflict graphs for convex complete drawings through six vertices and
additional abstract examples, exhausts their Boolean edge cubes, verifies
(4)--(10), computes the exact greedy law from (13), checks (11)--(12), and
tests the component formulas.

This module is classical independent-set probability in geometric clothing.
It makes no priority claim.  In particular:

- it does not handle arbitrary rerouting of a newly added edge;
- it does not handle the unlabeled re-embeddings in Clifton--Salia's
  definition;
- it does not improve their sharp \(1/16\) constant;
- it does not turn the partisan game of Hex into plane saturation.

The later maximal-planar study
[arXiv:2412.06068](https://arxiv.org/abs/2412.06068) also distinguishes
labeled and unlabeled notions.  The present prescribed-arc model is stricter
than both of those plane-graph notions.

For a reversible single-edge sampler of the plane states, see
[Hard-core heat-bath mixing on conflict hypercubes](HardCoreConflictWalkMixing.md).

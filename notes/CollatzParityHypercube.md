# Collatz parity cubes: maximal, Gaussian, and residue estimates

## Status and scope

This note gives exact finite-horizon statements for the accelerated Collatz
map, followed by classical probabilistic consequences.  The Boolean-cube
identification is classical; the maximal estimate comes directly from
Kolmogorov's inequality, and the normal approximation comes directly from
Berry--Esseen.  The finite-field calculation is an elementary transfer-chain
corollary.

None of these statements proves that every positive Collatz orbit reaches
one.  In particular, a theorem for a uniformly random residue class, Haar
random \(2\)-adic seed, or fixed-horizon prime sample must not be promoted to
a theorem about every fixed positive integer.

## 1. The exact parity hypercube

Use the accelerated map

```math
T(n)=
\begin{cases}
n/2,&n\equiv0\pmod2,\\
(3n+1)/2,&n\equiv1\pmod2.
\end{cases}
\tag{1}
```

For a horizon \(m\), define

```math
\epsilon_j(n)=T^j(n)\bmod2,\qquad
s_k(n)=\sum_{j=0}^{k-1}\epsilon_j(n).
\tag{2}
```

The parity-vector map

```math
n\bmod 2^m
\longmapsto
(\epsilon_0(n),\ldots,\epsilon_{m-1}(n))
\tag{3}
```

is a bijection from \(\mathbb Z/2^m\mathbb Z\) to
\(\{0,1\}^m\).  This is the finite quotient of the Bernstein--Lagarias
\(2\)-adic conjugacy between \(T\) and the shift.  Consequently, if the
starting residue is uniform modulo \(2^m\), its first \(m\) parity bits are
*exactly* independent fair Bernoulli variables.  No heuristic independence
assumption is needed.

There is also an exact affine iterate formula:

```math
2^kT^k(n)=3^{s_k(n)}n+C_k(\epsilon),
\tag{4}
```

where

```math
C_0=0,\qquad
C_{j+1}=3^{\epsilon_j}C_j+\epsilon_j2^j.
\tag{5}
```

Equation (4) follows by induction.  It cleanly separates the multiplicative
parity walk from the additive \(+1\) correction.

## 2. Sheffield's maximal inequality gives a path estimate

Put

```math
R_j=2\epsilon_j-1,\qquad
M_k=\sum_{j=0}^{k-1}R_j=2s_k-k.
\tag{6}
```

Under the uniform-residue model, the \(R_j\) are independent Rademacher
variables.  Kolmogorov's maximal inequality from
[Sheffield's Lecture 10](https://math.mit.edu/~sheffield/175/Lecture10.pdf)
therefore gives, for every \(r>0\),

```math
\boxed{
\Pr\!\left(\max_{1\leq k\leq m}|2s_k-k|\geq r\right)
\leq \frac{m}{r^2}.}
\tag{7}
```

The multiplicative log surrogate is

```math
G_k=s_k\log3-k\log2,\qquad
\mu=\frac{\log3}{2}-\log2=\log\frac{\sqrt3}{2}<0.
\tag{8}
```

Since \(G_k-k\mu=(\log3)M_k/2\), (7) is equivalent to

```math
\boxed{
\Pr\!\left(
\max_{1\leq k\leq m}|G_k-k\mu|\geq u
\right)
\leq
\frac{m(\log3)^2}{4u^2}.}
\tag{9}
```

This is a genuine improvement in *type* over an endpoint estimate: it
controls every prefix up to time \(m\) simultaneously.  It is not the
sharpest numerical tail bound for Rademachers; exponential martingale or
reflection bounds are stronger for large \(u\).

The same lecture supplies a short pathwise consequence.  Apply (7) on
dyadic horizons \(m=2^a\) with threshold \(r=\eta2^a\).  The resulting
probabilities are summable.  Borel--Cantelli gives

```math
\frac{M_k}{k}\longrightarrow0,\qquad
\frac{G_k}{k}\longrightarrow\mu<0
\quad\text{almost surely}
\tag{10}
```

in the independent parity model.  The additive term in (4), and the
difference between Haar-typical \(2\)-adic seeds and a fixed positive
integer, are precisely why (10) is not a proof of the Collatz conjecture.

## 3. Vershynin's Berry--Esseen proof gives the endpoint rate

[Vershynin's friendly proof](https://arxiv.org/pdf/2602.06234) states the
non-i.i.d. Berry--Esseen theorem in the normalized form

```math
\sup_a\left|\Pr(S\leq a)-\Phi(a)\right|
\leq C\sum_j\mathbb E|X_j|^3,
\quad
\sum_j\operatorname{Var}(X_j)=1.
\tag{11}
```

Take \(X_j=R_j/\sqrt m\).  Then

```math
\sum_{j=1}^m\mathbb E|X_j|^3=\frac1{\sqrt m},
```

so

```math
\boxed{
\sup_a
\left|
\Pr\!\left(\frac{2s_m-m}{\sqrt m}\leq a\right)-\Phi(a)
\right|
\leq\frac{C}{\sqrt m}.}
\tag{12}
```

The same statement holds for the standardized \(G_m\).  Here the exact law
is binomial, so exact binomial tails or a continuity-corrected approximation
can be sharper.  The value of (12) is that the proof extends immediately to
independent non-identically distributed parity surrogates.

If \(\epsilon_j\sim\operatorname{Bernoulli}(q_j)\) independently and
\(V=\sum_jq_j(1-q_j)>0\), Vershynin's theorem gives

```math
\sup_a\left|
\Pr\!\left(
\frac{\sum_j(\epsilon_j-q_j)}{\sqrt V}\leq a
\right)-\Phi(a)
\right|
\leq
C\,
\frac{
\sum_jq_j(1-q_j)\bigl(q_j^2+(1-q_j)^2\bigr)
}{V^{3/2}}.
\tag{13}
```

This is the appropriate Berry--Esseen upgrade for biased product measures
on a hypercube.  It still requires independence; measured correlations in
an orbit cannot simply be inserted into (13).

## 4. Berry--Esseen is not a modulo-one theorem

Vershynin's proof deliberately controls characteristic functions only near
frequency zero.  That is enough for a distribution-function approximation,
but not for equidistribution modulo one.  An integer-valued sum can satisfy a
central limit theorem while remaining identically zero modulo one.

For the Collatz multiplicative surrogate one can instead compute every fixed
circle Fourier coefficient exactly.  For a base \(B>1\), put

```math
\alpha=\log_B3,\qquad\beta=\log_B2.
```

Then for every integer \(\ell\),

```math
\mathbb E\exp\!\left(
2\pi i\ell\frac{G_m}{\log B}
\right)
=
e^{-2\pi i\ell m\beta}
\left(\frac{1+e^{2\pi i\ell\alpha}}2\right)^m,
\tag{14}
```

and hence its magnitude is

```math
|\cos(\pi\ell\alpha)|^m.
\tag{15}
```

If \(\alpha\) is irrational, every fixed nonzero Fourier mode decays
exponentially, so the laws of \(G_m/\log B\bmod1\) converge weakly to the
uniform circle law.  A uniform discrepancy rate requires Diophantine control
of \(\alpha\) as \(|\ell|\) grows.  This is a Bernoulli-shift problem, not a
Gauss--Kuzmin continued-fraction problem.

For the actual Collatz map, the correction \(C_m\) in (4) must still be
controlled.  Lagarias and Soundararajan proved a related finite-horizon
Benford theorem for most initial seeds; (14) should be viewed as the exact
multiplicative model behind that direction, not as a replacement for their
argument.

## 5. Prime starting values

There are two different meanings of “prime” here.

First, suppose the *starting seed* is prime.  Fix \(m\), let \(x\to\infty\),
and choose uniformly among primes at most \(x\).  The prime number theorem in
arithmetic progressions makes the odd residue classes modulo \(2^m\)
asymptotically equiprobable.  Combining this with (3) gives

```math
(\epsilon_0,\ldots,\epsilon_{m-1})
\ \Longrightarrow\
(1,\operatorname{Bernoulli}(1/2),\ldots,
\operatorname{Bernoulli}(1/2)).
\tag{16}
```

Thus

```math
s_m\Longrightarrow1+\operatorname{Binomial}(m-1,1/2)
\tag{17}
```

for fixed \(m\).  Equations (7) and (12) apply to the remaining \(m-1\)
bits in this limit.  If \(m\) grows with \(x\), a quantitative
primes-in-progressions theorem is additionally required; fixed-modulus
Dirichlet equidistribution is not enough.

Second, an odd prime can be used as a *modulus* for the affine correction.
This produces the arithmetic diagnostic that is suitable for the
supercongruence repository.

## 6. The affine correction modulo primes

Choose the parity vector uniformly and let

```math
\nu_{m,p}(a)=\Pr(C_m\equiv a\pmod p).
\tag{18}
```

For \(e_p(x)=e^{2\pi ix/p}\), define

```math
\widehat\nu_{m,p}(t)=\mathbb E e_p(tC_m).
```

The recurrence (5) gives the exact finite-field Fourier recursion

```math
\boxed{
\widehat\nu_{j+1,p}(t)
=
\frac12\widehat\nu_{j,p}(t)
+\frac12e_p(t2^j)\widehat\nu_{j,p}(3t).}
\tag{19}
```

For \(p\geq5\), the two maps \(x\mapsto x\) and
\(x\mapsto3x+2^j\) are permutations of \(\mathbb F_p\), so each transition
is doubly stochastic.  Let \(d=\operatorname{ord}_p(2)\) and group the
transitions into blocks of length \(d\).  A block has a self-loop and
contains the one-selection maps

```math
f_0(x)=3x+1,\qquad f_1(x)=3x+2.
```

The finite permutation group generated by them contains

```math
f_1f_0^{-1}(x)=x+1,
```

so it acts transitively on \(\mathbb F_p\).  The block transition is
therefore irreducible and aperiodic.  It follows that

```math
\boxed{\nu_{m,p}\longrightarrow\operatorname{Unif}(\mathbb F_p)
\quad\text{exponentially for every }p\geq5.}
\tag{20}
```

The first persistent exceptional prime is \(p=3\), exactly where the
Collatz multiplier ceases to be invertible.  If a selected bit occurs at
time \(j\), (5) resets the state modulo \(3\) to \(2^j=(-1)^j\).  Therefore
\(C_m=0\) only for the all-zero vector, and the other two residues are
determined by the parity of the last selected position.  For even \(m\),

```math
\boxed{
\nu_{m,3}(0)=2^{-m},\quad
\nu_{m,3}(1)=\frac{1-2^{-m}}3,\quad
\nu_{m,3}(2)=\frac{2(1-2^{-m})}3.}
\tag{21}
```

For odd \(m\),

```math
\boxed{
\nu_{m,3}(0)=2^{-m},\quad
\nu_{m,3}(1)=\frac23\left(1-2^{-(m+1)}\right),\quad
\nu_{m,3}(2)=\frac13\left(1-2^{-(m-1)}\right).}
\tag{22}
```

The limiting collision probability is \(5/9\), rather than the uniform
value \(1/3\).  Thus the “interesting prime” extracted by the affine-residue
experiment is not a mysterious rare prime: it is the singular prime \(3\),
and its anomaly has an exact algebraic explanation.

## 7. What was actually improved

The combined tools produce four distinct, rigorous upgrades:

1. exact hypercube independence for every finite parity horizon;
2. a simultaneous prefix bound from Kolmogorov's maximal inequality;
3. an \(O(m^{-1/2})\) endpoint normal approximation from Berry--Esseen; and
4. exact circle and finite-field Fourier recursions, including the complete
   exceptional-prime-\(3\) law.

Only the fourth item belongs near the finite-field/entropy portion of a
supercongruence library, and even there it is a residue-distribution theorem,
not an adjacent-prime-power congruence.

## References

- D. J. Bernstein and J. C. Lagarias,
  [The \(3x+1\) conjugacy map](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/3x-1-conjugacy-map/6975BB4A8C46CF6842217043AAF9EC13).
- J. C. Lagarias and K. Soundararajan,
  [Benford's law for the \(3x+1\) function](https://arxiv.org/abs/math/0509175).
- S. Sheffield,
  [18.175 Lecture 10: zero--one laws and maximal inequalities](https://math.mit.edu/~sheffield/175/Lecture10.pdf).
- R. Vershynin,
  [A friendly proof of the Berry--Esseen theorem](https://arxiv.org/abs/2602.06234).

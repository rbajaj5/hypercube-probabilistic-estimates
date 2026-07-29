# Sharp small-deviation bounds in dimensions at most four

## 1. Setup

Let \(X_1,\ldots,X_n\) be independent nonnegative random variables with

```math
\mu_i:=\mathbb E X_i\leq1,
\qquad
S:=\sum_{i=1}^n X_i.
```

The recent bound of Fu--Han--Wang--Yan--Zhang--Zhou
([arXiv:2607.23980](https://arxiv.org/abs/2607.23980)) states

```math
\mathbb P(S<\mathbb ES+\delta)\geq b_{n,\delta},
```

where

```math
b_{n,\delta}
=
\begin{cases}
\displaystyle
\delta\left(\frac{n}{n+\delta}\right)^n,
&0<\delta<1,\\[7pt]
\displaystyle
\left(1-\frac1{n+\delta}\right)^n,
&\delta\geq1.
\end{cases}
```

The second branch is sharp for every \(n\). Nie--Wei
([arXiv:2607.24528](https://arxiv.org/abs/2607.24528)) and Stander
([doi:10.5281/zenodo.21622951](https://doi.org/10.5281/zenodo.21622951))
concurrently obtained the same first branch from the same recent
Vlassis--Thomas input. See
[the concurrent-work comparison](../CONCURRENT_WORK.md). This note records
the exact low-dimensional refinement of that branch.

## 2. Sharp low-dimensional theorem

### Theorem

For \(1\leq n\leq4\) and every \(\delta>0\),

```math
\boxed{
\mathbb P(S<\mathbb ES+\delta)
\geq
c_{n,\delta}
:=
\min\left\{
\frac{\delta}{1+\delta},
\left(1-\frac1{n+\delta}\right)^n
\right\}.}
```

The bound is sharp for every \(n\leq4\) and every \(\delta>0\).

### Mean completion

Set

```math
Y_i=X_i+1-\mu_i.
```

The \(Y_i\) remain independent and nonnegative, and
\(\mathbb EY_i=1\).  Moreover,

```math
\sum_{i=1}^nY_i=S+n-\mathbb ES,
```

so the events agree exactly:

```math
\left\{S<\mathbb ES+\delta\right\}
=
\left\{\sum_{i=1}^nY_i<n+\delta\right\}.
```

It therefore suffices to work with independent nonnegative variables of mean
one.

### Samuels and Paulin

Samuels' extremal theorem, proved for at most four summands, says that the
minimum lower-tail probability at a threshold above the total mean is attained
by one of a finite family of two-point constructions.  With all means equal to
one and threshold \(n+\delta\), its candidate values are

```math
q_{k,\delta}
=
\left(1-\frac1{k+\delta}\right)^k,
\qquad 1\leq k\leq n.
```

Roland Paulin proved that Samuels' conjecture implies Feige's conjectured
dimension-free bound.  More precisely, the proof first obtains the stronger
finite-dimensional inequality

```math
\mathbb P(Z<1+\eta)
\geq
\min\left\{
\frac{\eta}{\eta+M},
\left(1-\frac{M}{1+\eta}\right)^{1/M}
\right\},
```

for weights summing to one and largest weight \(M\), before weakening the
second term to \(e^{-1}\).  Taking all weights equal to \(1/n\) and
\(\eta=\delta/n\) reduces the minimum of the Samuels candidates to

```math
\min_{1\leq k\leq n}q_{k,\delta}
=
\min\{q_{1,\delta},q_{n,\delta}\}
=
\min\left\{
\frac{\delta}{1+\delta},
\left(1-\frac1{n+\delta}\right)^n
\right\}.
```

Samuels' theorem through \(n=4\), followed by the finite-dimensional
inequality in Paulin's proof and the mean-completion step, proves the displayed
theorem.

### Sharpness

Both endpoint candidates are realized.

1. **One active coordinate.**  Let \(Y_2=\cdots=Y_n=1\), and let

   ```math
   Y_1=
   \begin{cases}
   1+\delta,&\text{with probability }1/(1+\delta),\\
   0,&\text{with probability }\delta/(1+\delta).
   \end{cases}
   ```

   Then \(\mathbb EY_i=1\), and

   ```math
   \mathbb P\left(\sum_iY_i<n+\delta\right)
   =\frac{\delta}{1+\delta}.
   ```

2. **All coordinates active.**  Independently for every \(i\), let

   ```math
   Y_i=
   \begin{cases}
   n+\delta,&\text{with probability }1/(n+\delta),\\
   0,&\text{otherwise}.
   \end{cases}
   ```

   The strict lower-tail event occurs exactly when every coordinate is zero,
   and hence has probability

   ```math
   \left(1-\frac1{n+\delta}\right)^n.
   ```

Whichever candidate is smaller supplies a matching extremizer.

## 3. Strict improvement over the small-slack branch

Fix \(n\geq2\) and \(0<\delta<1\).  Each endpoint candidate is strictly
larger than

```math
b_{n,\delta}
=
\delta\left(\frac{n}{n+\delta}\right)^n.
```

For the one-active candidate, the strict binomial inequality gives

```math
\left(1+\frac{\delta}{n}\right)^n>1+\delta,
```

and therefore

```math
\frac{\delta}{1+\delta}
>
\delta\left(\frac{n}{n+\delta}\right)^n.
```

For the all-active candidate, strict Bernoulli inequality gives

```math
\left(1-\frac{1-\delta}{n}\right)^n>\delta.
```

Multiplication by \((n/(n+\delta))^n\) yields

```math
\left(1-\frac1{n+\delta}\right)^n
>
\delta\left(\frac{n}{n+\delta}\right)^n.
```

The minimum of the two endpoint candidates is consequently also strictly
larger.  For \(2\leq n\leq4\), this is a proved sharp improvement.  For
\(n\geq5\), the algebra still compares the candidate formulas, but the
Samuels extremal theorem needed to certify their lower-bound status is not
available in this generality.

## 4. A modern bivariate proof

Nie--Wei explicitly record that the case \(n=2\) follows directly from the
admissible two-input merger of
Ming--Ramdas--Shen--Wang--Waudby-Smith
([arXiv:2607.18661](https://arxiv.org/abs/2607.18661)).  This gives a short
proof that does not invoke Samuels' case analysis. It is included here
self-contained for comparison, with no priority claim.

Let \(F=K_2^{\mathrm{ad}}\).  For independent nonnegative \(E_1,E_2\) with
means at most one, \(F\) is a valid merger:

```math
\mathbb P\{F(E_1,E_2)\leq\alpha\}\leq\alpha,
\qquad 0\leq\alpha\leq1.
```

The rule is symmetric and coordinatewise nonincreasing.  If
\(0\leq a<1<b\), its level sets have the form

```math
F(a,b)\leq\alpha
\quad\Longleftrightarrow\quad
b\geq B_\alpha(a)
:=
\frac1{(1-s)(1+as)},
\qquad
s=\sqrt{1-\alpha}.
```

On the upper quadrant, \(F(a,b)=1/(ab)\).

Put \(t=2+\delta\) and

```math
\alpha_\delta
=
\max\left\{
F(0,t),F(1,t-1)
\right\}
=
\max\left\{
\frac2t-\frac1{t^2},
\frac1{t-1}
\right\}.
```

We claim that

```math
x+y\geq t
\quad\Longrightarrow\quad
F(x,y)\leq\alpha_\delta.
```

By monotonicity it is enough to consider \(x+y=t\), and by symmetry write
\(a=\min(x,y)\), \(b=t-a\).

If \(a\geq1\), then

```math
F(a,t-a)=\frac1{a(t-a)}\leq\frac1{t-1}\leq\alpha_\delta,
```

because \(a(t-a)\) is minimized at \(a=1\) on
\([1,t/2]\).

If \(0\leq a\leq1\), define

```math
h(a)=t-a-B_{\alpha_\delta}(a).
```

The definition of \(\alpha_\delta\) and the endpoint formulas imply
\(h(0)\geq0\) and \(h(1)\geq0\).  Also,

```math
h''(a)
=
-\frac{2s^2}{(1-s)(1+as)^3}
\leq0.
```

Thus \(h\) is concave and lies above the chord joining its endpoint values, so
\(h(a)\geq0\) throughout \([0,1]\).  The level-set characterization gives
\(F(a,t-a)\leq\alpha_\delta\), proving the claim.

Apply the claim to the mean-completed variables \(Y_1,Y_2\).  Merger validity
then gives

```math
\mathbb P(Y_1+Y_2\geq2+\delta)
\leq\alpha_\delta.
```

Taking complements,

```math
\begin{aligned}
\mathbb P(Y_1+Y_2<2+\delta)
&\geq1-\alpha_\delta\\
&=
\min\left\{
\frac{\delta}{1+\delta},
\left(\frac{1+\delta}{2+\delta}\right)^2
\right\}.
\end{aligned}
```

The two sharpness constructions in Section 2 show that this bivariate
constant is exact.

## 5. What is and is not claimed

- The sharp \(n\leq4\) result is a corollary of classical work of Samuels plus
  the stronger finite-dimensional inequality inside Paulin's proof. It is not
  claimed as a new theorem.
- The bivariate argument is explicitly recorded by Nie--Wei from the 2026
  admissible merger. It is reproduced for comparison, with no
  literature-priority claim.
- Unlike the concurrent Dirichlet-calibration proofs, the Samuels--Paulin
  proof of the \(n\leq4\) theorem does not use the Vlassis--Thomas result.
- No improvement is claimed for the general \(n\geq5\), \(0<\delta<1\)
  problem.
- A separate 2025 preprint claims the full arbitrary-\(\delta\) conjecture,
  while the July 2026 papers cited here still describe that regime as open.
  This repository does not adjudicate that claim and does not use it.

## References

1. Weibo Fu, Yanjun Han, Guanyang Wang, Jun Yan, Peng Zhang, and Zhengqing
   Zhou, “Sharp small-deviation inequalities for sums of independent
   nonnegative random variables,”
   [arXiv:2607.23980](https://arxiv.org/abs/2607.23980), 2026.
2. S. M. Samuels, “On a Chebyshev-type inequality for sums of independent
   random variables,” *Annals of Mathematical Statistics* 37 (1966),
   248–259, [doi:10.1214/aoms/1177699614](https://doi.org/10.1214/aoms/1177699614).
3. S. M. Samuels, “More on a Chebyshev-type inequality for sums of independent
   random variables,” Technical Report 155, Purdue University, 1968.
4. Uriel Feige, “On sums of independent random variables with unbounded
   variance, and estimating the average degree in a graph,” STOC 2004,
   [author manuscript](https://www.wisdom.weizmann.ac.il/~feige/Others/newmarkov.pdf).
5. Roland Paulin, “On some conjectures of Samuels and Feige,”
   [arXiv:1703.05152](https://arxiv.org/abs/1703.05152), 2017.
6. Jiahao Ming, Aaditya Ramdas, Yi Shen, Ruodu Wang, and Ian Waudby-Smith,
   “Gaffke's confidence interval for the mean of bounded data is inadmissible
   but asymptotically efficient,”
   [arXiv:2607.18661](https://arxiv.org/abs/2607.18661), 2026.
7. Zipei Nie and Jiaye Wei, “On Feige's conjecture,”
   [arXiv:2607.24528](https://arxiv.org/abs/2607.24528), 2026.
8. Mark Stander, “A Conditional Proof of Feige's Conjecture with a Sharp
   Finite-Dimensional Bound,”
   [doi:10.5281/zenodo.21622951](https://doi.org/10.5281/zenodo.21622951),
   2026.

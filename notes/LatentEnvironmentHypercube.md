# Latent environments on a hypercube

## Scope and source boundary

The supplied Chinese transcript of a reading of Jared Diamond's preface
emphasizes a methodological question: how much of an observed outcome should
be attributed to an individual unit, and how much to the environment and path
that produced it?  This note extracts one precise probability model from that
question.  It neither tests nor endorses any historical claim in *Guns, Germs,
and Steel*.

The mathematical point is narrower and reusable.  Independent coordinates
conditional on a shared environment are generally dependent after that
environment is hidden.  This produces exact overdispersion and higher-order
Fourier coefficients that an independent Bernoulli null model misses.

## A two-environment model

Let \(Z\in\{0,1\}\), with

\[
\Pr(Z=1)=\alpha,\qquad 0<\alpha<1.
\]

Choose two rates \(0<u,v<1\), and set

\[
\Theta =
\begin{cases}
u,&Z=1,\\
v,&Z=0.
\end{cases}
\]

Conditional on \(\Theta\), let the coordinates of
\(X=(X_1,\ldots,X_d)\in\{0,1\}^d\) be independent Bernoulli variables with
common parameter \(\Theta\).  Thus

\[
\Pr(X=x)
=
\alpha u^{|x|}(1-u)^{d-|x|}
+(1-\alpha)v^{|x|}(1-v)^{d-|x|}.
\tag{1}
\]

This is a finite, exactly computable latent-environment law on the Boolean
cube.  It is exchangeable but, unless \(u=v\), it is not a product measure.

## The hidden environment creates correlation

Put

\[
p=\mathbb E\Theta=\alpha u+(1-\alpha)v,\qquad
\sigma_\Theta^2=\operatorname{Var}(\Theta)
=\alpha(1-\alpha)(u-v)^2.
\]

For every coordinate,

\[
\mathbb E X_i=p,\qquad \operatorname{Var}(X_i)=p(1-p).
\]

For distinct coordinates \(i\ne j\), conditional independence gives

\[
\operatorname{Cov}(X_i,X_j)
=\operatorname{Var}(\mathbb E[X_i\mid\Theta])
=\sigma_\Theta^2.
\tag{2}
\]

Consequently,

\[
\operatorname{Corr}(X_i,X_j)
=\frac{\alpha(1-\alpha)(u-v)^2}{p(1-p)}.
\tag{3}
\]

So units with the same one-coordinate marginal can be positively correlated
solely because they share an unobserved environment.

## Exact count law and overdispersion

Let \(K=|X|=\sum_i X_i\).  Its probability generating function is

\[
\mathbb E z^K
=
\alpha(1-u+uz)^d+(1-\alpha)(1-v+vz)^d,
\tag{4}
\]

and hence

\[
\Pr(K=k)
=
\binom dk\left[
\alpha u^k(1-u)^{d-k}
+(1-\alpha)v^k(1-v)^{d-k}
\right].
\tag{5}
\]

The mean agrees with the homogeneous binomial model
\(\operatorname{Bin}(d,p)\):

\[
\mathbb E K=dp.
\]

The variance does not:

\[
\operatorname{Var}(K)
=d\,p(1-p)+d(d-1)\sigma_\Theta^2.
\tag{6}
\]

The second term is the exact excess variance caused by the latent environment.
It grows quadratically in \(d\) when the environmental contrast remains fixed.
This is a concrete heavy-tail warning: a binomial fit can reproduce the mean
while severely understating dispersion.  Equation (6), not a universal
pointwise tail ordering, is the claim.

More generally, with \((a)_r=a(a-1)\cdots(a-r+1)\),

\[
\mathbb E (K)_r
=(d)_r\left[\alpha u^r+(1-\alpha)v^r\right].
\tag{7}
\]

These factorial moments identify all coefficients of (4).

## Centered Fourier interactions

The natural product-measure null with the same marginals is
\(\operatorname{Bernoulli}(p)^{\otimes d}\).  Under that null, every centered
mixed moment on a nonempty set of distinct coordinates vanishes.

In the environment mixture, for \(S\subseteq[d]\) and \(r=|S|\),

\[
\begin{aligned}
\mathbb E\prod_{i\in S}(X_i-p)
&=\mathbb E(\Theta-p)^r\\
&=(u-v)^r\left[
\alpha(1-\alpha)^r+(1-\alpha)(-\alpha)^r
\right].
\end{aligned}
\tag{8}
\]

The degree-one value is zero, the degree-two value is
\(\sigma_\Theta^2\), and higher degrees are generally nonzero.  Thus a shared
environment leaves a complete hierarchy of centered Boolean/Fourier
interactions, even though no interaction exists after conditioning on
\(\Theta\).  Checking only coordinate means cannot see this effect.

## Inferring the environment

The count \(K\) is sufficient for \(Z\).  Bayes' rule gives

\[
\Pr(Z=1\mid K=k)
=
\frac{\alpha u^k(1-u)^{d-k}}
{\alpha u^k(1-u)^{d-k}
+(1-\alpha)v^k(1-v)^{d-k}}.
\tag{9}
\]

If \(u>v\), the posterior odds are multiplied by

\[
\frac{u(1-v)}{v(1-u)}>1
\tag{10}
\]

whenever \(k\) increases by one.  A large count is therefore evidence for the
high-rate environment in this model, but it is not evidence for an intrinsic
difference between coordinates.

## Entropy accounting

Writing \(h\) for binary entropy,

\[
H(X)
=h(\alpha)
+d\bigl[\alpha h(u)+(1-\alpha)h(v)\bigr]
-H(Z\mid X).
\tag{11}
\]

It follows that

\[
d\bigl[\alpha h(u)+(1-\alpha)h(v)\bigr]
\le H(X)
\le
d\bigl[\alpha h(u)+(1-\alpha)h(v)\bigr]+h(\alpha).
\tag{12}
\]

The environment adds at most one binary variable's worth of uncertainty, but
it can change count variance by order \(d^2\).  Entropy and tail dispersion
therefore answer different questions.

## What is verified

`verification/verify_latent_environment_hypercube.py` uses exact rational
arithmetic to enumerate several finite cubes and checks:

1. normalization and exchangeability of (1);
2. every marginal, covariance, and centered interaction in (2) and (8);
3. the complete count distribution (5);
4. the mean, variance, and all factorial moments (6)--(7);
5. the posterior formula and monotone odds multiplier (9)--(10).

The checker proves the finite identities implemented there.  It does not test
historical causation, recover a real latent environment from observations, or
claim that every empirical heavy tail is an environment mixture.

## Relation to the other modules

- [Novelty perturbations](NoveltyRefreshHypercube.md) use independent,
  coordinate-specific stationary rates.  This note explains what changes when
  those rates share a hidden common cause.
- [Heavy-tail projector mixtures](ProjectorHeavyTailMixtures.md) mix a global
  scale.  Here the same mixture principle acts on Bernoulli coordinates and
  yields an exact beta-binomial-style overdispersion term.
- [Hypercube walks and matrix mixing](HypercubeWalkMatrixMixing.md) describe
  dynamical relaxation.  The present law is static; adding environment
  switching would require a larger Markov state space.

# Rare-event mass separation on a hypercube

## Scope and source boundary

[Li--Xia](https://arxiv.org/abs/2607.26549) construct a plurisubharmonic
singularity with unit residual Monge--Ampère mass and zero Lelong number.  At
finite stage, their normalization \(2^{-j}\) compensates for mapping degree
\(4^j\): the mass remains one while a first-order vanishing-order diagnostic
equals \(2^{-j}\).

That complex-geometric theorem does not transfer to a Boolean cube.  This
note records a deliberately elementary probability analogue: normalized mass
can remain one while its support and the Boolean boundary of its support both
vanish.  The analogy is only a warning about diagnostics; no Monge--Ampère,
Lelong-number, or pluripotential-theory claim is made here.  It is also
unrelated to the Bieberbach coefficient theorem.

## A normalized rare subcube

Let \(U\) be uniform on \(\Omega=\{0,1\}^d\).  Fix \(R\subseteq[d]\), with
\(|R|=r\), and a bit pattern \(a\in\{0,1\}^R\).  Define the codimension-\(r\)
subcube

\[
A=\{x:x_i=a_i\text{ for every }i\in R\},\qquad U(A)=2^{-r},
\tag{1}
\]

and the normalized density

\[
Z(x)=2^r\mathbf 1_A(x).
\tag{2}
\]

Then

\[
\mathbb E_U Z=1,\qquad
U(Z\ne0)=2^{-r},\qquad
\mathbb E_U Z^q=2^{r(q-1)}\quad(q>0).
\tag{3}
\]

In particular,

\[
\operatorname{Var}_U(Z)=2^r-1.
\tag{4}
\]

The total mass is fixed, but it is carried by an exponentially rare event
with an exponentially large density.  This is the finite-cube black-swan
mechanism in its simplest exact form.

Let \(\nu=Z\,U\).  It is uniform on \(A\).  Consequently,

\[
\operatorname{TV}(\nu,U)=1-2^{-r},
\qquad
H(\nu)=(d-r)\log2,
\tag{5}
\]

and every Rényi divergence of positive order, including relative entropy, is

\[
D_q(\nu\|U)=r\log2.
\tag{6}
\]

Thus keeping normalized mass fixed does not keep the tilted law close to the
reference law.

## A vanishing first-order boundary diagnostic

Put \(h=\mathbf1_A\).  With

\[
\operatorname{Inf}_i(h)
=\Pr_U\bigl(h(X)\ne h(X\oplus e_i)\bigr),
\]

direct enumeration gives

\[
\operatorname{Inf}_i(h)
=
\begin{cases}
2^{-(r-1)},&i\in R,\\
0,&i\notin R,
\end{cases}
\qquad
I(h)=r\,2^{1-r}.
\tag{7}
\]

Hence \(I(h)\to0\) although \(\mathbb E Z=1\).  There is no contradiction:
influence is being measured on the unscaled support indicator \(h\), whereas
the density height \(2^r\) diverges.  A small boundary for a rare event cannot
by itself control a normalized density supported on that event.

## Exact Fourier and noise laws

For \(\chi_S(x)=(-1)^{\sum_{i\in S}x_i}\), equation (2) factors as

\[
Z(x)=\prod_{i\in R}\left(1+(-1)^{a_i}\chi_{\{i\}}(x)\right).
\tag{8}
\]

Therefore

\[
\widehat Z(S)
=
\begin{cases}
(-1)^{\sum_{i\in S}a_i},&S\subseteq R,\\
0,&S\not\subseteq R.
\end{cases}
\tag{9}
\]

There are \(2^r\) nonzero coefficients, all of magnitude one, and Parseval
recovers \(\mathbb E Z^2=2^r\).

Let \(Y\) be obtained from \(X\sim U\) by flipping each bit independently
with probability \(\eta\).  The Walsh noise eigenvalue is
\(\theta=1-2\eta\), so

\[
\mathbb E[Z(X)Z(Y)]
=\sum_{S\subseteq R}\theta^{|S|}
=(1+\theta)^r
=\bigl(2(1-\eta)\bigr)^r.
\tag{10}
\]

Equivalently,

\[
\Pr(X\in A,Y\in A)
=\left(\frac{1-\eta}{2}\right)^r,
\qquad
\Pr(Y\in A\mid X\in A)=(1-\eta)^r.
\tag{11}
\]

The rare mass is extremely noise-sensitive once \(r\eta\) is appreciable.

## Sampling cost and the \(4^j\) dictionary

After \(N\) independent uniform samples,

\[
\Pr(\text{the support }A\text{ is never observed})
=(1-2^{-r})^N.
\tag{12}
\]

Thus observing the black-swan support requires the exponential scale \(2^r\).
This is the single-set specialization of the missed-mass estimates in
[State coverage and coordination](StateCoverageCoordinationHypercube.md).

If \(r=2j\), then

\[
|\,\{0,1\}^r\,|=4^j,\qquad
U(A)=4^{-j},\qquad
\sqrt{U(A)}=2^{-j}.
\tag{13}
\]

These are exactly the numerical degree and normalization scales appearing in
the finite stages of Li--Xia.  Equation (13) is only a scale dictionary:
collapsing cube atoms is not a model of their holomorphic maps, and Boolean
influence is not a Lelong number.

## What is verified

`verification/verify_rare_event_mass_separation.py` uses exact rational
enumeration to check:

1. normalization, support mass, moments, variance, and total variation;
2. support size and the entropy-effective cardinality;
3. every coordinate influence in (7);
4. every Walsh coefficient in (9) and Parseval's identity;
5. the joint and conditional noise laws (10)--(11); and
6. the sampling law and the \(r=2j\) scale dictionary.

The results are elementary and no literature-priority claim is made.

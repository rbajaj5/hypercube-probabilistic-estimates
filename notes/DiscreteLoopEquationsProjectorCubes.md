# Exact discrete loop equations for projector hypercubes

## 1. Source and boundary

Bourgade and Huang's
[*Loop Equations Characterize Random Matrix
Statistics*](https://arxiv.org/abs/2607.07617) derives microscopic loop
equations for Wigner matrices by cumulant expansion and for random regular
graphs by a discrete switching integration-by-parts formula followed by a
finite-rank Woodbury expansion.  The supplied excerpt contains precisely
those two mechanisms.

For the finite sign cubes already used in this repository, the switching
idea has a particularly clean form: pair the two vertices that differ in
one sign and retain the complete resolvent difference.  This gives an
**exact finite-\(d\)** loop equation with no expansion remainder.

That exactness does not imply Sine-\(\beta\) or Airy-\(\beta\) universality.
The Bourgade--Huang characterization assumes a particle-generated
Nevanlinna function satisfying the full closed microscopic hierarchy (and
currently rational \(\beta>0\)).  Its ensemble convergence applications
then verify approximate hierarchies using local laws and a large-matrix
limit.  A fixed projector cube generally has correlated matrix entries and
supplies none of those conclusions by itself.

## 2. Boolean integration by parts

Let \(\epsilon=(\epsilon_1,\ldots,\epsilon_d)\) be uniform on
\(\{-1,1\}^d\).  For any scalar- or matrix-valued function \(F\),

```math
\boxed{
\mathbb E[\epsilon_iF(\epsilon)]
=
\frac12\,
\mathbb E_{\epsilon_{-i}}
\left[
F(\epsilon_{-i},+1)-F(\epsilon_{-i},-1)
\right].}
\tag{1}
```

This is the exact Rademacher analogue of integration by parts.  Unlike a
cumulant expansion, (1) is just a pairing of cube vertices and has no
remainder.

Let \(P_1,\ldots,P_d\) be real symmetric matrices and define

```math
A(\epsilon)=\sum_{j=1}^d\epsilon_jP_j,
\qquad
R_\epsilon(z)=(A(\epsilon)-zI)^{-1}.
\tag{2}
```

Take \(z\in\mathbb C\) outside the spectra of all cube vertices.  After
deleting coordinate \(i\), write

```math
A_i=\sum_{j\ne i}\epsilon_jP_j,
\qquad
R_{i,\pm}=(A_i\pm P_i-zI)^{-1}.
\tag{3}
```

The resolvent identity, equivalently the full finite-rank Woodbury
difference, is

```math
R_{i,+}-R_{i,-}
=-2R_{i,+}P_iR_{i,-}.
\tag{4}
```

Combining (1) and (4) gives the coordinate identity

```math
\boxed{
\mathbb E[\epsilon_iR_\epsilon(z)]
=
-\mathbb E_{\epsilon_{-i}}
\left[R_{i,+}P_iR_{i,-}\right].}
\tag{5}
```

No smallness or orthogonality assumption is used.

## 3. Exact matrix loop equation

Since

```math
A(\epsilon)R_\epsilon(z)=I+zR_\epsilon(z),
\tag{6}
```

left-multiplying (5) by \(P_i\) and summing over \(i\) proves

```math
\boxed{
I+z\,\mathbb E R_\epsilon(z)
=
-\sum_{i=1}^d
P_i\,
\mathbb E_{\epsilon_{-i}}
\left[R_{i,+}P_iR_{i,-}\right].}
\tag{7}
```

For the normalized Stieltjes transform

```math
m_d(z)=\frac1D\operatorname{tr}R_\epsilon(z),
\tag{8}
```

where \(D\) is the matrix dimension, taking traces gives

```math
\boxed{
1+z\,\mathbb E m_d(z)
=
-\frac1D\sum_{i=1}^d
\mathbb E_{\epsilon_{-i}}
\operatorname{tr}
\left(P_iR_{i,+}P_iR_{i,-}\right).}
\tag{9}
```

Equation (7) is a finite projector-cube loop equation.  It remains valid
for strongly overlapping, perturbed, or almost-orthogonal matrices.  Those
geometries change the cavity resolvents on the right, not the identity.

## 4. Rank-one cavity closure

Now take the repository's projector model

```math
P_i=u_iu_i^{\mathsf T},\qquad \|u_i\|_2=1,
\tag{10}
```

and put

```math
R_i^{(0)}=(A_i-zI)^{-1},
\qquad
h_i=u_i^{\mathsf T}R_i^{(0)}u_i.
\tag{11}
```

The rank-one Sherman--Morrison formula gives

```math
u_i^{\mathsf T}R_{i,+}u_i=\frac{h_i}{1+h_i},
\qquad
u_i^{\mathsf T}R_{i,-}u_i=\frac{h_i}{1-h_i}.
\tag{12}
```

Consequently,

```math
\operatorname{tr}(P_iR_{i,+}P_iR_{i,-})
=
\frac{h_i^2}{1-h_i^2},
\tag{13}
```

and (9) becomes the exact scalar cavity equation

```math
\boxed{
1+z\,\mathbb E m_d(z)
=
-\frac1D\sum_{i=1}^d
\mathbb E_{\epsilon_{-i}}
\frac{h_i^2}{1-h_i^2}.}
\tag{14}
```

This is not a closed equation for \(\mathbb Em_d\): the cavity variables
\(h_i\) retain the frame geometry.  It is nevertheless useful because both
sides can be estimated independently in a simulation.  Their discrepancy
is a direct implementation or sampling diagnostic.

## 5. Perturbations and almost orthogonality

Let

```math
G_{ij}=\langle u_i,u_j\rangle,
\qquad
B_{ij}=G_{ij}^2=\operatorname{tr}(P_iP_j).
\tag{15}
```

The earlier matrix-mixing note controls the quadratic statistic
\(\operatorname{tr}A^2\) using \(B_{ij}\).  Equation (14) goes further:
each \(h_i\) resums walks that repeatedly revisit overlapping projectors.
Thus changing an almost-orthogonal frame perturbs the full resolvent
density, not merely its second moment.

There are two sound uses:

1. compare the empirical Stieltjes transform against the independent
   cavity side of (14); and
2. expand \(h_i\) for large \(|z|\) to organize higher closed-walk moments.

There is also a firm limitation.  Almost orthogonality alone does not make
the \(h_i\) independent or produce a GOE local law.  Any claimed local
universality still needs a separate large-\(D,d\) theorem controlling the
whole loop hierarchy.

## 6. The branch condition is part of the data

The second supplied excerpt constructs solutions attached to different
bulk and edge branch factors, then uses sectorial rays, Nevanlinna
properties, asymptotic normalization, and a Volterra uniqueness argument to
select the desired solution.  This matters here: an algebraic loop equation
without an analytic branch condition need not determine a probability law.

The projector cube has its own canonical branch without any asymptotic
analysis.  Average the empirical spectral measures of all vertices:

```math
\nu_d
=
\frac{1}{2^dD}
\sum_{\epsilon\in\{-1,1\}^d}
\sum_{k=1}^D
\delta_{\lambda_k(A(\epsilon))}.
\tag{16}
```

Then

```math
\overline m_d(z)
:=\mathbb E m_d(z)
=
\int_{\mathbb R}\frac{1}{x-z}\,d\nu_d(x).
\tag{17}
```

Because every \(A(\epsilon)\) is real symmetric, this transform satisfies

```math
\operatorname{Im}\overline m_d(z)
=
\operatorname{Im}z
\int_{\mathbb R}\frac{d\nu_d(x)}{|x-z|^2}>0
\quad(\operatorname{Im}z>0),
\tag{18}
```

and \(-z\overline m_d(z)\to1\) as \(|z|\to\infty\) away from the
real axis.  Thus it is the probability-normalized Nevanlinna branch.
Stieltjes inversion, or simply the poles and residues in this finite case,
recovers \(\nu_d\).

There is a technical but decisive distinction from the paper's hypotheses.
Its "particle-generated" functions have a configuration as representing
measure, hence integer atomic masses.  The averaged measure (16) usually
has fractional masses \(1/(2^dD)\).  Therefore \(\overline m_d\) is a
Nevanlinna transform of a probability measure, but it is not
particle-generated in that specialized sense and does not satisfy the
entry condition of the Sine/Airy characterization theorem.

For rank-one projectors, \(\|A(\epsilon)\|_{\mathrm{op}}\leq d\).  Moreover
\(A(-\epsilon)=-A(\epsilon)\), so \(\nu_d\) is symmetric and
\(\overline m_d(-z)=-\overline m_d(z)\).  For \(|z|>d\), the exact
convergent expansion is

```math
\overline m_d(z)
=
-\frac1z-\frac{\mu_2}{z^3}-\frac{\mu_4}{z^5}-\cdots,
\qquad
\mu_k=\int x^k\,d\nu_d(x),
\tag{19}
```

with

```math
\boxed{
\mu_2=\frac dD,\qquad
\mu_4=
\frac1D\left[
d+\sum_{i<j}\left(4B_{ij}+2B_{ij}^2\right)
\right].}
\tag{20}
```

The second moment is blind to the frame geometry, but the fourth moment is
not.  Equation (20) is therefore the first large-\(|z|\) density statistic
that detects an almost-orthogonal perturbation.  It is a concrete upgrade
over using only the quadratic-energy mean, while still being far weaker
than a microscopic local law.

## 7. Exact finite-difference loop hierarchy

The third supplied excerpt differentiates products of entrywise resolvents
and Stieltjes transforms.  Its factor \(1+\delta_{ij}\) accounts for the
chosen parametrization of a symmetric matrix entry.  A sign-cube coordinate
is already one independent variable, so its corresponding normalization is
the vertex-pair identity (1).

The full analogue is not limited to the one-point equation.  For a scalar
observable \(\Phi(\epsilon)\), define on the \(i\)-th vertex pair

```math
\overline R_i=\frac{R_{i,+}+R_{i,-}}2,\qquad
\Delta_i\Phi=\frac{\Phi_{i,+}-\Phi_{i,-}}2,\qquad
\overline\Phi_i=\frac{\Phi_{i,+}+\Phi_{i,-}}2.
\tag{21}
```

The exact discrete product rule and (4) give

```math
\boxed{
\mathbb E[\epsilon_iR_\epsilon\Phi]
=
\mathbb E_{\epsilon_{-i}}
\left[
\overline R_i\,\Delta_i\Phi
-R_{i,+}P_iR_{i,-}\,\overline\Phi_i
\right].}
\tag{22}
```

Consequently the exact matrix hierarchy is

```math
\boxed{
\mathbb E[(I+zR_\epsilon(z))\Phi]
=
\sum_{i=1}^d P_i\,
\mathbb E_{\epsilon_{-i}}
\left[
\overline R_i(z)\,\Delta_i\Phi
-R_{i,+}(z)P_iR_{i,-}(z)\,\overline\Phi_i
\right].}
\tag{23}
```

Taking

```math
\Phi(\epsilon)=\prod_{a=2}^p m_d(z_a)
\tag{24}
```

produces a finite \(p\)-point **cavity** hierarchy.  It is not the closed
\(\beta\)-dependent hierarchy in the paper: its right side still contains
coordinate-deleted resolvents.  It is important not to silently replace
the finite difference by a first derivative.  If
\(m_{a,\pm}=\overline m_a\pm\Delta m_a\), then

```math
\Delta_i\Phi
=
\sum_{\substack{S\subseteq\{2,\ldots,p\}\\|S|\ {\rm odd}}}
\prod_{a\in S}\Delta m_a
\prod_{a\notin S}\overline m_a.
\tag{25}
```

The terms with three, five, and more differences are genuine discrete
corrections.  They become negligible only after a separate small-switching
or large-system estimate.  Keeping them makes (23) exact and also preserves
the noncommutative ordering \(R_{i,+}P_iR_{i,-}\).

This hierarchy is the closest finite-cube counterpart of the attached
entrywise calculation.  What remains absent is precisely what drives the
paper's universality theorem: a microscopic scaling limit, local law, Ward
estimates, and control of the hierarchy as dimension grows.

## 8. What replaces the switching remainder

The fourth supplied excerpt bounds a replacement error of order
\(\Lambda^3/\sqrt d\) after a graph switching is linearized.  No analogous
analytic error is needed for (7) or (23), because the complete finite
difference is retained.  If the cube is sampled instead of enumerated,
however, there is a statistical error.

For one sampled vertex, define the paired loop residual

```math
\mathcal X_z(\epsilon)
=
I+zR_\epsilon(z)
+\sum_{i=1}^d
P_iR_{i,+}(\epsilon_{-i};z)P_iR_{i,-}(\epsilon_{-i};z).
\tag{26}
```

Equation (7) says exactly that

```math
\mathbb E\mathcal X_z=0.
\tag{27}
```

Let

```math
\gamma(z)
=
\min_\epsilon\operatorname{dist}
\left(z,\operatorname{spec}A(\epsilon)\right)>0.
```

Since \(\|A(\epsilon)\|_{\rm op}\leq d\) and every resolvent in (26) has
norm at most \(\gamma(z)^{-1}\),

```math
\|\mathcal X_z(\epsilon)\|_{\rm F}
\leq
\sqrt D\,d\left(\gamma(z)^{-1}+\gamma(z)^{-2}\right).
\tag{28}
```

For \(K\) independent cube samples and
\(\overline{\mathcal X}_{z,K}=K^{-1}\sum_{k=1}^K\mathcal X_z(\epsilon^{(k)})\),
the zero mean gives

```math
\mathbb E\|\overline{\mathcal X}_{z,K}\|_{\rm F}^2
=
\frac1K\mathbb E\|\mathcal X_z\|_{\rm F}^2
\leq
\frac{Dd^2}{K}
\left(\gamma(z)^{-1}+\gamma(z)^{-2}\right)^2.
\tag{29}
```

In particular, Markov's inequality supplies a confidence-\(1-\delta\)
certificate by multiplying the root-mean-square bound by
\(\delta^{-1/2}\).  For real \(|z|>d\), one may use the explicit lower
bound \(\gamma(z)\geq |z|-d\).

This bound is intentionally conservative.  It identifies the correct
remaining problem: sampling variance becomes poor close to the spectral
support.  Improving it would require variance reduction or a genuine local
law.  The exponential-observable and Calogero--Moser--Sutherland equations
later in the excerpt also require the specific limiting
Sine-\(\beta\)/Airy-\(\beta\) hierarchy; the finite projector identity alone
does not justify importing those differential equations.

## Verification

Run

```text
python verification/verify_discrete_loop_projector_cubes.py
```

The checker uses exact rational arithmetic.  It exhausts several one-,
two-, and three-dimensional rational projector frames; verifies the
coordinate integration-by-parts formula, the matrix and trace loop
equations, the rank-one cavity reduction, and the determinant form of
Sherman--Morrison at both positive and negative spectral parameters.  It
also verifies the transform symmetry and the second- and fourth-moment
coefficients in (19)--(20), plus the three-observable hierarchy (23)--(25)
for noncommuting rational frames.  Finally, it exhausts the paired
residuals in (26), checks their exact zero mean, the deterministic bound
(28), and the \(K=2\) mean-square identity behind (29).
